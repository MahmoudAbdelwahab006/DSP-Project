%% =========================================================
%  PROJECT II: Multi-Band Speech Equalizer for Podcast Enhancement
%  DSP Course Project – Zewail City of Science and Technology
% =========================================================
%  HOW TO RUN:
%  Just press Run (F5). The script will prompt you interactively.
%  A sample speech file is auto-generated if none is provided.
% =========================================================

clc; clear; close all;

%% =========================================================
%  STEP 0: USER INPUTS (Interactive prompts)
% =========================================================

fprintf('========================================\n');
fprintf(' Multi-Band Speech Equalizer\n');
fprintf('========================================\n\n');

% --- Audio file ---
audio_file = input('Enter audio filename (leave blank to use built-in test signal): ','s');
if isempty(audio_file)
    fprintf('No file provided. Generating synthetic speech signal...\n');
    use_synthetic = true;
else
    use_synthetic = false;
end

% --- Filter type ---
fprintf('\nFilter type:\n  1 = FIR\n  2 = IIR\n');
ft_choice = input('Select (1 or 2): ');
if ft_choice == 1
    filter_type = 'FIR';
    fprintf('\nFIR Window type:\n  1=Hamming  2=Hanning  3=Blackman\n');
    wt = input('Select window (1-3): ');
    switch wt
        case 1, win_type = 'hamming';
        case 2, win_type = 'hanning';
        otherwise, win_type = 'blackman';
    end
    order_in = input('FIR filter order (leave blank for default 100): ');
    if isempty(order_in), fir_order = 100; else, fir_order = order_in; end
    if mod(fir_order,2)~=0, fir_order = fir_order+1; end % must be even
else
    filter_type = 'IIR';
    fprintf('\nIIR type:\n  1=Butterworth  2=Chebyshev I  3=Chebyshev II\n');
    it = input('Select IIR type (1-3): ');
    switch it
        case 1, iir_type = 'butter';
        case 2, iir_type = 'cheby1';
        otherwise, iir_type = 'cheby2';
    end
    order_in = input('IIR filter order (leave blank for default 4): ');
    if isempty(order_in), iir_order = 4; else, iir_order = order_in; end
end

% --- Mode ---
fprintf('\nMode:\n  1 = Preset (7 speech-optimized bands)\n  2 = Custom\n');
mode_choice = input('Select (1 or 2): ');

if mode_choice == 1
    % Preset bands (Hz)
    band_edges = [0, 100, 300, 800, 2000, 4000, 7000, 7900];
    n_bands = 7;
    band_names = {'0-100','100-300','300-800','800-2k','2-5k','5-10k','10-20k'};
    fprintf('\nDefault gains are 0 dB. Enter custom gains per band:\n');
    gains_dB = zeros(1, n_bands);
    for k = 1:n_bands
        g = input(sprintf('  Band %d (%s Hz) gain [dB, default 0]: ', k, band_names{k}));
        if ~isempty(g), gains_dB(k) = g; end
    end
else
    n_bands = input('\nNumber of bands (5-10): ');
    n_bands = max(5, min(10, n_bands));
    fprintf('Enter %d band edge frequencies (must start at 0, end at 20000):\n', n_bands+1);
    band_edges = zeros(1, n_bands+1);
    band_edges(1) = 0; band_edges(end) = 20000;
    for k = 2:n_bands
        band_edges(k) = input(sprintf('  Edge %d (Hz): ', k));
    end
    band_edges = sort(band_edges);
    gains_dB = zeros(1, n_bands);
    for k = 1:n_bands
        g = input(sprintf('  Band %d (%.0f-%.0f Hz) gain [dB, default 0]: ', ...
            k, band_edges(k), band_edges(k+1)));
        if ~isempty(g), gains_dB(k) = g; end
    end
end

% --- Output sample rate ---
fprintf('\nOutput sample rates:\n');
fs_out_choice = input('Enter output sample rate multiplier (1=original, 4=4x, 0.5=half): ');

%% =========================================================
%  STEP 1: LOAD / GENERATE AUDIO SIGNAL
% =========================================================

if use_synthetic
    % Generate realistic synthetic speech (voiced + unvoiced components)
    fs = 16000;   % 16 kHz (standard speech)
    t_dur = 5;    % 5 seconds
    t = (0:1/fs:t_dur-1/fs)';
    % Voiced speech: sum of harmonics with fundamental ~120 Hz
    f0 = 120;
    speech = zeros(size(t));
    for h = 1:20
        speech = speech + (1/h) * sin(2*pi*h*f0*t + rand*2*pi);
    end
    % Unvoiced (fricative): bandpass noise 2-8 kHz
    rng(1);
    noise_uv = randn(size(t));
    high_cutoff = min(7900, fs/2 - 100);

    [b_uv,a_uv] = butter(4, [2000 high_cutoff]/(fs/2), 'bandpass');

    unvoiced = filter(b_uv,a_uv,noise_uv);
    % Modulate with envelope
    env = 0.5*(1+sin(2*pi*0.5*t));
    audio_in = speech.*env + 0.3*unvoiced.*env;
    audio_in = audio_in / max(abs(audio_in));
    fprintf('Synthetic speech generated: %.1f s at %d Hz\n', t_dur, fs);
else
    [audio_in, fs] = audioread(audio_file);
    if size(audio_in,2) > 1
        audio_in = mean(audio_in,2);  % Mono
    end
    fprintf('Loaded: %s  |  fs=%d Hz  |  %.2f s\n', audio_file, fs, length(audio_in)/fs);
end

% Ensure column vector
audio_in = audio_in(:);
N_audio = length(audio_in);
t_audio = (0:N_audio-1)/fs;

%% =========================================================
%  STEP 2 & 3: DESIGN BAND FILTERS + APPLY GAIN
% =========================================================

gains_lin = 10.^(gains_dB/20);  % Convert dB to linear
audio_out = zeros(N_audio, 1);
filter_responses = cell(n_bands, 1);

fprintf('\nDesigning and applying %d band filters (%s)...\n', n_bands, filter_type);

for k = 1:n_bands
    f_lo = band_edges(k);
    f_hi = band_edges(k+1);
    nyq  = fs/2;

    % Clamp to valid normalized range
    f_lo_n = max(f_lo, 1) / nyq;
    f_hi_n = min(f_hi, nyq*0.999) / nyq;

    % Extra protection
    f_lo_n = min(f_lo_n, 0.999);
    f_hi_n = min(f_hi_n, 0.999);

    % ---- DESIGN FILTER ----
    if strcmp(filter_type, 'FIR')
        switch win_type
            case 'hamming',  win_fn = hamming(fir_order+1);
            case 'hanning',  win_fn = hanning(fir_order+1);
            otherwise,       win_fn = blackman(fir_order+1);
        end
        if f_lo < 1  % First band: lowpass
            b_k = fir1(fir_order, f_hi_n, 'low', win_fn);
        elseif f_hi >= nyq*0.999  % Last band: highpass
            b_k = fir1(fir_order, f_lo_n, 'high', win_fn);
        else
            b_k = fir1(fir_order, [f_lo_n, f_hi_n], 'bandpass', win_fn);
        end
        a_k = 1;
        filter_responses{k} = struct('b',b_k,'a',a_k,'type',filter_type,...
            'order',fir_order,'window',win_type);
    else  % IIR
        Rp = 1; Rs = 40;
        try
            if f_lo < 1
                switch iir_type
                    case 'butter',  [b_k,a_k] = butter(iir_order, f_hi_n, 'low');
                    case 'cheby1',  [b_k,a_k] = cheby1(iir_order, Rp, f_hi_n, 'low');
                    case 'cheby2',  [b_k,a_k] = cheby2(iir_order, Rs, f_hi_n, 'low');
                end
            elseif f_hi >= nyq*0.999
                switch iir_type
                    case 'butter',  [b_k,a_k] = butter(iir_order, f_lo_n, 'high');
                    case 'cheby1',  [b_k,a_k] = cheby1(iir_order, Rp, f_lo_n, 'high');
                    case 'cheby2',  [b_k,a_k] = cheby2(iir_order, Rs, f_lo_n, 'high');
                end
            else
                switch iir_type
                    case 'butter',  [b_k,a_k] = butter(iir_order, [f_lo_n,f_hi_n], 'bandpass');
                    case 'cheby1',  [b_k,a_k] = cheby1(iir_order, Rp, [f_lo_n,f_hi_n], 'bandpass');
                    case 'cheby2',  [b_k,a_k] = cheby2(iir_order, Rs, [f_lo_n,f_hi_n], 'bandpass');
                end
            end
        catch
            warning('Band %d filter design failed. Using identity.',k);
            b_k = 1; a_k = 1;
        end
        filter_responses{k} = struct('b',b_k,'a',a_k,'type',filter_type,...
            'order',iir_order,'iir_type',iir_type);
    end

    % ---- FILTER AND APPLY GAIN ----
    band_signal = filtfilt(b_k, a_k, audio_in);
    audio_out   = audio_out + gains_lin(k) * band_signal;

    fprintf('  Band %d [%.0f-%.0f Hz] -> gain=%.1f dB  OK\n', ...
        k, f_lo, f_hi, gains_dB(k));
end

% Normalize output
audio_out = audio_out / max(abs(audio_out) + 1e-8);

%% =========================================================
%  STEP 4: FILTER ANALYSIS PLOTS (per band)
% =========================================================
Nfft = 2048;
figure('Name','Band Filter Responses','NumberTitle','off');
colors = lines(n_bands);
for k = 1:n_bands
    b_k = filter_responses{k}.b;
    a_k = filter_responses{k}.a;
    [H,f_resp] = freqz(b_k,a_k,Nfft,fs);
    plot(f_resp, 20*log10(abs(H)+1e-12), 'Color',colors(k,:), 'LineWidth',1.5);
    hold on;
end
xlabel('Frequency (Hz)'); ylabel('Magnitude (dB)');
title(sprintf('All Band Filter Responses (%s)',filter_type));
legend(arrayfun(@(k) sprintf('Band %d (%.0f-%.0f Hz)', k, band_edges(k), band_edges(k+1)), ...
    1:n_bands, 'UniformOutput',false), 'Location','best');
ylim([-80 5]); xlim([0 min(fs/2,20000)]); grid on;

% --- Per-band: magnitude, phase, impulse, step, pole-zero ---
for k = 1:n_bands
    b_k = filter_responses{k}.b;
    a_k = filter_responses{k}.a;
    [H,f_resp] = freqz(b_k,a_k,Nfft,fs);
    imp_in = [1; zeros(255,1)];
    ir = filter(b_k,a_k,imp_in);

    figure('Name',sprintf('Band %d Analysis',k),'NumberTitle','off');
    subplot(2,3,1);
    plot(f_resp,20*log10(abs(H)+1e-12),'b','LineWidth',1.5);
    title(sprintf('B%d Magnitude (%.0f-%.0f Hz)',k,band_edges(k),band_edges(k+1)));
    xlabel('Hz'); ylabel('dB'); grid on; xlim([0 min(fs/2,20000)]);

    subplot(2,3,2);
    plot(f_resp,angle(H)*180/pi,'r','LineWidth',1.5);
    title('Phase Response'); xlabel('Hz'); ylabel('°'); grid on;

    subplot(2,3,3);
    stem(ir(1:min(end,100)),'g','MarkerSize',3);
    title('Impulse Response'); xlabel('Sample'); grid on;

    subplot(2,3,4);
    plot(cumsum(ir),'m','LineWidth',1.5);
    title('Step Response'); xlabel('Sample'); grid on;

    subplot(2,3,5);
    zplane(b_k,a_k); title('Pole-Zero Plot'); grid on;

    subplot(2,3,6);
    bar([band_edges(k), band_edges(k+1)], [gains_dB(k), gains_dB(k)]);
    title(sprintf('Gain = %.1f dB', gains_dB(k)));
    xlabel('Hz'); ylabel('dB'); ylim([-30 30]); grid on;
end

%% =========================================================
%  STEP 5: TIME-DOMAIN COMPARISON
% =========================================================
seg = 1:min(5*fs, N_audio);
figure('Name','Time Domain Comparison','NumberTitle','off');
subplot(2,1,1);
plot(t_audio(seg), audio_in(seg),'b','LineWidth',0.8);
title('Original Speech Signal'); xlabel('Time (s)'); ylabel('Amplitude'); grid on;
subplot(2,1,2);
plot(t_audio(seg), audio_out(seg),'r','LineWidth',0.8);
title('Equalized Speech Signal'); xlabel('Time (s)'); ylabel('Amplitude'); grid on;
sgtitle('Time-Domain: Original vs Equalized');

%% =========================================================
%  STEP 6: POWER SPECTRAL DENSITY
% =========================================================
figure('Name','PSD Comparison','NumberTitle','off');
win_psd = hamming(512);
[pxx_in,  f_psd] = pwelch(audio_in,  win_psd, 256, 1024, fs);
[pxx_out, ~]     = pwelch(audio_out, win_psd, 256, 1024, fs);
plot(f_psd, 10*log10(pxx_in), 'b','LineWidth',1.5); hold on;
plot(f_psd, 10*log10(pxx_out),'r','LineWidth',1.5);
legend('Original','Equalized','Location','best');
xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
title('Power Spectral Density – Original vs Equalized');
xlim([0 min(fs/2,20000)]); grid on;

%% =========================================================
%  STEP 7: SPECTROGRAM
% =========================================================
figure('Name','Spectrogram Comparison','NumberTitle','off');
subplot(1,2,1);
spectrogram(audio_in,  256, 200, 512, fs, 'yaxis');
title('Original Signal Spectrogram'); colorbar;
ylim([0 min(fs/2000,10)]);

subplot(1,2,2);
spectrogram(audio_out, 256, 200, 512, fs, 'yaxis');
title('Equalized Signal Spectrogram'); colorbar;
ylim([0 min(fs/2000,10)]);
sgtitle('Spectrograms: Before vs After Equalization');

%% =========================================================
%  STEP 8: PLAY AND SAVE OUTPUT
% =========================================================
fprintf('\nPlaying equalized output...\n');
soundsc(audio_out, fs);

% Save at original sample rate
audiowrite('equalized_output.wav', audio_out, fs);
fprintf('Saved: equalized_output.wav at %d Hz\n', fs);

% Save at output sample rate if different
fs_out = round(fs * fs_out_choice);
if fs_out ~= fs
    % Resample
    [p_res, q_res] = rat(fs_out/fs, 1e-4);
    audio_resampled = resample(audio_out, p_res, q_res);
    filename_out = sprintf('equalized_%dHz.wav', fs_out);
    audiowrite(filename_out, audio_resampled, fs_out);
    fprintf('Saved resampled: %s at %d Hz\n', filename_out, fs_out);
end

%% =========================================================
%  STEP 9: DEMONSTRATE FIR vs IIR + SAMPLE RATE CHANGES
% =========================================================
fprintf('\n=== Demonstrating all required sample rates ===\n');
fs_variants = [fs*4, fs/2];
fs_labels   = {'4x','half'};
for vi = 1:length(fs_variants)
    fs_v = round(fs_variants(vi));
    [pv, qv] = rat(fs_v/fs, 1e-4);
    audio_v  = resample(audio_out, pv, qv);
    fname    = sprintf('output_%s_rate.wav', fs_labels{vi});
    audiowrite(fname, audio_v, fs_v);
    fprintf('Saved: %s  (%d Hz, %d samples)\n', fname, fs_v, length(audio_v));
end

%% =========================================================
%  PERFORMANCE SUMMARY
% =========================================================
total_power_in  = sum(audio_in.^2);
total_power_out = sum(audio_out.^2);
fprintf('\n========== PERFORMANCE SUMMARY ==========\n');
fprintf('Filter type      : %s\n', filter_type);
fprintf('Number of bands  : %d\n', n_bands);
fprintf('Input RMS power  : %.4f\n', sqrt(total_power_in/N_audio));
fprintf('Output RMS power : %.4f\n', sqrt(total_power_out/N_audio));
for k = 1:n_bands
    fprintf('  Band %d [%.0f-%.0f Hz]: gain = %.1f dB\n', ...
        k, band_edges(k), band_edges(k+1), gains_dB(k));
end
fprintf('==========================================\n');
fprintf('\nProject II complete! All plots generated.\n');