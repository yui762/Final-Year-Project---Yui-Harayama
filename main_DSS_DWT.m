% This script is to apply the Dual spectral subtraction followed by a Butterworth 
% bandpass filtering and wavelet transform (decomposition and denoising) followed
% again by a Butterworth bandpass filtering

% close all;
clc; clear 

cd 'C:\ReSpeaker\usb_4_mic_array'
addpath('C:\ReSpeaker\usb_4_mic_array\recordings\real world noise sources\heavy breathing')
% Read the heart sound + noise signal (signal1) and noise-only signal (signal2)
% super nice profile with this output_device2_2024-05-06_22-44-18
% output_device2_2024-05-10_15-06-49
% output_device2_2024-05-10_14-41-53
% output_device2_2024-05-06_23-03-17
% output_device2_2024-05-05_22-56-49
[signal1, fs1] = audioread("output_device1_2024-06-01_10-09-11.wav"); 
[signal2, fs2] = audioread("output_device2_2024-06-01_10-09-11.wav"); 

% output_device1_2024-05-24_11-43-37
% Combine channels together
[combined1, combined2] = combine_channels(signal1, signal2, fs1, fs2);

% Extract wanted segment
signal1 = combined1';
% signal1 = signal1(8040:17120);
signal2 = combined2'; % noise
% signal2 = signal2(8040:17120);

% Initialisation
alpha = 0.5; % smoothing factor for a priori SNR estimation, if do single SS
beta = 2; % power domain
lambda = 7; % subtraction factor controlling amount of suppression (initially set at 3)
NFFT = 128; % length of Hamming window

% Apply DSS
[output_signal] = DSS(alpha, beta, lambda, NFFT, signal1, signal2, fs1, fs2);

% Write the filtered signal to a WAV file
% audiowrite('DSS1.wav', output_signal, fs1);
%% Power spectrum

% PSD estimate of heart sound + noise
[Pxx_primary, f_primary] = pwelch(signal1, rectwin(length(signal1)), [], fs1*5, fs1);
% [Pxx_s1, f1] = pwelch(signal1, hamming(length(signal1)), 0, 1024, 'onesided'); % 4 = # overlapped segments used to compute average PSD
% [Pxx_s1, f1] = periodogram()

% Generate spectrogram
% window = hamming(512); % Window function
% noverlap = 256; % Overlap between windows
% nfft = 1024; % FFT length

% figure;
% spectrogram(signal1, window, noverlap, nfft, fs1, 'yaxis');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Spectrogram of Combined Filtered Signal');

% ------------------------------ Signal 2 -----------------------------------%
% Construct spectrogram for noise-only signal (signal2)
% figure;
% [Pxx, f] = pwelch(signal2, hamming(length(signal2)), 0, 1024, 4, 'onesided');
% plot(f, 10*log10(Pxx));

% figure;
% spectrogram(signal2, window, noverlap, nfft, fs1, 'yaxis');
% xlabel('Time (s)');
% ylabel('Frequency (Hz)');
% title('Spectrogram of Combined Filtered Signal');

%%
% Plot Original, Noise, and Filtered Signals in time domain
figure;

% Original Signal (Primary)
tile = tiledlayout('vertical');
nexttile
t_original = (0:length(signal1) - 1) / fs1;
plot(t_original, signal1);
ylabel('Amplitude', FontSize=11);
title('Primary Signal', FontSize=12);
% xlim([0 10])
axis tight
% ylim([-0.21 0.21])

% Noise-Only Signal (Reference)
nexttile
t_noise = (0:length(signal2) - 1) / fs1;
plot(t_noise, signal2);
ylabel('Amplitude', FontSize=11);
title('Reference Signal', FontSize=12);
% xlim([0 10])
axis tight
% ylim([-0.21 0.21])

% Processed Signal
nexttile
t_filtered = (0:length(output_signal) - 1) / fs1;
norm_output_signal = 2 * (output_signal - min(output_signal)) / (max(output_signal) - min(output_signal)) -1;
plot(t_filtered, norm_output_signal);
ylabel('N. Amplitude', FontSize=11);
title('Denoised Signal', FontSize=12);
% xlim([0 10])
axis tight

% -------------------------------------------------------------------%
%                          Bandpass Filtering                        %
% -------------------------------------------------------------------%
% Sample parameters
fc_low = 20; % Lower cutoff frequency
fc_high = 150; % Upper cutoff frequency
order = 4; % Filter order

% Creating a Butterworth bandpass filter
[b, a] = butter(order, [fc_low/(fs1/2) fc_high/(fs1/2)], 'bandpass');
filtered_signal = filter(b, a, output_signal);
shifted = filtfilt(b,a,output_signal); % IIR filter

nexttile 
% normalise
norm_filtered_signal = 2 * (filtered_signal - min(filtered_signal)) / (max(filtered_signal) - min(filtered_signal)) -1;
% norm_shifted = 2 * (shifted - min(shifted)) / (max(shifted) - min(shifted)) -1;
plot(t_filtered, norm_filtered_signal)
% hold on
% plot(t_filtered, norm_shifted)

xlabel(tile, 'Time (s)', FontSize=11);
ylabel('N. Amplitude', FontSize=11);
title('Denoised and Bandpass Filtered Signal', FontSize=12);
% xlim([0 10])
axis tight
xt = [0.1935 0.41975 0.921 1.1905 1.65225 1.92625];
yt = [0.6 0.6 0.6 0.6 0.6 0.6];
str = {'S1','S2', 'S1','S2','S1','S2'};
text(xt,yt,str)

set(gcf, 'Position',  [0, 0, 1200, 600])

% figure
% plot(t_filtered, norm_filtered_signal)
% set(gcf, 'Position',  [0, 0, 800, 300])
% xlabel('Time (s)', FontSize=11);
% ylabel('N. Amplitude', FontSize=11);
% title('Denoised and Bandpass Filtered Signal', FontSize=12);
% axis tight
% xt = [0.3345 0.538 1.110 1.38625];
% yt = [0.6 0.6 0.6 0.6];
% str = {'S1','S2', 'S1','S2'};
% text(xt,yt,str)

% SS + Bandpassed Signal Only
% figure;
% plot(t_filtered, norm_filtered_signal);
% title('Spectral Subtraction and Filtered Result');
% xlabel('Time (s)'); ylabel('Amplitude (AU)');
% set(gcf, 'Position',  [100, 100, 1000, 250]) %row, col width height

figure
plot(t_filtered, norm_filtered_signal)
xlabel('Time (s)', FontSize=11);
ylabel('N. Amplitude', FontSize=11);
title('Denoised and Bandpass Filtered Signal', FontSize=12);
axis tight
xlim([0 10])
set(gcf, 'Position',  [0, 0, 1200, 300])


% hold on;
% scatter(x_clicked, y_clicked, 'r', 'filled');
% hold off;

% Prompt the user to click on the plot
% disp('Click on the plot to select points. Press Enter when finished.');
% [x_clicked, y_clicked] = ginput;
% 
% % Display the coordinates of the clicked points
% disp('Coordinates of clicked points:');
% disp([x_clicked, y_clicked]);

%% Save the coordinates into an Excel file
coordinates = [x_clicked, y_clicked];
xlswrite('clicked_points.xlsx', coordinates, 'Sheet1');

% save as csv
% Specify the file path
% C = horzcat(filtered_signal, t_filtered');
% filePath = 'output_device1_2024-05-10_15-06-49.csv';
% 
% % Save the array into a CSV file
% writematrix(filtered_signal, filePath);
% writematrix(C, filePath);
%% To show how reference also has heart sound
figure
tile = tiledlayout('vertical');
nexttile
t_original = (0:length(signal1) - 1) / fs1;
plot(t_original, signal1);
ylabel('Amplitude', FontSize=11);
title('Primary Signal', FontSize=12);
axis tight

% Define the vertices of the square
x1 = [0.22 0.43925 0.43925 0.22];
y1 = [-1 -1 1 1];
x2 = [0.89 1.039 1.039 0.89];
y2 = [-1 -1 1 1];
x3 = [1.175 1.28 1.28 1.175];
y3 = [-1 -1 1 1];
x4 = [1.593 1.83 1.83 1.593];
y4 = [-1 -1 1 1];
x5 = [1.926 2.0177 2.0177 1.926];
y5 = [-1 -1 1 1];
x6 = [0.45 0.63925 0.63925 0.45];
y6 = [-1 -1 1 1];

% Define the opacity (alpha value)
opacity = 0.1;
% Plot the filled square
patch(x1, y1, 'm', 'FaceAlpha', opacity);
patch(x2, y2, 'm', 'FaceAlpha', opacity);
patch(x3, y3, 'm', 'FaceAlpha', opacity);
patch(x4, y4, 'm', 'FaceAlpha', opacity);
patch(x5, y5, 'm', 'FaceAlpha', opacity);
patch(x6, y6, 'm', 'FaceAlpha', opacity);
ylim([-0.6 0.6])

xt = [0.2335 0.45975 0.921 1.1905 1.62225 1.93625];
yt = [0.45 0.45 0.45 0.45 0.45 0.45];
str = {'S1','S2', 'S1','S2','S1','S2'};
text(xt,yt,str)

% Noise-Only Signal (Reference)
nexttile
t_noise = (0:length(signal2) - 1) / fs1;
plot(t_noise, signal2);
ylabel('Amplitude', FontSize=11);
title('Reference Signal', FontSize=12);
% xlim([0 10])
axis tight
x1 = [0.22 0.43925 0.43925 0.22];
y1 = [-1 -1 1 1];
x2 = [0.89 1.039 1.039 0.89];
y2 = [-1 -1 1 1];
x3 = [1.175 1.28 1.28 1.175];
y3 = [-1 -1 1 1];
x4 = [1.593 1.83 1.83 1.593];
y4 = [-1 -1 1 1];
x5 = [1.926 2.0177 2.0177 1.926];
y5 = [-1 -1 1 1];
x6 = [0.45 0.63925 0.63925 0.45];
y6 = [-1 -1 1 1];

% Define the opacity (alpha value)
opacity = 0.1;
% Plot the filled square
patch(x1, y1, 'm', 'FaceAlpha', opacity);
patch(x2, y2, 'm', 'FaceAlpha', opacity);
patch(x3, y3, 'm', 'FaceAlpha', opacity);
patch(x4, y4, 'm', 'FaceAlpha', opacity);
patch(x5, y5, 'm', 'FaceAlpha', opacity);
patch(x6, y6, 'm', 'FaceAlpha', opacity);
ylim([-0.08 0.08])

xt = [0.2335 0.45975 0.921 1.1905 1.62225 1.93625];
yt = [0.06 0.06 0.06 0.06 0.06 0.06];
str = {'S1','S2', 'S1','S2','S1','S2'};
text(xt,yt,str)

set(gcf, 'Position',  [0, 0, 1200, 600])
%% With Wavelet transform: only need one microphone
% Observations: signal seems fine but when zoom in can see that its less
% smooth. For coeff = 2, less filtered and if zoom in, can see some peaks
% in the waves but as coeff increases, waveform becomes more square

figure;

% Original Signal
tile = tiledlayout('vertical');
nexttile
t_original = (0:length(signal1) - 1) / fs1;
plot(t_original, signal1);
ylabel('Amplitude', FontSize=11);
title('Primary Signal', fontsize = 12);
axis tight

% Wavelet transform denoised
% https://github.com/Sephoro/HSSaC/blob/master/Code/Segmentation/GenerateFeatures.m
nexttile
[w,DENOISEDCFS,ORIGCFS] = wdenoise(signal1, 5, ... % level = 5 so d1 - d5
            'Wavelet', 'coif5', 'DenoisingMethod', 'UniversalThreshold', ...
            'ThresholdRule', 'Soft');
% db7 good for heart 
% threshold selection rule: rigrsure, heursure, minimaxi, sqtwolog
% instead of plotting w can also plot(DENOISEDCFS{6})
norm_dwt_signal = 2 * (w - min(w)) / (max(w) - min(w)) -1;
plot(t_original, -norm_dwt_signal); % invert 
ylabel('N. Amplitude', FontSize=11);   
title('Denoised Signal', FontSize=12)
axis tight

% Filtered Signal
nexttile
fc_low = 20; % Lower cutoff frequency
fc_high = 150; % Upper cutoff frequency
order = 4; % Filter order
[b, a] = butter(order, [fc_low/(fs1/2) fc_high/(fs1/2)], 'bandpass');
filtered_w = filter(b, a, w);
norm_filtered_w = 2 * (filtered_w - min(filtered_w)) / (max(filtered_w) - min(filtered_w)) -1;
plot(t_original, norm_filtered_w);
ylabel('N. Amplitude', FontSize=11);
title('Denoised and Bandpass Filtered Signal', FontSize=12);
axis tight
xlabel(tile, 'Time (s)', FontSize=11);

xt = [0.1435 0.41975 0.921 1.1905 1.7707 2.005];
yt = [0.6 0.6 0.6 0.6 0.6 0.6];
str = {'S1','S2', 'S1','S2','S1','S2'};
text(xt,yt,str)
set(gcf, 'Position',  [0, 0, 1200, 500])

SISDRin = SI_SDR_real(signal1, signal2);
SISDRout = SI_SDR(w, signal2);
gain = SISDRout/SISDRin;

% audiowrite('wavelet.wav', w, fs1);
%% DWT explicit decomposition into coefficients
% [c,l] = wavedec(x,n,wname) returns the wavelet decomposition of the 1-D signal x at level n using the wavelet wname
% wavelet decomposition vector c and the bookkeeping vector l, which is
% used to parse c (Mathworks)

% The wavelet decomposition results in levels of approximated and detailed coefficients. 
[c,l] = wavedec(signal1,5,"db7");
approx = appcoef(c,l,"db7");
[cd1,cd2,cd3,cd4,cd5] = detcoef(c,l,[1 2 3 4 5]); % 5 coeff because level of decomposition = 5

figure
tiledlayout('vertical')
% nexttile
% plot(t_original,w)
% title('Denoised Signal')
% axis tight

nexttile
plot(approx)
title("Approximation Coefficients")
axis tight

nexttile
plot(cd5)
title("Level 5 Detail Coefficients")
axis tight

nexttile
plot(cd4)
title("Level 4 Detail Coefficients")
axis tight

nexttile
plot(cd3)
title("Level 3 Detail Coefficients")
axis tight

nexttile
plot(cd2)
title("Level 2 Detail Coefficients")
axis tight

nexttile
plot(cd1)
title("Level 1 Detail Coefficients")
axis tight
set(gcf, 'Position',  [0, 0, 500, 700])

%% to quantify noise reduction, plot the power spectrum of original noisy signal and cleaned signal 
% difference in power of the noise floor = reduction in noise 

% don't forget to change the window length according to the length of the
% extracted segment
[Pxx_filtered, f2] = pwelch(filtered_signal, hamming(length(signal1)), [], 1024, 'onesided'); % onesided = real
% [Pxx_filtered, f2] = pwelch(filtered_w, hamming(length(signal1)), 0, 1024, 'onesided'); % onesided = real
[Pxx_filtered, f_filtered] = pwelch(filtered_signal, rectwin(length(filtered_signal)), [], fs1*5, fs1);

figure
hold all
plot(f_primary, 10*log10(Pxx_primary))
plot(f_filtered,10*log10(Pxx_filtered))
xlim([0 300])

% figure;
% hold all
% plot(f1, 10*log10(Pxx_s1));
% plot(f2, 10*log10(Pxx_filtered));
xlabel('Frequency (Hz)');
ylabel('PSD (dB)')
legend('Original Signal', 'Processed Signal')
title('Power Spectrum Density Estimate')
% xlim([0 3])

%% Power spectrum plotting using pspectrum - whole signal
p_noise = pspectrum(signal1, fs1);
p = pspectrum(filtered_signal,fs1); % can precise the type
figure
plot(10*log10(p_noise), LineWidth=1.5)
hold on
plot(10*log10(p)+10, LineWidth=1.5) % shift up
xlim([0 500])
xlabel('Frequency (Hz)', 'FontSize', 11);
ylabel('Power Spectrum (dB)', FontSize=11)
lgd = legend('Primary Signal', 'Denoised and Bandpassed Signal');
fontsize(lgd,11,'points')

title('Power Spectrum', FontSize=12)
set(gcf, 'Position',  [0, 0, 800, 500])

%% Power spectrum plotting using pspectrum - s1
p_noise = pspectrum(signal1, fs1);
p = pspectrum(filtered_signal,fs1); % can precise the type
figure
plot(10*log10(p_noise), LineWidth=1.5)
hold on
plot(10*log10(p), LineWidth=1.5)
xlim([0 500])
xlabel('Frequency (Hz)', 'FontSize', 11);
ylabel('Power Spectrum (dB)', FontSize=11)
lgd = legend('Primary Signal', 'Denoised and Bandpassed Signal');
fontsize(lgd,11,'points')

title('Power Spectrum - S2', FontSize=12)
set(gcf, 'Position',  [0, 0, 800, 500])
