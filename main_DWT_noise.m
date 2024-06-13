% Discrete Wavelet Transform

close all; clc; clear
cd 'C:\ReSpeaker\usb_4_mic_array'

set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

% Part 1: Simulation Data with LMS
addpath('Michigan_Heart_Sounds')

[signal, fs] = audioread("01 Apex, Normal S1 S2, Supine, Bell.mp3");

segment_start = 100;
% segment_end = 0.053*length(signal);
segment_end = 0.5*length(signal);
% if 0.05 i.e. a smaller segment is extracted, seems to be diverging. With
% 0.5, mu = e-4 works well
signal = signal(segment_start:segment_end)';
t = (0:length(signal)-1) / fs;

% Parameters
N = length(signal);

%% non stat noise, plot for different snr 
nostat_noise = 0.1*(1/2*t + 1/2*t.*randn(1,N)); % scaling factor in front for amplitude of noise to be similar to that of signal
snrs = [-2 -10 -20];

% Initialisation
MSPE_mean_all_LMS = zeros(length(snrs), 1);
SISDR_all_lms = MSPE_mean_all_LMS;

% xcorrs = MSPE_mean_all_LMS; % not us
fig = figure;
tiledlayout('vertical', 'TileSpacing','compact')

for i = 1:length(snrs)
    nexttile
    eta = AddNoise(signal,nostat_noise, snrs(i));
    a = 1;
    b = [1 0 0.5]; 
    eta = filter(b, a, eta);
    s = signal + eta; % primary signal = heart sound + noise 

    [w,~] = wdenoise(s, 7, ... % level = 5 so d1 - d5
            'Wavelet', 'db8', 'DenoisingMethod', 'UniversalThreshold', ...
            'ThresholdRule', 'Soft');
    
    plot(t, s,  LineWidth=1.1)
    hold on
    plot(t, w, LineWidth=1.1)
    hold on
    plot(t, signal, LineWidth=1.1)
    axis tight
    if i == 1
        ylim([-1 1])
        lgd = legend('Primary', 'Estimated', 'Original');
    else
        ylim([-1 1])
    end 
    xlim([0 15])
    fontsize(lgd,11,'points')

    % LMS Compute prediction error between original clean and estimated signals
    MSPE_LMS = (signal - w).^2;
    MSPE_mean_all_LMS(i,1) = mean(MSPE_LMS);
    SISDR_all_lms(i,1) = SI_SDR(signal, w);

    switch i 
        case 1
            title(sprintf('DWT: SNR = -2dB, MSPE = %3f', MSPE_mean_all_LMS(i)), FontSize=12)
        case 2
            title(sprintf('DWT: SNR = -10dB, MSPE = %3f', MSPE_mean_all_LMS(i)), FontSize=12)
        case 3 
            title(sprintf('DWT: SNR = -20dB, MSPE = %3f', MSPE_mean_all_LMS(i)), FontSize=12)
    end
end

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Amplitude', FontSize=11);
xlabel(han,'Time (s)', FontSize=11);
set(gcf, 'Position',  [0, 0, 1000, 600])

%% white noise = stationary noise
a = 1;
b = [1 0 0.5]; 
desired_snr = -2;

% 1. Compute the power of the clean signal
clean_power = sum(signal.^2) / length(signal);

% 2. Calculate the desired noise power based on the desired SNR
desired_snr_linear = 10^(desired_snr / 10); % Convert SNR from dB to linear scale
desired_noise_power = clean_power / desired_snr_linear;

% 3. Generate white Gaussian noise with the calculated power
eta = sqrt(desired_noise_power) * randn(1,N); % power = variance of WGN
eta = filter(b, a, eta);
s = signal + eta; % primary signal = heart sound + noise 


%% other noises
addpath('C:\ReSpeaker\usb_4_mic_array\recordings\simulated noise sources')
[vacuum, fs_vacuum] = audioread('vacuum cleaner-pixabay.mp3');

% vacuum = vacuum(1:0.5*length(vacuum))';
% sampling_factor = 2;
% vacuum = downsample(vacuum, sampling_factor);
signal = signal(1:length(vacuum));
vacuum = 2 * (vacuum - min(vacuum)) / (max(vacuum) - min(vacuum)) -1;
% fs_vacuum = fs_vacuum/sampling_factor;
t = (0:length(vacuum)-1) / fs_vacuum;
figure
plot(t, vacuum)

%%
desired_snr = -15;
eta = AddNoise(signal,vacuum, desired_snr);
s = signal' + eta; % primary signal = heart sound + noise 

%%
% plotting waveforms 
N = length(signal);

% With Wavelet transform: only need one microphone
% Observations: signal seems fine but when zoom in can see that its less
% smooth. For coeff = 2, less filtered and if zoom in, can see some peaks
% in the waves but as coeff increases, waveform becomes more square

[w,DENOISEDCFS,ORIGCFS] = wdenoise(s, 7, ... % level = 5 so d1 - d5
            'Wavelet', 'db8', 'DenoisingMethod', 'UniversalThreshold', ...
            'ThresholdRule', 'Soft');
% db7 good for heart 
% threshold selection rule: rigrsure, heursure, minimaxi, sqtwolog
MSPE = mean((signal' - w).^2);
SISDR = SI_SDR(signal', w);

figure;
plot(t, s-mean(s), LineWidth=1.1);
hold on
plot(t, w, LineWidth=1.1)
hold on
plot(t, signal-mean(signal), LineWidth=1.1);

lgd = legend('Primary', 'Estimated', 'Original');
fontsize(lgd,11,'points')
title(sprintf('DWT: SNR = %ddB, MSPE = %3f', desired_snr, MSPE), FontSize=12)

axis tight
xlim([0 5])
ylim([-1.2 1.2])


set(gcf, 'Position',  [0, 0, 1200, 400])
ylabel('Amplitude', FontSize=11)
xlabel('Time (s)', FontSize=11)

%% Effect of decom level on performance (order has v little effect)
% plotting waveforms 
N = length(signal);

[w,DENOISEDCFS,ORIGCFS] = wdenoise(s, 7, ... % level = 5 so d1 - d5
            'Wavelet', 'db8', 'DenoisingMethod', 'UniversalThreshold', ...
            'ThresholdRule', 'Soft');
% db7 good for heart 
% threshold selection rule: rigrsure, heursure, minimaxi, sqtwolog
MSPE = mean((signal - w).^2);
SISDR = SI_SDR(signal', w);

t = (1:length(s))/24000;

figure;
tile = tiledlayout('vertical', 'TileSpacing', 'compact');
nexttile
plot(t, s-mean(s), LineWidth=1.1);
hold on
plot(t, w, LineWidth=1.1)
hold on
plot(t, signal-mean(signal), LineWidth=1.1);

axis tight
xlim([0 5])
ylim([-1.2 1.2])
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');
lgd = legend('Primary', 'Estimated', 'Original');
fontsize(lgd,15,'points')
title(sprintf('DWT: order = 8, level = 7, \n MSPE = %3f, SI-SDR = %3f', MSPE, SISDR), FontSize=16)

nexttile
[w,DENOISEDCFS,ORIGCFS] = wdenoise(s, 1, ... % level = 5 so d1 - d5
            'Wavelet', 'db8', 'DenoisingMethod', 'UniversalThreshold', ...
            'ThresholdRule', 'Soft');
% db7 good for heart 
% threshold selection rule: rigrsure, heursure, minimaxi, sqtwolog
MSPE = mean((signal - w).^2);
SISDR = SI_SDR(signal', w);
plot(t, s-mean(s), LineWidth=1.1);
hold on
plot(t, w, LineWidth=1.1)
hold on
plot(t, signal-mean(signal), LineWidth=1.1);
ylabel(tile, 'Amplitude', 'Interpreter', 'latex', FontSize=18)
xlabel(tile, 'Time (s)', 'Interpreter', 'latex', FontSize=18)
xlim([0 5])
ylim([-1.2 1.2])
ax = gca;
ax.FontSize = 15;  % Font Size of 15

title(sprintf('DWT: order = 8, level = 1, \n MSPE = %3f, SI-SDR = %3f', MSPE, SISDR), 'Interpreter', 'latex', FontSize=16)
sgtitle('SNR = -2dB', 'Interpreter', 'latex', fontsize=16)

set(gcf, 'Position',  [0, 0, 1200, 600])


%% individual plots (selected ones in report)
close all; clc; clear
set(groot,'defaultAxesTickLabelInterpreter','latex'); 
set(groot,'defaulttextinterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaulttextinterpreter','latex');

cd 'C:\ReSpeaker\usb_4_mic_array\recordings\with ECG'
addpath('C:\ReSpeaker\usb_4_mic_array')

% Plotting for visualisation purpose
% Load the signals
[signal1, fs1] = audioread('output_device1_2024-05-10_15-43-14.wav'); % primary signal
[signal2, fs2] = audioread('output_device2_2024-05-10_15-43-14.wav'); % reference signal

[signal1, signal2] = combine_channels(signal1, signal2, fs1, fs2);
signal1 = signal1(2.8*fs1:5*fs1);
signal2 = signal2(2.8*fs1:5*fs1);
t = (0:length(signal1) - 1) / fs1;

N = length(signal1);
[w,DENOISEDCFS,ORIGCFS] = wdenoise(signal1, 5, ... % level = 5 so d1 - d5
            'Wavelet', 'coif5', 'DenoisingMethod', 'UniversalThreshold', ...
            'ThresholdRule', 'Soft');

% Sample parameters
fc_low = 20; % Lower cutoff frequency
fc_high = 200; % Upper cutoff frequency
order = 4; % Filter order

[b, a] = butter(order, [fc_low/(fs1/2) fc_high/(fs1/2)], 'bandpass');
filtered_signal = filter(b, a, w);
norm_filtered_signal = 2 * (filtered_signal - min(filtered_signal)) / (max(filtered_signal) - min(filtered_signal)) -1;
t_filtered = (0:length(norm_filtered_signal) - 1)/fs1;

figure
plot(t_filtered, filtered_signal)
xlabel('Time (s)', 'Interpreter', 'latex',  FontSize=16);
ylabel('N. Amplitude',  'Interpreter', 'latex',  FontSize=16);
title('DWT', 'Interpreter', 'latex',  FontSize=22);
axis tight
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');
set(gcf, 'Position',  [0, 0, 500, 200])

figure
plot(t, signal2)
xlabel('Time (s)', 'Interpreter', 'latex',  FontSize=16);
ylabel('N. Amplitude',  'Interpreter', 'latex',  FontSize=16);
title('Reference Signal', 'Interpreter', 'latex',  FontSize=22);
axis tight
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');
set(gcf, 'Position',  [0, 0, 500, 200])

figure
plot(t, signal1)
xlabel('Time (s)', 'Interpreter', 'latex',  FontSize=16);
ylabel('N. Amplitude',  'Interpreter', 'latex',  FontSize=16);
title('Primary Signal', 'Interpreter', 'latex',  FontSize=22);
axis tight
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');
set(gcf, 'Position',  [0, 0, 500, 200])

