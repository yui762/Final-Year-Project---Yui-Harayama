% Adaptive Noise Cancellation on Pre-recorded signals 

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
% Initialisation
alpha = 0.5; % smoothing factor for a priori SNR estimation, if do single SS
beta = 2; % power domain
lambda = 1; % subtraction factor controlling amount of suppression (initially set at 3)
NFFT = 128; % length of Hamming window

%% nonstationary noise Gaussian noise 
% (from https://uk.mathworks.com/matlabcentral/answers/311207-how-can-i-generate-a-non-stationary-gaussian-signal)

nostat_noise = 0.1*(1/2*t + 1/2*t.*randn(1,N)); % scaling factor in front for amplitude of noise to be similar to that of signal
desired_snr = -10;
eta = AddNoise(signal,nostat_noise, desired_snr);
a = 1;
b = [1 0 0.5]; 
eta = filter(b, a, eta);
s = signal + eta; % primary signal = heart sound + noise 
cs = [0.1 0.5 1];

% Initialisation
MSPE_mean_all_LMS = zeros(length(cs), 1);
MSPE_mean_all_NLMS = MSPE_end_all_LMS;
SISDR_all_lms = MSPE_mean_all_LMS;

% xcorrs = MSPE_mean_all_LMS; % not us
fig = figure;
tiledlayout('vertical', 'TileSpacing','compact')
for i = 1:length(cs)
    nexttile
    % Construct a reference noise correlated to noise in input signal for ANC
    c = cs(i);
    ref = 2*eta + c;
    
    [e_LMS] = DSS(alpha, beta, lambda, NFFT, s, ref, fs, fs);
    
    plot(t, s,  LineWidth=1.1)
    hold on
    plot(t, e_NLMS, LineWidth=1.1)
    hold on
    plot(t, signal, LineWidth=1.1)
    axis tight
    if i == 1
        ylim([-1.5 1.5])
        lgd = legend('Primary', 'Estimated', 'Original');
    else
        ylim([-1 1])
    end 
    xlim([0 15])
    fontsize(lgd,11,'points')

    % LMS Compute prediction error between original clean and estimated signals
    MSPE_LMS = (signal' - e_LMS).^2;
    MSPE_mean_all_LMS(i,1) = mean(MSPE_LMS);
    
    SISDR_all_lms(i,1) = SI_SDR(signal, e_LMS);

    switch i 
        case 1
            title(sprintf('DSS: Correlation Factor = 0.1, MSPE = %3f', MSPE_mean_all_LMS(i)), FontSize=12)
        case 2
            title(sprintf('DSS: Correlation Factor = 0.5, MSPE = %3f', MSPE_mean_all_LMS(i)), FontSize=12)
        case 3 
            title(sprintf('DSS: Correlation Factor = 1, MSPE = %3f', MSPE_mean_all_LMS(i)), FontSize=12)
        case 4
            title(sprintf('LMS: Correlation Factor = 5, MSPE = %3f', MSPE_mean_all_LMS(i)), FontSize=12)
    end
end

% Give common xlabel, ylabel and title to your figure
han=axes(fig,'visible','off'); 
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Amplitude', FontSize=11);
xlabel(han,'Time (s)', FontSize=11);
set(gcf, 'Position',  [0, 0, 1000, 600])

%% non stat noise, plot for different snr 
nostat_noise = 0.1*(1/2*t + 1/2*t.*randn(1,N)); % scaling factor in front for amplitude of noise to be similar to that of signal
snrs = [-2, -10, -20];

% Initialisation
MSPE_mean_all_LMS = zeros(length(snrs), 1);
MSPE_mean_all_NLMS = MSPE_mean_all_LMS;
SISDR_all_lms = MSPE_mean_all_LMS;

% xcorrs = MSPE_mean_all_LMS; % not us
fig = figure;
tiledlayout('vertical', 'TileSpacing','compact')

for i = 1:length(snrs)
    nexttile
    eta = AddNoise(signal,nostat_noise, snrs(i));
    c = 0.1;
    ref = 2*eta + c;
    a = 1;
    b = [1 0 0.5]; 
    eta = filter(b, a, eta);
    s = signal + eta; % primary signal = heart sound + noise 

    [e_LMS] = DSS(alpha, beta, lambda, NFFT, s, ref, fs, fs);
    t_DSS = (0:length(e_LMS)-1)/fs;
    signal_prime = signal(1:length(e_LMS));
    
    plot(t, s,  LineWidth=1.1)
    hold on
    plot(t, signal, LineWidth=1.1)
    hold on
    plot(t_DSS, e_LMS, LineWidth=1.1)
    axis tight
    if i == 1
        ylim([-1.5 1.5])
        lgd = legend('Primary', 'Original', 'Estimated');
    else
        ylim([-1 1])
    end 
    xlim([0 15])
    fontsize(lgd,11,'points')

    % LMS Compute prediction error between original clean and estimated signals
    MSPE_LMS = mean((signal_prime' - e_LMS).^2);
    SISDR_LMS = SI_SDR(signal_prime', e_LMS);
    MSPE_mean_all_LMS(i,1) = MSPE_LMS;
    SISDR_all_lms(i,1) = SI_SDR(signal, e_LMS);

    switch i 
        case 1
            title(sprintf('DSS: SNR = -2dB, MSPE = %3f', MSPE_mean_all_LMS(i)), FontSize=12)
        case 2
            title(sprintf('DSS: SNR = -10dB, MSPE = %3f', MSPE_mean_all_LMS(i)), FontSize=12)
        case 3 
            title(sprintf('DSS: SNR = -20dB, MSPE = %3f', MSPE_mean_all_LMS(i)), FontSize=12)
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
orders = 10;
desired_snr = -2;

% 1. Compute the power of the clean signal
clean_power = sum(signal.^2) / length(signal);

% 2. Calculate the desired noise power based on the desired SNR
desired_snr_linear = 10^(desired_snr / 10); % Convert SNR from dB to linear scale
desired_noise_power = clean_power / desired_snr_linear;

% 3. Generate white Gaussian noise with the calculated power
eta = sqrt(desired_noise_power) * randn(1,N); % power = variance of WGN

% eta = sqrt(var)*randn(1, N); % noise 
eta = filter(b, a, eta);
s = signal + eta; % primary signal = heart sound + noise 
c = 0.1;
ref = 2*eta + c;

%% other noises
addpath('C:\ReSpeaker\usb_4_mic_array\recordings\simulated noise sources')
[vacuum, fs_vacuum] = audioread('clothes-rustle-6817.mp3');

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
c = 0.1;
ref = 2*eta + c;

%% plotting waveforms for stat noise
N = length(signal);
% Initialisation
alpha = 0.5; % smoothing factor for a priori SNR estimation, if do single SS
beta = 2; % power domain
lambda = 8; % subtraction factor controlling amount of suppression (initially set at 3)
NFFT = 128; % length of Hamming window

[e_LMS] = DSS(alpha, beta, lambda, NFFT, s, ref, fs, fs);
t_DSS = (0:length(e_LMS)-1)/24000;
signal_prime = signal(1:length(e_LMS));
MSPE_lms = mean((signal_prime' - e_LMS).^2);
SISDR_lms = SI_SDR(signal_prime', e_LMS);

figure
tile = tiledlayout('vertical', 'TileSpacing','compact');

t = (1:length(s))/24000;
nexttile
plot(t, s-mean(s),  LineWidth=1.1)
hold on
plot(t, signal-mean(signal), LineWidth=1.1)
hold on
plot(t_DSS, e_LMS, LineWidth=1.1)
hold on
axis tight
set(groot,'defaultLegendInterpreter','latex');
lgd = legend('Primary', 'Original', 'Estimated');

title(sprintf('DSS:  $\\lambda$ = %.1f, MSPE = %3f, SI-SDR = %3f', lambda, MSPE_lms, SISDR_lms),'Interpreter','latex', FontSize=16)
axis tight
xlim([1 1.3])
ylim([-1.2 1.2])
fontsize(lgd,15,'points')
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
xt = 1.162;
yt = 0.81;
str = {'S2'};
text(xt,yt,str, FontSize=16)

% another lambda
lambda = 1;
NFFT = 128; % length of Hamming window

[e_LMS] = DSS(alpha, beta, lambda, NFFT, s, ref, fs, fs);
signal_prime = signal(1:length(e_LMS));
MSPE_lms = mean((signal_prime' - e_LMS).^2);
SISDR_lms = SI_SDR(signal_prime', e_LMS);

nexttile
plot(t, s-mean(s),  LineWidth=1.1)
hold on
plot(t, signal-mean(signal), LineWidth=1.1)
hold on
plot(t_DSS, e_LMS, LineWidth=1.1)
hold on
axis tight
title(sprintf('DSS:  $\\lambda$ = %.1f, MSPE = %3f, SI-SDR = %3f', lambda, MSPE_lms, SISDR_lms),'Interpreter','latex', FontSize=16)


% title(sprintf('DSS: %c = %d, MSPE = %3f, SI-SDR = %3f', 955, lambda, MSPE_lms, SISDR_lms), FontSize=16)
axis tight
xlim([1 1.3])
ylim([-1.2 1.2])
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
xt = 1.162;
yt = 0.81;
str = {'S2'};
text(xt,yt,str, FontSize=16)

% another lambda
lambda = 0.8; % subtraction factor controlling amount of suppression (initially set at 3)
NFFT = 128; % length of Hamming window

[e_LMS] = DSS(alpha, beta, lambda, NFFT, s, ref, fs, fs);
signal_prime = signal(1:length(e_LMS));
MSPE_lms = mean((signal_prime' - e_LMS).^2);
SISDR_lms = SI_SDR(signal_prime', e_LMS);

nexttile
plot(t, s-mean(s),  LineWidth=1.1)
hold on
plot(t, signal-mean(signal), LineWidth=1.1)
hold on
plot(t_DSS, e_LMS, LineWidth=1.1)
hold on
axis tight
title(sprintf('DSS:  $\\lambda$ = %.1f, MSPE = %3f, SI-SDR = %3f', lambda, MSPE_lms, SISDR_lms),'Interpreter','latex', FontSize=16)

axis tight
xlim([1 1.3])
ylim([-1.2 1.2])
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
xt = 1.162;
yt = 0.81;
str = {'S2'};
text(xt,yt,str, FontSize=16)

sgtitle('SNR = -2dB', 'Interpreter','latex', fontsize=16)
set(gcf, 'Position',  [0, 0, 700, 700])
ylabel(tile, 'Amplitude', 'Interpreter','latex', FontSize=18)
xlabel('Time (s)', 'Interpreter','latex', FontSize=18)

%% plotting for other noises
N = length(signal);
% Initialisation
alpha = 0.5; % smoothing factor for a priori SNR estimation, if do single SS
beta = 2; % power domain
lambda = 1; % subtraction factor controlling amount of suppression (initially set at 3)
NFFT = 128; % length of Hamming window

[e_LMS] = DSS(alpha, beta, lambda, NFFT, s, ref, fs, fs);
t_DSS = (0:length(e_LMS)-1)/fs;
signal_prime = signal(1:length(e_LMS));
MSPE_lms = mean((signal_prime' - e_LMS).^2);
SISDR_lms = SI_SDR(signal_prime', e_LMS);

figure
tile = tiledlayout('vertical', 'TileSpacing','compact');

nexttile
plot(t, s-mean(s),  LineWidth=1.1)
hold on
plot(t, signal-mean(signal), LineWidth=1.1)
hold on
plot(t_DSS, e_LMS, LineWidth=1.1)
hold on
axis tight
lgd = legend('Primary', 'Original', 'Estimated');
title(sprintf('DSS: SNR = %ddB, MSPE = %3f', desired_snr, MSPE_lms), FontSize=16)

axis tight
xlim([0 5])
ylim([-1.2 1.2])
fontsize(lgd,11,'points')

set(gcf, 'Position',  [0, 0, 1200, 600])
ylabel(tile, 'Amplitude', FontSize=16)
xlabel('Time (s)', FontSize=16)

