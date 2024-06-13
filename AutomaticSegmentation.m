close all; clc; clear

cd 'C:\ReSpeaker\usb_4_mic_array'
addpath('Michigan_Heart_Sounds')

[signal, fs] = audioread("01 Apex, Normal S1 S2, Supine, Bell.mp3");
var = 0.0001;
N = length(signal);
% signal = signal + sqrt(var).*randn(N, 1);
% 06 Apex, Early Sys Mur, Supine, Bell.mp3 is tricky because of the murmurs

% --------------------------------------------------------%
% cd 'C:\ReSpeaker\usb_4_mic_array'
% addpath ('phonocardiogram-heart-sound-analysis-main\data\training-a')
% [signal, fs] = audioread("a0012.wav"); % if want to use PhysioNet
% --------------------------------------------------------%

% done = apex S3, S4, holosystolic 

% --------------------------------------------------------%
% signal = filtered_signal; % if want to use MY data
% signal = signal - mean(signal); % remove mean 
% fs = fs1; % if want to use MY data
% --------------------------------------------------------%

segment_start = 33000;
% segment_start = 1;
segment_end = 0.06*length(signal);
% segment_end = 3.6*fs;
signal = signal(segment_start:segment_end)';
t = (0:length(signal)-1) / fs; 

figure
plot(t, signal)
axis tight

% autocorrelation
% [acf,lags] = xcorr(signal, 'biased');
% figure
% plot(acf)

%% Add WGN to clean signal 

desired_snr_db = 20;
[noisy_signal] = AddWGN(signal, desired_snr_db);

% Plot the clean signal and the noisy signal
tiledlayout('vertical')
plot(t, signal)
xlabel('Time (s)')
ylabel('Amplitude')
title('Clean Signal')
axis tight

nexttile
plot(t, noisy_signal)
xlabel('Time (s)')
ylabel('Amplitude')
title(['Noisy Signal (SNR = ', num2str(desired_snr_db), ' dB)'])
axis tight

%%  Filtering out the noise
fc_low = 20; % Lower cutoff frequency
fc_high = 150; % Upper cutoff frequency
order = 3; % Filter order, either 4 or 5 (both from literature and empirical observation)

% Creating a Butterworth bandpass filter
[b, a] = butter(order, [fc_low/(fs/2) fc_high/(fs/2)], 'bandpass');
filtered_signal = filter(b, a, real(noisy_signal));

subplot(3, 1, 3)
plot(t, filtered_signal)
xlabel('Time')
ylabel('Amplitude')
title('Filtered Signal')
axis tight

signal = filtered_signal;

%% Error between filtered and original signals
MSE = mean(abs(signal - filtered_signal).^2);
MSEdB = 10*log10(MSE); % or pow2db
error = signal - filtered_signal;
figure
plot(t, error)
axis tight

%% Envelope Extraction

% ---------------- Using filtered signal from DSS --------------------%
% in report: output_device1_2024-05-10_15-43-14(3.5*fs1:7.1*fs1) was used
% signal = filtered_signal;
% fs = fs1;
% t = (0:length(signal)-1) / fs1; 
% --------------------------------------------------------------------%
% signal = noisy_signal;
addpath('C:\ReSpeaker\usb_4_mic_array')
[absolute_n, energy_n, shannon_energy_n, shannon_entropy_n] = EnvelopeExtraction(signal);
figure
tile = tiledlayout('vertical');
set(gcf, 'Position',  [0,  0, 1200, 600])
nexttile
plot(t, absolute_n)
title('Absolute', 'Interpreter','latex', FontSize=18)
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')

nexttile
plot(t, energy_n)
title('Energy', 'Interpreter','latex', FontSize=18)
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')

nexttile
plot(t, -shannon_energy_n+1) % negative and add 1 to normalise
title('Shannon Energy', 'Interpreter','latex', FontSize=18)
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')

nexttile
plot(t, shannon_entropy_n)
title('Shannon Entropy', 'Interpreter','latex', FontSize=18)
ylabel(tile, 'N. Amplitude', 'Interpreter','latex', FontSize=16)
xlabel(tile, 'Time (s)', 'Interpreter','latex', FontSize=16)
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(gcf, 'Position',  [0, 0, 1000, 700])

% type_envelope = input('Which envelope do you want to compute the HT of?');
type_envelope = 3;
switch type_envelope
    case 1
        type_envelope = 'Absolute';
    case 2
        type_envelope = 'Energy';
    case 3
        type_envelope = 'Shannon Energy';
    case 4
        type_envelope = 'Shannon Entropy';
end

%% Shannon Energy Hilbert Transform (SEHT)

[selected_envelope,contour_env] = SelectedEnvelopeExtraction(signal, type_envelope);
% SEenv = envelope(selected_envelope, 100, 'peak');

% Initialisation
N = 550; % length of rectwin
M = 12000; % length of MA filter 

[smooth, zn] = SEHT(selected_envelope, N, M);

set(gcf, 'Position',  [0, 0, 1200, 500])
figure
tiledlayout('vertical', 'TileSpacing','compact')
nexttile
% normalise signal between -1 and +1
min_signal = min(signal);
max_signal = max(signal);
normalised_signal = ((signal - min_signal) / (max_signal- min_signal)) * 2 - 1;
plot(t, normalised_signal)
title('PCG Signal', 'Interpreter', 'latex', fontsize=18);
ylabel('N. Amplitude', 'Interpreter', 'latex', fontsize=16);
axis tight
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');

nexttile
if strcmp(type_envelope,'Absolute') || strcmp(type_envelope,'Energy')
    plot(t, selected_envelope)
    hold on
    plot(t, contour_env)
elseif strcmp(type_envelope,'Shannon Energy')
    plot(t, -selected_envelope+1)
    hold on
    plot(t, contour_env+1)
else
    plot(t, selected_envelope)
    hold on
    plot(t, contour_env)
end
title(sprintf('%s', type_envelope), fontsize=18);
ylabel('N. Amplitude', 'Interpreter', 'latex', fontsize=16);
lgd = legend('', 'Envelope');
fontsize(lgd,15,'points')
axis tight
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');

% nexttile
% plot(t, smooth)
% hold on
% title('Smoothed Envelope', FontSize=14) % HT of shannon energy
% yline(0)
% ylabel('Amplitude', FontSize=11)
% axis tight
% ylim([0.5 1.5])
% ax = gca;
% ax.FontSize = 15;  % Font Size of 15
% set(gca,'TickLabelInterpreter','latex')
% set(groot,'defaultLegendInterpreter','latex');

nexttile
plot(t, zn)
hold on
title('Hilbert Transform of Smoothed Envelope', FontSize=14) % HT of shannon energy
yline(0)
xlabel('Time (s)', FontSize=12)
ylabel('Amplitude', FontSize=11)
axis tight
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');
hold on
%% Find peaks from Hilbert Transform of Shannon Energy

% Step 1: Find peaks of the signal which correspond to the zero crossing
% from positive to negative of the HT of the Shannon energy z(n)
peak_locs_temp = FindPeaks(zn);
scatter(peak_locs_temp./fs, 0, 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5); % plot zero-crossing of the HT
lgd = legend('', '', 'Zero Crossing');
fontsize(lgd,15,'points')
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');
set(gcf, 'Position',  [0, 0, 1000, 700])

% Step 2: Normalise signal between -1 and +1
figure
min_signal = min(signal);
max_signal = max(signal);
normalised_signal = ((signal - min_signal) / (max_signal- min_signal)) * 2 - 1;
normalised_signal = normalised_signal - mean(normalised_signal); % remove mean
plot(t, normalised_signal)
hold on
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');

% Step 3: Take envelope of the normalised signal
env_normalised_signal = envelope(normalised_signal, 500, 'peak');
% plot(t, env_normalised_signal); % can plot or not the envelope itself.
% as want a rough envelope to get the amplitude at the peak_locs_temp
% input arg = 500

% Step 4: Take amplitude of the envelope at the peaks 
amplitudes = zeros(1, length(peak_locs_temp));

for i = 1:length(peak_locs_temp)
    amplitudes(i) = env_normalised_signal(peak_locs_temp(i)); 
    scatter(t(peak_locs_temp(i)), amplitudes(i), 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
    hold on
end
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');

title('Detected Peaks by using Zero Crossing of z(n)', 'Interpreter', 'latex', FontSize=18)
xlabel('Time (s)', 'Interpreter', 'latex', FontSize=16)
ylabel('Amplitude', 'Interpreter', 'latex', FontSize=16)
axis tight
ylim([-1.1 1.1])
lgd = legend('', 'Peaks');
fontsize(lgd,15,'points')
set(gcf, 'Position',  [0,  0, 1000, 400])
%% Endpoint determination = adaptive thresholding

% Step 1: sort peak vector in ascending order
amplitudes_ascending = sort(amplitudes);

% Step 2: detection threshold = 0.01 of average of few peak values extracted
% from the sorted peak vector
amplitudes_mean = mean(amplitudes_ascending(1:3));
threshold = 0.6*amplitudes_mean; % thresholding factor is super important!

% Step 3: adaptive thresholding rule
g = AdaptiveThresholding(signal, contour_env, threshold);

figure
set(gcf, 'Position',  [0, 0, 1000, 300])
plot(t, normalised_signal) % normalised between -1 and +1
hold on

% Mark by a circle the boundaries of each S1, S2, S3 peak
[zero_to_one_transitions, one_to_zero_transitions] = BoundaryDetection(g);

hold all
scatter(zero_to_one_transitions./fs, mean(normalised_signal), 'o', 'MarkerEdgeColor', 'b', 'LineWidth', 1.5); % divide by fs to get time 
scatter(one_to_zero_transitions./fs, mean(normalised_signal), 'o', 'MarkerEdgeColor', 'g', 'LineWidth', 1.5);
% add the peaks found with z(n)
scatter(peak_locs_temp./fs, amplitudes, 'o', 'MarkerEdgeColor', 'r', 'LineWidth', 1.5);
title('Detected Boundaries and Peaks', 'Interpreter', 'latex', FontSize=18)
axis tight 
ylim([-1.1 1.1])
ylabel('N. Amplitude', 'Interpreter', 'latex', FontSize=16)
xlabel('Time (s)', 'Interpreter', 'latex', FontSize=16)
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');

qw{1} = plot(nan, 'ro', 'LineWidth', 1.5);
qw{2} = plot(nan, 'bo', 'LineWidth', 1.5);
qw{3} = plot(nan, 'go', 'LineWidth', 1.5);

lgd = legend([qw{:}], {'Peak', 'Start', 'End'}, 'location', 'northeast', 'NumColumns', 3);
fontsize(lgd,15, 'points')
set(gcf, 'Position',  [0,  0, 1000, 400])
%% Assigning systole and diastole (systole shorter than diastole)
sys_dias = cell(length(one_to_zero_transitions), 1); % Initialize as a cell array

for i = 1:length(one_to_zero_transitions)
    % while i + 2 <length(sys_dias)
        if (abs(one_to_zero_transitions(i)-zero_to_one_transitions(i+1))) < abs((one_to_zero_transitions(i+1)-zero_to_one_transitions(i+2)))
            sys_dias{i} = 'systole'; % Use curly braces {} for cell array indexing
        else
            sys_dias{i} = 'diastole'; % Use curly braces {} for cell array indexing
        end
    % end
end

%%
figure
set(gcf, 'Position',  [0, 0, 1000, 400])
plot(t, normalised_signal)
hold on

plot(t, g, LineWidth=1.5)
title('PCG Signal with Gate Signal', 'Interpreter','latex', FontSize=18);
lgd = legend('', 'Gate Signal'); % adaptive thresholding rule
fontsize(lgd,15,'points')
axis tight
ylim([-1.1 1.1])
ylabel('N. Amplitude', 'Interpreter','latex', FontSize=16)
ax = gca;
ax.FontSize = 15;  % Font Size of 15
set(gca,'TickLabelInterpreter','latex')
set(groot,'defaultLegendInterpreter','latex');

%% Phase (doesnt really work)
% figure
% plot(t, signal)
% hold on
% phi = atan(zn./SEenv);
% plot(t, phi)
