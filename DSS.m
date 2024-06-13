function [yy] = DSS(alpha, beta, lambda, NFFT, signal1, signal2, fs1, fs2)

window = hanning(NFFT); % hanning or hamming
overlap = floor(0.5 * NFFT); % number of windows samples without overlapping, the larger the value 
% the worse (actl depends, magnitude became smaller tho) - 0.5 = best
% result, also matches conventional overlapping value

% --------------------------------------------------------------------
% Construct spectrogram for heart sound + noise (signal1) 
% window NFFT = to alleviate discontinuity at endpoints of each segment
% spectro = freq vs time
% --------------------------------------------------------------------
spectro_signal1 = spectrogram(signal1, window, overlap, NFFT, fs1); 
power_signal1 = abs(spectro_signal1).^2; % Power
N = size(spectro_signal1, 2);

% --------------------------------------------------------------------
%       Construct spectrogram of noise only (signal2)
% --------------------------------------------------------------------
spectro_signal2 = spectrogram(signal2, window, overlap, NFFT, fs2);
power_signal2 = abs(spectro_signal2).^2; % Power

% Estimate a priori SNR for heart sound + noise (signal1)
spectro_noise = max(0, (power_signal1 ./ power_signal2) - 1); % take max for when mic ref only collects a noise but not primary

% Smooth the a priori SNR estimate (low-pass)
spectro_noise = filter(1 - alpha, 1, spectro_noise); % filter(b, a, ) this doesn't affect much the results
% est_spectro_noise = lowpass(spectro_noise,alpha);

% attenuation factor g(omega) 
g_omega = max(0, (1 - lambda * (1 ./ (spectro_noise + 1)).^beta)); % (signal 1 power spectrum - noise power spect)/ signal1 power spec
Y_omega = spectro_signal1.*g_omega; 

% IFFT of denoised signal Y(omega) --> y(n)
yy = zeros((N - 1) * overlap + NFFT, 1);
for i = 1:N
    y_ifft = real(ifft(Y_omega(:, i), NFFT));
    yy(overlap*(i-1)+(1:NFFT)) = yy(overlap*(i-1)+(1:NFFT))+y_ifft.* window;
end

