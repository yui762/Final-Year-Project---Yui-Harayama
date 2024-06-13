function[noisy_signal] = AddWGN(x, desired_snr_db)

% 1. Compute the power of the clean signal
clean_power = sum(x.^2) / length(x);

% 2. Calculate the desired noise power based on the desired SNR
desired_snr_linear = 10^(desired_snr_db / 10); % Convert SNR from dB to linear scale
desired_noise_power = clean_power / desired_snr_linear;

% 3. Generate white Gaussian noise with the calculated power
noise = sqrt(desired_noise_power) * randn(size(x)); % power = variance of WGN

% 4. Add the noise to the clean signal
noisy_signal = x + noise;
