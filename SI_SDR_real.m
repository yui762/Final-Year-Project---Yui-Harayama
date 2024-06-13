function [SI_SDR] = SI_SDR_real(input_signal, noise_signal)

% input signal = either the primary signal x or denoised signal y

    % Ensure the signals are column vectors
    if isrow(input_signal)
        input_signal = input_signal';
    end
    if isrow(noise_signal)
        noise_signal = noise_signal';
    end

    % Ensure the signals have the same length
    len = min(length(input_signal), length(noise_signal));
    input_signal = input_signal(1:len);
    noise_signal = noise_signal(1:len);

    % Compute the optimal scaling factor
    alpha = (input_signal' * noise_signal) / (input_signal' * input_signal);

    % Decompose the estimated signal
    s_target = alpha * input_signal;
    e_noise = noise_signal;

    % Compute the power of the target and noise
    s_target_power = sum(s_target .^ 2);
    e_noise_power = sum(e_noise .^ 2);

    % Compute the SI-SDR
    SI_SDR = 10 * log10(s_target_power / e_noise_power);
end
