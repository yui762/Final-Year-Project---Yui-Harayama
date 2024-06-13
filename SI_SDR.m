function [SI_SDR] = SI_SDR(reference_signal, estimated_signal)

% reference in here is the clean simulated signal NOT the reference signal in the DSS
% context!

    % Ensure the signals are column vectors
    if isrow(reference_signal)
        reference_signal = reference_signal';
    end
    if isrow(estimated_signal)
        estimated_signal = estimated_signal';
    end

    % Ensure the signals have the same length
    len = min(length(reference_signal), length(estimated_signal));
    reference_signal = reference_signal(1:len);
    estimated_signal = estimated_signal(1:len);

    % Compute the optimal scaling factor
    alpha = (reference_signal' * estimated_signal) / (reference_signal' * reference_signal);

    % Decompose the estimated signal
    s_target = alpha * reference_signal;
    e_noise = estimated_signal - s_target;

    % Compute the power of the target and noise
    s_target_power = sum(s_target .^ 2);
    e_noise_power = sum(e_noise .^ 2);

    % Compute the SI-SDR
    SI_SDR = 10 * log10(s_target_power / e_noise_power);
end
