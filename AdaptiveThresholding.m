function [g] = AdaptiveThresholding(x, shannon_energy_envelope, threshold)

g = zeros(1, length(x));

for i = 1:length(shannon_energy_envelope)
    if shannon_energy_envelope(i)+1 <= threshold
        g(i) = 0;
    else
        g(i) = 1;
    end
end