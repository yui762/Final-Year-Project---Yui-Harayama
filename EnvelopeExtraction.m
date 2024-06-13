function [absolute_n, energy_n, shannon_energy_n, shannon_entropy_n] = EnvelopeExtraction(x)

    % ------------------------ Absolute Value -----------------------%
    absolute = abs(x);
    % Normalise
    min_val_abs = min(absolute);
    max_val_abs = max(absolute);
    absolute_n = (absolute - min_val_abs) / (max_val_abs - min_val_abs);

    % ---------------------------- Energy ---------------------------%
    energy = x.^2; % Square each sample to compute energy
    % normalise
    min_val_energy = min(energy);
    max_val_energy = max(energy);
    energy_n = (energy - min_val_energy) / (max_val_energy - min_val_energy);

    % ----------------------- Shannon Energy ------------------------%
    ShannonEnergy = @(x) x.^2 .* log(x.^2);
    SEdata = ShannonEnergy(x - mean(x));
    
    % normalise
    min_val_shen = min(SEdata);
    max_val_shen = max(SEdata);
    shannon_energy_n = (SEdata - min_val_shen) / (max_val_shen - min_val_shen);
    % shannon_energy_n = -shannon_energy_n + 1; % keep it commented out and
    % adjust in the main as scaling + biasing differ 
    % SEenv = envelope(-shen_n, 100, 'peak');
    
    % ----------------------- Shannon Entropy ------------------------%
    absolute = abs(x);
    shannon_entropy = -absolute.*log(absolute + eps);
    min_val_abs = min(shannon_entropy);
    max_val_abs = max(shannon_entropy);
    % normalise
    shannon_entropy_n = (shannon_entropy - min_val_abs) / (max_val_abs - min_val_abs);


  