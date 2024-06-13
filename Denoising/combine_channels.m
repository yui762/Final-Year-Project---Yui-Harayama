function [combined1, combined2] = combine_channels(signal1, signal2, fs1, fs2)

    % Combine channels together
    numChannels = size(signal1, 2) - 2; % Get the number of channels
    channels = cell(1, numChannels); % Each cell contains an array
    combined1 = zeros(length(signal1),1);
    combined2 = zeros(length(signal2),1);
    
    % Ensure both signals have the same sampling rate
    if fs1 ~= fs2
        error('Sampling rates of the two signals do not match.');
    end
    
    % Ensure both signals have the same length
    if length(signal1) ~= length(signal2)
        error('Sample number of the two signals do not match.');
    end
    
    fig = figure;
    for i = 2:5 % to get the raw mic data only, not the audio processed for voice recognition
        channels{i} = signal1(:, i);
        t = (0:numel(channels{i})-1) / fs1; 

        % Combine relevant channels together
        combined1 = channels{i} + combined1;

        % And plot individual channels (not normalised)
        subplot(4, 2, i-1);
        hold all;
        plot(t, channels{i}, 'color', [0 0.4470 0.7410]);
        title(sprintf('Channel %d', i-1), FontSize=12);
        lgd = legend('Primary');
        fontsize(lgd,9,'points')
    end
    
    for i = 2:5 % to get the raw mic data only, not the audio processed for voice recognition
        channels{i} = signal2(:, i);
        t = (0:numel(channels{i})-1) / fs2; 
        
        % Combine relevant channels together
        combined2 = channels{i} + combined2;

        % And plot individual channels (not normalised)
        subplot(4, 2, i+3);
        hold all;
        plot(t, channels{i}, 'color', [0.8500 0.3250 0.0980]);
        title(sprintf('Channel %d', i-1), FontSize=12);
        lgd = legend('Reference');
        fontsize(lgd,9,'points')
    end

    % Give common xlabel, ylabel and title to your figure
    han=axes(fig,'visible','off'); 
    han.XLabel.Visible='on';
    han.YLabel.Visible='on';
    ylabel(han,'Amplitude');
    xlabel(han,'Time (s)');
