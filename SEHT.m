function [smoothed2_SE, out] = SEHT(SE, N, M)

% -----------------------------------------------------------------------------%
% zero-phase filtering with rectangular impulse response that provides a smooth envelogram of the PCG signal. 
% The decision stage is based on Hilbert transform and moving average
% filter. N = length of the rectangular window 
% -----------------------------------------------------------------------------%

% two pass smoothing 
% rect_filtre = rectwin(N);
% smoothed1_SE = conv(SE, rect_filtre, 'same');
% smoothed1_SE = fliplr(smoothed1_SE);
% smoothed2_SE = conv(smoothed1_SE,rect_filtre,'same');
% smoothed2_SE = fliplr(smoothed2_SE);

% Define the filter coefficients for the rectangular window
rect_filter = rectwin(N);
% Normalize the filter coefficients
rect_filter = rect_filter / sum(rect_filter);

% Apply zero-phase filtering using filfilt
smoothed2_SE = filtfilt(rect_filter, 1, SE);


% Hilbert transform of the Shannon energy envelope s(n) but looking at code
% seems like its more the HT of the shannon energy not envelope
ht = imag(hilbert(smoothed2_SE));

% MA filtre, M = length of filtre: about 3 times sampling rate (drift)
MA_filtre = (1/M)*ones(1,M);
MA_out = filter(MA_filtre,1,ht); % num, denom, input
out = ht - MA_out;
% out = ht;