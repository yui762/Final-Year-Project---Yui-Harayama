function [peak_locs_temp] = FindPeaks(zn)
% peaks correspond to the zero crossing from positive to negative of the HT
% of s(n)
% odd-symmetry function, find the zero cross points

peak_locs_temp = [];
for i = 2:length(zn)-1
    if zn(i-1) >= 0 && zn(i+1) <= 0 % from positive to negative
        if ismember(i-1,peak_locs_temp) == 0
            peak_locs_temp = [peak_locs_temp, i];
        end
    end
end