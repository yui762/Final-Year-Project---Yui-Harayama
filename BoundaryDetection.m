function [zero_to_one_transitions, one_to_zero_transitions] = BoundaryDetection(g)

zero_to_one_transitions = [];
one_to_zero_transitions = [];

% Iterate through the signal to detect transitions
for i = 2:length(g)
    if g(i) == 1 && g(i-1) == 0
        zero_to_one_transitions = [zero_to_one_transitions, i];
    elseif g(i) == 0 && g(i-1) == 1
        one_to_zero_transitions = [one_to_zero_transitions, i];
    end
end