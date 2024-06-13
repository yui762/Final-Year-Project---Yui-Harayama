function [e, y, w] = ANC_LMS_FYP(primary, reference, step_size, order)
% x = reference 
% d = primary
% e = error (clean) 
% y = estimate of noise

L = length(primary);
x_shifted = zeros(order,1);
weights = zeros(order,1);
y = zeros(L,1); 
e = zeros(L,1);

for n = 1:L
    x_shifted = [x_shifted(2:order);reference(n)];

    % Estimate noise of primary signal using reference noise
    y(n) = weights' * x_shifted; 

    % Error signal 
    e(n) = primary(n) - y(n); % primary - estimate of noise = clean signal

    % Update weights
    weights = weights +  step_size * e(n) * x_shifted;
    w(:,n) = weights;
end
