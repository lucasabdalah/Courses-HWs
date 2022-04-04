function [y,e,w] = lms_regular(x, L, mu, start_state)
% lms_regular - Compute the classical LMS algorithm
%
% Syntax: [y,e,w] = lms_regular(x, L, mu, start_state)
%
% 2021/08/26 - Lucas Abdalah

N = length(x); 

if isequal(start_state, 'zeros')
    w = zeros(N - L, L);    
end

if isequal(start_state, 'random')
    w = randn(N - L, L);
end

d = x(L + 1:end);
e = zeros(N - L, 1);
y = zeros(N - L, 1);

% LMS - Update the weights 
for ii = 1:1:(N - L - 1)
    y(ii) = w(ii,:) * x(ii:ii+L-1).';
    e(ii) = d(ii) - y(ii); 
    w(ii + 1, :) = w(ii, :) + 2 * mu * e(ii) * x(ii:ii+L-1);
end 

end