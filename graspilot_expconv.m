function F = graspilot_expconv(u, t, k)
%GRASPILOT_EXPCONV Convolution with exp(-k t) of an input that is linear between samples.
%
%   F = graspilot_expconv(u, t, k)
%
%   u   input: [1 x T] (shared) or [n x T] (one per rate constant)
%   t   [1 x T] sample times (increasing, uneven allowed)
%   k   [n x 1] rate constants, in the inverse unit of t
%   F   [n x T], F(i, j) = int_t(1)^t(j) u_i(s) exp(-k_i (t(j) - s)) ds
%
%   Exact for u linear between the samples (recursive update per interval);
%   k = 0 gives the running integral of u.
k = k(:);
n = numel(k);
T = numel(t);
F = zeros(n, T);
if T < 2, return; end
d = diff(t(:)');
uni = all(abs(d - d(1)) <= 1e-6 * max(abs(d(1)), eps));
if uni, [E, g1, g2] = coefs(k * d(1)); end
for i = 1:T - 1
    if ~uni, [E, g1, g2] = coefs(k * d(i)); end
    F(:, i + 1) = E .* F(:, i) + d(i) * (u(:, i) .* g1 + (u(:, i + 1) - u(:, i)) .* g2);
end
end

function [E, g1, g2] = coefs(x)
% E = exp(-x), g1 = (1 - E)/x, g2 = (1 - g1)/x (series for small x)
E = exp(-x);
g1 = (1 - E) ./ x;
g2 = (1 - g1) ./ x;
s = abs(x) < 1e-3;
xs = x(s);
g1(s) = 1 - xs / 2 + xs .^ 2 / 6;
g2(s) = 0.5 - xs / 6 + xs .^ 2 / 24;
end
