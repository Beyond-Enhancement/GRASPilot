function [ca, info] = graspilot_vcai(A, N, s, t, bf, opts)
%GRASPILOT_VCAI Venous-calibrated arterial input (VCAI) for RATIO.
%
%   [ca, info] = graspilot_vcai(A, N, s, t, baseFrames, opts)
%
%   A     [1 x T] mean signal of an arterial ROI (MCA or distal ICA)
%   N     [1 x T] mean signal of a ring of neighbouring tissue (partial-volume
%         correction), or [] for a static background (N/N0 = 1)
%   s     [1 x T] sinus anchor: phase difference of the superior sagittal
%         sinus from its baseline (rad; linear in concentration), or the
%         sinus magnitude converted to concentration (mM)
%   t     [1 x T] frame times (s);  baseFrames: pre-contrast frames
%   opts  T10Blood (s, default 1.664), TR (s), FlipAngle (deg), r1
%         (1/(mM s)): blood linearization; WRange [0.05 1], TauRange
%         [0.5 30] s, DelayRange [-5 25] s
%
%   ca    [1 x T] arterial blood concentration (mM): shape calibrated
%         against the sinus; amplitude from the arterial curve linearized
%         with the blood T10 and the fitted blood fraction w (plasma =
%         ca / (1 - Hct))
%   info  w (signal-weighted blood fraction), Tau (s, venous dispersion),
%         Delta (s, artery -> sinus delay), Gain g, R2 of the sinus fit,
%         SinusModel [1 x T] = g [k_tau * ca](t - Delta), Ok
%
%   Model (RATIO, M. Lee, unpublished, 2026, Section 3):
%       A(t)/A0 = w Sb(ca(t))/Sb(0) + (1 - w) N(t)/N0
%       s(t) - s0 ~ g [k_tau * ca](t - Delta),  k_tau(t) = exp(-t/tau)/tau
%   Sb is the spoiled gradient echo blood signal. w, tau and Delta are found
%   on a grid and refined with fminsearch; g is solved linearly.

if nargin < 6, opts = struct(); end
def = struct('T10Blood', 1.664, 'TR', 4e-3, 'FlipAngle', 12, 'r1', 4.5, 'WRange', [0.05 1], ...
    'TauRange', [0.5 30], 'DelayRange', [-5 25]);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end
end
A = double(A(:)'); s = double(s(:)'); t = double(t(:)');
T = numel(t);
bf = bf(bf >= 1 & bf <= T);
if isempty(bf), bf = 1; end
rA = A / mean(A(bf));
if isempty(N)
    rN = ones(1, T);
else
    N = double(N(:)');
    rN = N / mean(N(bf));
end
s = s - mean(s(bf));
o = opts;
model = @(w, tau, dl) dispModel(rA, rN, t, w, tau, dl, o);

% grid
wg = linspace(max(o.WRange(1), 0.02), o.WRange(2), 20);
tg = logspace(log10(o.TauRange(1)), log10(o.TauRange(2)), 12);
dg = o.DelayRange(1):1:o.DelayRange(2);
best = Inf; x0 = [wg(end) tg(1) 0];
ss = sum(s .^ 2);
for w = wg
    ca = toConc(rA, rN, w, o);
    if any(~isfinite(ca)), continue; end
    M = graspilot_expconv(ca, t, 1 ./ tg') ./ repmat(tg', 1, T);     % [ntau x T]
    for dl = dg
        Ms = shiftRows(M, t, dl);
        g = (Ms * s') ./ sum(Ms .^ 2, 2);
        sse = ss - g .* (Ms * s');
        [m, i] = min(sse);
        if m < best, best = m; x0 = [w tg(i) dl]; end
    end
end
% refine (bounded by transforms)
lo = [max(o.WRange(1), 0.02) o.TauRange(1) o.DelayRange(1)];
hi = [o.WRange(2) o.TauRange(2) o.DelayRange(2)];
to = @(x) [log(x(1:2) - lo(1:2) + 1e-9) - log(hi(1:2) - x(1:2) + 1e-9), x(3)];
from = @(z) [lo(1:2) + (hi(1:2) - lo(1:2)) ./ (1 + exp(-z(1:2))), min(hi(3), max(lo(3), z(3)))];
obj = @(z) sinusSSE(model, from(z), s);
z = to(x0);
if isfinite(best)
    try
        z = fminsearch(obj, z, optimset('Display', 'off', 'MaxFunEvals', 600, 'TolX', 1e-4, 'TolFun', 1e-10));
    catch
    end
end
x = from(z);
if sinusSSE(model, x, s) > best, x = x0; end
[sse, g, Ms] = sinusSSE(model, x, s);
ca = toConc(rA, rN, x(1), o);
info = struct('w', x(1), 'Tau', x(2), 'Delta', x(3), 'Gain', g, 'R2', 1 - sse / ss, ...
    'SinusModel', g * Ms, 'Ok', isfinite(sse) && all(isfinite(ca)));
end

function [sse, g, Ms] = sinusSSE(model, x, s)
Ms = model(x(1), x(2), x(3));
if any(~isfinite(Ms)) || ~any(Ms)
    sse = Inf; g = NaN;
    return
end
g = (Ms * s') / sum(Ms .^ 2);
sse = sum((s - g * Ms) .^ 2);
end

function Ms = dispModel(rA, rN, t, w, tau, dl, o)
ca = toConc(rA, rN, w, o);
Ms = shiftRows(graspilot_expconv(ca, t, 1 / tau) / tau, t, dl);
end

function ca = toConc(rA, rN, w, o)
% blood signal ratio Sb(ca)/Sb(0) from the partial-volume model -> concentration
rb = (rA - (1 - w) * rN) / w;
ca = graspilot_signal2conc(rb, 1, o.T10Blood, o.TR, o.FlipAngle, o.r1);
end

function Y = shiftRows(X, t, d)
% rows of X delayed by d seconds (linear interpolation, edge values held)
T = numel(t);
tq = min(max(t - d, t(1)), t(end));
fi = interp1(t, 1:T, tq);
i0 = min(max(floor(fi), 1), T - 1);
w = fi - i0;
Y = X(:, i0) .* repmat(1 - w, size(X, 1), 1) + X(:, i0 + 1) .* repmat(w, size(X, 1), 1);
end
