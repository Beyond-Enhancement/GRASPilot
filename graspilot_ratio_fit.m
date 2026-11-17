function [R, Fit] = graspilot_ratio_fit(C, t, u, opts)
%GRASPILOT_RATIO_FIT RATIO kinetic core: scale-free two-compartment exchange with model selection.
%
%   R = graspilot_ratio_fit(C, t, u)
%   [R, Fit] = graspilot_ratio_fit(C, t, u, opts)
%
%   RATIO (Rate-based, Amplitude-invariant, Transit-resolved kinetics with
%   an Input anchored to the venous Outflow; M. Lee, unpublished, 2026).
%   Every curve has its own free gain G and delay d:
%       two-compartment (2C):  y(t) = G [h(.; kF, E, rho) * u](t - d)
%       one-compartment (1C):  y(t) = G [k exp(-k .) * u](t - d)
%   h is the two-compartment exchange impulse response per unit plasma
%   volume, written in the scale-free parameters kF = Fp/vp, E = PS/(Fp+PS)
%   and rho = kE/kF (kE = PS/ve), with rho <= RhoMax (interstitial exchange
%   slower than vascular transit: fixes which exponential is vascular).
%   The gain absorbs vp, T10, B1, relaxivity and the input amplitude, so
%   the primary outputs are gain-free rates. AICc selects 2C or 1C per curve.
%
%   C     [N x T] tissue curves (concentration with an assumed T10; any scale)
%   t     [1 x T] frame times (s)
%   u     [1 x T] input (plasma concentration shape; e.g. the VCAI)
%   opts  struct, all optional:
%     Basis      [T x n] orthonormal temporal basis: fit the basis
%                coefficients (n observations) instead of the frames;
%                [] (default) = frames. AICc uses n, so the 2C model
%                (p = 5) needs n >= 7 (rank rule)
%     DelayRange [-10 20] s (positive = tissue later than the input),
%     DelayStep  2 s (grid; refined continuously)
%     KRange     [1e-3 200] 1/min, rate constants searched
%     RhoMax     0.5;  EThresh 0.03 (leakage call);  KFast 4, KSlow 1.5
%                (1/min: one-compartment classes)
%     Refine     true (Levenberg-Marquardt after the grid), MaxIter 15
%
%   R fields, [N x 1] each (NaN where not reported):
%     Class   1 two-compartment, 2 one-compartment fast (k >= KFast),
%             3 one-compartment slow (k <= KSlow), 4 one-compartment
%             indeterminate (transit and leakage cannot be separated)
%     kep     rho kF (1 - E) (2C) or k (slow)                    1/min
%     kpe     kF E = E/MTT (2C)                                  1/min
%     MTT     1/kF (2C) or 1/k (fast): vascular mean transit time  s
%     E       extraction fraction (2C); VeVp = kpe/kep (2C)
%     k1      one-compartment rate (1C classes)                  1/min
%     Leak    1 = leakage (2C with E >= EThresh), 0 = none (2C below,
%             1C fast), NaN = no call (1C slow, indeterminate)
%     G       fitted gain; apparent values (exploratory, scale-dependent):
%     vp      G (2C, fast);  Ktrans = G kpe (2C), G k (slow)  1/min
%     Fp      100 G kF (2C)  ml/100 ml/min;  ve = Ktrans/kep (2C), G (slow)
%     Delay   s;  R2 of the selected model;  AICc1, AICc2;  NObs
%   Fit       [N x T] fitted curves of the selected model (frames)

if nargin < 4 || isempty(opts), opts = struct(); end
def = struct('Basis', [], 'DelayRange', [-10 20], 'DelayStep', 2, 'KRange', [1e-3 200], ...
    'RhoMax', 0.5, 'EThresh', 0.03, 'KFast', 4, 'KSlow', 1.5, 'Refine', true, 'MaxIter', 15, ...
    'NGrid1', 80, 'NGrid2', 36);
fn = fieldnames(def);
for i = 1:numel(fn)
    if ~isfield(opts, fn{i}) || isempty(opts.(fn{i})), opts.(fn{i}) = def.(fn{i}); end
end
if ~isfield(opts, 'Basis'), opts.Basis = []; end

C = double(C);
[N, T] = size(C);
t = double(t(:)');
u = double(u(:)');
tm = (t - t(1)) / 60;
ok = all(isfinite(C), 2);
C(~ok, :) = 0;
Q = double(opts.Basis);
if isempty(Q)
    P = @(X) X; nObs = T;
else
    P = @(X) X * Q; nObs = size(Q, 2);
end
Cq = P(C);
SS = sum(Cq .^ 2, 2);
SST = sum(P(C - repmat(mean(C, 2), 1, T)) .^ 2, 2);
SST(SST <= 0) = NaN;
lk = log(opts.KRange);
dr = opts.DelayRange;
dg = dr(1):opts.DelayStep:dr(2);
if isempty(dg), dg = 0; end

% ---------------------------------------------------------------- one-compartment
kg = logspace(log10(opts.KRange(1)), log10(opts.KRange(2)), opts.NGrid1)';
best = -Inf(N, 1); k1 = NaN(N, 1); d1 = zeros(N, 1);
for j = 1:numel(dg)
    B = P(graspilot_expconv(shiftIn(u, t, dg(j)), tm, kg));
    g = sum(B .^ 2, 2)';
    b = Cq * B';
    a = max(0, b ./ repmat(g, N, 1));
    red = a .* (2 * b - a .* repmat(g, N, 1));
    [mx, i] = max(red, [], 2);
    up = mx > best;
    best(up) = mx(up); k1(up) = kg(i(up)); d1(up) = dg(j);
end
th1 = [log(k1) d1];
f1 = @(th, rows) proj1(Cq(rows, :), SS(rows), t, tm, u, P, th);
if opts.Refine
    th1 = graspilot_varpro_lm(f1, th1, [lk(1) dr(1)], [lk(2) dr(2)], opts.MaxIter, [1e-4 0.02], [2 4]);
end
[sse1, a1] = f1(th1, (1:N)');
k1 = exp(th1(:, 1)); d1 = th1(:, 2);

% ---------------------------------------------------------------- two-compartment
kg2 = logspace(log10(opts.KRange(1)), log10(opts.KRange(2)), opts.NGrid2);
[J2, I2] = meshgrid(1:numel(kg2), 1:numel(kg2));
sel = I2 > J2;
I2 = I2(sel)'; J2 = J2(sel)';
np = numel(I2);
th2 = NaN(N, 3);
[~, jd] = min(abs(repmat(d1, 1, numel(dg)) - repmat(dg, N, 1)), [], 2);
for j = unique(jd)'
    rows = find(jd == j);
    nr = numel(rows);
    B = P(graspilot_expconv(shiftIn(u, t, dg(j)), tm, kg2'));
    G = B * B';
    ng = numel(kg2);
    bb = Cq(rows, :) * B';
    bst = -Inf(nr, 1); bi = zeros(nr, 1);
    blk = max(1, floor(3e6 / nr));
    for p0 = 1:blk:np
        q = p0:min(np, p0 + blk - 1);
        g11 = G(sub2ind([ng ng], I2(q), I2(q)));
        g22 = G(sub2ind([ng ng], J2(q), J2(q)));
        g12 = G(sub2ind([ng ng], I2(q), J2(q)));
        D = repmat(g11 .* g22 - g12 .^ 2, nr, 1);
        b1 = bb(:, I2(q)); b2 = bb(:, J2(q));
        a1x = (repmat(g22, nr, 1) .* b1 - repmat(g12, nr, 1) .* b2) ./ D;
        a2x = (repmat(g11, nr, 1) .* b2 - repmat(g12, nr, 1) .* b1) ./ D;
        red = a1x .* b1 + a2x .* b2;
        feas = feasible2(repmat(kg2(I2(q)), nr, 1), repmat(kg2(J2(q)), nr, 1), a1x, a2x, opts.RhoMax) & D > 0;
        red(~feas) = -Inf;
        [mx, i] = max(red, [], 2);
        up = mx > bst;
        bst(up) = mx(up); bi(up) = q(i(up));
    end
    okr = isfinite(bst);
    th2(rows(okr), :) = [log(kg2(I2(bi(okr))))' log(kg2(J2(bi(okr))))' repmat(dg(j), nnz(okr), 1)];
end
f2 = @(th, rows) proj2(Cq(rows, :), SS(rows), t, tm, u, P, th, opts.RhoMax);
if opts.Refine
    th2 = graspilot_varpro_lm(f2, th2, [lk(1) lk(1) dr(1)], [lk(2) lk(2) dr(2)], opts.MaxIter, ...
        [1e-4 1e-4 0.02], [2 2 4]);
end
[sse2, a2] = f2(th2, (1:N)');
sse2(~all(isfinite(th2), 2)) = Inf;

% ---------------------------------------------------------------- selection and outputs
ic1 = aicc(sse1, SS, nObs, 3);
ic2 = aicc(sse2, SS, nObs, 5);
two = isfinite(sse2) & ic2 < ic1;
nanv = NaN(N, 1);
R = struct('Class', nanv, 'kep', nanv, 'kpe', nanv, 'MTT', nanv, 'E', nanv, 'VeVp', nanv, 'k1', nanv, ...
    'Leak', nanv, 'G', nanv, 'vp', nanv, 'Ktrans', nanv, 'Fp', nanv, 've', nanv, 'Delay', nanv, ...
    'R2', nanv, 'AICc1', ic1, 'AICc2', ic2, 'NObs', nObs * ones(N, 1));
% two-compartment
Kp = exp(max(th2(:, 1:2), [], 2)); Km = exp(min(th2(:, 1:2), [], 2));
[kF, E, rho, G2] = map2(Kp, Km, a2(:, 1), a2(:, 2));
i2 = two;
R.Class(i2) = 1;
R.kep(i2) = rho(i2) .* kF(i2) .* (1 - E(i2));
R.kpe(i2) = kF(i2) .* E(i2);
R.MTT(i2) = 60 ./ kF(i2);
R.E(i2) = E(i2);
R.VeVp(i2) = R.kpe(i2) ./ R.kep(i2);
R.Leak(i2) = double(E(i2) >= opts.EThresh);
R.G(i2) = G2(i2);
R.vp(i2) = G2(i2);
R.Ktrans(i2) = G2(i2) .* R.kpe(i2);
R.Fp(i2) = 100 * G2(i2) .* kF(i2);
R.ve(i2) = R.Ktrans(i2) ./ R.kep(i2);
R.Delay(i2) = th2(i2, 3);
R.R2(i2) = 1 - sse2(i2) ./ SST(i2);
% one-compartment
i1 = ~two & isfinite(sse1) & a1(:, 1) > 0;
G1 = a1(:, 1) ./ k1;
fast = i1 & k1 >= opts.KFast;
slow = i1 & k1 <= opts.KSlow;
ind = i1 & ~fast & ~slow;
R.Class(fast) = 2; R.Class(slow) = 3; R.Class(ind) = 4;
R.k1(i1) = k1(i1);
R.G(i1) = G1(i1);
R.Delay(i1) = d1(i1);
R.R2(i1) = 1 - sse1(i1) ./ SST(i1);
R.MTT(fast) = 60 ./ k1(fast);
R.vp(fast) = G1(fast);
R.Leak(fast) = 0;
R.kep(slow) = k1(slow);
R.ve(slow) = G1(slow);
R.Ktrans(slow) = G1(slow) .* k1(slow);
fn = fieldnames(R);
for i = 1:numel(fn)
    R.(fn{i})(~ok) = NaN;
end
if nargout > 1
    Fit = NaN(N, T);
    for i = find(two & ok)'
        ud = shiftIn(u, t, th2(i, 3));
        Fit(i, :) = a2(i, 1) * graspilot_expconv(ud, tm, Kp(i)) + a2(i, 2) * graspilot_expconv(ud, tm, Km(i));
    end
    for i = find(i1 & ok)'
        Fit(i, :) = a1(i, 1) * graspilot_expconv(shiftIn(u, t, d1(i)), tm, k1(i));
    end
end
end

% =====================================================================
function v = aicc(sse, SS, n, p)
% corrected Akaike information criterion; Inf if n <= p + 1 (rank rule)
v = Inf(size(sse));
if n - p - 1 <= 0, return; end
s = max(sse, 1e-12 * max(SS, eps));
v = n * log(s / n) + 2 * p + 2 * p * (p + 1) / (n - p - 1);
v(~isfinite(sse)) = Inf;
end

function [sse, coef, r] = proj1(Cq, SS, t, tm, u, P, th)
% one-compartment: nonnegative coefficient a = G k for (log k, delay)
F = P(graspilot_expconv(shiftIn(u, t, th(:, 2)), tm, exp(th(:, 1))));
g = sum(F .^ 2, 2);
b = sum(Cq .* F, 2);
a = max(0, b ./ g);
a(~isfinite(a)) = 0;
sse = max(0, SS - a .* (2 * b - a .* g));
sse(~all(isfinite(th), 2)) = Inf;
coef = a;
if nargout > 2
    r = Cq - repmat(a, 1, size(Cq, 2)) .* F;
end
end

function [sse, coef, r] = proj2(Cq, SS, t, tm, u, P, th, rhoMax)
% two-compartment: coefficients of exp(-K+ t), exp(-K- t) for (log K+, log K-, delay)
Kp = exp(max(th(:, 1:2), [], 2)); Km = exp(min(th(:, 1:2), [], 2));
ud = shiftIn(u, t, th(:, 3));
F1 = P(graspilot_expconv(ud, tm, Kp));
F2 = P(graspilot_expconv(ud, tm, Km));
g11 = sum(F1 .^ 2, 2); g22 = sum(F2 .^ 2, 2); g12 = sum(F1 .* F2, 2);
b1 = sum(Cq .* F1, 2); b2 = sum(Cq .* F2, 2);
det = g11 .* g22 - g12 .^ 2;
a1 = (g22 .* b1 - g12 .* b2) ./ det;
a2 = (g11 .* b2 - g12 .* b1) ./ det;
sse = max(0, SS - a1 .* b1 - a2 .* b2);
bad = ~(det > 1e-12 * g11 .* g22) | ~feasible2(Kp, Km, a1, a2, rhoMax) | ~isfinite(sse);
sse(bad) = Inf;
coef = [a1 a2];
coef(bad, :) = NaN;
if nargout > 2
    r = Cq - repmat(a1, 1, size(Cq, 2)) .* F1 - repmat(a2, 1, size(Cq, 2)) .* F2;
end
end

function tf = feasible2(Kp, Km, a1, a2, rhoMax)
% valid two-compartment exchange response with rho = kE/kF <= rhoMax
[kF, E, rho] = map2(Kp, Km, a1, a2);
tf = a1 >= 0 & a2 >= 0 & a1 + a2 > 0 & kF > 0 & E > 0 & E < 1 & rho > 0 & rho <= rhoMax & Kp > Km;
end

function [kF, E, rho, G] = map2(Kp, Km, a1, a2)
% a1 exp(-K+ t) + a2 exp(-K- t) = G kF [A exp(-K+ t) + (1 - A) exp(-K- t)] -> G, kF, E, rho
Fg = a1 + a2;                                 % G kF
c = (a2 .* Kp + a1 .* Km) ./ Fg;              % kE + kPS
kF = Kp + Km - c;                             % Fp / vp
kE = Kp .* Km ./ kF;                          % PS / ve
kPS = c - kE;                                 % PS / vp
E = kPS ./ (kPS + kF);
rho = kE ./ kF;
G = Fg ./ kF;                                 % vp (times every scale factor)
end

function us = shiftIn(u, t, d)
% input delayed by d seconds (d scalar or [m x 1]), linear interpolation,
% held at the first / last value outside the sampled range -> [m x T]
T = numel(t);
d = d(:);
tq = repmat(t, numel(d), 1) - repmat(d, 1, T);
tq = min(max(tq, t(1)), t(end));
fi = interp1(t, 1:T, tq);
fi = reshape(fi, size(tq));
i0 = min(max(floor(fi), 1), T - 1);
w = fi - i0;
us = u(i0) .* (1 - w) + u(i0 + 1) .* w;
end
