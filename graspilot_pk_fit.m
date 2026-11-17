function [R, Fit] = graspilot_pk_fit(C, t, Cp, models, opts)
%GRASPILOT_PK_FIT Pharmacokinetic model fits: standard Tofts, extended Tofts, 2CXM.
%
%   R = graspilot_pk_fit(C, t, Cp)
%   [R, Fit] = graspilot_pk_fit(C, t, Cp, models, opts)
%
%   C       [N x T] tissue concentration curves (mM), one per row (voxels or
%           ROI means), e.g. from graspilot_signal2conc
%   t       [1 x T] frame times (s); the model time origin is t(1)
%   Cp      [1 x T] arterial (or venous) plasma concentration (mM): the AIF
%   models  any of 'tofts', 'etofts', '2cxm' (char or cell; default all)
%   opts    struct, all fields optional:
%     Refine   true (default): refine the grid solution by Levenberg-
%              Marquardt; false: grid search only
%     MaxIter  Levenberg-Marquardt iterations (default 15)
%     KRange   [min max] range of the rate constants searched (1/min);
%              default [1e-3 50] for kep, [1e-3 200] for the 2CXM
%     NGrid    grid points (log-spaced) for the rate search: default 120
%              (Tofts models), 45 per rate (2CXM)
%
%   R.tofts   Ktrans (1/min), kep (1/min), ve, R2, SSE           [N x 1] each
%   R.etofts  Ktrans, kep, ve, vp, R2, SSE
%   R.cxm     Fp (ml/100 ml/min), PS (ml/100 ml/min), vp, ve,
%             Ktrans = E Fp with E = PS/(PS + Fp) (1/min), E, R2, SSE
%   Fit.<model> [N x T] fitted curves (only if requested)
%
%   Models (Sourbron & Buckley, Phys Med Biol 2013;58:861):
%     Tofts           C(t) = Ktrans int Cp(u) exp(-kep (t - u)) du,  ve = Ktrans / kep
%     extended Tofts  C(t) = vp Cp(t) + Ktrans int Cp(u) exp(-kep (t - u)) du
%     2CXM            C(t) = Fp int Cp(u) [A exp(-K+ (t - u)) + (1 - A) exp(-K- (t - u))] du
%                     K+-: rate constants of the plasma (vp) and interstitial (ve)
%                     compartments exchanging with PS, plasma flow Fp.
%   Convolutions are exact for an AIF that is linear between frames (also
%   for uneven frame times). The rate constants are found on a log-spaced
%   grid with the linear coefficients (Ktrans, vp, Fp A, Fp (1 - A))
%   solved exactly at every grid point (variable projection; nonnegative
%   for the Tofts models, physically valid 2CXM parameters only), then
%   refined per curve by Levenberg-Marquardt. Curves with NaN values give
%   NaN parameters; ve and vp are fractions (ml/ml).

if nargin < 4 || isempty(models), models = {'tofts', 'etofts', '2cxm'}; end
if ischar(models), models = {models}; end
if nargin < 5 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'Refine') || isempty(opts.Refine), opts.Refine = true; end
if ~isfield(opts, 'MaxIter') || isempty(opts.MaxIter), opts.MaxIter = 15; end
if ~isfield(opts, 'KRange'), opts.KRange = []; end
if ~isfield(opts, 'NGrid'), opts.NGrid = []; end

C = double(C);
[N, T] = size(C);
t = double(t(:)');
tm = (t - t(1)) / 60;                       % minutes from the first frame
Cp = double(Cp(:)');
if numel(tm) ~= T || numel(Cp) ~= T
    error('graspilot_pk_fit:size', 'C, t and Cp must have the same number of frames.');
end
ok = all(isfinite(C), 2);
C(~ok, :) = 0;
SST = sum((C - repmat(mean(C, 2), 1, T)) .^ 2, 2);
SST(SST <= 0) = NaN;
SS = sum(C .^ 2, 2);
R = struct();
Fit = struct();
wantFit = nargout > 1;
for m = 1:numel(models)
    switch lower(models{m})
        case 'tofts'
            [rate, sse, coef] = fitTofts(C, SS, tm, Cp, false, opts);
            P.Ktrans = coef(:, 1); P.kep = rate; P.ve = P.Ktrans ./ P.kep;
            fld = 'tofts';
        case 'etofts'
            [rate, sse, coef] = fitTofts(C, SS, tm, Cp, true, opts);
            P.Ktrans = coef(:, 1); P.kep = rate; P.ve = P.Ktrans ./ P.kep; P.vp = coef(:, 2);
            fld = 'etofts';
        case {'2cxm', 'cxm'}
            [rate, sse, coef] = fit2CXM(C, tm, Cp, opts);
            P = cxmParams(rate(:, 1), rate(:, 2), coef(:, 1), coef(:, 2));
            fld = 'cxm';
        otherwise
            error('graspilot_pk_fit:model', 'Unknown model "%s".', models{m});
    end
    P.ve(P.Ktrans == 0) = 0;
    P.R2 = 1 - sse ./ SST;
    P.SSE = sse;
    fn = fieldnames(P);
    for i = 1:numel(fn)
        P.(fn{i})(~ok | ~isfinite(sse)) = NaN;
    end
    R.(fld) = P;
    if wantFit
        switch fld
            case 'tofts', Y = repmat(coef(:, 1), 1, T) .* expconv(Cp, tm, rate);
            case 'etofts', Y = repmat(coef(:, 1), 1, T) .* expconv(Cp, tm, rate) + coef(:, 2) * Cp;
            otherwise
                Y = repmat(coef(:, 1), 1, T) .* expconv(Cp, tm, rate(:, 1)) + ...
                    repmat(coef(:, 2), 1, T) .* expconv(Cp, tm, rate(:, 2));
        end
        Y(~ok | ~isfinite(sse), :) = NaN;
        Fit.(fld) = Y;
    end
    clear P
end
end

% =====================================================================
function [k, sse, coef] = fitTofts(C, SS, tm, Cp, ext, opts)
% Tofts (ext false) / extended Tofts (ext true): grid search over kep with
% the nonnegative linear coefficients solved exactly, then LM on log(kep)
N = size(C, 1);
kr = opts.KRange; if isempty(kr), kr = [1e-3 50]; end
ng = opts.NGrid; if isempty(ng), ng = 120; end
kg = logspace(log10(kr(1)), log10(kr(2)), ng);
B = expconv(Cp, tm, kg');                    % [ng x T]
g11 = sum(B .^ 2, 2)';                       % [1 x ng]
b1 = C * B';                                 % [N x ng]
if ext
    g22 = sum(Cp .^ 2);
    g12 = (B * Cp')';
    b2 = repmat(C * Cp', 1, ng);
    red = nnls2red(b1, b2, repmat(g11, N, 1), repmat(g12, N, 1), g22);
else
    a = max(0, b1 ./ repmat(g11, N, 1));
    red = a .* (2 * b1 - a .* repmat(g11, N, 1));
end
[~, j] = max(red, [], 2);
k = kg(j)';
if opts.Refine
    th = log(k);
    f = @(th, ii) projTofts(C(ii, :), SS(ii), tm, Cp, exp(th), ext);
    th = lm(f, th, log(kr), opts.MaxIter);
    k = exp(th);
end
[sse, coef] = projTofts(C, SS, tm, Cp, k, ext);
end

function [sse, coef, r] = projTofts(C, SS, tm, Cp, k, ext)
% SSE and nonnegative linear coefficients for per-curve kep (variable projection)
N = size(C, 1);
F = expconv(Cp, tm, k);                      % [N x T]
g11 = sum(F .^ 2, 2);
b1 = sum(C .* F, 2);
if ext
    g22 = sum(Cp .^ 2);
    g12 = F * Cp';
    b2 = C * Cp';
    [red, a1, a2] = nnls2red(b1, b2, g11, g12, g22);
    coef = [a1 a2];
else
    a1 = max(0, b1 ./ g11);
    a1(~isfinite(a1)) = 0;
    red = a1 .* (2 * b1 - a1 .* g11);
    coef = [a1 zeros(N, 1)];
end
sse = max(0, SS - red);
if nargout > 2
    r = C - repmat(coef(:, 1), 1, numel(tm)) .* F - coef(:, 2) * Cp;
end
end

function [red, a1, a2] = nnls2red(b1, b2, g11, g12, g22)
% Reduction of the SSE, ||C||^2 - ||C - a1 F1 - a2 F2||^2, for the
% nonnegative least-squares coefficients of two basis curves (element-wise
% for arrays of normal equations [g11 g12; g12 g22] [a1; a2] = [b1; b2])
g22 = g22 + zeros(size(b1));
det = g11 .* g22 - g12 .^ 2;
u1 = (g22 .* b1 - g12 .* b2) ./ det;
u2 = (g11 .* b2 - g12 .* b1) ./ det;
inside = u1 >= 0 & u2 >= 0 & det > 0;
s1 = max(0, b1 ./ g11); s1(~isfinite(s1)) = 0;
s2 = max(0, b2 ./ g22); s2(~isfinite(s2)) = 0;
r1 = s1 .* b1;                               % one-coefficient fits: a = b/g -> a b
r2 = s2 .* b2;
use1 = r1 >= r2;
a1 = s1 .* use1; a2 = s2 .* ~use1;
red = max(r1, r2);
a1(inside) = u1(inside); a2(inside) = u2(inside);
red(inside) = u1(inside) .* b1(inside) + u2(inside) .* b2(inside);
end

% =====================================================================
function [K, sse, coef] = fit2CXM(C, tm, Cp, opts)
% 2CXM: grid over the pair of rate constants K+ > K- with the coefficients
% (Fp A, Fp (1 - A)) solved exactly; only pairs that map to positive Fp,
% PS, vp, ve are accepted. Then LM on log(K+), log(K-).
[N, T] = size(C);
kr = opts.KRange; if isempty(kr), kr = [1e-3 200]; end
ng = opts.NGrid; if isempty(ng), ng = 45; end
kg = logspace(log10(kr(1)), log10(kr(2)), ng);
B = expconv(Cp, tm, kg');
G = B * B';
bb = C * B';                                  % [N x ng]
[J, I] = meshgrid(1:ng, 1:ng);
sel = I > J;                                  % K+ = kg(I) > K- = kg(J)
I = I(sel)'; J = J(sel)';
np = numel(I);
best = -Inf(N, 1); bi = ones(N, 1);
blk = max(1, floor(4e6 / max(1, N)));          % pairs per block (memory)
for p0 = 1:blk:np
    q = p0:min(np, p0 + blk - 1);
    g11 = G(sub2ind([ng ng], I(q), I(q)));
    g22 = G(sub2ind([ng ng], J(q), J(q)));
    g12 = G(sub2ind([ng ng], I(q), J(q)));
    det = g11 .* g22 - g12 .^ 2;
    b1 = bb(:, I(q)); b2 = bb(:, J(q));
    D = repmat(det, N, 1);
    a1 = (repmat(g22, N, 1) .* b1 - repmat(g12, N, 1) .* b2) ./ D;
    a2 = (repmat(g11, N, 1) .* b2 - repmat(g12, N, 1) .* b1) ./ D;
    red = a1 .* b1 + a2 .* b2;
    feas = cxmFeasible(repmat(kg(I(q)), N, 1), repmat(kg(J(q)), N, 1), a1, a2) & D > 0;
    red(~feas) = -Inf;
    [mx, j] = max(red, [], 2);
    up = mx > best;
    best(up) = mx(up);
    bi(up) = q(j(up));
end
K = [kg(I(bi))' kg(J(bi))'];
noFit = ~isfinite(best);
if opts.Refine
    th = log(K);
    SS = sum(C .^ 2, 2);
    f = @(th, ii) proj2CXM(C(ii, :), SS(ii), tm, Cp, exp(th));
    th(noFit, :) = NaN;
    th = lm(f, th, log(kr), opts.MaxIter);
    K = exp(th);
end
K = [max(K, [], 2) min(K, [], 2)];
[sse, coef] = proj2CXM(C, sum(C .^ 2, 2), tm, Cp, K);
sse(noFit | ~isfinite(sse)) = NaN;
coef(noFit, :) = NaN;
K(noFit, :) = NaN;
end

function [sse, coef, r] = proj2CXM(C, SS, tm, Cp, K)
% SSE and coefficients for per-curve (K+, K-); Inf where not a valid 2CXM
N = size(C, 1);
Kp = max(K, [], 2); Km = min(K, [], 2);
F1 = expconv(Cp, tm, Kp);
F2 = expconv(Cp, tm, Km);
g11 = sum(F1 .^ 2, 2); g22 = sum(F2 .^ 2, 2); g12 = sum(F1 .* F2, 2);
b1 = sum(C .* F1, 2); b2 = sum(C .* F2, 2);
det = g11 .* g22 - g12 .^ 2;
a1 = (g22 .* b1 - g12 .* b2) ./ det;
a2 = (g11 .* b2 - g12 .* b1) ./ det;
sse = max(0, SS - a1 .* b1 - a2 .* b2);
bad = ~(det > 1e-12 * g11 .* g22) | ~cxmFeasible(Kp, Km, a1, a2) | ~isfinite(sse);
sse(bad) = Inf;
coef = [a1 a2];
coef(bad, :) = NaN;
if nargout > 2
    r = C - repmat(a1, 1, size(C, 2)) .* F1 - repmat(a2, 1, size(C, 2)) .* F2;
end
if N == 0, sse = zeros(0, 1); end
end

function tf = cxmFeasible(Kp, Km, a1, a2)
% coefficients of exp(-K+ t), exp(-K- t) that correspond to Fp, PS, vp, ve > 0
Fp = a1 + a2;
c = (a2 .* Kp + a1 .* Km) ./ Fp;              % kep + kpe
kpo = Kp + Km - c;                            % Fp / vp
kep = Kp .* Km ./ kpo;                        % PS / ve
kpe = c - kep;                                % PS / vp
tf = Fp > 0 & a1 >= 0 & a2 >= 0 & kpo > 0 & kep > 0 & kpe > 0 & Kp > Km;
end

function P = cxmParams(Kp, Km, a1, a2)
Fp = a1 + a2;                                 % ml/ml/min
c = (a2 .* Kp + a1 .* Km) ./ Fp;
kpo = Kp + Km - c;
kep = Kp .* Km ./ kpo;
kpe = c - kep;
vp = Fp ./ kpo;
PS = kpe .* vp;
ve = PS ./ kep;
E = PS ./ (PS + Fp);
P.Fp = 100 * Fp;
P.PS = 100 * PS;
P.vp = vp;
P.ve = ve;
P.Ktrans = E .* Fp;
P.E = E;
end

% =====================================================================
function th = lm(f, th, bounds, maxIter)
% Levenberg-Marquardt on the nonlinear parameters th [N x M] (log rate
% constants) of the projected residual: the linear coefficients are
% re-solved at every evaluation. [sse, coef, r] = f(th, rows) evaluates the
% curves C(rows, :); the Jacobian of r is formed by forward differences.
[N, M] = size(th);
act = all(isfinite(th), 2);
if ~any(act) || maxIter < 1, return; end
lam = 1e-2 * ones(N, 1);
h = 1e-4;
for it = 1:maxIter
    idx = find(act);
    if isempty(idx), break; end
    [s0, ~, r0] = f(th(idx, :), idx);
    good = isfinite(s0);
    T = size(r0, 2);
    Jc = zeros(numel(idx), T, M);
    for q = 1:M
        tq = th(idx, :);
        tq(:, q) = tq(:, q) + h;
        [~, ~, rq] = f(tq, idx);
        Jc(:, :, q) = (rq - r0) / h;
    end
    Jc(~isfinite(Jc)) = 0;
    r0(~isfinite(r0)) = 0;
    L = lam(idx);
    if M == 1
        A11 = sum(Jc(:, :, 1) .^ 2, 2);
        g1 = sum(Jc(:, :, 1) .* r0, 2);
        d = -g1 ./ (A11 .* (1 + L) + 1e-12);
    else
        A11 = sum(Jc(:, :, 1) .^ 2, 2); A22 = sum(Jc(:, :, 2) .^ 2, 2);
        A12 = sum(Jc(:, :, 1) .* Jc(:, :, 2), 2);
        g1 = sum(Jc(:, :, 1) .* r0, 2); g2 = sum(Jc(:, :, 2) .* r0, 2);
        A11d = A11 .* (1 + L) + 1e-12; A22d = A22 .* (1 + L) + 1e-12;
        dt = A11d .* A22d - A12 .^ 2;
        d = [-(A22d .* g1 - A12 .* g2) ./ dt, -(A11d .* g2 - A12 .* g1) ./ dt];
    end
    d(~isfinite(d)) = 0;
    d = max(-2, min(2, d));                    % at most a factor e^2 per step
    tn = max(bounds(1), min(bounds(2), th(idx, :) + d));
    sn = f(tn, idx);
    better = sn < s0 & good;
    th(idx(better), :) = tn(better, :);
    lam(idx(better)) = max(1e-6, lam(idx(better)) / 3);
    lam(idx(~better)) = lam(idx(~better)) * 4;
    done = (better & (s0 - sn) <= 1e-8 * s0) | (~better & lam(idx) > 1e3) | ~good;
    act(idx(done)) = false;
end
end

% =====================================================================
function F = expconv(Cp, tm, k)
% F(k, t) = int_0^t Cp(u) exp(-k (t - u)) du for an input linear between
% the frames; k [n x 1] (1/min), tm (min) -> F [n x T]
k = k(:);
n = numel(k);
T = numel(tm);
F = zeros(n, T);
d = diff(tm);
if T < 2, return; end
uni = all(abs(d - d(1)) <= 1e-6 * max(abs(d(1)), eps));
if uni, [E, g1, g2] = coefs(k * d(1)); end
for i = 1:T - 1
    if ~uni, [E, g1, g2] = coefs(k * d(i)); end
    F(:, i + 1) = E .* F(:, i) + d(i) * (Cp(i) * g1 + (Cp(i + 1) - Cp(i)) * g2);
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
