function th = graspilot_varpro_lm(f, th, lo, hi, maxIter, h, maxStep)
%GRASPILOT_VARPRO_LM Levenberg-Marquardt for many curves at once (variable projection).
%
%   th = graspilot_varpro_lm(f, th, lo, hi, maxIter, h, maxStep)
%
%   f        [sse, coef, r] = f(th, rows): fit of curves `rows` for nonlinear
%            parameters th(rows, :) with the linear coefficients re-solved
%            (sse Inf where th is not admissible); r = residual vectors
%   th       [N x M] start values (rows with NaN are left alone)
%   lo, hi   [1 x M] bounds; h [1 x M] forward-difference steps;
%            maxStep [1 x M] largest change per iteration
%   Every curve has its own damping; a step is kept only if it lowers sse.
[N, M] = size(th);
act = all(isfinite(th), 2);
if ~any(act) || maxIter < 1, return; end
lam = 1e-2 * ones(N, 1);
for it = 1:maxIter
    idx = find(act);
    if isempty(idx), break; end
    n = numel(idx);
    [s0, ~, r0] = f(th(idx, :), idx);
    good = isfinite(s0);
    r0(~isfinite(r0)) = 0;
    T = size(r0, 2);
    J = zeros(n, T, M);
    for q = 1:M
        tq = th(idx, :);
        tq(:, q) = tq(:, q) + h(q);
        [~, ~, rq] = f(tq, idx);
        J(:, :, q) = (rq - r0) / h(q);
    end
    J(~isfinite(J)) = 0;
    A = zeros(n, M, M); g = zeros(n, M);
    for p = 1:M
        g(:, p) = sum(J(:, :, p) .* r0, 2);
        for q = p:M
            A(:, p, q) = sum(J(:, :, p) .* J(:, :, q), 2);
            A(:, q, p) = A(:, p, q);
        end
    end
    for p = 1:M
        A(:, p, p) = A(:, p, p) .* (1 + lam(idx)) + 1e-12;
    end
    d = batchSolve(A, -g);
    d(~isfinite(d)) = 0;
    d = bsxfun(@max, -maxStep, bsxfun(@min, maxStep, d));
    tn = bsxfun(@max, lo, bsxfun(@min, hi, th(idx, :) + d));
    sn = f(tn, idx);
    better = sn < s0 & good;
    th(idx(better), :) = tn(better, :);
    lam(idx(better)) = max(1e-6, lam(idx(better)) / 3);
    lam(idx(~better)) = lam(idx(~better)) * 4;
    done = (better & (s0 - sn) <= 1e-8 * s0) | (~better & lam(idx) > 1e3) | ~good;
    act(idx(done)) = false;
end
end

function x = batchSolve(A, b)
% Gaussian elimination for N small systems A(i, :, :) x = b(i, :) (SPD, damped)
M = size(b, 2);
for k = 1:M
    piv = A(:, k, k);
    piv(abs(piv) < 1e-300) = 1e-300;
    for i = k + 1:M
        fct = A(:, i, k) ./ piv;
        for j = k:M
            A(:, i, j) = A(:, i, j) - fct .* A(:, k, j);
        end
        b(:, i) = b(:, i) - fct .* b(:, k);
    end
end
x = zeros(size(b));
for k = M:-1:1
    s = b(:, k);
    for j = k + 1:M
        s = s - A(:, k, j) .* x(:, j);
    end
    x(:, k) = s ./ A(:, k, k);
end
end
