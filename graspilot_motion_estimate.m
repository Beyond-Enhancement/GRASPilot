function M = graspilot_motion_estimate(D, varargin)
%GRASPILOT_MOTION_ESTIMATE Frame-by-frame rigid head motion of a dynamic series.
%
%   M = graspilot_motion_estimate(D)
%   M = graspilot_motion_estimate(D, 'Name', Value, ...)
%
%   Every reconstructed (raw, unbinned) frame is registered rigidly (3
%   translations, 3 rotations) to a reference volume by maximizing the
%   normalized mutual information (Studholme et al. 1999), which tolerates
%   the signal changes caused by contrast enhancement. Registration is
%   multi-resolution (Gaussian smoothing 3 mm, then 1 mm) with a simplex
%   search (fminsearch) and starts from the estimate of the previous frame.
%   The cost uses up to MaxPoints voxels of the head (Region): nothing
%   outside the head, where the streak artifacts of radial (stack-of-stars)
%   sampling rotate from frame to frame. Optionally strongly enhancing voxels
%   (vessels, scalp, dura, enhancing lesions) are left out as well, whose
%   enhancement can be read as motion. The estimates are then smoothed over
%   time with a running median (SmoothFrames), which removes frame-to-frame
%   estimation noise but keeps sudden (step) motion. No toolboxes are needed.
%
%   Geometry: canonical voxel order of GRASPilotData (axis 1 R->L, 2 A->P,
%   3 I->S) in mm. A head point x of the reference is found at
%       y = Rot * (x - Centre) + Centre + [tL tP tS]
%   in frame f, Rot = Rz(rS) * Ry(rP) * Rx(rL) (degrees, rotation about the
%   L-R, A-P and I-S axes through the head centre). Params(f, :) =
%   [tL tP tS rL rP rS]. graspilot_rigid_matrix(Params(f, :)) gives Rot.
%
%   Options (defaults)
%     'Reference'   'baseline' (mean of RefFrames) | raw frame number ('baseline')
%     'RefFrames'   raw frames averaged for 'baseline' (1:min(5, NTRaw))
%     'Frames'      raw frames to register (1:NTRaw)
%     'MaxPoints'   sample voxels used by the cost function (30000)
%     'Bins'        histogram bins for mutual information (32)
%     'Region'      'head' (default: head voxels only), 'brain' (head minus
%                   MarginMM along its outer surface and the first and last
%                   slice: no scalp) or 'outline' (head and a 4 mm band
%                   outside it, the former behaviour)
%     'MarginMM'    8
%     'ExcludeEnhancing'  false; true: leave out voxels whose median signal
%                   in the frames after the reference rises by more than
%                   EnhFrac (0.15) of the median head signal (dilated 1 voxel)
%     'SmoothFrames' running-median length (5; 1 = none)
%     'WarnAbsMM'   warning threshold, mean head displacement from the
%                   reference (one voxel: geometric mean voxel size, >= 1 mm)
%     'WarnRelMM'   warning threshold, mean head displacement between
%                   consecutive frames (half of WarnAbsMM)
%     'Progress'    @(fraction, message) progress report
%     'Cancel'      @() returning true to stop (frames done so far are kept)
%
%   Output M (struct)
%     Params [NTRaw x 6], NaN for frames not registered (smoothed)
%     ParamsRaw   [NTRaw x 6] per-frame estimates before smoothing
%     AbsMean, AbsMax  [NTRaw x 1] mean / max displacement (mm) of the head
%                      voxels relative to the reference
%     RelMean, RelMax  [NTRaw x 1] same, relative to the previous frame
%     Centre, VoxelMM, RefFrames, Reference (text), Frames, NMI
%     WarnAbsMM, WarnRelMM, ExcessFrames, Excess (true when any frame
%     exceeds a threshold), Complete (false if cancelled), Message
%
%   Head motion that happened during the acquisition is also built into a
%   subspace (low-rank) reconstruction itself (blurring, ghosting that the
%   rigid realignment of frames cannot remove). Excess motion should be
%   noted even when the frames are corrected.

p = struct('Reference', 'baseline', 'RefFrames', [], 'Frames', [], 'MaxPoints', 30000, 'Bins', 32, ...
    'Region', 'head', 'MarginMM', 8, 'SmoothFrames', 5, 'ExcludeEnhancing', false, 'EnhFrac', 0.15, ...
    'WarnAbsMM', [], 'WarnRelMM', [], 'Progress', [], 'Cancel', []);
for i = 1:2:numel(varargin)
    fn = fieldnames(p);
    m = find(strcmpi(fn, varargin{i}), 1);
    if isempty(m), error('graspilot_motion_estimate:opt', 'Unknown option %s.', varargin{i}); end
    p.(fn{m}) = varargin{i + 1};
end
T = D.NTRaw;
n = D.canonSize();
vox = D.canonVox();
if isempty(p.RefFrames), p.RefFrames = 1:min(5, T); end
p.RefFrames = p.RefFrames(p.RefFrames >= 1 & p.RefFrames <= T);
if isempty(p.Frames), p.Frames = 1:T; end
vmm = max(1, prod(vox) ^ (1 / 3));
if isempty(p.WarnAbsMM), p.WarnAbsMM = vmm; end
if isempty(p.WarnRelMM), p.WarnRelMM = p.WarnAbsMM / 2; end
report = @(f, msg) 0;
if ~isempty(p.Progress), report = p.Progress; end

% reference volume
report(0, 'Motion: reference volume...');
if ischar(p.Reference)
    rf = p.RefFrames;
    R = zeros(n, 'single');
    for f = rf
        R = R + D.rawCanonVolume(f);
    end
    R = R / numel(rf);
    if numel(rf) > 1
        reftxt = sprintf('mean of raw frames %d-%d', rf(1), rf(end));
    else
        reftxt = sprintf('raw frame %d', rf);
    end
else
    rf = p.Reference;
    R = D.rawCanonVolume(rf);
    reftxt = sprintf('raw frame %d', rf);
end
R = double(R);

% head mask, sample region (head + 4 mm outline), head centre
head = R > otsu(R(1:max(1, floor(numel(R) / 300000)):end));
if nnz(head) < 100
    error('graspilot_motion_estimate:head', 'Could not find the head in the reference volume.');
end
if strcmpi(p.Region, 'outline')
    region = dilate(head, max(1, round(4 ./ vox)));
elseif strcmpi(p.Region, 'brain')
    region = interior(head, vox, p.MarginMM);
    region(:, :, [1 end]) = false;
    if nnz(region) < 2000                       % small head or slab: fall back
        region = head;
    end
else
    region = head;                              % head voxels only: nothing outside the head
end
nexcl = 0;
if p.ExcludeEnhancing && T > numel(rf) + 1
    % median signal change after the reference (robust to transient motion)
    fs = unique(round(linspace(max(rf) + 1, T, min(9, T - max(rf)))));
    Rs = smooth3d(R, 2 ./ vox);
    Ch = zeros([n numel(fs)], 'single');
    for i = 1:numel(fs)
        Ch(:, :, :, i) = single(smooth3d(double(D.rawCanonVolume(fs(i))), 2 ./ vox) - Rs);
    end
    enh = median(Ch, 4) > p.EnhFrac * median(R(head));
    enh = dilate(enh, [1 1 1]);
    if nnz(region & ~enh) >= 2000
        nexcl = nnz(region & enh);
        region = region & ~enh;
    end
    clear Ch Rs
end
hl = find(head);
[a1, a2, a3] = ind2sub(n, hl);
Cmm = [mean(a1 - 1) mean(a2 - 1) mean(a3 - 1)] .* vox;
% displacement evaluation points: up to 20000 head voxels
hs = hl(round(linspace(1, numel(hl), min(numel(hl), 20000))));
[e1, e2, e3] = ind2sub(n, hs);
Xd = ([e1 e2 e3] - 1) .* repmat(vox, numel(hs), 1);

% registration levels: smoothing (mm), points, simplex scale, tolerance
lev = struct('sig', {3, 1}, 'np', {min(8000, p.MaxPoints), p.MaxPoints}, 's', {30, 6}, ...
    'tol', {0.01, 0.004}, 'maxf', {250, 250});
rl = find(region);
for L = 1:numel(lev)
    Rs = smooth3d(R, lev(L).sig ./ vox);
    pick = rl(round(linspace(1, numel(rl), min(numel(rl), lev(L).np))));
    [s1, s2, s3] = ind2sub(n, pick);
    lev(L).X = ([s1 s2 s3] - 1) .* repmat(vox, numel(pick), 1);
    lev(L).pick = pick;
    rv = Rs(pick);
    q = prctile1(rv, [0.5 99.5]);
    lev(L).rb = min(p.Bins, max(1, floor((rv - q(1)) / max(eps, q(2) - q(1)) * p.Bins) + 1));
end
clear Rs

% register frames
P = NaN(T, 6);
nmi = NaN(T, 1);
fr = p.Frames(:)';
pp = zeros(1, 6);
done = false(1, T);
complete = true;
opt = optimset('Display', 'off');
for ii = 1:numel(fr)
    f = fr(ii);
    if ~isempty(p.Cancel) && p.Cancel()
        complete = false;
        break
    end
    report((ii - 1) / numel(fr), sprintf('Motion: registering frame %d of %d...', f, T));
    V = double(D.rawCanonVolume(f));
    for L = 1:numel(lev)
        Vs = smooth3d(V, lev(L).sig ./ vox);
        % intensity range of the moving frame inside the sample region
        mq = prctile1(Vs(lev(L).pick), [0.5 99.5]);
        s = lev(L).s;
        cost = @(z) -nmiCost(pp + s * (z(:)' - 1), Vs, lev(L).X, lev(L).rb, Cmm, vox, n, mq, p.Bins);
        o = optimset(opt, 'TolX', lev(L).tol, 'TolFun', 1e-6, 'MaxFunEvals', lev(L).maxf, ...
            'MaxIter', lev(L).maxf);
        [z, fv] = fminsearch(cost, ones(1, 6), o);
        pp = pp + s * (z(:)' - 1);
    end
    P(f, :) = pp;
    nmi(f) = -fv;
    done(f) = true;
end
report(1, 'Motion: done');

% temporal smoothing: running median per parameter (keeps steps)
Praw = P;
w = max(1, round(p.SmoothFrames));
if w > 1
    h = floor(w / 2);
    dn = find(done);
    for f = dn
        nb = dn(abs(dn - f) <= h);
        P(f, :) = median(Praw(nb, :), 1);
    end
end

% displacements
absM = NaN(T, 1); absX = NaN(T, 1); relM = NaN(T, 1); relX = NaN(T, 1);
Y = cell(1, T);
for f = find(done)
    Y{f} = transformPts(P(f, :), Xd, Cmm);
    d = sqrt(sum((Y{f} - Xd) .^ 2, 2));
    absM(f) = mean(d); absX(f) = max(d);
    if f > 1 && done(f - 1)
        d = sqrt(sum((Y{f} - Y{f - 1}) .^ 2, 2));
        relM(f) = mean(d); relX(f) = max(d);
    end
end
ex = find(absM > p.WarnAbsMM | relM > p.WarnRelMM)';
M = struct('Params', P, 'ParamsRaw', Praw, 'SmoothFrames', w, 'Region', lower(p.Region), ...
    'ExcludedEnhancing', nexcl, 'AbsMean', absM, 'AbsMax', absX, 'RelMean', relM, 'RelMax', relX, ...
    'Centre', Cmm, 'VoxelMM', vox, 'RefFrames', rf, 'Reference', reftxt, 'Frames', find(done), ...
    'NMI', nmi, 'WarnAbsMM', p.WarnAbsMM, 'WarnRelMM', p.WarnRelMM, 'ExcessFrames', ex, ...
    'Excess', ~isempty(ex), 'Complete', complete, 'Message', '');
M.Message = graspilot_motion_summary(M);
end

% -------------------------------------------------------------------------
function v = nmiCost(p, V, X, rb, c, vox, n, mq, nb)
Y = transformPts(p, X, c);
q = Y ./ repmat(vox, size(Y, 1), 1) + 1;           % voxel coordinates (1-based)
in = all(q >= 1 & q <= repmat(n, size(q, 1), 1), 2);
if nnz(in) < 0.5 * numel(in)
    v = 0;
    return
end
m = interpn(V, q(in, 1), q(in, 2), q(in, 3), 'linear');
mb = min(nb, max(1, floor((m - mq(1)) / max(eps, mq(2) - mq(1)) * nb) + 1));
J = accumarray([rb(in) mb], 1, [nb nb]);
J = J / sum(J(:));
pr = sum(J, 2); pm = sum(J, 1);
h = @(x) -sum(x(x > 0) .* log(x(x > 0)));
v = (h(pr) + h(pm)) / h(J);
end

function Y = transformPts(p, X, c)
Rm = graspilot_rigid_matrix(p);
C = repmat(c, size(X, 1), 1);
Y = (X - C) * Rm' + C + repmat(p(1:3), size(X, 1), 1);
end

function B = smooth3d(A, sig)
% separable Gaussian smoothing, sigma in voxels per axis
B = A;
for a = 1:3
    s = sig(a);
    if s < 0.3, continue; end
    r = ceil(2.5 * s);
    g = exp(-(-r:r) .^ 2 / (2 * s ^ 2));
    g = g / sum(g);
    sz = ones(1, 3); sz(a) = numel(g);
    B = convn(B, reshape(g, sz), 'same');
    % renormalize at the borders
    w = convn(ones(size(B, a), 1), g(:), 'same');
    sw = ones(1, 3); sw(a) = numel(w);
    B = B ./ repmat(reshape(w, sw), size(B) ./ sw);
end
end

function R = interior(head, vox, mm)
% head voxels at least mm (in-plane) from the outer head surface: distance to
% the first and last head voxel of the voxel's row and column (dark skull or
% CSF inside the head do not count as outside)
n = size(head); n(end + 1:3) = 1;
d = Inf(n, 'single');
for dim = 1:2
    m = any(head, dim);
    [~, f] = max(head, [], dim);
    [~, l] = max(flip(head, dim), [], dim);
    l = n(dim) - l + 1;
    sz = ones(1, 3); sz(dim) = n(dim);
    idx = reshape(single(1:n(dim)), sz);
    rep = n; rep(dim) = 1;
    dd = min(repmat(idx, rep) - repmat(single(f), sz), repmat(single(l), sz) - repmat(idx, rep)) * vox(dim);
    dd(~repmat(m, sz)) = -Inf;
    d = min(d, dd);
end
R = d >= mm;
end

function M = dilate(M, r)
for a = 1:3
    sz = ones(1, 3); sz(a) = 2 * r(a) + 1;
    M = convn(double(M), ones(sz), 'same') > 0;
end
end

function q = prctile1(v, pc)
v = sort(v(isfinite(v)));
if isempty(v), q = [0 1]; return; end
q = v(max(1, min(numel(v), round(pc / 100 * numel(v)))));
q = q(:)';
if q(2) <= q(1), q(2) = q(1) + 1; end
end

function th = otsu(vals)
vals = double(vals(isfinite(vals)));
lo = min(vals); hi = max(vals);
if hi <= lo, th = lo; return; end
nb = 256;
idx = min(nb, floor((vals - lo) / (hi - lo) * nb) + 1);
c = accumarray(idx(:), 1, [nb 1]);
pp = c / sum(c);
ctr = lo + ((1:nb)' - 0.5) * (hi - lo) / nb;
om = cumsum(pp);
mu = cumsum(pp .* ctr);
sb = (mu(end) * om - mu) .^ 2 ./ (om .* (1 - om));
sb(~isfinite(sb)) = 0;
[~, k] = max(sb);
th = lo + k * (hi - lo) / nb;
end
