function [lin, info] = gdce_detect_sss(D, varargin)
%GDCE_DETECT_SSS Find the posterior superior sagittal sinus just above the torcula.
%
%   [lin, info] = gdce_detect_sss(D)
%   [lin, info] = gdce_detect_sss(D, 'Name', Value, ...)
%
%   Returns canonical linear indices of an ROI in the thickest part of the
%   posterior superior sagittal sinus (SSS) above the torcula (confluence of
%   sinuses), excluding the transverse sinuses and the straight sinus.
%   Orientation must be correct (canonical axes R->L, A->P, I->S; see
%   GraspData.setAxisCodes): the method relies on "posterior", "midline"
%   and "superior".
%
%   Method, on the peak enhancement map (max over time - baseline):
%   1. head mask (Otsu on the baseline image) and its centre
%   2. venous structure: vessel mask = enhancement >= RelThr x the maximum
%      within +/-MidlineMM of the head centre in the posterior half of the
%      head; from up to Seeds strongly enhancing voxels there (>= 60% of
%      that maximum) grow 3D 6-connected regions and keep the largest (the
%      SSS-torcula-transverse sinus complex; arteries and scalp veins give
%      small regions)
%   3. per axial slice: full left-right width of that structure and the
%      2D cross-section of its most posterior part (the SSS against the
%      occipital inner table), no deeper front-to-back than it is wide so
%      that the straight sinus joining from anterior is excluded
%   4. torcula level = most superior slice where the structure is wider
%      than TorculaWidthMM (transverse sinuses leaving the midline)
%   5. candidate slices GapMM .. GapMM+SearchMM above the torcula where
%      the structure is no wider than MaxWidthMM (midline SSS only);
%      thickest = largest posterior cross-section (3-slice average)
%   6. ROI = SlabMM of slices centred on the thickest slice; in each slice
%      the cross-section voxels enhancing >= Core x their maximum
%      (excludes partial-volume edge voxels)
%
%   Options (defaults): 'BaseFrames' (1:min(5,T)), 'MidlineMM' (10),
%   'LateralMM' (50), 'RelThr' (0.35), 'PosteriorMM' (12),
%   'TorculaWidthMM' (30), 'MaxWidthMM' (20), 'GapMM' (5), 'SearchMM' (35),
%   'SlabMM' (5), 'Core' (0.5), 'Seeds' (6), 'Progress' ([])
%
%   info: Seed, TorculaSlice (NaN if not found), Slices, NVoxels, VolumeML,
%   Centroid [i j k], MeanPeakEnh, Width and Area per slice, Message.

p = struct('BaseFrames', [], 'MidlineMM', 10, 'LateralMM', 50, 'RelThr', 0.35, ...
    'PosteriorMM', 12, 'TorculaWidthMM', 30, 'MaxWidthMM', 20, 'GapMM', 5, 'SearchMM', 35, ...
    'SlabMM', 5, 'Core', 0.5, 'Seeds', 6, 'Progress', []);
for i = 1:2:numel(varargin)
    fn = fieldnames(p);
    m = find(strcmpi(fn, varargin{i}), 1);
    if isempty(m), error('gdce_detect_sss:opt', 'Unknown option %s.', varargin{i}); end
    p.(fn{m}) = varargin{i + 1};
end
lin = zeros(0, 1);
info = struct('Seed', [], 'TorculaSlice', NaN, 'Slices', [], 'NVoxels', 0, 'VolumeML', 0, ...
    'Centroid', [], 'MeanPeakEnh', NaN, 'Width', [], 'Area', [], 'Message', '');

n = D.canonSize();
vox = D.canonVox();
T = D.NT;
bf = p.BaseFrames;
if isempty(bf), bf = 1:max(1, min(5, T - 1)); end
if ~D.hasMaps(bf), D.computeMaps(bf, p.Progress); end
E = double(D.canonMap('PeakEnh'));
B = double(D.canonMap('Baseline'));
E(~isfinite(E)) = 0;

% 1. head mask and centre
step = max(1, floor(numel(B) / 300000));
head = B > otsu(B(1:step:end));
if nnz(head) < 50
    info.Message = 'Could not find the head in the baseline image.';
    return
end
[hx, hy, hz] = ind2sub(n, find(head)); %#ok<ASGLU>
xc = round(median(hx));
yc = round(median(hy));

% 2. venous structure: largest region grown from the brightest posterior midline voxels
xr = max(1, xc - round(p.MidlineMM / vox(1))):min(n(1), xc + round(p.MidlineMM / vox(1)));
yr = yc:n(2);
lx = max(1, xc - round(p.LateralMM / vox(1))):min(n(1), xc + round(p.LateralMM / vox(1)));
region = false(n);
region(lx, yr, :) = true;
search = false(n);
search(xr, yr, :) = true;
search = search & head;
best = []; bestN = 0; bestSeed = [];
Es = E;
Es(~search) = -Inf;
emax1 = max(Es(:));
vessel = (E >= p.RelThr * emax1) & region;     % one vessel threshold for all seeds
for s = 1:p.Seeds
    Es = E;
    Es(~search) = -Inf;
    [emax, k] = max(Es(:));
    if ~(emax > 0) || emax < 0.6 * emax1, break; end   % strongly enhancing seeds only
    [a, b, c] = ind2sub(n, k);
    C = GraspDCEViewer.floodFill3D(vessel, [a b c]);
    if nnz(C) > bestN
        best = C; bestN = nnz(C); bestSeed = [a b c];
    end
    search(C) = false;
end
if isempty(best)
    info.Message = 'No enhancing structure found near the posterior midline.';
    return
end
C = best;
info.Seed = bestSeed;

% 3. per-slice width and posterior cross-section
W = zeros(1, n(3)); A = zeros(1, n(3));
blob = cell(1, n(3));
pr = round(p.PosteriorMM / vox(2));
for k = 1:n(3)
    Ck = C(:, :, k);
    if ~any(Ck(:)), continue; end
    [ix, iy] = find(Ck);
    W(k) = (max(ix) - min(ix) + 1) * vox(1);
    yp = max(iy);
    cand = find(iy == yp);
    [~, j] = min(abs(ix(cand) - xc));
    seed2 = [ix(cand(j)) yp 1];
    M2 = Ck;
    M2(:, 1:max(0, yp - pr - 1)) = false;
    R2 = GraspDCEViewer.floodFill3D(M2, seed2);
    % depth (A-P) limited to the left-right width of the cross-section:
    % the SSS is about as deep as it is wide, the straight sinus joins
    % from anterior
    [bx, by] = find(R2);
    wd = max(bx) - min(bx) + 1;
    dlim = min(pr, round(wd * vox(1) / vox(2)));
    if any(by < yp - dlim + 1)
        M2(:, 1:max(0, yp - dlim)) = false;
        R2 = GraspDCEViewer.floodFill3D(M2, seed2);
    end
    blob{k} = find(R2);
    A(k) = numel(blob{k}) * vox(1) * vox(2);
end
info.Width = W;
info.Area = A;
present = find(A > 0);

% 4. torcula level
wide = find(W >= p.TorculaWidthMM);
if isempty(wide)
    k0 = min(present);
    k1 = k0 + floor(p.SearchMM / vox(3));
    tmsg = 'Torcula not identified (no transverse sinuses in the slab); used the lowest posterior SSS.';
else
    info.TorculaSlice = max(wide);
    k0 = info.TorculaSlice + ceil(p.GapMM / vox(3));
    k1 = k0 + floor(p.SearchMM / vox(3));
    tmsg = sprintf('Torcula at slice %d.', info.TorculaSlice);
end

% 5. thickest midline segment above the torcula
ok = present(present >= k0 & present <= k1 & W(present) <= p.MaxWidthMM);
if isempty(ok)
    info.Message = ['No midline sinus segment found above the torcula. Check the orientation ' ...
        '(View > Geometry) and that the slab covers the posterior SSS.'];
    return
end
As = zeros(1, n(3));
As(ok) = A(ok);
Asm = conv(As, ones(1, 3) / 3, 'same');
[~, i] = max(Asm(ok));
kb = ok(i);
h = max(0, round((p.SlabMM / vox(3) - 1) / 2));
ks = ok(ok >= kb - h & ok <= kb + h);

% 6. core voxels of the cross-section in those slices
np = n(1) * n(2);
for k = ks
    l2 = blob{k};
    e = E(l2 + (k - 1) * np);
    keep = e >= p.Core * max(e);
    M2 = false(n(1), n(2));
    M2(l2(keep)) = true;
    [~, im] = max(e);
    [sx, sy] = ind2sub([n(1) n(2)], l2(im));
    R2 = GraspDCEViewer.floodFill3D(M2, [sx sy 1]);
    lin = [lin; find(R2) + (k - 1) * np]; %#ok<AGROW>
end
[ii, jj, kk] = ind2sub(n, lin);
info.Slices = ks;
info.NVoxels = numel(lin);
info.VolumeML = numel(lin) * prod(vox) / 1000;
info.Centroid = [mean(ii) mean(jj) mean(kk)];
info.MeanPeakEnh = mean(E(lin));
info.Message = sprintf(['SSS ROI: %d voxels (%.3f mL) on slices %d-%d, thickest at slice %d ' ...
    '(%.0f mm2). %s'], numel(lin), info.VolumeML, min(ks), max(ks), kb, A(kb), tmsg);
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
