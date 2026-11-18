function [lin, info] = graspilot_detect_sss(D, varargin)
%GRASPILOT_DETECT_SSS Find the posterior superior sagittal sinus just above the torcula.
%
%   [lin, info] = graspilot_detect_sss(D)
%   [lin, info] = graspilot_detect_sss(D, 'Name', Value, ...)
%
%   Returns canonical linear indices of an ROI in the posterior superior
%   sagittal sinus (SSS): from its thickest part above the torcula
%   (confluence of sinuses) up along the sinus toward the vertex, excluding
%   the transverse sinuses and the straight sinus.
%   Orientation must be correct (canonical axes R->L, A->P, I->S; see
%   GRASPilotData.setAxisCodes): the method relies on "posterior", "midline"
%   and "superior".
%
%   Method, on the peak enhancement map (max after the baseline - baseline):
%   1. head mask (Otsu on the baseline image), its centre, and a band of
%      BandMM along the outer head surface (in-plane distance to the first
%      and last head voxel of each row and column, so the dark skull does
%      not count as outside). The sinuses lie in this band; the straight
%      sinus, deep arteries and veins and most tumours do not.
%   2. venous reference = maximum enhancement within +/-MidlineMM of the
%      head centre, in the posterior half and in the band; vessel mask =
%      enhancement >= RelThr x reference in the band (excludes the
%      moderately enhancing dura, falx and tentorium). From up to Seeds
%      strongly enhancing voxels there grow 3D 6-connected regions and keep
%      the largest (scalp vessels give small regions)
%   3. per axial slice: the most posterior vessel voxel near the midline;
%      the left-right width of its 2D-connected structure (transverse
%      sinuses at the torcula, cortical veins higher up), and the SSS
%      cross-section: the 2D region around it within MaxWidthMM/2 of it and
%      PosteriorMM of the posterior edge, no deeper front-to-back than it is
%      wide (excludes the straight sinus joining from anterior)
%   4. torcula level = top of the lowest run of slices whose posterior
%      structure is wider than TorculaWidthMM (transverse sinuses). Wide
%      slices higher up (bridging cortical veins joining the SSS) are
%      ignored
%   5. candidate slices GapMM .. GapMM+SearchMM above the torcula, preferring
%      those whose posterior structure is no wider than MaxWidthMM (no vein
%      joining); thickest = largest cross-section (3-slice average)
%   6. ROI slices: from InferiorMM below the thickest slice (not closer to
%      the torcula than GapMM) up to SuperiorMM above it, following the
%      sinus slice by slice. The extension stops early where the sinus
%      cross-section is missing, no longer overlaps the cross-section of the
%      slice below, or where the sinus turns forward over the vertex (its
%      posterior edge moving anteriorly by more than MaxSlope mm per mm of
%      height over >= 4 mm: the axial cut becomes oblique and partial-volumed).
%      All distances are in mm: a wrong voxel size shortens or clips the ROI
%   7. in each ROI slice the cross-section voxels enhancing >= Core x their
%      maximum (excludes partial-volume edge voxels), no wider than 1.5 x
%      the sinus at the thickest slice; voxels peaking more than 3 s (at
%      least 1.5 frames) before the ROI's median time to peak (arteries
%      touching the sinus) are dropped
%
%   Venous blood washes out after the first pass; enhancing scalp, dura,
%   falx and tumour keep rising or plateau. Only voxels whose late
%   enhancement (mean of the last 3 analysed frames) is below
%   1 - MinWashout of their peak enhancement count as vessels. If nothing is
%   found that way (e.g. an acquisition that ends before the sinus washes
%   out), the search is repeated without this test and the message says so.
%
%   'Seed' [i j k] (canonical voxel on the sinus, e.g. the crosshair) skips
%   the search: the region is grown from the most enhancing voxel within
%   3 voxels of it, the reference is its enhancement and the ROI starts at
%   its slice. Its curve must wash out after the peak (late enhancement below
%   1 - MinWashout (0.15) of the peak), as venous blood does and enhancing
%   tumour, dura and falx do not.
%
%   Options (defaults): 'BaseFrames' (automatic, GRASPilotData.detectBaseline), 'MidlineMM' (10),
%   'LateralMM' (50), 'BandMM' (25), 'RelThr' (0.5), 'PosteriorMM' (15),
%   'TorculaWidthMM' (30), 'MaxWidthMM' (20), 'GapMM' (5), 'SearchMM' (35),
%   'InferiorMM' (2), 'SuperiorMM' (50), 'MaxSlope' (1.5), 'Core' (0.5),
%   'Seeds' (6), 'Seed' ([]), 'MinWashout' (0.15), 'UseWashout' (true),
%   'Progress' ([]). 'SlabMM' (empty) gives the
%   former ROI: SlabMM of slices centred on the thickest slice (sets
%   InferiorMM and SuperiorMM to SlabMM/2).
%
%   info: Seed, TorculaSlice (NaN if not found), ThickestSlice, Slices,
%   ExtentMM, StopReason, NVoxels, VolumeML, Centroid [i j k], MeanPeakEnh,
%   Width (posterior structure) and Area (cross-section) per slice, Message.

[lin, info] = detectCore(D, varargin{:});
given = any(strcmpi(varargin(1:2:end), 'UseWashout'));
if isempty(lin) && ~given
    [lin2, info2] = detectCore(D, varargin{:}, 'UseWashout', false);
    if ~isempty(lin2)
        lin = lin2; info = info2;
        info.Message = [info.Message ' No washing-out venous structure was found, so the washout test was ' ...
            'not applied: check that the ROI is not in the scalp.'];
    end
end
end

function [lin, info] = detectCore(D, varargin)
p = struct('BaseFrames', [], 'MidlineMM', 10, 'LateralMM', 50, 'BandMM', 25, 'RelThr', 0.5, ...
    'PosteriorMM', 15, 'TorculaWidthMM', 30, 'MaxWidthMM', 20, 'GapMM', 5, 'SearchMM', 35, ...
    'InferiorMM', 2, 'SuperiorMM', 50, 'MaxSlope', 1.5, 'SlabMM', [], 'Core', 0.5, 'Seeds', 6, ...
    'Seed', [], 'MinWashout', 0.15, 'UseWashout', true, 'Progress', []);
for i = 1:2:numel(varargin)
    fn = fieldnames(p);
    m = find(strcmpi(fn, varargin{i}), 1);
    if isempty(m), error('graspilot_detect_sss:opt', 'Unknown option %s.', varargin{i}); end
    p.(fn{m}) = varargin{i + 1};
end
if ~isempty(p.SlabMM)
    p.InferiorMM = p.SlabMM / 2;
    p.SuperiorMM = p.SlabMM / 2;
end
lin = zeros(0, 1);
info = struct('Seed', [], 'TorculaSlice', NaN, 'ThickestSlice', NaN, 'Slices', [], 'ExtentMM', 0, ...
    'StopReason', '', 'NVoxels', 0, 'VolumeML', 0, 'Centroid', [], 'MeanPeakEnh', NaN, 'Width', [], ...
    'Area', [], 'Message', '');
hint = [' If the sinus is visible: put the crosshair on it above the torcula and use Edit > SSS ROI ' ...
    'from crosshair.'];

n = D.canonSize();
vox = D.canonVox();
T = D.NT;
bf = p.BaseFrames;
if isempty(bf), bf = D.detectBaseline(); end
if ~D.hasMaps(bf), D.computeMaps(bf, p.Progress); end
E = double(D.canonMap('PeakEnh'));
B = double(D.canonMap('Baseline'));
E(~isfinite(E)) = 0;
if p.UseWashout                                  % washing-out (venous) voxels only
    Pl = double(D.canonMap('Plateau'));
    Wm = Pl < (1 - p.MinWashout) * E;
else
    Wm = true(n);
end

% 1. head mask, centre and the band along the head surface
step = max(1, floor(numel(B) / 300000));
head = B > otsu(B(1:step:end));
if nnz(head) < 50
    info.Message = 'Could not find the head in the baseline image.';
    return
end
[hx, hy, hz] = ind2sub(n, find(head)); %#ok<ASGLU>
xc = round(median(hx));
yc = round(median(hy));
band = surfaceBand(head, vox, p.BandMM);
[XX, YY] = ndgrid(1:n(1), 1:n(2));
seedGiven = ~isempty(p.Seed);
if seedGiven
    s0 = round(p.Seed(:)');
    s0 = min(max(s0, 1), n);
    xc = s0(1);
end

% 2. vessel mask and the venous structure
xr = max(1, xc - round(p.MidlineMM / vox(1))):min(n(1), xc + round(p.MidlineMM / vox(1)));
lx = max(1, xc - round(p.LateralMM / vox(1))):min(n(1), xc + round(p.LateralMM / vox(1)));
region = false(n);
region(lx, :, :) = true;
if seedGiven
    lo = max(1, s0 - [3 3 1]); hi = min(n, s0 + [3 3 1]);
    sub = E(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));
    [emax1, k] = max(sub(:));
    [a, b, c] = ind2sub(size(sub), k);
    sd = [a b c] + lo - 1;
    if ~(emax1 > 0)
        info.Message = 'The crosshair is not on an enhancing vessel.';
        return
    end
    % sinus blood washes out after the first pass; tumour, dura and falx keep
    % rising or plateau
    Pl = double(D.canonMap('Plateau'));
    wf = 1 - Pl(sd(1), sd(2), sd(3)) / emax1;
    if ~(wf >= p.MinWashout)
        info.Message = sprintf(['The crosshair voxel does not wash out like venous blood (late signal %.0f%% ' ...
            'of its peak enhancement): put the crosshair on the sinus.'], 100 * (1 - wf));
        return
    end
    inBand = band(sd(1), sd(2), sd(3));
    vessel = (E >= p.RelThr * emax1) & region & Wm;
    if inBand, vessel = vessel & band; end
    best = GRASPilot.floodFill3D(vessel, sd);
    bestSeed = sd;
else
    search = false(n);
    search(xr, yc:n(2), :) = true;
    search = search & head & band & Wm;
    Es = E;
    Es(~search) = -Inf;
    emax1 = max(Es(:));
    vessel = (E >= p.RelThr * emax1) & region & band & Wm;
    best = []; bestN = 0; bestSeed = [];
    for s = 1:p.Seeds
        Es = E;
        Es(~search) = -Inf;
        [emax, k] = max(Es(:));
        if ~(emax > 0) || emax < 0.6 * emax1, break; end   % strongly enhancing seeds only
        [a, b, c] = ind2sub(n, k);
        C = GRASPilot.floodFill3D(vessel, [a b c]);
        if nnz(C) > bestN
            best = C; bestN = nnz(C); bestSeed = [a b c];
        end
        search(C) = false;
    end
end
if isempty(best) || ~any(best(:))
    info.Message = ['No enhancing structure found near the posterior midline.' hint];
    return
end
C = best;
info.Seed = bestSeed;

% 3. per slice: width of the posterior structure and the SSS cross-section
W = zeros(1, n(3)); A = zeros(1, n(3));
Wb = zeros(1, n(3)); Yp = NaN(1, n(3));      % cross-section width (mm), posterior edge (voxel)
blob = cell(1, n(3));
pr = round(p.PosteriorMM / vox(2));
hw = round(p.MaxWidthMM / 2 / vox(1));
strip = abs(XX - xc) <= round(p.MidlineMM / vox(1));
for k = 1:n(3)
    Ck = C(:, :, k);
    Cm = Ck & strip;
    if ~any(Cm(:)), continue; end
    [ix, iy] = find(Cm);
    yp = max(iy);
    cand = find(iy == yp);
    [~, j] = min(abs(ix(cand) - xc));
    seed2 = [ix(cand(j)) yp 1];
    Rt = GRASPilot.floodFill3D(Ck, seed2);
    [tx, ~] = find(Rt);
    W(k) = (max(tx) - min(tx) + 1) * vox(1);
    M2 = Ck & abs(XX - seed2(1)) <= hw;
    M2(:, 1:max(0, yp - pr - 1)) = false;
    R2 = GRASPilot.floodFill3D(M2, seed2);
    % depth (A-P) limited to the left-right width of the cross-section:
    % the SSS is about as deep as it is wide, the straight sinus joins
    % from anterior
    [bx, by] = find(R2);
    wd = max(bx) - min(bx) + 1;
    dlim = min(pr, round(wd * vox(1) / vox(2)));
    if any(by < yp - dlim + 1)
        M2(:, 1:max(0, yp - dlim)) = false;
        R2 = GRASPilot.floodFill3D(M2, seed2);
    end
    blob{k} = find(R2);
    A(k) = numel(blob{k}) * vox(1) * vox(2);
    Wb(k) = wd * vox(1);
    Yp(k) = yp;
end
info.Width = W;
info.Area = A;
present = find(A > 0);
if isempty(present)
    info.Message = ['No sinus cross-section found near the posterior midline.' hint];
    return
end

% 4. torcula level: top of the lowest run of wide slices
wide = W >= p.TorculaWidthMM;
k1w = find(wide, 1);
if seedGiven && ~isempty(k1w) && k1w > s0(3), k1w = []; end
if isempty(k1w)
    k0 = min(present);
    tmsg = 'Torcula not identified (no transverse sinuses in the slab); used the lowest posterior SSS.';
else
    kt = k1w;
    while kt < n(3) && (wide(kt + 1) || (kt + 2 <= n(3) && A(kt + 1) == 0 && wide(kt + 2)))
        kt = kt + 1;
    end
    info.TorculaSlice = kt;
    k0 = kt + ceil(p.GapMM / vox(3));
    tmsg = sprintf('Torcula at slice %d.', info.TorculaSlice);
end
k1 = k0 + floor(p.SearchMM / vox(3));

% 5. thickest midline segment above the torcula (or the crosshair slice)
if seedGiven
    [~, i] = min(abs(present - s0(3)));
    kb = present(i);
    ok = kb;
    k0 = min(k0, kb);
    tmsg = sprintf('%s Started at the crosshair slice %d.', tmsg, kb);
else
    inr = present(present >= k0 & present <= k1);
    ok = inr(W(inr) <= p.MaxWidthMM);
    if isempty(ok), ok = inr; end              % veins join on every slice: use them all
    if isempty(ok)
        info.Message = ['No midline sinus segment found above the torcula. Check the orientation ' ...
            '(View > Geometry) and that the slab covers the posterior SSS.' hint];
        return
    end
    As = zeros(1, n(3));
    As(ok) = A(ok);
    Asm = conv(As, ones(1, 3) / 3, 'same');
    [~, i] = max(Asm(ok));
    kb = ok(i);
end

info.ThickestSlice = kb;
kl = max(k0, kb - round(p.InferiorMM / vox(3)));
ks = kl:kb;
ks = ks(ismember(ks, present));

% 6. follow the sinus superiorly toward the vertex
np = n(1) * n(2);
kw = kb + floor(p.SuperiorMM / vox(3) + 1e-9);   % wanted top slice
kt = min(n(3), kw);
span = max(3, ceil(4 / vox(3)));                 % slices over which the forward turn is measured
why = '';
for k = kb + 1:kt
    if A(k) == 0
        why = 'sinus not found higher up'; break
    end
    if ~overlaps(blob{k}, blob{k - 1}, [n(1) n(2)])
        why = 'cross-section not continuous with the slice below'; break
    end
    j = k - span;
    if j >= ks(1) && (Yp(j) - Yp(k)) * vox(2) / (span * vox(3)) > p.MaxSlope
        why = 'sinus turning forward over the vertex'; break
    end
    ks(end + 1) = k; %#ok<AGROW>
end
if isempty(why) && kw > kb
    if kt < kw, why = 'top of the imaged volume';
    else, why = sprintf('%g mm above the thickest slice (SuperiorMM)', p.SuperiorMM); end
end
info.StopReason = why;

% 7. core voxels of the cross-section in those slices, no wider than 1.5 x
%    the sinus at the thickest slice (veins joining are cut off)
hwk = Inf;
for pass = 1:2
    lin = zeros(0, 1);
    for k = ks
        l2 = blob{k};
        e = E(l2 + (k - 1) * np);
        keep = e >= p.Core * max(e);
        [bx, ~] = ind2sub([n(1) n(2)], l2);
        [~, im] = max(e);
        sx = bx(im);
        xm = round(median(bx(keep)));             % centre of the sinus core
        keep = keep & abs(bx - xm) <= hwk;
        M2 = false(n(1), n(2));
        M2(l2(keep)) = true;
        if ~keep(im)                              % brightest voxel trimmed: start from the brightest kept one
            ek = e; ek(~keep) = -Inf;
            [~, im] = max(ek);
            sx = bx(im);
        end
        [~, sy] = ind2sub([n(1) n(2)], l2(im));
        R2 = GRASPilot.floodFill3D(M2, [sx sy 1]);
        lin = [lin; find(R2) + (k - 1) * np]; %#ok<AGROW>
        if pass == 1 && k == kb
            [rx, ~] = find(R2);
            hwk = max(1.5 * (max(rx) - min(rx) + 1) / 2, 4 / vox(1));
        end
    end
    if ~isfinite(hwk), break; end
end
% arterial voxels touching the sinus peak earlier than venous blood: drop
% voxels peaking more than 1.5 frames before the ROI's median time to peak
if numel(lin) > 10
    TTP = double(D.canonMap('TTP'));
    tt = D.paramTimes();
    lag = 1.5 * median(diff(tt));
    if D.hasTiming(), lag = max(lag, 3); end    % arteries peak several seconds before the sinus
    tp = TTP(lin);
    lin = lin(~(tp < median(tp) - lag));
end
[ii, jj, kk] = ind2sub(n, lin);
info.Slices = ks;
info.ExtentMM = numel(ks) * vox(3);
info.NVoxels = numel(lin);
info.VolumeML = numel(lin) * prod(vox) / 1000;
info.Centroid = [mean(ii) mean(jj) mean(kk)];
info.MeanPeakEnh = mean(E(lin));
info.Message = sprintf(['SSS ROI: %d voxels (%.3f mL) on slices %d-%d (%.0f mm), from the thickest ' ...
    'slice %d (%.0f mm2) up toward the vertex. %s'], numel(lin), info.VolumeML, min(ks), max(ks), ...
    info.ExtentMM, kb, A(kb), tmsg);
if ~isempty(why), info.Message = sprintf('%s Upper end: %s.', info.Message, why); end
end

function band = surfaceBand(head, vox, mm)
% voxels within mm of the outer head surface, in-plane: distance to the first
% and last head voxel of the voxel's row and column (holes such as the dark
% skull do not count as outside)
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
    F = repmat(single(f), sz); L = repmat(single(l), sz); I = repmat(idx, rep);
    dd = min(I - F, L - I) * vox(dim);
    dd(~repmat(m, sz) | dd < 0) = Inf;
    d = min(d, dd);
end
band = d <= mm;
end

function tf = overlaps(a, b, n2)
% 2D cross-sections (linear indices in a slice) touch or overlap
if isempty(a) || isempty(b), tf = false; return; end
M = false(n2);
M(b) = true;
M = M | circshift(M, [1 0]) | circshift(M, [-1 0]) | circshift(M, [0 1]) | circshift(M, [0 -1]);
tf = any(M(a));
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
