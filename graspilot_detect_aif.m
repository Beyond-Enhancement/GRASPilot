function [lin, info] = graspilot_detect_aif(D, varargin)
%GRASPILOT_DETECT_AIF Automatic arterial input function: arterial voxels of a DCE series.
%
%   [lin, info] = graspilot_detect_aif(D)
%   [lin, info] = graspilot_detect_aif(D, 'Name', Value, ...)
%
%   Canonical linear indices of NVoxels voxels with an arterial curve: among
%   the voxels inside the head (at least MarginMM from its outer surface,
%   in-plane, so scalp arteries are left out, and away from the first and
%   last Slices of the slab) that wash out after the peak (late enhancement
%   < (1 - MinWashout) x peak enhancement), the TopPct percent with the
%   highest peak enhancement are candidates; of those, the ones peaking
%   earliest (time to peak within the first EarlyPct percent of the
%   candidates) and, among them, the NVoxels with the highest peak
%   enhancement are kept. Arteries peak before veins and tissue.
%
%   Small arteries are partial-volumed: the AIF amplitude is then too low
%   (Ktrans, vp too high). Check the voxels and the curve; the VCAI corrects
%   the amplitude with a neighbour ring and calibrates the shape against the
%   superior sagittal sinus.
%
%   Options (defaults): 'BaseFrames' (automatic), 'NVoxels' (20),
%   'MarginMM' (12), 'Slices' (2), 'MinWashout' (0.3), 'TopPct' (5),
%   'EarlyPct' (20)
%
%   info: NVoxels, MeanPeakEnh, TTP (median, s or frames), Centroid, Message

p = struct('BaseFrames', [], 'NVoxels', 20, 'MarginMM', 12, 'Slices', 2, 'MinWashout', 0.3, ...
    'TopPct', 5, 'EarlyPct', 20);
for i = 1:2:numel(varargin)
    fn = fieldnames(p);
    m = find(strcmpi(fn, varargin{i}), 1);
    if isempty(m), error('graspilot_detect_aif:opt', 'Unknown option %s.', varargin{i}); end
    p.(fn{m}) = varargin{i + 1};
end
lin = zeros(0, 1);
info = struct('NVoxels', 0, 'MeanPeakEnh', NaN, 'TTP', NaN, 'Centroid', [], 'Message', '');
bf = p.BaseFrames;
if isempty(bf), bf = D.detectBaseline(); end
if ~D.hasMaps(bf), D.computeMaps(bf); end
E = double(D.canonMap('PeakEnh'));
B = double(D.canonMap('Baseline'));
P = double(D.canonMap('Plateau'));
TTP = double(D.canonMap('TTP'));
n = D.canonSize(); n(end + 1:3) = 1;
vox = D.canonVox();

head = B > GRASPilotData.otsu(B(1:max(1, floor(numel(B) / 300000)):end));
inner = interior(head, vox, p.MarginMM);
s = min(p.Slices, floor((n(3) - 1) / 2));
if s > 0, inner(:, :, [1:s, n(3) - s + 1:n(3)]) = false; end
ok = inner & isfinite(E) & E > 0 & P < (1 - p.MinWashout) * E & isfinite(TTP);
if nnz(ok) < p.NVoxels
    info.Message = ['No arterial voxels found (strongly enhancing, washing out, inside the head). ' ...
        'Draw an arterial ROI (+ Artery) instead.'];
    return
end
e = sort(E(ok), 'descend');
thr = e(max(1, min(numel(e), round(p.TopPct / 100 * nnz(ok)))));
cand = find(ok & E >= thr);
if numel(cand) < p.NVoxels
    [~, o] = sort(E(ok), 'descend');
    okl = find(ok);
    cand = okl(o(1:p.NVoxels));
end
tc = sort(TTP(cand));
tEarly = tc(max(1, round(p.EarlyPct / 100 * numel(tc))));
early = cand(TTP(cand) <= tEarly);
if numel(early) < p.NVoxels                      % not enough at the earliest times: next earliest
    [~, o] = sort(TTP(cand) - 1e-6 * E(cand));
    early = cand(o(1:min(numel(cand), p.NVoxels)));
end
[~, o] = sort(E(early), 'descend');
lin = early(o(1:min(numel(early), p.NVoxels)));
[ii, jj, kk] = ind2sub(n, lin);
info.NVoxels = numel(lin);
info.MeanPeakEnh = mean(E(lin));
info.TTP = median(TTP(lin));
info.Centroid = [mean(ii) mean(jj) mean(kk)];
unit = 's'; if ~D.hasTiming(), unit = 'frames'; end
info.Message = sprintf(['Auto AIF: %d arterial voxels (earliest peak among the most enhancing voxels inside ' ...
    'the head; time to peak %.3g %s). Check them: small arteries are partial-volumed (AIF amplitude ' ...
    'too low).'], numel(lin), info.TTP, unit);
end

function R = interior(head, vox, mm)
% head voxels at least mm (in-plane) from the outer head surface, measured to
% the first and last head voxel of the voxel's row and column
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
R = head & d >= mm;
end
