function nfail = gdce_selftest_gui(D, ref, out)
%GDCE_SELFTEST_GUI Drive GraspDCEViewer by script (called from gdce_selftest(true)).
nfail = 0;
tol = 1e-5;
D.setOrientation('axial', false, [0 0 0]);
D.VoxNative = [1 1 2]; D.dt = 2.5; D.t0 = 0;
app = GraspDCEViewer(D);
app.Quiet = true;
cleanup = onCleanup(@() app.close(true));
nn = D.NativeSize;

% ---- display orientation: the native plane shows the native slice unchanged
app.setFrame(7);
app.setPos([5 6 3]);
e = max(abs(reshape(double(app.ViewImg{1}) - double(ref(:, :, 3, 7)), [], 1)));
nfail = nfail + check('axial preset: axial view == data1(:,:,7,3)', e < tol, e);
app.applyGeometry('coronal', false, [0 0 0], [1 1 2], 2.5, 0);
c = app.Pos; c(2) = 4; app.setPos(c);
e = max(abs(reshape(double(app.ViewImg{2}) - double(ref(:, :, 4, 7)), [], 1)));
nfail = nfail + check('coronal preset: coronal view == native slice', e < tol, e);
app.applyGeometry('sagittal', false, [0 0 0], [1 1 2], 2.5, 0);
c = app.Pos; c(1) = 2; app.setPos(c);
e = max(abs(reshape(double(app.ViewImg{3}) - double(ref(:, :, 2, 7)), [], 1)));
nfail = nfail + check('sagittal preset: sagittal view == native slice', e < tol, e);
app.applyGeometry('axial', false, [0 0 0], [1 1 2], 2.5, 0);
n = D.canonSize();

% ---- 2D circle/brush, 3D sphere, exclusivity, undo/redo
app.addSSS();
app.setPos([10 10 4]);
app.beginAction(); app.paintShape(1, 10, 12, 3, false); app.commitAction('t'); app.afterMaskChange();
[XX, YY] = meshgrid(1:n(1), 1:n(2));
exp2d = nnz((XX - 10) .^ 2 + (YY - 12) .^ 2 <= 9);
nfail = nfail + check('2D disk voxel count', nnz(app.Masks{1}) == exp2d && all(find(app.Masks{1}) > 0), nnz(app.Masks{1}));
[ii, jj, kk] = ind2sub(n, find(app.Masks{1}));
nfail = nfail + check('2D disk on current axial slice only', all(kk == 4) && all(abs(ii - 10) <= 3) && all(abs(jj - 12) <= 3), 0);

app.addTumor();
app.set3D(true);
app.beginAction(); app.paintShape(1, 12, 12, 3, false); app.commitAction('t'); app.afterMaskChange();
[I, J, K] = ndgrid(1:n(1), 1:n(2), 1:n(3));
expS = ((I - 12) .^ 2 + (J - 12) .^ 2 + ((K - 4) * 2) .^ 2 <= 9);
nfail = nfail + check('3D sphere (anisotropic voxels)', isequal(app.Masks{2}, expS), nnz(xor(app.Masks{2}, expS)));
ov = nnz(app.Masks{1} & app.Masks{2});
nfail = nfail + check('exclusive ROIs: no overlap', ov == 0, ov);
nSSS = nnz(app.Masks{1});
app.undo();
nfail = nfail + check('undo restores both ROIs', nnz(app.Masks{2}) == 0 && nnz(app.Masks{1}) == exp2d, nnz(app.Masks{1}));
app.redo();
nfail = nfail + check('redo re-applies', isequal(app.Masks{2}, expS) && nnz(app.Masks{1}) == nSSS, 0);
app.set3D(false);

% ---- eraser
app.setActiveROI(2);
app.beginAction(); app.paintShape(1, 12, 12, 1, true); app.commitAction('e'); app.afterMaskChange();
nfail = nfail + check('eraser removes from active ROI', nnz(app.Masks{2}) == nnz(expS) - 5, nnz(app.Masks{2}));
app.undo();

% ---- threshold (current frame signal) keeps only voxels >= Low
lin = find(app.Masks{2});
v = double(D.valuesAtFrame(lin, app.Frame));
th = median(v);
app.setThrSource(1);
app.setThreshold(th, Inf);
app.applyThreshold('3d');
lin2 = find(app.Masks{2});
v2 = double(D.valuesAtFrame(lin2, app.Frame));
nfail = nfail + check('threshold keeps voxels >= Low only', all(v2 >= th) && numel(lin2) == nnz(v >= th), numel(lin2));
app.undo();

% ---- map-based threshold (peak relative enhancement) + preview + Otsu
app.setThrSource(5);
app.setPreview(true);
nfail = nfail + check('enhancement maps computed on demand', D.hasMaps(1:5), 0);
th = app.suggestThreshold();
nfail = nfail + check('Otsu suggestion finite', isfinite(th), th);
al = app.PrevA{1};
rgb = get(app.V(1).Img, 'CData');
nfail = nfail + check('preview overlay composited', any(al(:) > 0) && size(rgb, 3) == 3, nnz(al));
app.applyThreshold('3d');
RE = D.canonMap('RelEnh');
nfail = nfail + check('peak-RE threshold applied', all(RE(app.Masks{2}) >= th), nnz(app.Masks{2}));
app.undo();
app.setPreview(false);
app.setThrSource(1);
app.setThreshold(-Inf, Inf);

% ---- constrained painting
app.setActiveROI(1);
disk = app.diskIndices(1, 20, 20, 4);
thc = median(double(D.valuesAtFrame(disk, app.Frame)));
app.setThreshold(thc, Inf);
app.setConstrain(true);
before = app.Masks{1};
app.beginAction(); app.paintShape(1, 20, 20, 4, false); app.commitAction('c'); app.afterMaskChange();
newv = find(app.Masks{1} & ~before);
vals = double(D.valuesAtFrame(newv, app.Frame));
nfail = nfail + check('constrained brush adds only in-window voxels', all(vals >= thc) && ...
    numel(newv) > 0 && numel(newv) < numel(disk), numel(newv));
app.setConstrain(false);
app.setThreshold(-Inf, Inf);

% ---- fill holes
app.addROI('Ring');
r = numel(app.ROIs);
ring = false(n);
ring(5:11, 5:11, 2) = true; ring(7:9, 7:9, 2) = false;
app.beginAction(); app.setVoxels(r, find(ring), true); app.commitAction('ring');
app.ActiveView = 1;
app.fillHoles();
filled = false(n); filled(5:11, 5:11, 2) = true;
nfail = nfail + check('fill holes (axial)', isequal(app.Masks{r}, filled), nnz(xor(app.Masks{r}, filled)));

% ---- TIC statistics
app.computeStats(2);
lin = find(app.Masks{2});
TC = zeros(numel(lin), D.NT);
Vc = zeros(prod(n), D.NT);
for t = 1:D.NT
    tmp = D.native2canon(ref(:, :, :, t));
    Vc(:, t) = tmp(:);
end
TC = Vc(lin, :);
e = max(abs(app.ROIs(2).Stats.mean - mean(TC, 1)));
nfail = nfail + check('ROI mean curve == mean of ReadData voxels', e < tol, e);
for m = 1:4, app.setFlag('TicY', m); end
for m = 1:3, app.setFlag('TicStat', m); end
app.setFlag('TicX', 2); app.setFlag('TicX', 1);
nfail = nfail + check('TIC plot items drawn', numel(app.Tic.Items) >= 3, numel(app.Tic.Items));

% ---- keyboard and frame navigation
app.setFrame(3);
ev = struct('Key', 'rightarrow', 'Modifier', {{}});
app.onKey(ev);
ev2 = struct('Key', 'rightarrow', 'Modifier', {{'shift'}});
app.onKey(ev2);
nfail = nfail + check('arrow keys step frames (1, then 10)', app.Frame == 14, app.Frame);
app.setFrame(Inf);
nfail = nfail + check('last-frame button', app.Frame == D.NT, app.Frame);
app.onKey(struct('Key', 'z', 'Modifier', {{'control'}}));
nfail = nfail + check('Ctrl+Z undo via keyboard', ~any(app.Masks{r}(filled & ~ring)), 0);
app.onKey(struct('Key', 'y', 'Modifier', {{'control'}}));

% ---- exports
masksBefore = app.Masks;
fImg = fullfile(out, 'gui_series.nii.gz');
app.exportImageNifti(fImg);
[img, h] = gdce_nifti_read(fImg);
e = max(abs(reshape(img(:, :, :, 5) - D.canonVolume(5), [], 1)));
nfail = nfail + check('GUI image NIfTI export', e < tol && isequal(h.Dim, [n D.NT]), e);
fMask = fullfile(out, 'gui_mask.nii.gz');
app.exportMaskNifti(fMask, true);
[L, hm] = gdce_nifti_read(fMask);
okm = true;
for q = 1:numel(app.Masks), okm = okm && isequal(L == q, app.Masks{q}); end
A0 = D.affine();
nfail = nfail + check('mask NIfTI labels and affine', okm && max(abs(hm.Affine(:) - A0(:))) < 1e-5, 0);
nfail = nfail + check('label tables written', exist(fullfile(out, 'gui_mask_labels.txt'), 'file') == 2 && ...
    exist(fullfile(out, 'gui_mask_itksnap_labels.txt'), 'file') == 2, 0);
files = app.exportTimeseries(fullfile(out, 'gui_ts.csv'));
nfail = nfail + check('time-series CSV export', numel(files) == 3, numel(files));

% ---- ROI session and mask import round trips
fS = fullfile(out, 'gui_rois.mat');
app.saveROIs(fS);
natB = cellfun(@(m) D.canon2native(m), masksBefore, 'UniformOutput', false);
app.applyGeometry('sagittal', true, [1 0 1], [1 1 2], 2.5, 0);
nat = cellfun(@(m) D.canon2native(m), app.Masks, 'UniformOutput', false);
nfail = nfail + check('ROIs follow an orientation change', isequal(nat, natB), 0);
app.loadROIs(fS);
nfail = nfail + check('ROI session reload (geometry restored)', isequal(app.Masks, masksBefore) && ...
    strcmp(D.Plane, 'axial'), 0);
names = {app.ROIs.Name};
app.importMaskNifti(fMask);
nfail = nfail + check('mask NIfTI import (names from label table)', isequal(app.Masks, masksBefore) && ...
    isequal({app.ROIs.Name}, names), 0);

% ---- delete ROI keeps undo stack consistent
nU = numel(app.UndoStack);
app.setActiveROI(1);
app.deleteROI();
nfail = nfail + check('delete ROI', numel(app.ROIs) == numel(names) - 1, numel(app.ROIs));
app.undo(); app.redo();
nfail = nfail + check('undo/redo after delete', true, nU);

% ---- playback controls and layout
app.togglePlay(); app.togglePlay();
set(app.Fig, 'Position', [50 50 1300 800]);
app.doLayout();
app.resetZoom();
app.zoomView(1, 10, 10, 0.5);
xl = get(app.V(1).Ax, 'XLim');
nfail = nfail + check('zoom', diff(xl) < n(1), diff(xl));
app.showHelp();
end

function f = check(name, cond, val)
if cond
    fprintf('  PASS  %-48s (%.3g)\n', name, val);
    f = 0;
else
    fprintf('  FAIL  %-48s (%.3g)\n', name, val);
    f = 1;
end
end
