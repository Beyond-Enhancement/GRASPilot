function nfail = gdce_selftest_gui_sss(out)
%GDCE_SELFTEST_GUI_SSS Auto SSS button and orientation check on the venous phantom.
nfail = 0;
[D, tr] = gdce_sss_phantom();
app = GraspDCEViewer(D);
app.Quiet = true;
cleanup = onCleanup(@() app.close(true));
nfail = nfail + check('orientation check passes (LPS phantom)', app.checkOrientation(), 0);
app.addTumor();
k = app.autoSSS();
lin = find(app.Masks{k});
[ref, info] = gdce_detect_sss(D, 'BaseFrames', app.BaseFrames(1):app.BaseFrames(2));
nfail = nfail + check('Auto SSS button creates "SSS (auto)" ROI', k > 0 && strcmp(app.ROIs(k).Name, 'SSS (auto)') ...
    && isequal(sort(lin), sort(ref)) && all(tr.SSS(lin)), numel(lin));
nfail = nfail + check('crosshair moved to the SSS ROI', isequal(app.Pos, round(info.Centroid)), 0);
k2 = app.autoSSS();
nfail = nfail + check('second run replaces, does not duplicate', k2 == k && numel(app.ROIs) == 2 && ...
    nnz(app.Masks{k}) == numel(ref), numel(app.ROIs));
app.undo(); app.undo();
nfail = nfail + check('Auto SSS undoable', ~any(app.Masks{k}(:)), 0);
app.computeStats(k);
app.redo();
app.computeStats(k);
st = app.ROIs(k).Stats;
[~, pk] = max(st.mean);
nfail = nfail + check('SSS curve peaks in the venous phase (frame 6)', pk == 6, pk);
app.applyGeometry('PLS', [1 1 1], 4, 0);
nfail = nfail + check('rotated view flagged by orientation check', ~app.checkOrientation(), 0);
nfail = nfail + check('ROI follows axis-code change', isequal(D.canon2native(app.Masks{k}) , tr.SSS & ...
    D.canon2native(app.Masks{k})), 0);
app.applyGeometry('LPS', [1 1 1], 4, 0);
nfail = nfail + check('ROI back in place after restoring LPS', isequal(sort(find(app.Masks{k})), sort(ref)), 0);

% ---- semi-quantitative parameters, normalization to the SSS, statistics, export
n = D.canonSize();
t = 1;                                               % the tumor ROI added above
M = false(n); M(40:46, 40:46, 50:54) = true;         % brain tissue box
app.beginAction(); app.setVoxels(t, find(M), true); app.commitAction('box'); app.afterMaskChange();
ok = app.computeParams();
lin = find(app.Masks{t});
bf = app.BaseFrames(1):app.BaseFrames(2);
P = gdce_curve_params(double(D.timecourses(lin)), D.paramTimes(), bf, D.ParamOpts);
app.computeStats(k);
R = gdce_curve_params(app.ROIs(k).Stats.mean, D.paramTimes(), bf, D.ParamOpts);
e = 0;
pairs = {'PE', 'PE', 1; 'RE', 'RE', 1; 'WIS', 'WashIn', 1; 'WOS', 'WashOut', 1; 'AUC', 'AUC', 1; ...
    'TTP', 'TTP', 1; 'PE_SSS', 'PE', R.PE; 'WIS_SSS', 'WashIn', R.WashIn; 'WOS_SSS', 'WashOut', R.WashOut; ...
    'AUC_SSS', 'AUC', R.AUC};
for q = 1:size(pairs, 1)
    v = app.paramValues(pairs{q, 1}, lin);
    w = P.(pairs{q, 2}) / pairs{q, 3};
    e = max(e, max(abs(v - w)) / max(abs(w)));
end
nfail = nfail + check('parameter values (absolute and /SSS)', ok && e < 1e-4, e);
nfail = nfail + check('SSS reference: PE > 0, washout < 0', R.PE > 0.5 && R.WashOut < 0, R.PE);
[hdr, rows] = app.paramStatsTable();
ir = find(strcmp(rows(:, 2), 'Peak enhancement / SSS'));
d = GraspDCEViewer.describeVals(P.PE / R.PE);
app.computeStats(t);
cm = gdce_curve_params(app.ROIs(t).Stats.mean, D.paramTimes(), bf, D.ParamOpts);
okT = size(rows, 1) == numel(app.ParamKeys) && all(strcmp(rows(:, 1), app.ROIs(t).Name)) && numel(ir) == 1 && ...
    abs(rows{ir, 5} - mean(P.PE / R.PE)) < 1e-5 && abs(rows{ir, 7} - median(P.PE / R.PE)) < 1e-5 && ...
    rows{ir, 4} == numel(lin) && abs(rows{ir, 6} - std(P.PE / R.PE)) < 1e-5 && ...
    abs(rows{ir, end} - cm.PE / R.PE) < 1e-6 && abs(d.p90 - rows{ir, 11}) < 1e-6 && strcmp(hdr{end}, 'ROI mean curve');
nfail = nfail + check('tumor statistics table (SSS excluded)', okT, size(rows, 1));
app.IncludeRefStats = true;
[~, rows2] = app.paramStatsTable();
app.IncludeRefStats = false;
nfail = nfail + check('statistics include the SSS when asked', size(rows2, 1) == 2 * numel(app.ParamKeys), size(rows2, 1));

% colour map overlay and colour bar
app.setPos(round([43 43 52]));
app.setParamKey('PE_SSS');
app.setParamDisplay(true);
app.setParamCLim([0 0.2]);
okM = app.ParamOn && ~isempty(app.ParamA{1}) && any(app.ParamA{1}(:) > 0) && ...
    strcmp(get(app.CB.Img, 'Visible'), 'on') && isequal(app.ParamCLims.PE_SSS, [0 0.2]);
cd1 = get(app.V(1).Img, 'CData');
sp = app.canon2screen(1, app.Pos);
px = squeeze(cd1(round(sp(2)), round(sp(1)), :));
nfail = nfail + check('parameter map overlay + colour bar', okM && std(px) > 0.05, std(px));
c = app.autoParamCLim('rois');
vv = sort(P.PE / R.PE);
nfail = nfail + check('auto colour range from tumor ROIs (2-98%)', abs(c(1) - vv(max(1, round(0.02 * numel(vv))))) < 1e-5 ...
    && abs(c(2) - vv(round(0.98 * numel(vv)))) < 1e-5, c(2));
app.onKey(struct('Key', 'm', 'Modifier', {{}}, 'Character', 'm'));
nfail = nfail + check('key M toggles the map', ~app.ParamOn && isempty(app.ParamA{1}), 0);
app.onKey(struct('Key', 'm', 'Modifier', {{}}, 'Character', 'm'));
app.openParamWindow();
tb = get(app.PW.table, 'Data');
okW = ishghandle(app.PW.Fig) && size(tb, 1) == numel(app.ParamKeys) && ...
    strcmp(get(app.PW.cmin, 'String'), sprintf('%.4g', c(1)));
set(app.PW.cmap, 'Value', 2); app.onParamCmap(2);
set(app.PW.washIn, 'Value', 2); app.onParamOptions(); app.onParamCompute();
okW = okW && strcmp(app.ParamCmap, 'turbo') && strcmp(D.ParamOpts.WashIn, 'max') && D.hasMaps(bf);
nfail = nfail + check('parameter window (table, colour map, options)', okW, size(tb, 1));
set(app.PW.washIn, 'Value', 1); app.onParamOptions(); app.onParamCompute();
if nargin > 0 && ~isempty(out)
    shot = fullfile(out, 'param_window.png');
    try, print(app.PW.Fig, '-dpng', '-r60', shot); catch, end
end

% export
fcsv = fullfile(out, 'params.csv');
files = app.exportParams(fcsv);
fid = fopen(fullfile(out, 'params_param_voxels.csv')); h1 = fgetl(fid); l1 = fgetl(fid); fclose(fid);
cols = regexp(l1, ',', 'split');
ip = find(strcmp(regexp(h1, ',', 'split'), '"PE_SSS [ratio]"'));
okE = numel(files) == 3 && all(cellfun(@(f) exist(f, 'file') == 2, files)) && numel(ip) == 1 && ...
    abs(str2double(cols{ip}) - P.PE(1) / R.PE) < 1e-5;
fid = fopen(fullfile(out, 'params_param_stats.csv')); nl = 0;
while ischar(fgetl(fid)), nl = nl + 1; end
fclose(fid);
nfail = nfail + check('parameter export (stats + voxel values csv)', okE && nl == 1 + numel(app.ParamKeys), nl);
if exist('writecell', 'file') == 2
    fx = app.exportParams(fullfile(out, 'params.xlsx'));
    nfail = nfail + check('parameter export (xlsx)', exist(fx{1}, 'file') == 2, 0);
end
app.exportMaps(out, 'pmaps');
[V, hv] = gdce_nifti_read(fullfile(out, 'pmaps_WOS_SSS.nii.gz'));
Wm = double(D.canonMap('WashOut')) / R.WashOut;
e = max(abs(double(V(lin)) - Wm(lin)));
nfail = nfail + check('parameter maps NIfTI (WOS / SSS)', e < 1e-4 && isequal(hv.Dim(1:3), n), e);

% normalization needs the reference ROI
app.ROIs(k).Name = 'vein';
app.afterMaskChange();
v = app.paramValues('PE_SSS', lin);
nfail = nfail + check('no SSS ROI: normalized values unavailable', all(isnan(v)) && isempty(app.paramVolume('PE_SSS')), 0);
app.RefROIName = 'vein';
app.afterMaskChange();
v = app.paramValues('PE_SSS', lin);
nfail = nfail + check('reference ROI chosen by name', max(abs(v - P.PE / R.PE)) < 1e-4, 0);
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
