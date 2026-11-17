function nfail = graspilot_selftest_pk(out, doGUI)
%GRASPILOT_SELFTEST_PK Pharmacokinetic models: RATIO conversion, Tofts, extended Tofts, 2CXM.
%
%   nfail = graspilot_selftest_pk(out, doGUI)   (called from graspilot_selftest)
%
%   Uses graspilot_pk_phantom: spoiled gradient echo signals of regions with
%   known Tofts, extended Tofts and 2CXM parameters and a vessel (SSS) whose
%   blood concentration is the population AIF.
if nargin < 1 || isempty(out), out = tempname; mkdir(out); end
if nargin < 2, doGUI = false; end
nfail = 0;

% ---- signal ratio (RATIO) conversion is the exact inverse of the SPGR equation
TR = 4e-3; FA = 12; r1 = 4.5; T10 = 1.4;
spgr = @(R1) sind(FA) * (1 - exp(-TR * R1)) ./ (1 - cosd(FA) * exp(-TR * R1));
C = [0 0.01 0.1 0.5 1 3 8];
S = 0.37 * spgr(1 / T10 + r1 * C);
[Cc, valid] = graspilot_signal2conc([S; 2 * S], [S(1); 2 * S(1)], T10, TR, FA, r1);
e = max(max(abs(Cc - [C; C])));
nfail = nfail + check('RATIO conversion inverts the SPGR equation', e < 1e-9 && all(valid(:)), e);
[Cc, valid] = graspilot_signal2conc([1 1 100], 1, T10, TR, FA, r1);
nfail = nfail + check('signal ratio above the SPGR limit -> NaN', ~valid(3) && isnan(Cc(3)) && abs(Cc(1)) < 1e-12, 0);
Cg = graspilot_signal2conc(S(1:5), S(1), 1.0, TR, FA, r1);
nfail = nfail + check('shorter global T10 -> higher concentration', all(Cg(2:5) > C(2:5)), Cg(5) - C(5));

% ---- population AIF
tt = -10:0.5:300;
Cb = graspilot_pop_aif(tt);
[pk, ip] = max(Cb);
nfail = nfail + check('population AIF: 0 before injection, first pass ~10 s', all(Cb(tt < 0) == 0) && ...
    abs(tt(ip) - 10.2) < 1 && pk > 5 && pk < 7, tt(ip));

% ---- fits of curves on uneven frame times (Tofts, extended Tofts, 2CXM)
t = [0 3 6 9 12 14 16 18 20 23 26 30 35 40 50 60 75 90 110 130 150 180 210 240];
Cp = graspilot_pop_aif(t - 13) / 0.55;
tf = 0:0.01:t(end); Cpf = interp1(t, Cp, tf); d = 0.01 / 60;
ex = @(k) interp1(tf, filter(d * [0.5 0.5], [1 -exp(-k * d)], Cpf), t);
Ct = [0.3 * ex(1.5); 0.12 * ex(0.4) + 0.04 * Cp];
R = graspilot_pk_fit(Ct, t, Cp, {'tofts', 'etofts'});
e = max(abs([R.tofts.Ktrans(1) - 0.3, R.tofts.ve(1) - 0.2, R.etofts.Ktrans(2) - 0.12, ...
    R.etofts.ve(2) - 0.3, R.etofts.vp(2) - 0.04]));
nfail = nfail + check('Tofts / extended Tofts, uneven frame times', e < 1e-3, e);
nfail = nfail + check('extended Tofts on a Tofts curve: vp = 0', R.etofts.vp(1) < 1e-3 && ...
    abs(R.etofts.Ktrans(1) - 0.3) < 1e-3, R.etofts.vp(1));
[R, Fit] = graspilot_pk_fit([Ct; NaN(1, numel(t))], t, Cp, {'tofts'});
nfail = nfail + check('NaN curve -> NaN parameters; fitted curves returned', isnan(R.tofts.Ktrans(3)) && ...
    max(abs(Fit.tofts(1, :) - Ct(1, :))) < 1e-3 * max(Ct(1, :)), 0);

% ---- data model: maps of the phantom with the SSS as AIF
[D, tr] = graspilot_pk_phantom();
aif = mean(double(D.timecourses(find(tr.SSS))), 1);
bf = tr.BaseFrames;
D.computePK(bf, aif, tr.Head, []);
truth = {'Tofts', 'KT_TM', 0.25; 'Tofts', 'VE_TM', 0.25; 'Tofts', 'KEP_TM', 1; ...
    'ETofts', 'KT_ETM', 0.10; 'ETofts', 'VE_ETM', 0.20; 'ETofts', 'VP_ETM', 0.05; ...
    'CXM', 'FP_2CXM', 50; 'CXM', 'PS_2CXM', 12; 'CXM', 'VP_2CXM', 0.08; 'CXM', 'VE_2CXM', 0.25; ...
    'CXM', 'KT_2CXM', tr.P.CXM.Ktrans};
e = pkErr(D, tr, truth);
nfail = nfail + check('maps (SSS AIF, RATIO): Tofts, ext. Tofts, 2CXM', e < 0.01, e);
A = D.canonPK('KT_TM');
nfail = nfail + check('voxels outside the region are not fitted (NaN)', all(isnan(A(~tr.Head))) && ...
    all(isfinite(A(tr.Head))), 0);
[Cp, tq] = D.pkAIF(bf, aif);
e = max(abs(Cp - tr.Cp)) / max(tr.Cp);
nfail = nfail + check('AIF: SSS signal -> plasma concentration', e < 1e-4 && isequal(tq, tr.t), e);
nfail = nfail + check('maps valid for these settings', D.hasPK(bf, aif, GRASPilotData.maskHash(tr.Head)), 0);
D.setPKOptions('T10', 1200);
nfail = nfail + check('changing the global T10 invalidates the maps', ~D.hasPK(bf, aif, ...
    GRASPilotData.maskHash(tr.Head)), 0);
R = D.fitPKCurves(mean(double(D.timecourses(find(tr.Tofts))), 1), bf, aif);
nfail = nfail + check('wrong (shorter) tissue T10 biases Ktrans upward', R.tofts.Ktrans > 0.27, R.tofts.Ktrans);
D.setPKOptions('T10', 1400, 'AIFShift', -3);
R = D.fitPKCurves(mean(double(D.timecourses(find(tr.Tofts))), 1), bf, aif);
nfail = nfail + check('AIF shift changes the fit', abs(R.tofts.Ktrans - 0.25) > 1e-3, R.tofts.Ktrans);
D.setPKOptions('AIFShift', 0, 'AIF', 'population', 'Models', {'tofts'});
R = D.fitPKCurves(mean(double(D.timecourses(find(tr.Tofts))), 1), bf, []);
nfail = nfail + check('population AIF (injection at the end of the baseline)', ...
    abs(R.tofts.Ktrans - 0.25) < 0.02 && isequal(fieldnames(R), {'tofts'}), R.tofts.Ktrans);
D.setPKOptions('AIF', 'ref', 'Models', {'tofts', 'etofts', '2cxm'});
try
    D.setPKOptions('Hct', 1.2);
    bad = false;
catch
    bad = true;
end
nfail = nfail + check('invalid PK option rejected', bad && D.PKOpts.Hct == 0.45, 0);

% ---- noise
[Dn, tr] = graspilot_pk_phantom(0.02);
aif = mean(double(Dn.timecourses(find(tr.SSS))), 1);
Dn.computePK(bf, aif, tr.Head, []);
e = pkErr(Dn, tr, truth);
nfail = nfail + check('maps with noise (2% of S0): medians unbiased', e < 0.03, e);

if ~doGUI, return; end
% ---- viewer: fit, statistics, maps, exports, invalidation, sessions
[D, tr] = graspilot_pk_phantom(0.01);
app = GRASPilot(D);
app.Quiet = true;
cleanup = onCleanup(@() app.close(true));
app.BaseFrames = [1 5];
rs = app.addROI('SSS', [1 0 0]); app.setVoxels(rs, find(tr.SSS), true);
rt = app.addROI('Tofts box', [1 1 0]); app.setVoxels(rt, find(tr.Tofts), true);
rc = app.addROI('2CXM box', [0 1 1]); app.setVoxels(rc, find(tr.CXM), true);
app.afterMaskChange();
app.openPKWindow();
nfail = nfail + check('PK keys hidden before the fit', ~any(strcmp(app.paramKeyList(), 'KT_TM')), 0);
ok = app.computePK('Models', {'tofts', '2cxm'}, 'Region', 'rois');
nfail = nfail + check('viewer: models fitted (ROIs only)', ok && app.pkReady() && ...
    D.PKMaps.Info.NVoxels == nnz(tr.SSS | tr.Tofts | tr.CXM) && ...
    any(strcmp(app.paramKeyList(), 'KT_TM')) && ~any(strcmp(app.paramKeyList(), 'KT_ETM')), D.PKMaps.Info.NVoxels);
[~, rows] = app.paramStatsTable();
v1 = statVal(rows, 'Tofts box', 'Tofts Ktrans', 7);
v2 = statVal(rows, 'Tofts box', 'Tofts Ktrans', 14);
v3 = statVal(rows, '2CXM box', '2CXM plasma flow Fp', 7);
nfail = nfail + check('statistics: Tofts Ktrans median / ROI curve, 2CXM Fp', abs(v1 - 0.25) < 0.005 && ...
    abs(v2 - 0.25) < 0.005 && abs(v3 - 50) < 1.5, v1);
app.setParamKey('KT_TM');
app.setParamDisplay(true);
nfail = nfail + check('Ktrans colour map shown', app.ParamOn && ~isempty(app.paramVolume('KT_TM')), 0);
files = app.exportMaps(out, 'pk', {'KT_TM', 'VP_2CXM'});
img = graspilot_nifti_read(files{1});
A = double(app.paramVolume('KT_TM'));
ok = isequaln(single(img), single(A));
nfail = nfail + check('PK map -> NIfTI', numel(files) == 2 && ok, 0);
files = app.exportParams(fullfile(out, 'pk.csv'));
fid = fopen(files{2}); hl = fgetl(fid); fclose(fid);
txt = fileread(files{3});
nfail = nfail + check('parameter export: PK columns and settings', ~isempty(strfind(hl, 'KT_TM [1/min]')) && ...
    ~isempty(strfind(hl, 'FP_2CXM [ml/100ml/min]')) && ~isempty(strfind(txt, 'RATIO method')) && ...
    ~isempty(strfind(txt, 'PK tissue T10')), 0);
fs = fullfile(out, 'pk_session.mat');
D.setPKOptions('T10', 1300);
nfail = nfail + check('settings change -> maps stale, map hidden', ~app.pkReady() && ...
    isempty(app.paramVolume('KT_TM')) && ~any(strcmp(app.paramKeyList(), 'KT_TM')), 0);
app.saveROIs(fs);
D.setPKOptions('T10', 1400);
app.loadROIs(fs);
nfail = nfail + check('ROI session restores the PK settings', D.PKOpts.T10 == 1300, D.PKOpts.T10);
D.setPKOptions('T10', 1400, 'Region', 'head');
app.computePK();
nfail = nfail + check('refit inside the head', app.pkReady() && D.PKMaps.Info.NVoxels == nnz(tr.Head), ...
    D.PKMaps.Info.NVoxels);
app.setVoxels(rs, find(tr.SSS, 4), false); app.afterMaskChange();
nfail = nfail + check('editing the SSS (AIF) invalidates the maps', ~app.pkReady(), 0);
app.refreshPKWindow();
nfail = nfail + check('PK window refreshes', ishghandle(app.KW.Fig), 0);
end

function e = pkErr(D, tr, truth)
% largest relative error of the median map value in each region
e = 0;
for i = 1:size(truth, 1)
    A = double(D.canonPK(truth{i, 2}));
    v = A(tr.(truth{i, 1}));
    e = max(e, abs(median(v) - truth{i, 3}) / truth{i, 3});
end
end

function v = statVal(rows, roi, par, col)
i = find(strcmp(rows(:, 1), roi) & strcmp(rows(:, 2), par), 1);
v = NaN;
if ~isempty(i), v = rows{i, col}; end
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
