function nfail = graspilot_selftest_pk(out, doGUI)
%GRASPILOT_SELFTEST_PK Pharmacokinetic models: signal-ratio conversion, Tofts, extended Tofts, 2CXM, RATIO.
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
nfail = nfail + check('signal-ratio conversion inverts the SPGR equation', e < 1e-9 && all(valid(:)), e);
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

% ---- automatic baseline: bolus arrival of the earliest strongly enhancing voxels
bfs = {};
[D1, t1] = graspilot_pk_phantom(0.02); bfs{1} = D1.detectBaseline();
[D2, t2] = graspilot_sss_phantom(); bfs{2} = D2.detectBaseline();
[D3, t3] = graspilot_sss_phantom([], 'realistic'); bfs{3} = D3.detectBaseline();
[D4, t4] = graspilot_ratio_phantom(0.005); bfs{4} = D4.detectBaseline();
D4.setBinning(2); bfs{5} = D4.detectBaseline();
ok = isequal(bfs{1}, 1:5) && isequal(bfs{2}, 1:3) && isequal(bfs{3}, 1:3) && isequal(bfs{4}, 1:5) && ...
    isequal(bfs{5}, 1:2);
nfail = nfail + check('automatic baseline (4 phantoms; binned x2)', ok, numel(bfs{4}));
Dn = GRASPilotData.fromDense(rand(8, 8, 2, 12, 'single'), 'xyzt');
[bfn, in] = Dn.detectBaseline();
nfail = nfail + check('no bolus: baseline falls back to frames 1-5', numel(bfn) <= 5 && bfn(1) == 1 && ...
    ~isempty(in.Message), numel(bfn));

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
nfail = nfail + check('maps (SSS AIF, signal ratio): Tofts, ETofts, 2CXM', e < 0.01, e);
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

% ---- RATIO kinetic core: classes and gain-free rates
t = 0:1.5:240;
u = graspilot_pop_aif(t - 15) / 0.55;
tf = 0:0.01:t(end); d = 0.01 / 60;
ex = @(k, del) interp1(tf, filter(d * [0.5 0.5], [1 -exp(-k * d)], interp1(t, u, max(tf - del, 0), ...
    'linear', 0)), t);
Fp = 0.5; PS = 0.12; vp = 0.08; ve = 0.25;
kpo = Fp / vp; kpe = PS / vp; kep = PS / ve; a = kpo + kpe; b = kep;
sq = sqrt((a + b) ^ 2 - 4 * kpo * b);
Kp = (a + b + sq) / 2; Km = (a + b - sq) / 2; A = (Kp - (kep + kpe)) / (Kp - Km);
E = PS / (PS + Fp);
Y = [7 * Fp * (A * ex(Kp, 3) + (1 - A) * ex(Km, 3)); 0.9 * ex(0.5, 4); 0.6 * 12 * ex(12, 2); 0.5 * ex(2.5, 2)];
randn('state', 5); %#ok<RAND>
Y = Y + 0.003 * repmat(max(Y, [], 2), 1, numel(t)) .* randn(size(Y));    % AICc needs noise, not round-off
R = graspilot_ratio_fit(Y, t, u);
e = max(abs([R.kep(1) / (kep * (1 - E)), R.kpe(1) / (kpo * E), R.MTT(1) / (60 / kpo), R.E(1) / E] - 1));
nfail = nfail + check('RATIO 2C: kep, kpe, MTT, E (arbitrary gain, delay)', R.Class(1) == 1 && e < 0.03 && ...
    abs(R.Delay(1) - 3) < 0.1 && R.Leak(1) == 1, e);
nfail = nfail + check('RATIO 1C classes: slow (kep = k), fast (MTT), indeterminate', isequal(R.Class(2:4)', [3 2 4]) && ...
    abs(R.kep(2) - 0.5) < 0.01 && abs(R.MTT(3) - 5) < 0.15 && R.Leak(3) == 0 && isnan(R.Leak(4)) && ...
    abs(R.k1(4) - 2.5) < 0.05 && isnan(R.kpe(4)), R.k1(4));
Q = orth([ones(numel(t), 1), t(:), t(:) .^ 2, exp(-t(:) / 30), exp(-t(:) / 100)]);
R = graspilot_ratio_fit(Y, t, u, struct('Basis', Q));
nfail = nfail + check('RATIO rank rule: 5 basis coefficients -> no 2C class', ~any(R.Class == 1) && ...
    all(R.NObs == 5), 0);

% ---- VCAI and RATIO on a complex subspace phantom
[Dr, tr] = graspilot_ratio_phantom();
bf = tr.BaseFrames;
nfail = nfail + check('complex subspace data: phase and temporal basis', Dr.hasPhase() && ~D.hasPhase() && ...
    size(Dr.pkBasis(1:8), 2) == 8, size(Dr.pkBasis(1:8), 2));
mtc = @(M) mean(double(Dr.timecourses(find(M))), 1);
ringM = false(size(tr.Artery)); ringM(13:18, 4:9, 2:7) = true; ringM(tr.Artery) = false;
aif = struct('Ref', mtc(tr.SSS), 'Art', mtc(tr.Artery), 'Ring', mtc(ringM), ...
    'Phase', Dr.phaseCurve(find(tr.SSS), bf));
Dr.setPKOptions('AIF', 'vcai', 'Models', {'ratio'});
[~, ~, Cb, vi] = Dr.pkAIF(bf, aif);
e = max(abs(Cb - tr.ca)) / max(tr.ca);
nfail = nfail + check('VCAI (sinus phase): w, tau, delay, arterial curve', strcmp(vi.Anchor, 'phase') && ...
    abs(vi.w - tr.WSignal) < 0.005 && abs(vi.Tau - tr.Tau) < 0.1 && abs(vi.Delta - tr.Delta) < 0.1 && e < 1e-3, e);
Dr.setPKOptions('VCAIAnchor', 'magnitude');
[~, ~, Cb, vi] = Dr.pkAIF(bf, aif);
nfail = nfail + check('VCAI (sinus magnitude anchor)', strcmp(vi.Anchor, 'magnitude') && ...
    max(abs(Cb - tr.ca)) / max(tr.ca) < 1e-3, vi.w);
Dr.setPKOptions('VCAIAnchor', 'phase');
cls = [1 3 2]; nm = {'Lesion2C', 'LesionSlow', 'LesionFast'};
okc = true;
for i = 1:3
    R = Dr.fitPKCurves(mtc(tr.(nm{i})), bf, aif, find(tr.(nm{i})));
    okc = okc && R.ratio.Class == cls(i);
    if i == 1
        P = tr.P2C;
        e = max(abs([R.ratio.kep / P.kep, R.ratio.kpe / P.kpe, R.ratio.MTT / P.MTT, R.ratio.E / P.E, ...
            R.ratio.vp / P.vp, R.ratio.Ktrans / P.Ktrans] - 1));
        nobs = R.ratio.NObs;
    end
end
nfail = nfail + check('RATIO on the phantom: classes, 2C rates and apparent vp, Ktrans', okc && e < 0.06 && ...
    nobs == 8, e);
Dr.computePK(bf, aif, tr.Lesion2C | tr.LesionSlow, []);
A = Dr.canonPK('CLASS_RATIO'); B = Dr.canonPK('KPE_RATIO');
nfail = nfail + check('RATIO maps (basis coefficients per slice)', all(A(tr.Lesion2C) == 1) && ...
    all(A(tr.LesionSlow) == 3) && abs(median(B(tr.Lesion2C)) / tr.P2C.kpe - 1) < 0.05 && ...
    all(isnan(B(tr.LesionSlow))) && Dr.PKMaps.Info.RatioNObs == 8, median(B(tr.Lesion2C)));

if ~doGUI, return; end
% ---- viewer: fit, statistics, maps, exports, invalidation, sessions
[D, tr] = graspilot_pk_phantom(0.01);
app = GRASPilot(D);
app.Quiet = true;
cleanup = onCleanup(@() app.close(true));
nfail = nfail + check('viewer: baseline detected on load (Auto on)', app.AutoBaseline && ...
    isequal(app.BaseFrames, [1 5]) && strcmp(get(app.H.base2, 'String'), '5'), app.BaseFrames(2));
set(app.H.base2, 'String', '3'); app.onBaseline();
ok = ~app.AutoBaseline && isequal(app.BaseFrames, [1 3]) && get(app.H.baseAuto, 'Value') == 0;
app.setAutoBaseline(true);
nfail = nfail + check('typed baseline turns Auto off; Auto re-detects', ok && isequal(app.BaseFrames, [1 5]) && ...
    get(app.H.baseAuto, 'Value') == 1, app.BaseFrames(2));
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
    ~isempty(strfind(hl, 'FP_2CXM [ml/100ml/min]')) && ~isempty(strfind(txt, 'signal-ratio conversion')) && ...
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
clear cleanup

% ---- viewer: RATIO with the VCAI (artery ROI, neighbour ring, SSS phase)
[Dr, tr] = graspilot_ratio_phantom(0.003);
app = GRASPilot(Dr);
app.Quiet = true;
cleanup = onCleanup(@() app.close(true));
app.BaseFrames = [1 5];
r = app.addROI('SSS', [1 0 0]); app.setVoxels(r, find(tr.SSS), true);
ra = app.addROI('Artery', [1 0.5 0]); app.setVoxels(ra, find(tr.Artery), true);
r = app.addROI('Lesion 2C', [1 1 0]); app.setVoxels(r, find(tr.Lesion2C), true);
app.afterMaskChange();
app.openPKWindow();
ring = app.ringMask(ra);
nfail = nfail + check('arterial ROI found by name; neighbour ring', app.arteryIndex() == ra && any(ring(:)) && ...
    ~any(ring(:) & (tr.Artery(:) | tr.SSS(:))), nnz(ring));
ok = app.computePK('Models', {'etofts', 'ratio'}, 'AIF', 'vcai', 'Region', 'rois');
[~, rows] = app.paramStatsTable();
v1 = statVal(rows, 'Lesion 2C', 'RATIO kpe = E/MTT', 7);
v2 = statVal(rows, 'Lesion 2C', 'RATIO class: two-compartment fraction', 5);
nfail = nfail + check('viewer: RATIO kpe and class fraction; vessels not in tumor stats', ok && ...
    abs(v1 / tr.P2C.kpe - 1) < 0.06 && v2 > 95 && ~any(strcmp(rows(:, 1), 'Artery')), v1);
S = app.pkSettings();
nfail = nfail + check('settings: RATIO model, VCAI fit', any(strcmp(S(:, 1), 'RATIO model')) && ...
    any(strcmp(S(:, 1), 'VCAI fit')), 0);
app.setVoxels(ra, find(tr.Artery, 2), false); app.afterMaskChange();
nfail = nfail + check('editing the arterial ROI invalidates the maps', ~app.pkReady(), 0);
app.refreshPKWindow();
clear cleanup

% ---- viewer: Auto AIF and the arterial ROI as AIF (head phantom)
[Dh, th] = graspilot_head_phantom();
[linA, infoA] = graspilot_detect_aif(Dh);
nfail = nfail + check('Auto AIF finds the artery (not the sinus or tumour)', numel(linA) >= 10 && ...
    mean(th.Art(linA)) > 0.9 && ~any(th.SSS(linA)) && ~any(th.Ring(linA)), mean(th.Art(linA)));
app = GRASPilot(Dh);
app.Quiet = true;
cleanup = onCleanup(@() app.close(true));
rt = app.addROI('Tumor 1', [1 1 0]); app.setVoxels(rt, find(th.Ring), true);
app.afterMaskChange();
app.openPKWindow();
ok = app.computePK('Models', {'etofts'}, 'AIF', 'auto', 'Region', 'rois');
ka = app.autoAIFIndex();
[~, rows] = app.paramStatsTable();
okA = ok && ka > 0 && mean(th.Art(app.Masks{ka})) > 0.9 && app.pkReady() && ...
    ~any(strcmp(rows(:, 1), 'AIF (auto)')) && any(strcmp(rows(:, 1), 'Tumor 1'));
nfail = nfail + check('AIF = Auto AIF: ROI "AIF (auto)" made and used, not in tumor stats', okA, ka);
nfail = nfail + check('Auto AIF puts the crosshair on one of its voxels', ka > 0 && ...
    app.Masks{ka}(app.Pos(1), app.Pos(2), app.Pos(3)), 0);
ra = app.addROI('Artery', [1 0.5 0]); app.setVoxels(ra, find(th.Art), true); app.afterMaskChange();
ok = app.computePK('AIF', 'artery');
a = app.aifInput();
nfail = nfail + check('AIF = Arterial ROI: its mean signal is the input', ok && app.pkReady() && ...
    max(abs(a.Art - app.ROIs(app.arteryIndex()).Stats.mean)) < 1e-9, 0);
app.KWRoi = 'Tumor 1'; app.refreshPKWindow();
tt = get(get(app.KW.Ax(2), 'Title'), 'String');
nfail = nfail + check('concentration plot titled with the ROI shown', ~isempty(strfind(tt, 'Tumor 1')), 0);
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
