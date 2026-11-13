function ok = graspilot_selftest(doGUI)
%GRASPILOT_SELFTEST Synthetic end-to-end test of the GRASPilot toolbox.
%
%   graspilot_selftest          core tests (data model vs ReadData.m, orientation,
%                          NIfTI round trip, mask reorientation, exports)
%   graspilot_selftest(true)    also opens the viewer and drives it by script
%
%   The reference is the original ReadData.m logic using the TempPCASub
%   class in legacy/@TempPCASub, run on random complex data.

if nargin < 1, doGUI = false; end
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(fullfile(root, 'legacy'));
out = tempname; mkdir(out);
fprintf('GRASPilot self-test (output in %s)\n', out);
nfail = 0;

% ---------------------------------------------------------------- data
rand('state', 1); randn('state', 2); %#ok<RAND>
nz = 7; nxF = 40; nyF = 40; K = 4; T = 23;
vk_all = complex(randn(nz, nxF, nyF, K, 'single'), randn(nz, nxF, nyF, K, 'single'));
basis_all = complex(randn(nz, T, K, 'single'), randn(nz, T, K, 'single'));

% ---------------------------------------------------------------- reference (ReadData.m verbatim)
clear data1
Recon_Allslice = permute(vk_all, [2 3 4 1]);
for ii = 1:size(basis_all, 1)
    Phi = single(squeeze(squeeze(basis_all(ii, :, :))));
    PCA = TempPCASub(Phi);
    data1(:, :, :, ii) = abs(PCA' * Recon_Allslice(:, :, :, ii)); %#ok<AGROW>
end
data1 = CropImg(data1 / max(data1(:)), size(vk_all, 2) / 2, size(vk_all, 2) / 2);
ref = permute(data1, [1 2 4 3]);                     % [nx ny nz T]
tol = 1e-5;

% ---------------------------------------------------------------- 1. subspace model
D = GRASPilotData.fromSubspace(vk_all, basis_all);
e = 0;
for t = 1:T
    for z = 1:nz
        e = max(e, max(max(abs(D.nativeSlice(3, z, t) - ref(:, :, z, t)))));
    end
    for x = [1 7 20]
        e = max(e, max(max(abs(D.nativeSlice(1, x, t) - squeeze(ref(x, :, :, t))))));
        e = max(e, max(max(abs(D.nativeSlice(2, x, t) - squeeze(ref(:, x, :, t))))));
    end
    e = max(e, max(abs(reshape(D.nativeVolume(t) - ref(:, :, :, t), [], 1))));
end
nfail = nfail + check('subspace frames == ReadData.m', e < tol, e);
lin = [1; 77; 500; numel(ref(:, :, :, 1))];
R2 = reshape(ref, [], T);
e = max(max(abs(D.nativeTC(lin) - R2(lin, :))));
nfail = nfail + check('subspace voxel time courses', e < tol, e);
e = max(abs(D.nativeValuesAtFrame(lin, 9) - R2(lin, 9)));
nfail = nfail + check('subspace values at frame', e < tol, e);
e = max(max(abs(D.nativeTC(lin(2)) - R2(lin(2), :))));
nfail = nfail + check('single-voxel time course', e < tol, e);

% shared basis [T K] and basis stored as [nz K T]
Dp = GRASPilotData.fromSubspace(vk_all, permute(basis_all, [1 3 2]));
e = max(abs(reshape(Dp.nativeVolume(5) - ref(:, :, :, 5), [], 1)));
nfail = nfail + check('basis stored as [slices K T]', e < tol, e);

% ---------------------------------------------------------------- 2. dense model (ReadData data1 order)
Dd = GRASPilotData.fromDense(data1, 'xytz');
e = max(abs(reshape(Dd.nativeVolume(11) - ref(:, :, :, 11), [], 1)));
e = max(e, max(max(abs(Dd.nativeTC(lin) - R2(lin, :)))));
nfail = nfail + check('dense model (xytz)', e < tol, e);

% ---------------------------------------------------------------- 3. orientation algebra
planes = {'axial', 'coronal', 'sagittal'};
e = 0; ec = 0;
for ip = 1:3
    for sw = 0:1
        for fl = 0:7
            f = logical(bitget(fl, 1:3));
            D.setOrientation(planes{ip}, sw, f);
            t = 1 + mod(ip + sw + fl, T);
            Vc = D.native2canon(ref(:, :, :, t));
            e = max(e, max(abs(reshape(D.canonVolume(t) - Vc, [], 1))));
            n = D.canonSize();
            for a = 1:3
                q = max(1, round(n(a) / 3));
                switch a
                    case 1, S = reshape(Vc(q, :, :), n(2), n(3));
                    case 2, S = reshape(Vc(:, q, :), n(1), n(3));
                    case 3, S = Vc(:, :, q);
                end
                e = max(e, max(max(abs(D.canonSlice(a, q, t) - S))));
            end
            cl = [1; 5; 333; prod(n)];
            Vall = zeros(prod(n), T, 'single');
            for tt = 1:T
                tmp = D.native2canon(ref(:, :, :, tt));
                Vall(:, tt) = tmp(:);
            end
            e = max(e, max(max(abs(D.timecourses(cl) - Vall(cl, :)))));
            M = rand(n) > 0.7;
            ec = ec + ~isequal(D.native2canon(D.canon2native(M)), M);
        end
    end
end
nfail = nfail + check('canonical slices/volumes/TCs, 48 orientations', e < tol, e);
nfail = nfail + check('canon2native/native2canon inverse', ec == 0, ec);

% ---------------------------------------------------------------- 3b. axis codes
e = 0;
for ip = 1:3
    for sw = 0:1
        for fl = 0:7
            f = logical(bitget(fl, 1:3));
            D.setOrientation(planes{ip}, sw, f);
            P0 = D.Perm; F0 = D.Flip; c0 = D.AxisCodes;
            D.setAxisCodes(c0);
            e = e + ~isequal(D.Perm, P0) + ~isequal(D.Flip, F0) + ~strcmp(D.AxisCodes, c0);
        end
    end
end
nfail = nfail + check('axis codes <-> plane/swap/flip (48 combos)', e == 0, e);
D.setAxisCodes('LPS');
e = max(abs(reshape(D.canonVolume(3) - ref(:, :, :, 3), [], 1)));
nfail = nfail + check('LPS = stored array order (identity)', isequal(D.Perm, [1 2 3]) && ~any(D.Flip) && e < tol, e);
e = 0;
cs = {'LPS', 'PLS', 'RAS', 'ALI', 'ILP', 'SRA'};
for ic = 1:numel(cs)
    D.setAxisCodes(cs{ic});
    a0 = D.canonSlice(3, 2, 4).';                     % axial screen image
    D.setAxisCodes(GRASPilotData.rotateCodesCW(cs{ic}));
    a1 = D.canonSlice(3, 2, 4).';
    e = e + ~isequal(a1, rot90(a0, -1));
end
nfail = nfail + check('Rotate axial 90 CW == rot90(screen, -1)', e == 0, e);
D.setAxisCodes('PLS');
nfail = nfail + check('flip L-R code', strcmp(GRASPilotData.flipCodes('LPS', 'LR'), 'RPS'), 0);
bad = 0;
try, D.setAxisCodes('LLS'); catch, bad = 1; end
nfail = nfail + check('invalid axis codes rejected', bad == 1, 0);

% ---------------------------------------------------------------- 3c. SSS detection (venous phantom)
[Dp, tr] = graspilot_sss_phantom();
[linS, infoS] = graspilot_detect_sss(Dp);
np = Dp.canonSize();
M = false(np); M(linS) = true;
[~, ~, kk] = ind2sub(np, linS);
ok = numel(linS) > 50 && all(tr.SSS(linS)) && ~any(M(:) & (tr.TS(:) | tr.SS(:) | tr.SV(:) | tr.AR(:)));
nfail = nfail + check('SSS ROI only in SSS (not TS/straight/scalp/artery)', ok, numel(linS));
ok = abs(infoS.TorculaSlice - tr.TorculaZ) <= 2 && min(kk) >= infoS.TorculaSlice + 5 && min(kk) <= tr.TorculaZ + 10 ...
    && abs(infoS.Centroid(1) - 50) < 2 && isequal(unique(kk)', min(kk):max(kk));
nfail = nfail + check('SSS ROI starts just above torcula, midline, contiguous', ok, infoS.TorculaSlice);
% extends up the posterior SSS toward the vertex: 50 mm above the thickest slice, or until the
% sinus turns forward (phantom SSS: posterior edge moves forward > 1 mm per mm above z ~ 70)
ok = max(kk) >= infoS.ThickestSlice + 40 && max(kk) <= 75 && numel(linS) > 2000 && ...
    ~isempty(strfind(infoS.StopReason, 'forward'));
nfail = nfail + check('SSS ROI extends superiorly toward the vertex', ok, max(kk));
[linV, infoV] = graspilot_detect_sss(Dp, 'SuperiorMM', 20);
[~, ~, kv] = ind2sub(np, linV);
[linB, infoB] = graspilot_detect_sss(Dp, 'SlabMM', 5);
[~, ~, kB] = ind2sub(np, linB);
ok = max(kv) == infoV.ThickestSlice + 20 && all(tr.SSS(linV)) && max(kB) - min(kB) <= 4 && ...
    all(ismember(linB, linS)) && ~isempty(strfind(infoV.StopReason, 'SuperiorMM'));
nfail = nfail + check('SSS ROI extent options (SuperiorMM, SlabMM)', ok, max(kv));
[Dp2, tr2] = graspilot_sss_phantom(25);
[lin2, info2] = graspilot_detect_sss(Dp2);
ok = numel(lin2) > 50 && all(tr2.SSS(lin2)) && isnan(info2.TorculaSlice);
nfail = nfail + check('SSS detection without torcula in slab', ok, numel(lin2));
Dp2.setAxisCodes('PLS');
[lin3, info3] = graspilot_detect_sss(Dp2);
nfail = nfail + check('SSS detection reports failure, wrong orientation', ...
    isempty(lin3) || ~all(tr2.SSS(Dp2.canon2nativeLin(lin3))), numel(lin3));

% ---------------------------------------------------------------- 4. enhancement maps
D.setOrientation('axial', false, [0 0 0]);
D.computeMaps(1:3);
b = mean(ref(:, :, :, 1:3), 4); mx = max(ref(:, :, :, 4:end), [], 4);   % peak after the baseline
e = max(abs(reshape(D.Maps.Baseline - b, [], 1)));
e = max(e, max(abs(reshape(D.Maps.Max - mx, [], 1))));
nfail = nfail + check('enhancement maps', e < tol, e);

% ---------------------------------------------------------------- 4a. curve parameters (analytic)
% baseline 1 (frames 1-5, t = 0-4 s), linear rise 0.2/s to 3 at 14 s,
% linear washout -0.05/s to 2.2 at 30 s
tt = 0:30;
Sa = ones(size(tt));
up = tt > 4 & tt <= 14; Sa(up) = 1 + 0.2 * (tt(up) - 4);
dn = tt > 14; Sa(dn) = 3 - 0.05 * (tt(dn) - 14);
Sp = Sa; Sp(dn) = 3 + 0.05 * (tt(dn) - 14);            % persistent: keeps rising
P = graspilot_curve_params([Sa; Sp], tt, 1:5, struct('OnsetFrac', 0.05));   % onset = end of baseline
e = max(abs([P.PE(1) - 2, P.RE(1) - 200, P.TTP(1) - 10, P.Onset(1), P.WashIn(1) - 0.2, ...
    P.WashOut(1) + 0.05, P.AUC(1) - 35.6]));
nfail = nfail + check('curve parameters, washout curve (analytic)', e < 1e-9, e);
nfail = nfail + check('curve parameters, persistent curve (analytic)', ...
    abs(P.WashOut(2) - 0.05) < 1e-9 && abs(P.TTP(2) - 26) < 1e-9, P.WashOut(2));
P2 = graspilot_curve_params(Sa, tt, 1:5, struct('WashIn', 'max', 'WashOut', 'endpoints', 'AUCWindow', 10));
e = max(abs([P2.WashIn - 0.2, P2.WashOut + 0.05, P2.AUC - 10]));
nfail = nfail + check('curve parameters, max slope / endpoints / IAUC10', e < 1e-9, e);
% T90 at 4 + 9 s; MITR 2/10; early frame 10 s after baseline (t = 14, the peak), late = last:
% SER = 2/1.2, RE early 200%, washout % = 100 (3 - 2.2)/3, plateau = mean(1.3 1.25 1.2), type 3
P3 = graspilot_curve_params([Sa; Sp], tt, 1:5, struct('EarlyTime', 10));
e = max(abs([P3.T90(1) - 9, P3.MITR(1) - 0.2, P3.SER(1) - 2 / 1.2, P3.REearly(1) - 200, ...
    P3.WashoutPct(1) - 100 * 0.8 / 3, P3.Plateau(1) - 1.25, P3.CurveType(1) - 3, P3.CurveType(2) - 1]));
nfail = nfail + check('T90, MITR, SER, early RE, washout %, plateau, type', e < 1e-9, e);
P4 = graspilot_curve_params(Sa, tt, 1:5, struct('EndTime', 20));       % first 20 s: frames 1-21
P5 = graspilot_curve_params(Sa(1:21), tt(1:21), 1:5);
P6 = graspilot_curve_params(Sa, tt, 1:5, struct('EndFrame', 21));
nfail = nfail + check('analysis endpoint (EndTime, EndFrame)', P4.EndFrame == 21 && isequal(P4.AUC, P5.AUC) && ...
    isequal(P4.WashOut, P5.WashOut) && isequal(P6.AUC, P5.AUC) && abs(P4.AUC - (10 + 6 * (2 + 1.7) / 2)) < 1e-9, P4.AUC);
% reference rising 0.1/s (t 4-6), 0.5/s (6-10), 0.3/s (10-12), 0.05/s (12-14), then falling
Sr = ones(size(tt));
Sr(tt > 4) = 1 + 0.1 * (tt(tt > 4) - 4); Sr(tt > 6) = 1.2 + 0.5 * (tt(tt > 6) - 6);
Sr(tt > 10) = 3.2 + 0.3 * (tt(tt > 10) - 10); Sr(tt > 12) = 3.8 + 0.05 * (tt(tt > 12) - 12);
Sr(tt > 14) = 3.9 - 0.02 * (tt(tt > 14) - 14);
[wr, ir] = graspilot_ref_window(Sr, tt, 1:5);
Pw = graspilot_curve_params([Sa; Sr], tt, 1:5, struct('WashIn', 'sss', 'WashInFrames', wr));
nfail = nfail + check('SSS wash-in period and slopes within it', isequal(wr, [7 13]) && abs(ir.MaxSlope - 0.5) < 1e-9 && ...
    abs(Pw.WashIn(1) - 0.2) < 1e-9 && abs(Pw.WashIn(2) - 0.5) < 1e-9, wr(1));
Sn = [Sa; 2 * Sa; Sp] + 0.01 * randn(3, numel(tt));
Pv = graspilot_curve_params(Sn, tt, 1:5);
e = 0;
for r = 1:3
    Pr = graspilot_curve_params(Sn(r, :), tt, 1:5);
    for fn = {'S0', 'PE', 'TTP', 'Onset', 'WashIn', 'WashOut', 'AUC'}
        e = max(e, abs(Pr.(fn{1}) - Pv.(fn{1})(r)));
    end
end
nfail = nfail + check('curve parameters, vectorized == per curve', e < 1e-12, e);

% ---------------------------------------------------------------- 4a2. parameter maps == curve parameters
D.dt = 2.5; D.t0 = 0;
D.computeMaps(1:3);
nc = prod(D.canonSize());
Pm = graspilot_curve_params(double(D.timecourses((1:nc)')), D.paramTimes(), 1:3, D.ParamOpts);
e = mapErr(D, Pm, {'Baseline', 'S0'; 'PeakEnh', 'PE'; 'RelEnh', 'RE'; 'TTP', 'TTP'; 'Onset', 'Onset'; ...
    'WashIn', 'WashIn'; 'WashOut', 'WashOut'; 'AUC', 'AUC'});
nfail = nfail + check('parameter maps == graspilot_curve_params (s)', e < 1e-4, e);
D.setParamOptions('WashIn', 'max', 'WashOut', 'endpoints', 'AUCWindow', 20);
stale = ~D.hasMaps(1:3);
D.computeMaps(1:3);
Pm = graspilot_curve_params(double(D.timecourses((1:nc)')), D.paramTimes(), 1:3, D.ParamOpts);
e = mapErr(D, Pm, {'WashIn', 'WashIn'; 'WashOut', 'WashOut'; 'AUC', 'AUC'});
nfail = nfail + check('parameter options invalidate and recompute maps', stale && e < 1e-4, e);
D.setParamOptions('WashIn', 'onset', 'WashOut', 'regression', 'AUCWindow', Inf);
D.dt = NaN;
nfail = nfail + check('timing change invalidates maps', ~D.hasMaps(1:3) && isequal(D.paramTimes(), 1:T), 0);
D.dt = 2.5;

% ---------------------------------------------------------------- 4b. temporal binning
% each binned frame = mean magnitude of n consecutive ReadData.m frames
for mode = 1:2
    if mode == 1, Db = GRASPilotData.fromSubspace(vk_all, basis_all); lab = 'subspace';
    else, Db = GRASPilotData.fromDense(data1, 'xytz'); lab = 'dense'; end
    Db.setOrientation('coronal', true, [1 0 0]);
    Db.dt = 2; Db.t0 = 1;
    Db.setBinning(3);                                  % 23 frames -> 7 bins, frames 22-23 dropped
    e = 0;
    okn = Db.NT == 7 && isequal(Db.droppedFrames(), [22 23]);
    refB = zeros([size(ref, 1) size(ref, 2) size(ref, 3) Db.NT], 'single');
    for b = 1:Db.NT
        refB(:, :, :, b) = mean(ref(:, :, :, (b - 1) * 3 + (1:3)), 4);
    end
    for b = [1 4 7]
        e = max(e, max(abs(reshape(Db.nativeVolume(b) - refB(:, :, :, b), [], 1))));
        for d = 1:3
            q = 3;
            switch d
                case 1, R = squeeze(refB(q, :, :, b));
                case 2, R = squeeze(refB(:, q, :, b));
                case 3, R = refB(:, :, q, b);
            end
            e = max(e, max(max(abs(Db.nativeSlice(d, q, b) - R))));
        end
        e = max(e, max(abs(Db.nativeValuesAtFrame(lin, b) - reshape(refB(lin + (b - 1) * numel(refB(:, :, :, 1))), [], 1))));
    end
    RB = reshape(refB, [], Db.NT);
    e = max(e, max(max(abs(Db.nativeTC(lin) - RB(lin, :)))));
    ft = Db.frameTimes();
    okt = abs(ft(1) - 3) < 1e-9 && abs(ft(2) - 9) < 1e-9 && abs(Db.effectiveDt() - 6) < 1e-9;
    nfail = nfail + check(sprintf('binning x3 == mean of ReadData frames (%s)', lab), ...
        e < tol && okn && okt, e);
    Db.computeMaps(1:2);
    e = max(abs(reshape(Db.Maps.Baseline - mean(refB(:, :, :, 1:2), 4), [], 1)));
    e = max(e, max(abs(reshape(Db.Maps.Max - max(refB(:, :, :, 3:end), [], 4), [], 1))));
    nfail = nfail + check(sprintf('maps on binned series (%s)', lab), e < tol, e);
end
Db.setBinning(5, true);                                % keep partial last bin (frames 21-23)
e = max(abs(reshape(Db.nativeVolume(5) - mean(ref(:, :, :, 21:23), 4), [], 1)));
nfail = nfail + check('binning keeps partial last bin', Db.NT == 5 && e < tol, e);
Db.setBinning(1);
e = max(abs(reshape(Db.nativeVolume(T) - ref(:, :, :, T), [], 1)));
nfail = nfail + check('binning reset to 1', Db.NT == T && e < tol, e);

% ---------------------------------------------------------------- 4c. subtraction map (last - first)
D.setOrientation('sagittal', false, [0 1 0]);
S = D.canonSubtraction(D.NT, 1);
e = max(abs(reshape(S - D.native2canon(ref(:, :, :, T) - ref(:, :, :, 1)), [], 1)));
S2 = D.canonSubtraction(D.NT, 1);
nfail = nfail + check('subtraction map == data1(last) - data1(first)', e < tol && isequal(S, S2), e);
D.setOrientation('axial', false, [0 0 0]);

% ---------------------------------------------------------------- 4d. motion estimate and correction
Pm = [1.5 -2 1 2 -1.5 3];
[Dm, trm] = graspilot_motion_phantom(Pm, 7:12);
Mm = graspilot_motion_estimate(Dm);
np = Dm.canonSize(); vxm = Dm.canonVox();
R0 = double(trm.Still.rawCanonVolume(1));
[a1, a2, a3] = ind2sub(np, find(R0 > 0.2));
Xh = ([a1 a2 a3] - 1) .* repmat(vxm, numel(a1), 1);
fe = zeros(1, Dm.NTRaw);
for f = 1:Dm.NTRaw
    Rt = graspilot_rigid_matrix(trm.Params(f, :)); Re = graspilot_rigid_matrix(Mm.Params(f, :));
    Yt = (Xh - repmat(trm.Centre, size(Xh, 1), 1)) * Rt' + repmat(trm.Centre + trm.Params(f, 1:3), size(Xh, 1), 1);
    Ye = (Xh - repmat(Mm.Centre, size(Xh, 1), 1)) * Re' + repmat(Mm.Centre + Mm.Params(f, 1:3), size(Xh, 1), 1);
    fe(f) = mean(sqrt(sum((Yt - Ye) .^ 2, 2)));
end
nfail = nfail + check('motion estimate: 2-3 mm / 2-3 deg motion recovered (mean field error)', ...
    max(fe) < 0.2 && max(fe(1:6)) < 0.05, max(fe));
nfail = nfail + check('excess motion flagged in the moved frames only', Mm.Excess && isequal(Mm.ExcessFrames, 7:12) ...
    && ~isempty(strfind(Mm.Message, 'WARNING')) && Mm.RelMean(7) > 2 && all(Mm.RelMean(8:12) < 0.2), Mm.RelMean(7));
Ms = graspilot_motion_estimate(trm.Still, 'Frames', [1 6 9 12]);
nfail = nfail + check('no motion: no warning', ~Ms.Excess && max(Ms.AbsMean) < 0.05 && isequal(Ms.Frames, [1 6 9 12]), ...
    max(Ms.AbsMean));
Dm.applyMotionCorrection(Mm);
hm = R0 > 0.2;
e1 = zeros(1, 2); e0 = e1;
for q = 1:2
    f = [9 12]; f = f(q);
    B = double(trm.Still.rawCanonVolume(f));
    A1 = double(Dm.rawCanonVolume(f)) - B; A0 = double(Dm.native2canon(Dm.Orig.D(:, :, :, f))) - B;
    e1(q) = sqrt(mean(A1(hm) .^ 2));
    e0(q) = sqrt(mean(A0(hm) .^ 2));
end
same = isequal(Dm.rawCanonVolume(3), Dm.native2canon(Dm.Orig.D(:, :, :, 3)));
nfail = nfail + check('motion correction realigns moved frames, keeps still ones', Dm.MotionCorrected && ...
    all(e1 < 0.4 * e0) && same && strcmp(Dm.Mode, 'dense'), max(e1 ./ e0));
Dm.setBinning(3);
tcm = Dm.timecourses(find(hm, 5));
okB = Dm.NT == 4 && abs(tcm(1, 3) - mean(Dm.D(Dm.canon2nativeLin(find(hm, 1)) + (6:8) * prod(np)))) < 1e-6;
Dm.setBinning(1);
Dm.revertMotionCorrection();
nfail = nfail + check('binning on corrected frames; revert restores the original', okB && ~Dm.MotionCorrected && ...
    isempty(Dm.Orig) && isequal(Dm.rawCanonVolume(9), Dm.native2canon(Dm.D(:, :, :, 9))) && ...
    ~isempty(strfind(Dm.motionDescription(), 'not applied')), 0);

% ---------------------------------------------------------------- 5. NIfTI round trip
D.setOrientation('coronal', true, [1 0 1]);
D.VoxNative = [0.9 1.1 2.5]; D.dt = 3.2; D.t0 = 0;
f4 = fullfile(out, 'series.nii.gz');
graspilot_write_dce_nifti(D, f4);
[img, h] = graspilot_nifti_read(f4);
e = 0;
for t = [1 9 T]
    e = max(e, max(abs(reshape(img(:, :, :, t) - D.canonVolume(t), [], 1))));
end
nfail = nfail + check('4D NIfTI round trip (data)', e < tol && isequal(h.Dim, [D.canonSize() T]), e);
A0 = D.affine();
e = max(abs(h.Affine(:) - A0(:)));
nfail = nfail + check('4D NIfTI affine', e < 1e-5 && abs(h.dt - 3.2) < 1e-5, e);
nfail = nfail + check('JSON sidecar written', exist(fullfile(out, 'series.json'), 'file') == 2, 0);

% qform-only reading path
hq = struct('Dim', D.canonSize(), 'Vox', D.canonVox(), 'Affine', D.affine(), 'Datatype', 'uint8');
fq = fullfile(out, 'q.nii');
graspilot_nifti_write(fq, uint8(D.native2canon(ref(:, :, :, 1)) > 0.3), hq);
fid = fopen(fq, 'r+', 'ieee-le'); fseek(fid, 254, 'bof'); fwrite(fid, 0, 'int16'); fclose(fid);
[~, h2] = graspilot_nifti_read(fq);
e = max(abs(h2.Affine(:) - A0(:)));
nfail = nfail + check('qform -> affine', e < 1e-4, e);

% ---------------------------------------------------------------- 6. mask reorientation
n = D.canonSize();
M = false(n); M(3:9, 2:min(5, n(2)), 2:5) = true; M(2, 2, 1) = true;
A = D.affine();
arr2 = permute(flip(M, 1), [2 1 3]);     % another tool's voxel order
Mx = [0 -1 0 n(1) - 1; 1 0 0 0; 0 0 1 0; 0 0 0 1];
fm = fullfile(out, 'reoriented_mask.nii.gz');
graspilot_nifti_write(fm, uint8(arr2), struct('Dim', size(arr2), 'Vox', D.canonVox(), ...
    'Affine', A * Mx, 'Datatype', 'uint8'));
[im2, h3] = graspilot_nifti_read(fm);
back = graspilot_reorient_to_ref(im2, h3.Affine, A, n) > 0;
nfail = nfail + check('mask reorientation on import', isequal(back, M), nnz(xor(back, M)));

% ---------------------------------------------------------------- 7. time-series export (csv)
masks = {M, false(n)};
masks{2}(10:12, n(2) - 1:n(2), 6:7) = true;
fcsv = fullfile(out, 'ts.csv');
graspilot_export_timeseries(D, masks, {'SSS (control)', 'Tumor 1'}, fcsv, 'BaseFrames', 1:3);
fid = fopen(fullfile(out, 'ts_voxels.csv')); hdrl = fgetl(fid); l1 = fgetl(fid); fclose(fid);
vals = str2double(regexp(l1, ',', 'split'));
lin1 = find(M, 1);
tc1 = D.timecourses(lin1);
e = max(abs(vals(10:end) - double(tc1)));
nfail = nfail + check('CSV voxel export values', e < 1e-5 && numel(vals) == 9 + T, e);
if exist('writecell', 'file') == 2
    fx = fullfile(out, 'ts.xlsx');
    graspilot_export_timeseries(D, masks, {'SSS (control)', 'Tumor 1'}, fx, 'BaseFrames', 1:3);
    nfail = nfail + check('XLSX export', exist(fx, 'file') == 2, 0);
end

% ---------------------------------------------------------------- 8. GUI
if doGUI
    nfail = nfail + graspilot_selftest_gui(D, ref, out);
    nfail = nfail + graspilot_selftest_gui_sss(out);
    nfail = nfail + graspilot_selftest_gui_motion(out);
end

ok = nfail == 0;
if ok
    fprintf('ALL TESTS PASSED\n');
else
    fprintf('%d TEST(S) FAILED\n', nfail);
end
end

function e = mapErr(D, P, pairs)
% max relative error between maps and curve parameters (NaN must match)
n = D.canonSize();
e = 0;
for k = 1:size(pairs, 1)
    A = double(D.canonMap(pairs{k, 1}));
    B = reshape(P.(pairs{k, 2}), n);
    if ~isequal(isnan(A), isnan(B)), e = Inf; return; end
    ok = ~isnan(A);
    sc = max(1e-6, max(abs(B(ok))));
    e = max(e, max(abs(A(ok) - B(ok))) / sc);
end
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
