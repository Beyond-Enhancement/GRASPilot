function nfail = graspilot_selftest_gui_io(out)
%GRASPILOT_SELFTEST_GUI_IO NIfTI / DICOM input in the viewer and NIfTI parameter-map output.
nfail = 0;
[S, tr] = graspilot_sss_phantom();
n0 = S.canonSize();
T = S.NTRaw;
% 2 mm phantom stored in RAS order (as most NIfTI files), scanner units ~ 0-1500, off-centre origin
Vc = zeros([ceil(n0 / 2) T], 'single');
for t = 1:T
    v = S.rawCanonVolume(t);
    Vc(:, :, :, t) = 1000 * v(1:2:end, 1:2:end, 1:2:end);
end
n = size(Vc(:, :, :, 1));
Vn = flip(flip(Vc, 1), 2);                                  % LPS -> RAS voxel order
A = [2 0 0 -47; 0 2 0 -61; 0 0 2 12; 0 0 0 1];
fin = fullfile(out, 'ras_dyn.nii.gz');
graspilot_nifti_write(fin, Vn, struct('Dim', size(Vn), 'Vox', [2 2 2], 'Affine', A, 'Datatype', 'float32', ...
    'Descrip', 'test', 'dt', 4));
app = GRASPilot();
app.Quiet = true;
cleanup = onCleanup(@() app.close(true));
app.openFile(fin);
D = app.Data;
okI = ~isempty(D) && strcmp(D.AxisCodes, 'RAS') && strcmp(D.SourceType, 'nifti') && abs(D.dt - 4) < 1e-6 && ...
    isequal(D.canonVolume(6), Vc(:, :, :, 6)) && max(abs(D.canonVox() - 2)) < 1e-6;
nfail = nfail + check('viewer opens 4D NIfTI: orientation, voxel size, timing from header', okI, D.dt);
% anatomy-aware tools work on it: Auto SSS finds the sinus
k = app.autoSSS();
ss = tr.SSS(1:2:end, 1:2:end, 1:2:end);
okS = k > 0 && nnz(app.Masks{k} & ss) >= 0.9 * nnz(app.Masks{k});
nfail = nfail + check('Auto SSS on NIfTI input (orientation from the affine)', okS, nnz(app.Masks{max(k, 1)}));
% parameter maps -> NIfTI, selected keys, in the scanner geometry of the input
app.computeParams();
files = app.exportMaps(out, 'io', {'PE', 'AUC_SSS', 'Baseline'});
[Pm, hp] = graspilot_nifti_read(files{1});
[Am, ~] = graspilot_nifti_read(files{2});
PE = app.paramVolume('PE');
% world position of canonical voxel c equals that of the source voxel it came from
c = [7 11 5];
wc = hp.Affine * [c - 1, 1]';
ws = A * [n(1) - c(1), n(2) - c(2), c(3) - 1, 1]';
okM = numel(files) == 3 && isequal(size(Pm), n) && max(abs(double(Pm(:)) - double(PE(:)))) < 1e-3 * max(PE(:)) && ...
    max(abs(wc - ws)) < 1e-4 && all(isfinite(Am(app.Masks{k}))) && ~isempty(strfind(files{2}, 'io_AUC_SSS.nii.gz'));
nfail = nfail + check('parameter maps -> NIfTI (chosen maps, source geometry)', okM, max(abs(wc - ws)));
fm = fullfile(out, 'io_mask.nii.gz');
app.exportMaskNifti(fm, false);
[Mk, hm] = graspilot_nifti_read(fm);
nfail = nfail + check('mask NIfTI in the source geometry', max(abs(hm.Affine(:) - hp.Affine(:))) < 1e-5 && ...
    isequal(Mk > 0, app.Masks{k}), 0);
% DICOM folder with two series
gd = struct('X', [1 0 0], 'Y', [0 1 0], 'ipp0', [-47 -61 12], 'ps', [2 2], 'dz', 2, 'dt', 4, 'Implicit', false, ...
    'SeriesNumber', 3, 'Description', 'dyn');
dd = fullfile(out, 'dcm_gui');
Vu = uint16(Vc(:, :, 1:4:end, 1:6));                         % LPS order = DICOM axial
graspilot_dicom_write(fullfile(dd, 's3'), Vu, gd);
g2 = gd; g2.SeriesNumber = 7; g2.Description = 'other';
graspilot_dicom_write(fullfile(dd, 's7'), Vu(:, :, :, 1:2), g2);
app.openFile(dd);                                            % quiet: the largest series
D = app.Data;
okD = strcmp(D.SourceType, 'dicom') && D.NTRaw == 6 && strcmp(D.AxisCodes, 'LPS') && abs(D.dt - 4) < 1e-6 && ...
    isequal(D.canonVolume(5), single(Vu(:, :, :, 5)));
app.openImages(dd, 'dicom', struct('Series', 2));
okD = okD && app.Data.NTRaw == 2 && ~isempty(strfind(app.Data.SourceVars, 'other'));
f1 = dir(fullfile(dd, 's3')); f1 = f1(~[f1.isdir]);
app.openFile(fullfile(dd, 's3', f1(1).name));                % one image of a series
okD = okD && app.Data.NTRaw == 6;
nfail = nfail + check('viewer opens DICOM (folder, series choice, single file)', okD, app.Data.NTRaw);
% uneven frame timing: display, parameters, export / re-import, typed s/frame
tu = [0 4 8 12 20 28];
gu = gd; gu.SeriesNumber = 5; gu.Description = 'uneven'; gu.Times = tu;
du = fullfile(out, 'dcm_uneven');
graspilot_dicom_write(du, Vu, gu);
app.openFile(du);
D = app.Data;
info = get(app.H.info, 'String'); if iscell(info), info = strjoin(info, ' '); end
okT = D.unevenTiming() && isequal(D.frameTimes(), tu) && ~isempty(strfind(info, 'uneven')) && ...
    isequal(D.paramTimes(), tu);
fe = fullfile(out, 'uneven_dce.nii.gz');
app.exportImageNifti(fe, 1:D.NT);
De = GRASPilotData.fromNifti(fe);
okT = okT && max(abs(De.rawTimes() - tu)) < 1e-4;
app.setTemporalResolution(5);
okT = okT && ~D.unevenTiming() && isequal(D.frameTimes(), (0:5) * 5);
nfail = nfail + check('uneven frame times: shown, used, exported, overridable', okT, D.dt);
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
