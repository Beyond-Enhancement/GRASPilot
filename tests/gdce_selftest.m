function ok = gdce_selftest(doGUI)
%GDCE_SELFTEST Synthetic end-to-end test of the GRASP DCE Viewer toolbox.
%
%   gdce_selftest          core tests (data model vs ReadData.m, orientation,
%                          NIfTI round trip, mask reorientation, exports)
%   gdce_selftest(true)    also opens the viewer and drives it by script
%
%   The reference is the original ReadData.m logic using the TempPCASub
%   class in legacy/@TempPCASub, run on random complex data.

if nargin < 1, doGUI = false; end
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(fullfile(root, 'legacy'));
out = tempname; mkdir(out);
fprintf('GRASP DCE Viewer self-test (output in %s)\n', out);
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
D = GraspData.fromSubspace(vk_all, basis_all);
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
Dp = GraspData.fromSubspace(vk_all, permute(basis_all, [1 3 2]));
e = max(abs(reshape(Dp.nativeVolume(5) - ref(:, :, :, 5), [], 1)));
nfail = nfail + check('basis stored as [slices K T]', e < tol, e);

% ---------------------------------------------------------------- 2. dense model (ReadData data1 order)
Dd = GraspData.fromDense(data1, 'xytz');
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

% ---------------------------------------------------------------- 4. enhancement maps
D.setOrientation('axial', false, [0 0 0]);
D.computeMaps(1:3);
b = mean(ref(:, :, :, 1:3), 4); mx = max(ref, [], 4);
e = max(abs(reshape(D.Maps.Baseline - b, [], 1)));
e = max(e, max(abs(reshape(D.Maps.Max - mx, [], 1))));
nfail = nfail + check('enhancement maps', e < tol, e);

% ---------------------------------------------------------------- 4b. temporal binning
% each binned frame = mean magnitude of n consecutive ReadData.m frames
for mode = 1:2
    if mode == 1, Db = GraspData.fromSubspace(vk_all, basis_all); lab = 'subspace';
    else, Db = GraspData.fromDense(data1, 'xytz'); lab = 'dense'; end
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
    e = max(e, max(abs(reshape(Db.Maps.Max - max(refB, [], 4), [], 1))));
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

% ---------------------------------------------------------------- 5. NIfTI round trip
D.setOrientation('coronal', true, [1 0 1]);
D.VoxNative = [0.9 1.1 2.5]; D.dt = 3.2; D.t0 = 0;
f4 = fullfile(out, 'series.nii.gz');
gdce_write_dce_nifti(D, f4);
[img, h] = gdce_nifti_read(f4);
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
gdce_nifti_write(fq, uint8(D.native2canon(ref(:, :, :, 1)) > 0.3), hq);
fid = fopen(fq, 'r+', 'ieee-le'); fseek(fid, 254, 'bof'); fwrite(fid, 0, 'int16'); fclose(fid);
[~, h2] = gdce_nifti_read(fq);
e = max(abs(h2.Affine(:) - A0(:)));
nfail = nfail + check('qform -> affine', e < 1e-4, e);

% ---------------------------------------------------------------- 6. mask reorientation
n = D.canonSize();
M = false(n); M(3:9, 2:min(5, n(2)), 2:5) = true; M(2, 2, 1) = true;
A = D.affine();
arr2 = permute(flip(M, 1), [2 1 3]);     % another tool's voxel order
Mx = [0 -1 0 n(1) - 1; 1 0 0 0; 0 0 1 0; 0 0 0 1];
fm = fullfile(out, 'reoriented_mask.nii.gz');
gdce_nifti_write(fm, uint8(arr2), struct('Dim', size(arr2), 'Vox', D.canonVox(), ...
    'Affine', A * Mx, 'Datatype', 'uint8'));
[im2, h3] = gdce_nifti_read(fm);
back = gdce_reorient_to_ref(im2, h3.Affine, A, n) > 0;
nfail = nfail + check('mask reorientation on import', isequal(back, M), nnz(xor(back, M)));

% ---------------------------------------------------------------- 7. time-series export (csv)
masks = {M, false(n)};
masks{2}(10:12, n(2) - 1:n(2), 6:7) = true;
fcsv = fullfile(out, 'ts.csv');
gdce_export_timeseries(D, masks, {'SSS (control)', 'Tumor 1'}, fcsv, 'BaseFrames', 1:3);
fid = fopen(fullfile(out, 'ts_voxels.csv')); hdrl = fgetl(fid); l1 = fgetl(fid); fclose(fid);
vals = str2double(regexp(l1, ',', 'split'));
lin1 = find(M, 1);
tc1 = D.timecourses(lin1);
e = max(abs(vals(10:end) - double(tc1)));
nfail = nfail + check('CSV voxel export values', e < 1e-5 && numel(vals) == 9 + T, e);
if exist('writecell', 'file') == 2
    fx = fullfile(out, 'ts.xlsx');
    gdce_export_timeseries(D, masks, {'SSS (control)', 'Tumor 1'}, fx, 'BaseFrames', 1:3);
    nfail = nfail + check('XLSX export', exist(fx, 'file') == 2, 0);
end

% ---------------------------------------------------------------- 8. GUI
if doGUI
    nfail = nfail + gdce_selftest_gui(D, ref, out);
end

ok = nfail == 0;
if ok
    fprintf('ALL TESTS PASSED\n');
else
    fprintf('%d TEST(S) FAILED\n', nfail);
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
