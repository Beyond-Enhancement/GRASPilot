function r = graspilot_validate_readdata(matFile, slices, frames, coefVar, basisVar)
%GRASPILOT_VALIDATE_READDATA Check the viewer's frames against ReadData.m on a real file.
%
%   graspilot_validate_readdata('recon.mat')                     3 slices, all frames
%   graspilot_validate_readdata('recon.mat', [10 40 70], 1:10:150)
%   r = graspilot_validate_readdata(..., coefVar, basisVar)      other variable names
%
%   Runs the original ReadData.m arithmetic (TempPCASub adjoint, global
%   maximum over all slices and frames before CropImg) slice by slice, so
%   the full 4D array is never built, and compares it with
%   GRASPilotData.nativeSlice for the chosen slices/frames. Differences should
%   be at single-precision rounding level (~1e-7 relative to the maximum).

if nargin < 4 || isempty(coefVar), coefVar = 'vk_all'; end
if nargin < 5 || isempty(basisVar), basisVar = 'basis_all'; end
here = fileparts(mfilename('fullpath'));
root = fileparts(here);
addpath(root);
addpath(fullfile(root, 'legacy'));

fprintf('Loading %s ...\n', matFile);
S = load(matFile, coefVar, basisVar);
vk_all = S.(coefVar);
basis_all = S.(basisVar);
S = []; %#ok<NASGU>
nz = size(vk_all, 1);
if nargin < 2 || isempty(slices)
    slices = unique(round(linspace(1, nz, 5)));
    slices = slices(2:end - 1);
end

% reference: ReadData.m, one slice at a time
nxnew = size(vk_all, 2) / 2;
gmax = 0;
ref = cell(1, numel(slices));
for ii = 1:nz
    Phi = single(squeeze(squeeze(basis_all(ii, :, :))));
    PCA = TempPCASub(Phi);
    Recon = permute(vk_all(ii, :, :, :), [2 3 4 1]);   % == Recon_Allslice(:,:,:,ii)
    d = abs(PCA' * Recon);
    gmax = max(gmax, max(d(:)));
    k = find(slices == ii);
    if ~isempty(k), ref{k} = CropImg(d, nxnew, nxnew); end
end

fprintf('Building GRASPilotData ...\n');
D = GRASPilotData.fromSubspace(vk_all, basis_all);
if nargin < 3 || isempty(frames), frames = 1:D.NT; end
emax = 0;
for k = 1:numel(slices)
    R = ref{k} / gmax;
    for t = frames
        a = D.nativeSlice(3, slices(k), t);
        emax = max(emax, max(abs(a(:) - reshape(R(:, :, t), [], 1))));
    end
end
r = struct('maxAbsDiff', emax, 'slices', slices, 'frames', frames, 'globalMax', gmax, ...
    'size', D.NativeSize, 'K', D.NK, 'T', D.NT);
fprintf('Slices %s, %d frames: max |viewer - ReadData.m| = %.3g (data normalized to max 1)\n', ...
    mat2str(slices), numel(frames), emax);
if emax < 1e-4
    fprintf('OK: identical to ReadData.m within single-precision rounding.\n');
else
    fprintf('WARNING: difference larger than expected.\n');
end
end
