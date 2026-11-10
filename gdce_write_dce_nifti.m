function gdce_write_dce_nifti(D, outFile, frames, progress)
%GDCE_WRITE_DCE_NIFTI Write frames of a GraspData series to NIfTI + JSON sidecar.
%
%   gdce_write_dce_nifti(D, outFile)                all frames, 4D float32
%   gdce_write_dce_nifti(D, outFile, frames, progress)

if nargin < 3 || isempty(frames), frames = 1:D.NT; end
if nargin < 4, progress = []; end
n = D.canonSize();
hdr.Dim = [n numel(frames)];
hdr.Vox = D.canonVox();
ft = D.frameTimes();
if D.hasTiming() && numel(frames) > 1
    hdr.dt = median(diff(ft(frames)));
else
    hdr.dt = 0;
end
hdr.Affine = D.affine();
hdr.Datatype = 'float32';
hdr.Descrip = 'GRASP DCE Viewer';
gdce_nifti_write(outFile, @(i) D.canonVolume(frames(i)), hdr, progress);

% JSON sidecar
J = struct();
J.CreatedBy = 'GRASP DCE Viewer';
J.Created = datestr(now, 'yyyy-mm-ddTHH:MM:SS');
J.SourceFile = D.SourceFile;
J.SourceVariables = D.SourceVars;
J.DataMode = D.Mode;
J.Frames = frames;
if D.hasTiming()
    J.FrameTimes_s = ft(frames);
    J.TemporalResolution_s = D.effectiveDt();
    J.RawTemporalResolution_s = D.dt;
end
J.RawFrames = D.NTRaw;
J.BinSize = D.BinSize;
J.BinningNote = 'each frame = mean magnitude of BinSize consecutive reconstructed frames';
J.RawFramesPerFrame = D.BinIdx(frames);
J.DroppedRawFrames = D.droppedFrames();
J.VoxelSize_mm = D.canonVox();
J.NativeSize = D.NativeSize;
J.NativeAxisCodes = D.AxisCodes;
J.NativePlane = D.Plane;
J.SwapInPlane = D.SwapInPlane;
J.UserFlip_LR_AP_SI = D.UserFlip;
J.CropFactor = D.CropFactor;
J.NormalizedToGlobalMax = D.Normalized;
J.ScaleFactor = D.ScaleFactor;
J.VoxelOrder = 'LPS (i: R->L, j: A->P, k: I->S); affine maps to RAS mm, volume centred at 0';
jf = regexprep(outFile, '(\.nii\.gz|\.nii)$', '.json', 'ignorecase');
if strcmp(jf, outFile), jf = [outFile '.json']; end
try
    fn = fieldnames(J);
    for i = 1:numel(fn)
        if isa(J.(fn{i}), 'single'), J.(fn{i}) = double(J.(fn{i})); end
    end
    txt = jsonencode(J);
    fid = fopen(jf, 'w');
    fprintf(fid, '%s', txt);
    fclose(fid);
catch err
    warning('gdce:json', 'JSON sidecar not written: %s', err.message);
end
end
