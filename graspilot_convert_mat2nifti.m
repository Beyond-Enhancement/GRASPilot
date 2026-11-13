function outFile = graspilot_convert_mat2nifti(matFile, outFile, varargin)
%GRASPILOT_CONVERT_MAT2NIFTI Batch-convert a GRASP .mat reconstruction to 4D NIfTI.
%
%   graspilot_convert_mat2nifti('recon.mat')                 -> recon_dce.nii.gz
%   graspilot_convert_mat2nifti('recon.mat', 'out.nii.gz', 'Name', Value, ...)
%
%   Options (in addition to all GRASPilotData.fromFile options):
%   'AxisCodes'  direction of increasing rows, cols, slices ('LPS' =
%                scanner order of the GRASP reconstructions; 'PLS' = data1
%                slices already radiological; see GRASPilotData)    ('LPS')
%   'Plane', 'SwapInPlane', 'Flip'  older alternative to AxisCodes
%   'Voxel'      [row col slice] mm                    ([1 1 1])
%   'dt'         s per frame                           (NaN)
%   't0'         s, time of frame 1                    (0)
%   'BinSize'    average n consecutive frames          (1)
%   'Frames'     frame indices to write, after binning (all)
%   'MotionCorrect'  estimate rigid head motion and realign the raw frames
%                before binning and writing (false); the motion summary is
%                printed and stored in the JSON sidecar, the parameters in
%                <out>_motion.csv
%   'MotionCheck' estimate motion and warn on excess motion without
%                correcting (false)
%   'BaseFrames' raw frames averaged as the motion reference (1:5)
%
%   Uses exactly the same frame synthesis, orientation and affine as the
%   viewer, so masks exported from GRASPilot overlay this file. A JSON
%   sidecar with timing and provenance is written next to the image.
%
%   Example (loop over 8/16/32-spoke reconstructions):
%     f = dir('*.mat');
%     for k = 1:numel(f)
%         graspilot_convert_mat2nifti(f(k).name, '', 'Voxel', [0.94 0.94 2], 'dt', 3.3);
%     end

if nargin < 2 || isempty(outFile)
    [p, b] = fileparts(matFile);
    outFile = fullfile(p, [b '_dce.nii.gz']);
end
g = struct('AxisCodes', 'LPS', 'Plane', '', 'SwapInPlane', false, 'Flip', [false false false], ...
    'Voxel', [1 1 1], 'dt', NaN, 't0', 0, 'Frames', [], 'BinSize', 1, 'MotionCorrect', false, ...
    'MotionCheck', false, 'BaseFrames', []);
pass = {};
for i = 1:2:numel(varargin)
    fn = fieldnames(g);
    m = find(strcmpi(fn, varargin{i}), 1);
    if isempty(m)
        pass(end + 1:end + 2) = varargin(i:i + 1); %#ok<AGROW>
    else
        g.(fn{m}) = varargin{i + 1};
    end
end
prog = @(f, msg) fprintf('  %5.1f%%  %s\n', 100 * f, msg);
fprintf('Loading %s\n', matFile);
D = GRASPilotData.fromFile(matFile, pass{:});
if ~isempty(g.Plane)
    D.setOrientation(g.Plane, g.SwapInPlane, g.Flip);
else
    D.setAxisCodes(g.AxisCodes);
end
D.VoxNative = g.Voxel;
D.dt = g.dt;
D.setBinning(g.BinSize);
D.t0 = g.t0;
if g.MotionCorrect || g.MotionCheck
    M = graspilot_motion_estimate(D, 'RefFrames', g.BaseFrames, 'Progress', prog);
    D.setMotion(M);
    fprintf('%s\n', M.Message);
    if M.Excess, warning('graspilot:motion', 'Excess motion in %s', matFile); end
    if g.MotionCorrect
        D.applyMotionCorrection(M, prog);
    end
    mf = regexprep(outFile, '(\.nii\.gz|\.nii)$', '_motion.csv', 'ignorecase');
    writeMotionCsv(mf, D, M);
end
fr = g.Frames;
if isempty(fr), fr = 1:D.NT; end
graspilot_write_dce_nifti(D, outFile, fr, prog);
fprintf('Wrote %s\n', outFile);
end

function writeMotionCsv(file, D, M)
fid = fopen(file, 'w');
if fid < 0, warning('graspilot:motion', 'Cannot write %s', file); return; end
fprintf(fid, 'RawFrame,TransL_mm,TransP_mm,TransS_mm,RotLR_deg,RotAP_deg,RotIS_deg,MeanDispRef_mm,MeanDispPrev_mm,Excess\n');
for f = 1:size(M.Params, 1)
    fprintf(fid, '%d,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%d\n', f, M.Params(f, :), M.AbsMean(f), ...
        M.RelMean(f), any(M.ExcessFrames == f));
end
fprintf(fid, '# %s\n# Motion correction: %s\n', M.Message, D.motionDescription());
fclose(fid);
end
