function outFile = gdce_convert_mat2nifti(matFile, outFile, varargin)
%GDCE_CONVERT_MAT2NIFTI Batch-convert a GRASP .mat reconstruction to 4D NIfTI.
%
%   gdce_convert_mat2nifti('recon.mat')                 -> recon_dce.nii.gz
%   gdce_convert_mat2nifti('recon.mat', 'out.nii.gz', 'Name', Value, ...)
%
%   Options (in addition to all GraspData.fromFile options):
%   'AxisCodes'  direction of increasing rows, cols, slices ('LPS' =
%                scanner order of the GRASP reconstructions; 'PLS' = data1
%                slices already radiological; see GraspData)    ('LPS')
%   'Plane', 'SwapInPlane', 'Flip'  older alternative to AxisCodes
%   'Voxel'      [row col slice] mm                    ([1 1 1])
%   'dt'         s per frame                           (NaN)
%   't0'         s, time of frame 1                    (0)
%   'BinSize'    average n consecutive frames          (1)
%   'Frames'     frame indices to write, after binning (all)
%
%   Uses exactly the same frame synthesis, orientation and affine as the
%   viewer, so masks exported from GraspDCEViewer overlay this file. A JSON
%   sidecar with timing and provenance is written next to the image.
%
%   Example (loop over 8/16/32-spoke reconstructions):
%     f = dir('*.mat');
%     for k = 1:numel(f)
%         gdce_convert_mat2nifti(f(k).name, '', 'Voxel', [0.94 0.94 2], 'dt', 3.3);
%     end

if nargin < 2 || isempty(outFile)
    [p, b] = fileparts(matFile);
    outFile = fullfile(p, [b '_dce.nii.gz']);
end
g = struct('AxisCodes', 'LPS', 'Plane', '', 'SwapInPlane', false, 'Flip', [false false false], ...
    'Voxel', [1 1 1], 'dt', NaN, 't0', 0, 'Frames', [], 'BinSize', 1);
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
D = GraspData.fromFile(matFile, pass{:});
if ~isempty(g.Plane)
    D.setOrientation(g.Plane, g.SwapInPlane, g.Flip);
else
    D.setAxisCodes(g.AxisCodes);
end
D.VoxNative = g.Voxel;
D.dt = g.dt;
D.setBinning(g.BinSize);
D.t0 = g.t0;
fr = g.Frames;
if isempty(fr), fr = 1:D.NT; end
gdce_write_dce_nifti(D, outFile, fr, prog);
fprintf('Wrote %s\n', outFile);
end
