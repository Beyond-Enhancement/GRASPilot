function c = graspilot_folder_contents(folder)
%GRASPILOT_FOLDER_CONTENTS What a folder holds that GRASPilot can open.
%
%   c = graspilot_folder_contents(folder)
%
%   c.Dicom      true if the folder or one of its subfolders holds DICOM
%                files (any name, with or without extension; checked by
%                content, graspilot_is_dicom)
%   c.DicomFile  the first DICOM file found
%   c.Nifti      NIfTI files (.nii, .nii.gz) in the folder, sorted by name
%   c.Mat        .mat files in the folder, sorted by name
%   NIfTI and .mat files are looked for in the folder itself only.
if exist(folder, 'dir') ~= 7
    error('graspilot_folder_contents:dir', 'No such folder: %s', folder);
end
c = struct('Dicom', false, 'DicomFile', '', 'Nifti', {{}}, 'Mat', {{}});
d = dir(folder);
d = d(~cellfun(@(n) isempty(n) || n(1) == '.', {d.name}));
names = {d(~[d.isdir]).name};
low = lower(names);
isNii = ~cellfun(@isempty, regexp(low, '(\.nii|\.nii\.gz)$', 'once'));
isMat = ~cellfun(@isempty, regexp(low, '\.mat$', 'once'));
c.Nifti = sort(cellfun(@(n) fullfile(folder, n), names(isNii), 'UniformOutput', false));
c.Mat = sort(cellfun(@(n) fullfile(folder, n), names(isMat), 'UniformOutput', false));
c.DicomFile = findDicom(folder, 0);
c.Dicom = ~isempty(c.DicomFile);
end

function f = findDicom(folder, depth)
% first DICOM file in the folder, else in its subfolders (depth-first)
f = '';
skip = {'.nii', '.gz', '.mat', '.json', '.txt', '.csv', '.xlsx', '.xls', '.png', '.jpg', '.jpeg', '.pdf', ...
    '.xml', '.zip', '.m', '.md', '.bval', '.bvec', '.html', '.log', '.tar', '.docx', '.pptx'};
d = dir(folder);
d = d(~cellfun(@(n) isempty(n) || n(1) == '.', {d.name}));
files = d(~[d.isdir]);
nchk = 0;
for k = 1:numel(files)
    [~, ~, e] = fileparts(files(k).name);
    if any(strcmpi(e, skip)) || strcmpi(files(k).name, 'DICOMDIR') || files(k).bytes < 136, continue; end
    full = fullfile(folder, files(k).name);
    if graspilot_is_dicom(full), f = full; return; end
    nchk = nchk + 1;
    if nchk >= 20, break; end                % not DICOM: stop checking this folder
end
if depth >= 8, return; end
sub = d([d.isdir]);
for k = 1:numel(sub)
    f = findDicom(fullfile(folder, sub(k).name), depth + 1);
    if ~isempty(f), return; end
end
end
