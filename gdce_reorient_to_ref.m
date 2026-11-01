function [out, msg] = gdce_reorient_to_ref(img, Afile, Aref, refSize)
%GDCE_REORIENT_TO_REF Bring a NIfTI array onto the viewer's voxel grid.
%
%   [out, msg] = gdce_reorient_to_ref(img, Afile, Aref, refSize)
%
%   img      3D array read from a NIfTI file
%   Afile    its 4x4 affine (voxel 0-based -> RAS mm)
%   Aref     affine of the viewer grid (GraspData.affine)
%   refSize  size of the viewer grid
%
%   Handles axis permutations and flips (e.g. a mask re-saved by another
%   program in a different voxel order). Resampling is not performed: if
%   the grids are not related by a signed permutation, an error is thrown.
%   msg describes what was done (or a warning on a voxel-size mismatch).

msg = '';
T = Aref \ Afile;                       % file voxel -> reference voxel
R = T(1:3, 1:3);
[~, o] = max(abs(R), [], 1);            % file axis f -> reference axis o(f)
if numel(unique(o)) < 3
    error('gdce_reorient_to_ref:grid', 'The mask grid is oblique to the image grid; cannot reorient.');
end
lin = sub2ind([3 3], o, 1:3);
s = sign(R(lin));
mag = abs(R(lin));
off = abs(R);
off(lin) = 0;
if any(off(:) > 0.05)
    error('gdce_reorient_to_ref:grid', 'The mask grid is rotated relative to the image grid; cannot reorient.');
end
if any(abs(mag - 1) > 0.02)
    msg = 'Voxel sizes in the mask header differ from the current geometry; matched by index.';
end
P = zeros(1, 3);
P(o) = 1:3;                             % out dim a <- img dim P(a)
out = permute(img, P);
for a = 1:3
    if s(P(a)) < 0, out = flip(out, a); end
end
sz3 = @(x) [size(x, 1) size(x, 2) size(x, 3)];
if ~isequal(sz3(out), refSize(:)')
    if isequal(sz3(img), refSize(:)')
        out = img;
        msg = 'Header orientation could not be matched; the array was used as stored.';
    else
        error('gdce_reorient_to_ref:size', 'Mask size [%s] does not match the image [%s].', ...
            num2str(sz3(out)), num2str(refSize(:)'));
    end
elseif isempty(msg) && ~isequal(P, 1:3) || any(s < 0)
    if isempty(msg), msg = 'Mask was reoriented to the image voxel order.'; end
end
end
