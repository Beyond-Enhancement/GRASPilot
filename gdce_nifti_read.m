function [img, hdr] = gdce_nifti_read(fname)
%GDCE_NIFTI_READ Read a NIfTI-1 single-file volume (.nii or .nii.gz).
%
%   [img, hdr] = gdce_nifti_read(fname)
%
%   hdr fields: Dim, Vox, dt, Affine (voxel 0-based -> RAS mm, from sform
%   if set, else qform, else pixdim), Datatype, SformCode, QformCode.
%   Scaling (scl_slope/scl_inter) is applied when present. No toolbox needed.

tmpDir = '';
if numel(fname) > 3 && strcmpi(fname(end - 2:end), '.gz')
    tmpDir = tempname;
    mkdir(tmpDir);
    files = gunzip(fname, tmpDir);
    f = files{1};
else
    f = fname;
end
cleaner = onCleanup(@() removeTmp(tmpDir));

fid = fopen(f, 'r', 'ieee-le');
if fid < 0, error('gdce_nifti_read:open', 'Cannot open %s.', fname); end
sz = fread(fid, 1, 'int32');
if sz ~= 348
    fclose(fid);
    fid = fopen(f, 'r', 'ieee-be');
    sz = fread(fid, 1, 'int32');
    if sz ~= 348
        fclose(fid);
        error('gdce_nifti_read:format', '%s is not a NIfTI-1 file.', fname);
    end
end
fseek(fid, 40, 'bof');  dim = fread(fid, 8, 'int16')';
fseek(fid, 70, 'bof');  dtype = fread(fid, 1, 'int16');
fseek(fid, 76, 'bof');  pixdim = fread(fid, 8, 'float32')';
voxOffset = fread(fid, 1, 'float32');
slope = fread(fid, 1, 'float32');
inter = fread(fid, 1, 'float32');
fseek(fid, 252, 'bof');
qcode = fread(fid, 1, 'int16');
scode = fread(fid, 1, 'int16');
quat = fread(fid, 3, 'float32')';
qoff = fread(fid, 3, 'float32')';
srow = reshape(fread(fid, 12, 'float32'), 4, 3)';
fseek(fid, 344, 'bof');
magic = fread(fid, 4, 'uint8=>char')';
if ~strncmp(magic, 'n+1', 3)
    fclose(fid);
    error('gdce_nifti_read:format', 'Only single-file NIfTI-1 (.nii) is supported.');
end

nd = max(1, min(7, dim(1)));
dims = dim(2:1 + nd);
prec = precOf(dtype);
fseek(fid, round(voxOffset), 'bof');
img = fread(fid, prod(dims), [prec '=>' prec]);
fclose(fid);
if numel(img) ~= prod(dims)
    error('gdce_nifti_read:size', 'File %s is truncated.', fname);
end
img = reshape(img, [dims 1 1]);
if slope ~= 0 && ~(slope == 1 && inter == 0)
    img = double(img) * slope + inter;
end

hdr.Dim = dims;
hdr.Vox = pixdim(2:4);
hdr.dt = 0;
if nd >= 4, hdr.dt = pixdim(5); end
hdr.Datatype = prec;
hdr.SformCode = scode;
hdr.QformCode = qcode;
if scode > 0
    hdr.Affine = [srow; 0 0 0 1];
elseif qcode > 0
    hdr.Affine = quatern2mat(quat, qoff, pixdim);
else
    hdr.Affine = diag([pixdim(2:4) 1]);
end
end

% -------------------------------------------------------------------------
function p = precOf(code)
switch code
    case 2,   p = 'uint8';
    case 4,   p = 'int16';
    case 8,   p = 'int32';
    case 16,  p = 'single';
    case 64,  p = 'double';
    case 256, p = 'int8';
    case 512, p = 'uint16';
    case 768, p = 'uint32';
    otherwise
        error('gdce_nifti_read:type', 'Unsupported NIfTI datatype %d.', code);
end
end

function A = quatern2mat(q, off, pixdim)
b = q(1); c = q(2); d = q(3);
a = 1 - (b * b + c * c + d * d);
if a < 1e-7
    a = 1 / sqrt(b * b + c * c + d * d);
    b = b * a; c = c * a; d = d * a;
    a = 0;
else
    a = sqrt(a);
end
qfac = pixdim(1);
if qfac == 0, qfac = 1; end
R = [a*a+b*b-c*c-d*d, 2*b*c-2*a*d,     2*b*d+2*a*c;
     2*b*c+2*a*d,     a*a+c*c-b*b-d*d, 2*c*d-2*a*b;
     2*b*d-2*a*c,     2*c*d+2*a*b,     a*a+d*d-c*c-b*b];
S = diag([pixdim(2) pixdim(3) pixdim(4) * qfac]);
A = [R * S off(:); 0 0 0 1];
end

function removeTmp(d)
if ~isempty(d) && exist(d, 'dir') == 7
    try, rmdir(d, 's'); catch, end
end
end
