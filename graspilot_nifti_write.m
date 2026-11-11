function graspilot_nifti_write(fname, src, hdr, progress)
%GRASPILOT_NIFTI_WRITE Write a NIfTI-1 single-file volume (.nii or .nii.gz).
%
%   graspilot_nifti_write(fname, src, hdr)
%   graspilot_nifti_write(fname, src, hdr, progress)
%
%   src  numeric 3D/4D array, or a function handle f(t) that returns frame
%        t (3D array) for t = 1..nt. Frames are streamed to disk one at a
%        time, so a 4D series never has to exist in memory.
%   hdr  struct with fields
%        Dim       [n1 n2 n3] or [n1 n2 n3 nt]
%        Vox       [dx dy dz] mm
%        dt        seconds per frame (4D only; 0 if unknown)        [0]
%        Affine    4x4, voxel (0-based) -> RAS mm (sform and qform)  [diag(Vox)]
%        Datatype  'float32' | 'float64' | 'uint8' | 'int16' | 'uint16' | 'int32'
%        Descrip   text, up to 79 characters                         ['']
%        CalRange  [cal_min cal_max] display range                   [0 0]
%   progress  optional function handle progress(fraction, message)
%
%   No toolbox required. Little-endian output.

if nargin < 4, progress = []; end
if ~isfield(hdr, 'dt') || isempty(hdr.dt) || ~isfinite(hdr.dt), hdr.dt = 0; end
if ~isfield(hdr, 'Descrip'), hdr.Descrip = ''; end
if ~isfield(hdr, 'CalRange'), hdr.CalRange = [0 0]; end
if ~isfield(hdr, 'Affine') || isempty(hdr.Affine)
    hdr.Affine = diag([hdr.Vox(:)' 1]);
end

dim = hdr.Dim(:)';
if numel(dim) < 4, dim(4) = 1; end
nt = dim(4);
ndim = 3 + (nt > 1);
[code, bitpix, prec] = dtypeInfo(hdr.Datatype);

A = double(hdr.Affine);
[qb, qc, qd, qfac] = mat2quatern(A);

H = zeros(352, 1, 'uint8');
H = put(H, 0, 348, 'int32');                       % sizeof_hdr
H(39) = uint8('r');                                % regular
d = ones(1, 8); d(1) = ndim; d(2:4) = dim(1:3); d(5) = nt;
H = put(H, 40, d, 'int16');                        % dim
H = put(H, 70, code, 'int16');                     % datatype
H = put(H, 72, bitpix, 'int16');                   % bitpix
pix = zeros(1, 8); pix(1) = qfac; pix(2:4) = hdr.Vox(1:3); pix(5) = hdr.dt;
H = put(H, 76, pix, 'single');                     % pixdim
H = put(H, 108, 352, 'single');                    % vox_offset
H = put(H, 112, 1, 'single');                      % scl_slope
H = put(H, 116, 0, 'single');                      % scl_inter
H(124) = uint8(2 + 8);                             % xyzt_units: mm + s
H = put(H, 124, hdr.CalRange(2), 'single');        % cal_max
H = put(H, 128, hdr.CalRange(1), 'single');        % cal_min
ds = uint8(hdr.Descrip(1:min(79, numel(hdr.Descrip))));
H(149:148 + numel(ds)) = ds;                       % descrip
H = put(H, 252, 1, 'int16');                       % qform_code (scanner)
H = put(H, 254, 1, 'int16');                       % sform_code (scanner)
H = put(H, 256, [qb qc qd], 'single');             % quatern_b/c/d
H = put(H, 268, A(1:3, 4)', 'single');             % qoffset_x/y/z
H = put(H, 280, A(1, :), 'single');                % srow_x
H = put(H, 296, A(2, :), 'single');                % srow_y
H = put(H, 312, A(3, :), 'single');                % srow_z
H(345:348) = uint8([double('n+1') 0]);                     % magic

gz = numel(fname) > 3 && strcmpi(fname(end - 2:end), '.gz');
if gz
    folder = fileparts(fname);
    if isempty(folder), folder = pwd; end
    target = [tempname(folder) '.nii'];
else
    target = fname;
end

fid = fopen(target, 'w', 'ieee-le');
if fid < 0
    error('graspilot_nifti_write:open', 'Cannot open %s for writing.', target);
end
cleanup = onCleanup(@() closeIfOpen(fid));
fwrite(fid, H, 'uint8');
for t = 1:nt
    if isa(src, 'function_handle')
        vol = src(t);
    elseif nt > 1
        vol = src(:, :, :, t);
    else
        vol = src;
    end
    if numel(vol) ~= prod(dim(1:3))
        error('graspilot_nifti_write:size', 'Frame %d has %d voxels, expected %d.', t, numel(vol), prod(dim(1:3)));
    end
    fwrite(fid, cast(vol, prec), prec);
    if ~isempty(progress)
        try, progress(t / nt, sprintf('Writing frame %d/%d', t, nt)); catch, end
    end
end
fclose(fid);
fid = -1; %#ok<NASGU>
clear cleanup

if gz
    if ~isempty(progress)
        try, progress(1, 'Compressing...'); catch, end
    end
    gzip(target);
    if exist(fname, 'file') == 2, delete(fname); end
    movefile([target '.gz'], fname);
    delete(target);
end
end

% -------------------------------------------------------------------------
function H = put(H, offset, val, type)
b = typecast(cast(val(:)', type), 'uint8');
H(offset + 1:offset + numel(b)) = b;
end

function closeIfOpen(fid)
try
    if fid > 0 && ~isempty(fopen(fid)), fclose(fid); end
catch
end
end

function [code, bitpix, prec] = dtypeInfo(name)
switch lower(name)
    case 'uint8',   code = 2;   bitpix = 8;  prec = 'uint8';
    case 'int16',   code = 4;   bitpix = 16; prec = 'int16';
    case 'int32',   code = 8;   bitpix = 32; prec = 'int32';
    case 'float32', code = 16;  bitpix = 32; prec = 'single';
    case 'single',  code = 16;  bitpix = 32; prec = 'single';
    case 'float64', code = 64;  bitpix = 64; prec = 'double';
    case 'double',  code = 64;  bitpix = 64; prec = 'double';
    case 'uint16',  code = 512; bitpix = 16; prec = 'uint16';
    otherwise
        error('graspilot_nifti_write:type', 'Unsupported datatype %s.', name);
end
end

function [b, c, d, qfac] = mat2quatern(A)
% nifti_mat44_to_quatern (nifti1_io.c), rotation part only
R = A(1:3, 1:3);
v = sqrt(sum(R .^ 2, 1));
v(v == 0) = 1;
R = R ./ repmat(v, 3, 1);
if det(R) > 0
    qfac = 1;
else
    qfac = -1;
    R(:, 3) = -R(:, 3);
end
a = R(1, 1) + R(2, 2) + R(3, 3) + 1;
if a > 0.5
    a = 0.5 * sqrt(a);
    b = 0.25 * (R(3, 2) - R(2, 3)) / a;
    c = 0.25 * (R(1, 3) - R(3, 1)) / a;
    d = 0.25 * (R(2, 1) - R(1, 2)) / a;
else
    xd = 1 + R(1, 1) - (R(2, 2) + R(3, 3));
    yd = 1 + R(2, 2) - (R(1, 1) + R(3, 3));
    zd = 1 + R(3, 3) - (R(1, 1) + R(2, 2));
    if xd > 1
        b = 0.5 * sqrt(xd);
        c = 0.25 * (R(1, 2) + R(2, 1)) / b;
        d = 0.25 * (R(1, 3) + R(3, 1)) / b;
        a = 0.25 * (R(3, 2) - R(2, 3)) / b;
    elseif yd > 1
        c = 0.5 * sqrt(yd);
        b = 0.25 * (R(1, 2) + R(2, 1)) / c;
        d = 0.25 * (R(2, 3) + R(3, 2)) / c;
        a = 0.25 * (R(1, 3) - R(3, 1)) / c;
    else
        d = 0.5 * sqrt(zd);
        b = 0.25 * (R(1, 3) + R(3, 1)) / d;
        c = 0.25 * (R(2, 3) + R(3, 2)) / d;
        a = 0.25 * (R(2, 1) - R(1, 2)) / d;
    end
    if a < 0
        b = -b; c = -c; d = -d;
    end
end
end
