function graspilot_dicom_write(folder, V, g)
%GRASPILOT_DICOM_WRITE Minimal classic MR DICOM writer for the self-tests.
%
%   graspilot_dicom_write(folder, V, g)
%
%   V  uint16 [nx ny nz T] (dim 1 along the image rows / increasing column,
%      dim 2 down the columns), one file per slice and frame, written in a
%      scrambled order with unrelated file names.
%   g  struct: X, Y (row and column direction cosines, LPS), ipp0 (position
%      of the first voxel of slice 1, LPS mm), ps [row col] spacing, dz
%      (slice spacing along cross(X, Y)), dt (s between frames), Implicit
%      (true: implicit VR little endian), SeriesNumber, Description.
%   Test data only: not a complete DICOM object.
if exist(folder, 'dir') ~= 7, mkdir(folder); end
[nx, ny, nz, T] = size(V);
N = cross(g.X, g.Y);
uid = sprintf('1.2.826.0.1.3680043.8.498.%d', 1000 + g.SeriesNumber);
order = reshape(reshape(1:nz * T, nz, T)', 1, []);       % frame-interleaved
order = order(mod((0:numel(order) - 1) * 7, numel(order)) + 1);   % scrambled (7 coprime to count)
if numel(unique(order)) ~= nz * T, order = 1:nz * T; end
for c = 1:numel(order)
    [z, t] = ind2sub([nz T], order(c));
    sec = 36300 + g.dt * (t - 1);
    tm = sprintf('%02d%02d%09.6f', floor(sec / 3600), floor(mod(sec, 3600) / 60), mod(sec, 60));
    ipp = g.ipp0(:)' + (z - 1) * g.dz * N(:)';
    E = { ...
        '0008', '0032', 'TM', tm; ...
        '0008', '0060', 'CS', 'MR'; ...
        '0008', '103E', 'LO', g.Description; ...
        '0018', '0050', 'DS', sprintf('%g', g.dz); ...
        '0020', '000E', 'UI', uid; ...
        '0020', '0011', 'IS', sprintf('%d', g.SeriesNumber); ...
        '0020', '0013', 'IS', sprintf('%d', c); ...
        '0020', '0032', 'DS', sprintf('%.6f\\%.6f\\%.6f', ipp); ...
        '0020', '0037', 'DS', sprintf('%.8f\\%.8f\\%.8f\\%.8f\\%.8f\\%.8f', [g.X(:)' g.Y(:)']); ...
        '0028', '0002', 'US', 1; ...
        '0028', '0010', 'US', ny; ...
        '0028', '0011', 'US', nx; ...
        '0028', '0030', 'DS', sprintf('%.6f\\%.6f', g.ps); ...
        '0028', '0100', 'US', 16; ...
        '0028', '0103', 'US', 0};
    fid = fopen(fullfile(folder, sprintf('I%05d', mod(c * 7919, 99991))), 'w', 'ieee-le');
    fwrite(fid, zeros(1, 128, 'uint8'));
    fwrite(fid, 'DICM');
    if g.Implicit, ts = '1.2.840.10008.1.2'; else, ts = '1.2.840.10008.1.2.1'; end
    tsv = pad(ts, 'UI');
    fwrite(fid, [2 0 0 0], 'uint8'); fwrite(fid, 'UL'); fwrite(fid, 4, 'uint16'); fwrite(fid, 8 + numel(tsv), 'uint32');
    fwrite(fid, [2 0 16 0], 'uint8'); fwrite(fid, 'UI'); fwrite(fid, numel(tsv), 'uint16'); fwrite(fid, tsv);
    for k = 1:size(E, 1)
        fwrite(fid, hex2dec(E{k, 1}), 'uint16'); fwrite(fid, hex2dec(E{k, 2}), 'uint16');
        if strcmp(E{k, 3}, 'US')
            val = typecast(uint16(E{k, 4}), 'uint8');
        else
            val = pad(E{k, 4}, E{k, 3});
        end
        if g.Implicit
            fwrite(fid, numel(val), 'uint32');
        else
            fwrite(fid, E{k, 3}); fwrite(fid, numel(val), 'uint16');
        end
        fwrite(fid, val, 'uint8');
    end
    px = V(:, :, z, t);
    fwrite(fid, hex2dec('7FE0'), 'uint16'); fwrite(fid, hex2dec('0010'), 'uint16');
    if ~g.Implicit, fwrite(fid, 'OW'); fwrite(fid, 0, 'uint16'); end
    fwrite(fid, 2 * numel(px), 'uint32');
    fwrite(fid, px(:), 'uint16');
    fclose(fid);
end
end

function v = pad(s, vr)
v = uint8(s);
if mod(numel(v), 2) == 1
    if strcmp(vr, 'UI'), v(end + 1) = 0; else, v(end + 1) = 32; end
end
end
