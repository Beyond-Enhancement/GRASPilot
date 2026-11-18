function tf = graspilot_is_dicom(file)
%GRASPILOT_IS_DICOM True if a file looks like DICOM, whatever its name or extension.
%
%   tf = graspilot_is_dicom(file)
%
%   Part 10 files: 'DICM' after the 128-byte preamble. Bare data sets
%   (no preamble, as written by some older exporters): the first element is
%   in group 0008 (little endian), with an explicit VR or a plausible
%   implicit length.
tf = false;
fid = fopen(file, 'r');
if fid < 0, return; end
b = fread(fid, 132, '*uint8')';
fclose(fid);
if numel(b) == 132 && isequal(b(129:132), uint8('DICM'))
    tf = true;
elseif numel(b) >= 8 && b(1) == 8 && b(2) == 0 && b(4) == 0
    vr = char(b(5:6));
    len = double(b(5)) + 256 * double(b(6)) + 65536 * double(b(7)) + 16777216 * double(b(8));
    tf = all(vr >= 'A' & vr <= 'Z') || len < 1024;
end
end
