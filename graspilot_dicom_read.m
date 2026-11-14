function [V, info] = graspilot_dicom_read(src, varargin)
%GRASPILOT_DICOM_READ Read a dynamic (4D) or 3D MR series from DICOM files.
%
%   [V, info] = graspilot_dicom_read(folder)
%   [V, info] = graspilot_dicom_read(folder, 'Series', k)
%   scan      = graspilot_dicom_read(folder, 'ScanOnly', true)
%   [V, info] = graspilot_dicom_read(scan, 'Series', k)      reuse a scan
%
%   src       folder (searched recursively), one DICOM file (its folder is
%             scanned and the series of that file is read), a cell array of
%             files, or a scan struct from a previous call
%   V         single [nx ny nz T]: dim 1 along the image rows (increasing
%             column index), dim 2 down the columns, dim 3 slices sorted
%             along the slice normal, dim 4 time frames. Rescale slope and
%             intercept are applied.
%   info      Affine (4x4, voxel 0-based -> RAS mm, for V as returned), Vox
%             (mm), dt (s between frames; NaN if the headers give no usable
%             frame times), Times (s from the first frame), Series
%             (description), SeriesNumber, Modality, NFiles, Folder, Message
%
%   Options
%     'Series'    index into scan.Series or a SeriesInstanceUID. Several
%                 indices: the series are concatenated as time frames (e.g.
%                 one 3D series per frame); they must share the same grid.
%                 Default: the series with the most images.
%     'ScanOnly'  true: only list the series (scan.Series(k).Label etc.)
%     'Recursive' search subfolders (true)
%     'Progress'  @(fraction, message)
%
%   Supported without any toolbox: classic single-frame and enhanced
%   multi-frame MR images, 8/16/32-bit integer pixels, uncompressed transfer
%   syntaxes (implicit/explicit VR little endian, explicit VR big endian).
%   Compressed pixel data (JPEG, JPEG 2000, RLE) are read through dicomread
%   when the Image Processing Toolbox is installed; otherwise convert the
%   series first (e.g. dcm2niix to NIfTI, which GRASPilot also opens).
%   Frames are ordered by temporal position, acquisition time, trigger time,
%   acquisition number, content time or instance number (the first of these
%   that separates the frames of every slice). No patient-identifying
%   elements are read.

p = struct('Series', [], 'ScanOnly', false, 'Recursive', true, 'Progress', []);
for i = 1:2:numel(varargin)
    fn = fieldnames(p);
    m = find(strcmpi(fn, varargin{i}), 1);
    if isempty(m), error('graspilot_dicom_read:opt', 'Unknown option %s.', varargin{i}); end
    p.(fn{m}) = varargin{i + 1};
end
report = @(f, msg) 0;
if ~isempty(p.Progress), report = p.Progress; end

pickFile = '';
if isstruct(src) && isfield(src, 'Rec')
    scan = src;
else
    if iscell(src)
        files = src(:)';
        folder = fileparts(files{1});
    elseif exist(src, 'dir') == 7
        folder = src;
        files = listFiles(folder, p.Recursive);
    elseif exist(src, 'file') == 2
        folder = fileparts(src);
        if isempty(folder), folder = pwd; end
        files = listFiles(folder, false);
        pickFile = src;
    else
        error('graspilot_dicom_read:src', 'No such file or folder: %s', src);
    end
    scan = scanFiles(files, folder, report);
end
if isempty(scan.Series)
    error('graspilot_dicom_read:none', 'No readable DICOM images found in %s.', scan.Folder);
end
if p.ScanOnly
    V = scan; info = [];
    return
end

% choose the series
sel = p.Series;
if isempty(sel)
    if ~isempty(pickFile)
        [~, a] = fileparts(pickFile);
        fi = find(strcmp(scan.Files, pickFile), 1);
        if isempty(fi)
            for k = 1:numel(scan.Files)
                [~, b] = fileparts(scan.Files{k});
                if strcmp(a, b), fi = k; break; end
            end
        end
        if ~isempty(fi)
            r = find(scan.Rec.File == fi, 1);
            sel = scan.Rec.Series(r);
        end
    end
    if isempty(sel), [~, sel] = max([scan.Series.NImages]); end
elseif ischar(sel)
    k = find(strcmp({scan.Series.UID}, sel), 1);
    if isempty(k), error('graspilot_dicom_read:series', 'Series %s not found.', sel); end
    sel = k;
end
if any(sel < 1 | sel > numel(scan.Series))
    error('graspilot_dicom_read:series', 'Series index out of range (1-%d).', numel(scan.Series));
end

V = []; info = [];
tStart = [];
for q = 1:numel(sel)
    [Vq, iq] = readSeries(scan, sel(q), @(f, msg) report((q - 1 + f) / numel(sel), msg));
    if q == 1
        V = Vq; info = iq;
    else
        if ~isequal(size(Vq(:, :, :, 1)), size(V(:, :, :, 1))) || max(abs(iq.Affine(:) - info.Affine(:))) > 1e-2
            error('graspilot_dicom_read:grid', 'Series %d and %d do not share the same image grid.', sel(1), sel(q));
        end
        V = cat(4, V, Vq);
    end
    tStart(end + 1) = iq.StartTime; %#ok<AGROW>
end
if numel(sel) > 1
    info.Series = sprintf('%d series combined (%s ...)', numel(sel), info.Series);
    info.dt = NaN; info.Times = [];
    if all(isfinite(tStart)) && all(diff(tStart) > 0) && size(V, 4) == numel(sel)
        info.Times = tStart - tStart(1);
        info.dt = median(diff(tStart));
    end
    info.Message = sprintf('%d series combined as %d frames.', numel(sel), size(V, 4));
end
report(1, 'DICOM: done');
end

% =========================================================================
function files = listFiles(folder, recursive)
files = {};
d = dir(folder);
skip = {'.nii', '.gz', '.mat', '.json', '.txt', '.csv', '.xlsx', '.xls', '.png', '.jpg', '.jpeg', '.pdf', ...
    '.xml', '.zip', '.m', '.md', '.bval', '.bvec', '.html', '.log', '.tar', '.docx', '.pptx'};
for k = 1:numel(d)
    nm = d(k).name;
    if isempty(nm) || nm(1) == '.', continue; end
    full = fullfile(folder, nm);
    if d(k).isdir
        if recursive, files = [files listFiles(full, true)]; end %#ok<AGROW>
        continue
    end
    [~, ~, e] = fileparts(nm);
    if any(strcmpi(e, skip)) || strcmpi(nm, 'DICOMDIR'), continue; end
    files{end + 1} = full; %#ok<AGROW>
end
end

% =========================================================================
function scan = scanFiles(files, folder, report)
D = tagTable();
n = numel(files);
H = cell(1, n);
ok = false(1, n);
for k = 1:n
    if mod(k, 50) == 1
        report(0.6 * (k - 1) / max(1, n), sprintf('DICOM: reading headers (%d of %d files)...', k, n));
    end
    try
        h = readHeader(files{k}, D);
    catch
        h = [];
    end
    if isempty(h) || ~isfield(h, 'PixelOffset') || ~isfield(h, 'Rows'), continue; end
    H{k} = h;
    ok(k) = true;
end
idx = find(ok);
% one record per image (frame of a multi-frame file)
R = struct('File', [], 'Frame', [], 'Series', [], 'IPP', zeros(0, 3), 'IOP', zeros(0, 6), 'PS', zeros(0, 2), ...
    'Thick', [], 'TPos', [], 'AcqT', [], 'ContT', [], 'Trig', [], 'AcqNum', [], 'Inst', [], 'Slope', [], 'Inter', []);
uids = {};
S = struct('UID', {}, 'Number', {}, 'Description', {}, 'Modality', {}, 'NImages', {}, 'Rows', {}, ...
    'Columns', {}, 'NSlices', {}, 'NFrames', {}, 'Label', {});
FH = struct('Rows', {}, 'Cols', {}, 'Bits', {}, 'Signed', {}, 'BE', {}, 'Encaps', {}, 'Offset', {}, 'NF', {}, ...
    'Spp', {});
for k = idx
    h = H{k};
    uid = getf(h, 'SeriesUID', '(no series UID)');
    s = find(strcmp(uids, uid), 1);
    if isempty(s)
        uids{end + 1} = uid; %#ok<AGROW>
        s = numel(uids);
        S(s).UID = uid;
        S(s).Number = getf(h, 'SeriesNumber', NaN);
        S(s).Description = getf(h, 'SeriesDescription', '');
        S(s).Modality = getf(h, 'Modality', '');
        S(s).NImages = 0;
        S(s).Rows = h.Rows; S(s).Columns = h.Cols;
    end
    nf = max(1, round(getf(h, 'NumberOfFrames', 1)));
    FH(k).Rows = h.Rows; FH(k).Cols = h.Cols;
    FH(k).Bits = getf(h, 'BitsAllocated', 16);
    FH(k).Signed = getf(h, 'PixelRepresentation', 0) == 1;
    FH(k).BE = h.BE; FH(k).Encaps = h.PixelLength == 4294967295;
    FH(k).Offset = h.PixelOffset; FH(k).NF = nf;
    FH(k).Spp = getf(h, 'SamplesPerPixel', 1);
    sh = struct();
    if isfield(h, 'Shared'), sh = h.Shared; end
    for f = 1:nf
        fr = struct();
        if isfield(h, 'Frames') && numel(h.Frames) >= f, fr = h.Frames{f}; end
        g = @(name, def) pick3(fr, sh, h, name, def);
        R.File(end + 1) = k;
        R.Frame(end + 1) = f;
        R.Series(end + 1) = s;
        R.IPP(end + 1, :) = vec(g('IPP', [NaN NaN NaN]), 3);
        R.IOP(end + 1, :) = vec(g('IOP', [1 0 0 0 1 0]), 6);
        R.PS(end + 1, :) = vec(g('PixelSpacing', [1 1]), 2);
        R.Thick(end + 1) = first(g('SpacingBetweenSlices', g('SliceThickness', 1)));
        R.TPos(end + 1) = first(g('TemporalPositionIndex', g('TemporalPosition', NaN)));
        R.AcqT(end + 1) = first(g('FrameAcqDateTime', g('AcqDateTime', g('AcqTime', NaN))));
        R.ContT(end + 1) = first(g('ContentTime', NaN));
        R.Trig(end + 1) = first(g('TriggerTime', NaN)) / 1000;
        R.AcqNum(end + 1) = first(g('AcquisitionNumber', NaN));
        R.Inst(end + 1) = first(g('InstanceNumber', NaN)) + (f - 1) / (nf + 1);
        R.Slope(end + 1) = first(g('RescaleSlope', 1));
        R.Inter(end + 1) = first(g('RescaleIntercept', 0));
    end
    S(s).NImages = S(s).NImages + nf;
end
% slices x frames per series
for s = 1:numel(S)
    r = find(R.Series == s);
    [nz, nt] = layoutOf(R, r);
    S(s).NSlices = nz; S(s).NFrames = nt;
    ds = S(s).Description;
    if isempty(ds), ds = '(no description)'; end
    if isnan(nt)
        lay = sprintf('%d images', numel(r));
    else
        lay = sprintf('%d slices x %d frames', nz, nt);
    end
    S(s).Label = sprintf('Series %g: %s - %s, %dx%d', S(s).Number, ds, lay, S(s).Columns, S(s).Rows);
end
[~, o] = sort([S.Number]);
if ~isempty(o)
    S = S(o);
    newIdx = zeros(1, numel(o));
    newIdx(o) = 1:numel(o);
    R.Series = newIdx(R.Series);
end
scan = struct('Folder', folder, 'Files', {files}, 'Series', S, 'Rec', R, 'FH', FH);
end

function [nz, nt, zi, pos] = layoutOf(R, r)
% slice index of every record of a series; nt = NaN if slices have unequal counts
iop = R.IOP(r(1), :);
nrm = cross(iop(1:3), iop(4:6));
pos = R.IPP(r, :) * nrm(:);
pos(~isfinite(pos)) = 0;
[sp, o] = sort(pos);
tol = max(1e-3, 0.1 * min([R.Thick(r(1)), 1e6]));
grp = cumsum([1; diff(sp) > tol]);
zi = zeros(numel(r), 1);
zi(o) = grp;
nz = max(grp);
cnt = accumarray(zi, 1);
if all(cnt == cnt(1)), nt = cnt(1); else, nt = NaN; end
end

% =========================================================================
function [V, info] = readSeries(scan, s, report)
R = scan.Rec;
r = find(R.Series == s);
S = scan.Series(s);
[nz, nt, zi] = layoutOf(R, r);
if isnan(nt)
    error('graspilot_dicom_read:layout', ['Series %g (%s): the slices do not all have the same number of ' ...
        'images (%d images, %d slice positions). Is the series complete?'], S.Number, S.Description, numel(r), nz);
end
% frame order: first key that separates the frames of every slice
keys = {R.TPos(r), R.AcqT(r), R.Trig(r), R.AcqNum(r), R.ContT(r), R.Inst(r), 1:numel(r)};
names = {'temporal position', 'acquisition time', 'trigger time', 'acquisition number', 'content time', ...
    'instance number', 'file order'};
ti = zeros(numel(r), 1);
used = '';
for c = 1:numel(keys)
    kv = keys{c}(:);
    if any(~isfinite(kv)), continue; end
    good = true;
    for z = 1:nz
        m = find(zi == z);
        [u, ~, j] = unique(kv(m));
        if numel(u) ~= nt, good = false; break; end
        ti(m) = j;
    end
    if good, used = names{c}; break; end
end
if isempty(used)
    error('graspilot_dicom_read:order', 'Series %g: cannot order the images into time frames.', S.Number);
end
% frame times (slice 1)
m1 = find(zi == 1);
[~, o] = sort(ti(m1));
m1 = m1(o);
dt = NaN; times = []; t0 = NaN;
cand = {R.AcqT(r), R.Trig(r), R.ContT(r)};
for c = 1:numel(cand)
    tv = cand{c}(m1);
    if all(isfinite(tv)), t0 = tv(1); end
    if nt > 1 && all(isfinite(tv)) && all(diff(tv) > 0)
        times = tv(:)' - tv(1);
        dt = median(diff(tv));
        break
    end
end
msg = sprintf('Frames ordered by %s.', used);
if isnan(dt) && nt > 1
    msg = [msg ' No usable frame times in the headers: enter s/frame.'];
elseif nt > 1 && max(abs(diff(times) - dt)) > 0.1 * dt
    msg = sprintf('%s Frame times are not uniform (%.2f-%.2f s); the median %.2f s is used.', msg, ...
        min(diff(times)), max(diff(times)), dt);
end

% geometry
k1 = r(find(zi == 1, 1)); kN = r(find(zi == nz, 1));
X = R.IOP(k1, 1:3); Y = R.IOP(k1, 4:6);
ps = R.PS(k1, :);
if nz > 1
    sv = (R.IPP(kN, :) - R.IPP(k1, :)) / (nz - 1);
else
    sv = cross(X, Y) * R.Thick(k1);
end
A = [X(:) * ps(2), Y(:) * ps(1), sv(:), R.IPP(k1, :)'; 0 0 0 1];
A = diag([-1 -1 1 1]) * A;                       % LPS -> RAS
vox = [ps(2) ps(1) norm(sv)];

% pixels
fh = scan.FH(R.File(r(1)));
nx = fh.Cols; ny = fh.Rows;
V = zeros(nx, ny, nz, nt, 'single');
[fl, ~, fj] = unique(R.File(r));
for a = 1:numel(fl)
    if mod(a, 25) == 1
        report(0.6 + 0.4 * (a - 1) / numel(fl), sprintf('DICOM: reading images (%d of %d files)...', a, numel(fl)));
    end
    k = fl(a);
    f = scan.FH(k);
    if f.Cols ~= nx || f.Rows ~= ny
        error('graspilot_dicom_read:size', 'Series %g: images differ in matrix size.', S.Number);
    end
    if f.Spp ~= 1
        error('graspilot_dicom_read:rgb', 'Series %g: colour images are not supported.', S.Number);
    end
    mine = find(fj == a);
    P = readPixels(scan.Files{k}, f);            % [nx ny NF]
    for q = mine(:)'
        img = single(P(:, :, R.Frame(r(q))));
        sl = R.Slope(r(q)); ic = R.Inter(r(q));
        if sl ~= 1 || ic ~= 0, img = img * single(sl) + single(ic); end
        V(:, :, zi(q), ti(q)) = img;
    end
end
info = struct('Affine', A, 'Vox', vox, 'dt', dt, 'Times', times, 'StartTime', t0, 'Series', S.Description, ...
    'SeriesNumber', S.Number, 'Modality', S.Modality, 'NFiles', numel(fl), 'Folder', scan.Folder, 'Message', msg);
end

function P = readPixels(file, f)
n = f.Cols * f.Rows * f.NF;
if f.Encaps
    if exist('dicomread', 'file') ~= 2
        error('graspilot_dicom_read:compressed', ['%s holds compressed pixel data (JPEG, JPEG 2000 or RLE). ' ...
            'Reading it needs dicomread (Image Processing Toolbox). Without it, convert the series to ' ...
            'uncompressed DICOM or to NIfTI (e.g. dcm2niix) and open that.'], file);
    end
    P = squeeze(dicomread(file));                % [rows cols NF]
    P = permute(P, [2 1 3]);
    return
end
switch f.Bits
    case 8,  t = 'int8';
    case 16, t = 'int16';
    case 32, t = 'int32';
    otherwise
        error('graspilot_dicom_read:bits', '%d bits per pixel are not supported (%s).', f.Bits, file);
end
if ~f.Signed, t = ['u' t]; end
mf = 'ieee-le';
if f.BE, mf = 'ieee-be'; end
fid = fopen(file, 'r', mf);
if fid < 0, error('graspilot_dicom_read:open', 'Cannot open %s.', file); end
fseek(fid, f.Offset - 1, 'bof');
P = fread(fid, n, [t '=>' t]);
fclose(fid);
if numel(P) ~= n, error('graspilot_dicom_read:short', 'Pixel data of %s are truncated.', file); end
P = reshape(P, f.Cols, f.Rows, f.NF);
end

% =========================================================================
function h = readHeader(file, D)
h = [];
d = dir(file);
if isempty(d) || d.bytes < 136, return; end
chunk = 2 ^ 20;
fid = fopen(file, 'r');
if fid < 0, return; end
while true
    fseek(fid, 0, 'bof');
    b = fread(fid, min(chunk, d.bytes), '*uint8')';
    h = parseFile(b, D, numel(b) < d.bytes);
    if isempty(h) || ~isfield(h, 'Truncated') || numel(b) >= d.bytes, break; end
    chunk = chunk * 8;
end
fclose(fid);
end

function h = parseFile(b, D, partial)
h = [];
expl = false; be = false;
if numel(b) >= 132 && isequal(b(129:132), uint8('DICM'))
    % file meta information (group 0002): always explicit VR little endian
    pos = 133;
    ts = '';
    while pos + 8 <= numel(b)
        g = double(b(pos)) + 256 * double(b(pos + 1));
        if g ~= 2, break; end
        e = double(b(pos + 2)) + 256 * double(b(pos + 3));
        vr = char(b(pos + 4:pos + 5));
        if any(strcmp(vr, {'OB', 'OW', 'OF', 'SQ', 'UT', 'UN'}))
            len = u32(b, pos + 8, false); vp = pos + 12;
        else
            len = u16(b, pos + 6, false); vp = pos + 8;
        end
        if e == 16, ts = strtrim(char(b(vp:min(end, vp + len - 1)))); ts(ts == 0) = []; end
        pos = vp + len;
    end
    switch ts
        case '1.2.840.10008.1.2',       expl = false;
        case '1.2.840.10008.1.2.2',     expl = true; be = true;
        case '1.2.840.10008.1.2.1.99',  return        % deflated: not supported
        otherwise,                      expl = true;  % explicit LE and the encapsulated syntaxes
    end
elseif numel(b) >= 8 && b(1) == 8 && b(2) == 0
    % no preamble: assume a bare data set, implicit unless a VR follows
    pos = 1;
    vr = char(b(5:6));
    expl = all(vr >= 'A' & vr <= 'Z');
else
    return
end
h = struct('BE', be);
[h, ~] = parseElems(b, pos, numel(b), expl, be, h, D, true, true);
if ~isfield(h, 'PixelOffset') && partial, h.Truncated = true; end
if ~isfield(h, 'PixelOffset') && ~partial, h = []; end
end

function [h, pos] = parseElems(b, pos, endPos, expl, be, h, D, collect, top)
% data elements from pos; stops at endPos, an item/sequence delimiter or the pixel data
UNDEF = 4294967295;
while pos + 7 <= endPos
    g = u16(b, pos, be); e = u16(b, pos + 2, be);
    if g == 65534                                  % FFFE: delimiters
        pos = pos + 8;
        return
    end
    sq = false;
    if expl
        vr = char(b(pos + 4:pos + 5));
        if vr(1) == 'O' || vr(1) == 'U' && vr(2) ~= 'I' && vr(2) ~= 'L' && vr(2) ~= 'S' || ...
                strcmp(vr, 'SQ') || strcmp(vr, 'SV')
            len = u32(b, pos + 8, be); vp = pos + 12;
        else
            len = u16(b, pos + 6, be); vp = pos + 8;
        end
        sq = strcmp(vr, 'SQ');
    else
        len = u32(b, pos + 4, be); vp = pos + 8;
    end
    tag = g * 65536 + e;
    if tag == 2145386512 && top                    % (7FE0,0010) pixel data
        h.PixelOffset = vp; h.PixelLength = len;
        pos = endPos + 1;
        return
    end
    isSeqTag = any(tag == D.SeqTags);
    if len == UNDEF || sq || isSeqTag
        % sequence: descend when it is (inside) a functional group, else skip
        want = collect && (isSeqTag || ~top);
        if len == UNDEF, se = endPos; else, se = min(endPos, vp + len - 1); end
        [items, p2] = parseSeq(b, vp, se, expl, be, D, want, len == UNDEF);
        if len == UNDEF, pos = p2; else, pos = vp + len; end
        if want && ~isempty(items)
            if tag == 1375769136                   % (5200,9230) per-frame functional groups
                h.Frames = items;
            elseif tag == 1375769129               % (5200,9229) shared functional groups
                h.Shared = items{1};
            else                                   % nested macro: flatten into this item
                fn = fieldnames(items{1});
                for q = 1:numel(fn), h.(fn{q}) = items{1}.(fn{q}); end
            end
        end
        continue
    end
    if vp + len - 1 > endPos                       % value beyond the buffer
        h.Truncated = true;
        pos = endPos + 1;
        return
    end
    if collect && any(g == D.Groups)
        k = find(D.Tags == tag, 1);
        if ~isempty(k) && len > 0
            raw = b(vp:vp + len - 1);
            switch D.Types{k}
                case 's'
                    v = char(raw); v(v == 0) = [];
                    v = strtrim(v);
                case 'n'
                    v = char(raw); v(v == 0) = ' '; v(v == '\') = ' ';
                    v = sscanf(v, '%f')';
                case 'us'
                    v = u16(b, vp, be);
                case 'ul'
                    v = u32(b, vp, be);
                case 'tm'
                    v = tm2sec(char(raw));
                case 'dt'
                    v = dt2sec(char(raw));
            end
            if ~isempty(v), h.(D.Names{k}) = v; end
        end
    end
    pos = vp + len;
end
end

function [items, pos] = parseSeq(b, pos, endPos, expl, be, D, collect, undef)
items = {};
UNDEF = 4294967295;
while pos + 7 <= endPos
    g = u16(b, pos, be); e = u16(b, pos + 2, be); len = u32(b, pos + 4, be);
    pos = pos + 8;
    if g ~= 65534, break; end
    if e == 57565                                  % (FFFE,E0DD) sequence delimiter
        return
    end
    if e ~= 57344, continue; end                   % not an item
    if len == UNDEF
        [hi, pos] = parseElems(b, pos, endPos, expl, be, struct(), D, collect, false);
    else
        [hi, ~] = parseElems(b, pos, min(endPos, pos + len - 1), expl, be, struct(), D, collect, false);
        pos = pos + len;
    end
    if collect, items{end + 1} = hi; end %#ok<AGROW>
end
if ~undef, pos = endPos + 1; end
end

% =========================================================================
function D = tagTable()
T = { ...
    '0008', '0060', 'Modality', 's'; ...
    '0008', '103E', 'SeriesDescription', 's'; ...
    '0008', '0032', 'AcqTime', 'tm'; ...
    '0008', '0033', 'ContentTime', 'tm'; ...
    '0008', '002A', 'AcqDateTime', 'dt'; ...
    '0018', '0050', 'SliceThickness', 'n'; ...
    '0018', '0088', 'SpacingBetweenSlices', 'n'; ...
    '0018', '1060', 'TriggerTime', 'n'; ...
    '0018', '9074', 'FrameAcqDateTime', 'dt'; ...
    '0020', '000E', 'SeriesUID', 's'; ...
    '0020', '0011', 'SeriesNumber', 'n'; ...
    '0020', '0012', 'AcquisitionNumber', 'n'; ...
    '0020', '0013', 'InstanceNumber', 'n'; ...
    '0020', '0032', 'IPP', 'n'; ...
    '0020', '0037', 'IOP', 'n'; ...
    '0020', '0100', 'TemporalPosition', 'n'; ...
    '0020', '9128', 'TemporalPositionIndex', 'ul'; ...
    '0028', '0002', 'SamplesPerPixel', 'us'; ...
    '0028', '0008', 'NumberOfFrames', 'n'; ...
    '0028', '0010', 'Rows', 'us'; ...
    '0028', '0011', 'Cols', 'us'; ...
    '0028', '0030', 'PixelSpacing', 'n'; ...
    '0028', '0100', 'BitsAllocated', 'us'; ...
    '0028', '0103', 'PixelRepresentation', 'us'; ...
    '0028', '1052', 'RescaleIntercept', 'n'; ...
    '0028', '1053', 'RescaleSlope', 'n'};
n = size(T, 1);
D.Tags = zeros(1, n);
for k = 1:n
    D.Tags(k) = hex2dec(T{k, 1}) * 65536 + hex2dec(T{k, 2});
end
D.Names = T(:, 3)';
D.Types = T(:, 4)';
D.Groups = unique(floor(D.Tags / 65536));
% sequences to descend into: shared / per-frame functional groups and their macros
sq = {'5200', '9229'; '5200', '9230'; '0020', '9113'; '0020', '9116'; '0028', '9110'; '0020', '9111'; ...
    '0028', '9145'};
D.SeqTags = zeros(1, size(sq, 1));
for k = 1:size(sq, 1)
    D.SeqTags(k) = hex2dec(sq{k, 1}) * 65536 + hex2dec(sq{k, 2});
end
end

function v = u16(b, pos, be)
if be
    v = 256 * double(b(pos)) + double(b(pos + 1));
else
    v = double(b(pos)) + 256 * double(b(pos + 1));
end
end

function v = u32(b, pos, be)
if be
    v = 16777216 * double(b(pos)) + 65536 * double(b(pos + 1)) + 256 * double(b(pos + 2)) + double(b(pos + 3));
else
    v = double(b(pos)) + 256 * double(b(pos + 1)) + 65536 * double(b(pos + 2)) + 16777216 * double(b(pos + 3));
end
end

function s = tm2sec(t)
% DICOM TM (HHMMSS.FFFFFF) -> seconds
t = strtrim(t); t(t == 0 | t == ':') = [];
s = NaN;
if numel(t) < 2, return; end
x = [0 0 0];
x(1) = str2double(t(1:2));
if numel(t) >= 4, x(2) = str2double(t(3:4)); end
if numel(t) >= 6, x(3) = str2double(t(5:end)); end
s = x(1) * 3600 + x(2) * 60 + x(3);
end

function s = dt2sec(t)
% DICOM DT (YYYYMMDDHHMMSS.FFFFFF) -> seconds (day number x 86400 + time of day)
t = strtrim(t); t(t == 0) = [];
k = find(t == '+' | t == '-', 1);
if ~isempty(k), t = t(1:k - 1); end
s = NaN;
if numel(t) < 10, return; end
day = datenum(str2double(t(1:4)), str2double(t(5:6)), str2double(t(7:8)));
s = (day - 730000) * 86400 + tm2sec(t(9:end));
end

function v = getf(h, name, def)
if isfield(h, name) && ~isempty(h.(name)), v = h.(name); else, v = def; end
end

function v = pick3(a, b, c, name, def)
if isfield(a, name) && ~isempty(a.(name)), v = a.(name);
elseif isfield(b, name) && ~isempty(b.(name)), v = b.(name);
elseif isfield(c, name) && ~isempty(c.(name)), v = c.(name);
else, v = def;
end
end

function v = vec(x, n)
v = x(:)';
if numel(v) < n, v(end + 1:n) = NaN; end
v = v(1:n);
end

function v = first(x)
if isempty(x), v = NaN; else, v = x(1); end
end
