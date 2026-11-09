function files = gdce_export_timeseries(D, masks, names, file, varargin)
%GDCE_EXPORT_TIMESERIES Write voxel-wise ROI time series to .xlsx or .csv.
%
%   files = gdce_export_timeseries(D, masks, names, file, 'Name', Value, ...)
%
%   D       GraspData object
%   masks   cell array of canonical logical masks (same size as D.canonSize())
%   names   cell array of ROI names
%   file    output .xlsx or .csv
%
%   Options
%   'BaseFrames'  baseline frame indices for relative enhancement (default 1)
%   'Progress'    function handle progress(fraction, message)
%
%   .xlsx : sheets Info, Summary (per-ROI mean/SD/median/relative
%           enhancement per frame) and one sheet per ROI with one row per
%           voxel: Voxel, i, j, k, x_mm, y_mm, z_mm, F001..Fnnn.
%   .csv  : <base>_voxels.csv (all ROIs, ROI_ID and ROI_Name columns),
%           <base>_summary.csv and <base>_info.txt.
%
%   i, j, k are 1-based indices into the arrays exported to NIfTI by the
%   viewer (canonical LPS voxel order); x, y, z are RAS mm from the NIfTI
%   affine. Signal values are in the units of the loaded data (normalized
%   to the global maximum when ReadData.m normalization is on).

p = struct('BaseFrames', 1, 'Progress', []);
for i = 1:2:numel(varargin)
    fn = fieldnames(p);
    m = find(strcmpi(fn, varargin{i}), 1);
    if isempty(m), error('gdce_export_timeseries:opt', 'Unknown option %s.', varargin{i}); end
    p.(fn{m}) = varargin{i + 1};
end

[folder, base, ext] = fileparts(file);
if isempty(folder), folder = pwd; end
ext = lower(ext);
T = D.NT;
ft = D.frameTimes();
n = D.canonSize();
A = D.affine();
voxmL = prod(D.canonVox()) / 1000;
bf = p.BaseFrames(p.BaseFrames >= 1 & p.BaseFrames <= T);
if isempty(bf), bf = 1; end

keep = find(cellfun(@(m) any(m(:)), masks));
if isempty(keep)
    error('gdce_export_timeseries:empty', 'All ROIs are empty.');
end
R = numel(keep);

% gather data
roi = struct('id', {}, 'name', {}, 'ijk', {}, 'xyz', {}, 'TC', {});
for q = 1:R
    r = keep(q);
    lin = find(masks{r});
    [i1, i2, i3] = ind2sub(n, lin);
    xyz = A * [i1' - 1; i2' - 1; i3' - 1; ones(1, numel(lin))];
    roi(q).id = r;
    roi(q).name = names{r};
    roi(q).ijk = [i1 i2 i3];
    roi(q).xyz = xyz(1:3, :)';
    roi(q).TC = double(D.timecourses(lin));
    report(p.Progress, 0.5 * q / R, sprintf('Collecting %s', names{r}));
end

% summary table
sumHdr = {'Frame', 'Time_s'};
sumDat = [(1:T)' ft(:)];
for q = 1:R
    TC = roi(q).TC;
    mu = mean(TC, 1);
    if size(TC, 1) > 1, sd = std(TC, 0, 1); else, sd = zeros(1, T); end
    md = median(TC, 1);
    s0 = mean(mu(bf));
    re = 100 * (mu - s0) / s0;
    nm = roi(q).name;
    sumHdr = [sumHdr, {[nm '_mean'], [nm '_sd'], [nm '_median'], [nm '_RelEnh_pct'], [nm '_nvox']}]; %#ok<AGROW>
    sumDat = [sumDat, mu(:), sd(:), md(:), re(:), repmat(size(TC, 1), T, 1)]; %#ok<AGROW>
end

% info block
info = {
    'Created by', 'GRASP DCE Viewer (gdce_export_timeseries)';
    'Date', datestr(now, 'yyyy-mm-dd HH:MM:SS');
    'Source file', D.SourceFile;
    'Source variables', D.SourceVars;
    'Data mode', D.Mode;
    'Native size (rows cols partitions)', num2str(D.NativeSize);
    'Exported array size (i j k)', num2str(n);
    'Voxel size i j k (mm)', num2str(D.canonVox());
    'Native acquisition plane', D.Plane;
    'Frames', T;
    'Temporal binning (raw frames averaged per frame)', D.BinSize;
    'Raw frames', D.NTRaw;
    'Dropped raw frames (incomplete last bin)', num2str(D.droppedFrames());
    'Temporal resolution (s/frame, after binning)', D.effectiveDt();
    'Time of frame 1 (s)', D.t0;
    'Crop factor', D.CropFactor;
    'Normalized to global max (ReadData.m)', double(D.Normalized);
    'Baseline frames (relative enhancement)', sprintf('%d-%d', min(bf), max(bf));
    'Coordinates', 'i,j,k 1-based indices into the exported NIfTI arrays (LPS voxel order); x,y,z RAS mm from the NIfTI affine';
    'Relative enhancement', '100*(mean(t)-S0)/S0 of the ROI mean curve, S0 = mean of baseline frames';
    };
for q = 1:R
    info(end + 1, :) = {sprintf('ROI %d', roi(q).id), sprintf('%s: %d voxels, %.4f mL', ...
        roi(q).name, size(roi(q).TC, 1), size(roi(q).TC, 1) * voxmL)}; %#ok<AGROW>
end

fhdr = arrayfun(@(t) sprintf('F%03d', t), 1:T, 'UniformOutput', false);
vhdr = [{'Voxel', 'i', 'j', 'k', 'x_mm', 'y_mm', 'z_mm'}, fhdr];

switch ext
    case {'.xlsx', '.xls'}
        if exist('writecell', 'file') ~= 2 || exist('writematrix', 'file') ~= 2
            error('gdce_export_timeseries:xlsx', ['Writing .xlsx needs MATLAB R2019a or newer ' ...
                '(writecell/writematrix). Export to .csv instead.']);
        end
        for q = 1:R
            if size(roi(q).TC, 1) + 1 > 1048576
                error('gdce_export_timeseries:rows', ['ROI %s has %d voxels, more than an Excel ' ...
                    'sheet can hold. Export to .csv instead.'], roi(q).name, size(roi(q).TC, 1));
            end
        end
        if 7 + T > 16384
            error('gdce_export_timeseries:cols', 'Too many frames for an Excel sheet; use .csv.');
        end
        if exist(file, 'file') == 2, delete(file); end
        writecell(info, file, 'Sheet', 'Info');
        writecell(sumHdr, file, 'Sheet', 'Summary', 'Range', 'A1');
        writematrix(sumDat, file, 'Sheet', 'Summary', 'Range', 'A2');
        used = {'Info', 'Summary'};
        for q = 1:R
            sh = sheetName(roi(q).id, roi(q).name, used);
            used{end + 1} = sh; %#ok<AGROW>
            nv = size(roi(q).TC, 1);
            M = [(1:nv)' roi(q).ijk roi(q).xyz roi(q).TC];
            writecell(vhdr, file, 'Sheet', sh, 'Range', 'A1');
            writematrix(M, file, 'Sheet', sh, 'Range', 'A2');
            report(p.Progress, 0.5 + 0.5 * q / R, sprintf('Writing sheet %s', sh));
        end
        files = {file};

    case '.csv'
        fv = fullfile(folder, [base '_voxels.csv']);
        fs = fullfile(folder, [base '_summary.csv']);
        fi = fullfile(folder, [base '_info.txt']);
        % voxel file
        fid = fopen(fv, 'w');
        if fid < 0, error('gdce_export_timeseries:open', 'Cannot write %s.', fv); end
        fprintf(fid, '%s\n', strjoin([{'ROI_ID', 'ROI_Name'}, vhdr], ','));
        for q = 1:R
            nv = size(roi(q).TC, 1);
            nm = strrep(strrep(roi(q).name, '"', '""'), '%', '%%');
            fmt = ['%d,"' nm '",%d,%d,%d,%d,%.4f,%.4f,%.4f' repmat(',%.7g', 1, T) '\n'];
            M = [repmat(roi(q).id, nv, 1) (1:nv)' roi(q).ijk roi(q).xyz roi(q).TC];
            fprintf(fid, fmt, M');
            report(p.Progress, 0.5 + 0.5 * q / R, sprintf('Writing %s', roi(q).name));
        end
        fclose(fid);
        % summary
        fid = fopen(fs, 'w');
        hq = cellfun(@(s) ['"' strrep(s, '"', '""') '"'], sumHdr, 'UniformOutput', false);
        fprintf(fid, '%s\n', strjoin(hq, ','));
        fprintf(fid, ['%d' repmat(',%.7g', 1, size(sumDat, 2) - 1) '\n'], sumDat');
        fclose(fid);
        % info
        fid = fopen(fi, 'w');
        for k = 1:size(info, 1)
            v = info{k, 2};
            if isnumeric(v), v = num2str(v); end
            fprintf(fid, '%s: %s\n', info{k, 1}, v);
        end
        fprintf(fid, 'Frame times (s): %s\n', num2str(ft, '%.3f '));
        fclose(fid);
        files = {fv, fs, fi};

    otherwise
        error('gdce_export_timeseries:ext', 'Use a .xlsx or .csv file name.');
end
report(p.Progress, 1, 'Done');
end

% -------------------------------------------------------------------------
function sh = sheetName(id, name, used)
s = regexprep(name, '[\[\]:\*\?/\\]', '_');
sh = sprintf('ROI%d_%s', id, s);
sh = sh(1:min(31, numel(sh)));
k = 1;
base = sh;
while any(strcmpi(used, sh))
    suf = sprintf('_%d', k);
    sh = [base(1:min(numel(base), 31 - numel(suf))) suf];
    k = k + 1;
end
end

function report(fcn, frac, msg)
if ~isempty(fcn)
    try, fcn(frac, msg); catch, end
end
end
