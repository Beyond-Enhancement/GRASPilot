function files = gdce_export_params(file, D, hdr, rows, blocks, keys, units, settings)
%GDCE_EXPORT_PARAMS Write semi-quantitative parameter statistics and voxel values.
%
%   files = gdce_export_params(file, D, hdr, rows, blocks, keys, units, settings)
%
%   file      output .xlsx or .csv
%   D         GraspData object (geometry, binning)
%   hdr       statistics header (cell row), rows: statistics (cell array),
%             one row per ROI and parameter (GraspDCEViewer.paramStatsTable)
%   blocks    struct array per ROI: id, name, lin (canonical linear indices),
%             values [nvox x numel(keys)]
%   keys      parameter keys (e.g. 'PE', 'PE_SSS', ...), units: their units
%   settings  {name, value} cell array describing how the values were computed
%
%   .xlsx : sheets Stats, Settings and one sheet per ROI with one row per
%           voxel: Voxel, i, j, k, x_mm, y_mm, z_mm and one column per
%           parameter (header "<key> [unit]").
%   .csv  : <base>_param_stats.csv, <base>_param_voxels.csv (all ROIs, with
%           ROI_ID and ROI_Name columns) and <base>_param_settings.txt.
%   After temporal binning, _<N>frames is appended to the base name (unless
%   already present).
%
%   i, j, k are 1-based indices into the arrays exported to NIfTI by the
%   viewer (canonical LPS voxel order); x, y, z are RAS mm from the affine.

[folder, base, ext] = fileparts(file);
if isempty(folder), folder = pwd; end
ext = lower(ext);
if D.BinSize > 1
    tag = sprintf('_%dframes', D.NT);
    if isempty(strfind(base, tag)), base = [base tag]; end
    file = fullfile(folder, [base ext]);
end
n = D.canonSize();
A = D.affine();
phdr = cell(1, numel(keys));
for k = 1:numel(keys)
    phdr{k} = sprintf('%s [%s]', keys{k}, units{k});
end
vhdr = [{'Voxel', 'i', 'j', 'k', 'x_mm', 'y_mm', 'z_mm'}, phdr];
info = [{'Created by', 'GRASPilot (gdce_export_params)';
    'Date', datestr(now, 'yyyy-mm-dd HH:MM:SS');
    'Source file', D.SourceFile;
    'Frames', D.NT;
    'Temporal binning (raw frames averaged per frame)', D.BinSize;
    'Temporal resolution (s/frame, after binning)', D.effectiveDt();
    'Voxel size i j k (mm)', num2str(D.canonVox())}; settings];
for k = 1:numel(keys)
    info(end + 1, :) = {sprintf('Column %s', keys{k}), sprintf('unit %s', units{k})}; %#ok<AGROW>
end
for q = 1:numel(blocks)
    info(end + 1, :) = {sprintf('ROI %d', blocks(q).id), sprintf('%s: %d voxels', blocks(q).name, ...
        numel(blocks(q).lin))}; %#ok<AGROW>
end
for q = 1:numel(blocks)
    [i1, i2, i3] = ind2sub(n, blocks(q).lin(:));
    xyz = A * [i1' - 1; i2' - 1; i3' - 1; ones(1, numel(i1))];
    blocks(q).ijk = [i1 i2 i3];
    blocks(q).xyz = xyz(1:3, :)';
end

switch ext
    case {'.xlsx', '.xls'}
        if exist('writecell', 'file') ~= 2 || exist('writematrix', 'file') ~= 2
            error('gdce_export_params:xlsx', ['Writing .xlsx needs MATLAB R2019a or newer ' ...
                '(writecell/writematrix). Export to .csv instead.']);
        end
        for q = 1:numel(blocks)
            if numel(blocks(q).lin) + 1 > 1048576
                error('gdce_export_params:rows', ['ROI %s has more voxels than an Excel sheet ' ...
                    'can hold. Export to .csv instead.'], blocks(q).name);
            end
        end
        if exist(file, 'file') == 2, delete(file); end
        writecell([hdr; nanToEmpty(rows)], file, 'Sheet', 'Stats');
        writecell(stringify(info), file, 'Sheet', 'Settings');
        used = {'Stats', 'Settings'};
        for q = 1:numel(blocks)
            sh = sheetName(blocks(q).id, blocks(q).name, used);
            used{end + 1} = sh; %#ok<AGROW>
            nv = numel(blocks(q).lin);
            M = [(1:nv)' blocks(q).ijk blocks(q).xyz blocks(q).values];
            writecell(vhdr, file, 'Sheet', sh, 'Range', 'A1');
            writematrix(M, file, 'Sheet', sh, 'Range', 'A2');
        end
        files = {file};

    case '.csv'
        fs = fullfile(folder, [base '_param_stats.csv']);
        fv = fullfile(folder, [base '_param_voxels.csv']);
        fi = fullfile(folder, [base '_param_settings.txt']);
        % statistics
        fid = fopen(fs, 'w');
        if fid < 0, error('gdce_export_params:open', 'Cannot write %s.', fs); end
        fprintf(fid, '%s\n', strjoin(cellfun(@quote, hdr, 'UniformOutput', false), ','));
        for r = 1:size(rows, 1)
            c = cell(1, size(rows, 2));
            for j = 1:size(rows, 2)
                v = rows{r, j};
                if ischar(v), c{j} = quote(v);
                elseif isfinite(v), c{j} = sprintf('%.7g', v);
                else, c{j} = ''; end
            end
            fprintf(fid, '%s\n', strjoin(c, ','));
        end
        fclose(fid);
        % voxel values
        fid = fopen(fv, 'w');
        if fid < 0, error('gdce_export_params:open', 'Cannot write %s.', fv); end
        fprintf(fid, '%s\n', strjoin(cellfun(@quote, [{'ROI_ID', 'ROI_Name'}, vhdr], ...
            'UniformOutput', false), ','));
        for q = 1:numel(blocks)
            nv = numel(blocks(q).lin);
            nm = quote(blocks(q).name);
            for v = 1:nv
                x = blocks(q).values(v, :);
                xs = arrayfun(@num, x, 'UniformOutput', false);
                fprintf(fid, '%d,%s,%d,%d,%d,%d,%.4f,%.4f,%.4f,%s\n', blocks(q).id, nm, v, ...
                    blocks(q).ijk(v, :), blocks(q).xyz(v, :), strjoin(xs, ','));
            end
        end
        fclose(fid);
        % settings
        fid = fopen(fi, 'w');
        s = stringify(info);
        for k = 1:size(s, 1)
            fprintf(fid, '%s: %s\n', s{k, 1}, s{k, 2});
        end
        fclose(fid);
        files = {fs, fv, fi};

    otherwise
        error('gdce_export_params:ext', 'Use a .xlsx or .csv file name.');
end
end

% -------------------------------------------------------------------------
function s = num(x)
if isfinite(x), s = sprintf('%.7g', x); else, s = ''; end
end

function s = quote(s)
s = ['"' strrep(s, '"', '""') '"'];
end

function c = nanToEmpty(c)
for k = 1:numel(c)
    if isnumeric(c{k}) && isscalar(c{k}) && ~isfinite(c{k}), c{k} = []; end
end
end

function c = stringify(c)
for k = 1:numel(c)
    if isnumeric(c{k}) || islogical(c{k}), c{k} = num2str(c{k}); end
end
end

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
