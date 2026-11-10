function out = gdce_load_dialog(matFile, opts)
%GDCE_LOAD_DIALOG Choose variables, format and geometry for a GRASP .mat file.
%
%   opts = gdce_load_dialog(matFile, defaults)
%   opts = gdce_load_dialog('', defaults)      geometry only (data in memory)
%
%   Returns [] if cancelled. Fields of opts: Format ('subspace'|'dense'),
%   CoefVar, BasisVar, ImageVar, DimOrder, CropFactor, Normalize, Plane,
%   Voxel [row col slice] mm, dt (s/frame), t0 (s).

out = [];
geomOnly = isempty(matFile);
names = {}; info = [];
if ~geomOnly
    info = whos('-file', matFile);
    names = {info.name};
end
numIdx = [];
for i = 1:numel(info)
    if any(strcmp(info(i).class, {'double', 'single', 'int16', 'uint16', 'int32', 'uint8', 'uint32', 'int8'}))
        numIdx(end + 1) = i; %#ok<AGROW>
    end
end
numNames = names(numIdx);
if isempty(numNames), numNames = {'(none)'}; end

bg = [0.94 0.94 0.94];
if geomOnly, H = 290; else, H = 590; end
W = 640;
ss = get(0, 'ScreenSize');
d = figure('Name', 'GRASP DCE Viewer - load', 'NumberTitle', 'off', 'MenuBar', 'none', ...
    'ToolBar', 'none', 'Resize', 'off', 'Color', bg, 'Units', 'pixels', ...
    'Position', [max(1, (ss(3) - W) / 2), max(40, (ss(4) - H) / 2), W, H], ...
    'WindowStyle', 'modal', 'CloseRequestFcn', @onCancel);

y = H;
hFmt = []; hCoef = []; hBasis = []; hImg = []; hOrder = []; hCrop = []; hNorm = []; hInfo = [];
if ~geomOnly
    y = y - 30;
    [~, fb, fe] = fileparts(matFile);
    txt(12, y, W - 24, ['File: ' fb fe], 'FontWeight', 'bold');
    y = y - 96;
    desc = cell(1, numel(info));
    for i = 1:numel(info)
        c = '';
        if info(i).complex, c = ' complex'; end
        desc{i} = sprintf('%-16s [%s]  %s%s  (%.0f MB)', info(i).name, ...
            strjoin(arrayfun(@num2str, info(i).size, 'UniformOutput', false), 'x'), ...
            info(i).class, c, info(i).bytes / 2^20);
    end
    uicontrol(d, 'Style', 'listbox', 'Units', 'pixels', 'Position', [12 y W - 24 90], ...
        'String', desc, 'Enable', 'inactive', 'FontName', 'FixedWidth', 'BackgroundColor', 'w');

    hasSub = any(strcmp(names, opts.CoefVar)) && any(strcmp(names, opts.BasisVar));
    fmtIdx = 2;
    if hasSub || strcmpi(opts.Format, 'subspace'), fmtIdx = 1; end
    if strcmpi(opts.Format, 'dense'), fmtIdx = 2; end
    y = y - 32;
    txt(12, y, 110, 'Format');
    hFmt = pop(125, y + 2, W - 137, {'GRASP subspace: coefficients + temporal basis (ReadData.m)', ...
        '4D magnitude image array'}, fmtIdx);
    set(hFmt, 'Callback', @(~, ~) onFormat());

    y = y - 30;
    txt(12, y, 110, 'Coefficients');
    hCoef = pop(125, y + 2, 160, numNames, pickIdx(numNames, opts.CoefVar, 1));
    txt(300, y, 110, 'Temporal basis');
    hBasis = pop(415, y + 2, 213, numNames, pickIdx(numNames, opts.BasisVar, min(2, numel(numNames))));

    y = y - 30;
    txt(12, y, 110, 'Image array');
    hImg = pop(125, y + 2, 160, numNames, pickIdx(numNames, opts.ImageVar, 1));
    txt(300, y, 110, 'Dimension order');
    hOrder = edt(415, y + 2, 60, opts.DimOrder);
    set(hOrder, 'Callback', @(~, ~) updateInfo());
    txt(482, y, 140, 'e.g. xytz, xyzt', 'ForegroundColor', [0.4 0.4 0.4]);

    y = y - 30;
    txt(12, y, 110, 'Crop factor');
    hCrop = edt(125, y + 2, 50, '');
    set(hCrop, 'Callback', @(~, ~) updateInfo());
    txt(180, y, 120, '(2 = CropImg.m)', 'ForegroundColor', [0.4 0.4 0.4]);
    hNorm = uicontrol(d, 'Style', 'checkbox', 'Units', 'pixels', 'Position', [300 y + 2 328 20], ...
        'String', 'Normalize to global max (ReadData.m)', 'BackgroundColor', bg);

    y = y - 62;
    hInfo = uicontrol(d, 'Style', 'text', 'Units', 'pixels', 'Position', [12 y W - 24 56], ...
        'String', '', 'HorizontalAlignment', 'left', 'BackgroundColor', [1 1 0.92], ...
        'FontName', 'FixedWidth');
end

y = y - 34;
txt(12, y, 300, 'Acquisition geometry', 'FontWeight', 'bold');
y = y - 30;
txt(12, y, 160, 'Voxel orientation');
codes = {'LPS', 'PLS', 'ILP', 'IPL'};
olabels = {'LPS: axial, rows R->L, cols A->P (scanner order, GRASP recon)', ...
    'PLS: axial, rows A->P, cols R->L (data1 already radiological)', ...
    'ILP: coronal slices, rows S->I, cols R->L', ...
    'IPL: sagittal slices, rows S->I, cols A->P'};
cur = '';
if isfield(opts, 'AxisCodes'), cur = upper(opts.AxisCodes); end
if isempty(cur)
    [pm, fl] = GraspData.presetOrientation(opts.Plane, opts.SwapInPlane, opts.Flip);
    cur = GraspData.permFlipToCodes(pm, fl);
end
if ~any(strcmp(codes, cur))
    codes{end + 1} = cur;
    olabels{end + 1} = sprintf('%s: as set last time (View > Geometry)', cur);
end
hPlane = pop(178, y + 2, W - 190, olabels, find(strcmp(codes, cur), 1));
set(hPlane, 'Callback', []);
y = y - 30;
txt(12, y, 160, 'Pixel spacing (mm)');
txt(178, y, 34, 'rows');
hVr = edt(212, y + 2, 55, num2str(opts.Voxel(1)));
txt(277, y, 32, 'cols');
hVc = edt(309, y + 2, 55, num2str(opts.Voxel(2)));
txt(376, y, 100, 'slice spacing');
hVs = edt(478, y + 2, 55, num2str(opts.Voxel(3)));
y = y - 30;
txt(12, y, 166, 'Temporal resolution (s)');
dts = '';
if isfinite(opts.dt), dts = num2str(opts.dt); end
hDt = edt(178, y + 2, 70, dts);
txt(262, y, 140, 'Time of frame 1 (s)');
hT0 = edt(400, y + 2, 60, num2str(opts.t0));
y = y - 30;
txt(12, y, W - 24, 'Temporal resolution = seconds per reconstructed frame; leave blank if unknown.', ...
    'ForegroundColor', [0.4 0.4 0.4]);
txt(12, y - 18, W - 24, 'Geometry can be changed later under View > Geometry.', ...
    'ForegroundColor', [0.4 0.4 0.4]);

uicontrol(d, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [W - 210 12 95 30], ...
    'String', 'OK', 'FontWeight', 'bold', 'Callback', @onOK);
uicontrol(d, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [W - 108 12 95 30], ...
    'String', 'Cancel', 'Callback', @onCancel);

if ~geomOnly
    onFormat(true);
end
uiwait(d);

% ------------------------------------------------------------------ nested
    function h = txt(x, yy, w, s, varargin)
        h = uicontrol(d, 'Style', 'text', 'Units', 'pixels', 'Position', [x yy w 18], ...
            'String', s, 'HorizontalAlignment', 'left', 'BackgroundColor', bg, varargin{:});
    end
    function h = edt(x, yy, w, s)
        h = uicontrol(d, 'Style', 'edit', 'Units', 'pixels', 'Position', [x yy w 22], ...
            'String', s, 'BackgroundColor', 'w', 'HorizontalAlignment', 'left');
    end
    function h = pop(x, yy, w, s, v)
        h = uicontrol(d, 'Style', 'popupmenu', 'Units', 'pixels', 'Position', [x yy w 22], ...
            'String', s, 'Value', v, 'BackgroundColor', 'w', 'Callback', @(~, ~) updateInfo());
    end

    function onFormat(initial)
        if nargin < 1, initial = false; end
        sub = get(hFmt, 'Value') == 1;
        onoff = {'off', 'on'};
        set([hCoef hBasis], 'Enable', onoff{sub + 1});
        set([hImg hOrder], 'Enable', onoff{~sub + 1});
        cf = opts.CropFactor; nm = opts.Normalize;
        if ~initial || isempty(cf)
            if sub, cf = 2; else, cf = 1; end
        end
        if ~initial || isempty(nm)
            nm = sub;
        end
        set(hCrop, 'String', num2str(cf));
        set(hNorm, 'Value', double(logical(nm)));
        updateInfo();
    end

    function updateInfo()
        if geomOnly || isempty(hInfo), return; end
        sub = get(hFmt, 'Value') == 1;
        cf = str2double(get(hCrop, 'String'));
        if ~isfinite(cf) || cf < 1, cf = 1; end
        try
            if sub
                ci = info(numIdx(get(hCoef, 'Value')));
                bi = info(numIdx(get(hBasis, 'Value')));
                s = ci.size; b = bi.size;
                if numel(s) < 4, s(end + 1:4) = 1; end
                K = s(4);
                if numel(b) == 3 && b(3) == K, T = b(2);
                elseif numel(b) == 3 && b(2) == K, T = b(3);
                elseif numel(b) == 2 && b(2) == K, T = b(1);
                else, T = NaN;
                end
                cplx = 1 + ci.complex;
                memC = prod(s) / cf ^ 2 * 4 * cplx;
                msg = sprintf(['%d slices, %dx%d -> crop %dx%d, K = %d coefficients, T = %d frames\n' ...
                    'Basis [%s] must be [slices T K] as in ReadData.m%s\n' ...
                    'Memory: ~%.2f GB after loading (peak ~%.2f GB while reading)'], ...
                    s(1), s(2), s(3), floor(s(2) / cf), floor(s(3) / cf), K, T, num2str(b), ...
                    ternary(isnan(T) || b(1) ~= s(1), '  <-- size mismatch', ''), ...
                    memC / 2^30, (ci.bytes + memC) / 2^30);
            else
                ii = info(numIdx(get(hImg, 'Value')));
                s = ii.size;
                o = lower(strrep(get(hOrder, 'String'), ' ', ''));
                lab = '';
                for k = 1:min(numel(o), numel(s))
                    lab = [lab sprintf('%s=%d ', o(k), s(k))]; %#ok<AGROW>
                end
                msg = sprintf('Array [%s] interpreted as %s\nMemory: ~%.2f GB as single (peak ~%.2f GB)', ...
                    num2str(s), lab, prod(s) / cf ^ 2 * 4 / 2^30, (ii.bytes + 2 * prod(s) * 4) / 2^30);
            end
        catch
            msg = 'Select the variables to use.';
        end
        set(hInfo, 'String', msg);
    end

    function onOK(~, ~)
        o = opts;
        if ~geomOnly
            sub = get(hFmt, 'Value') == 1;
            if sub
                o.Format = 'subspace';
                o.CoefVar = numNames{get(hCoef, 'Value')};
                o.BasisVar = numNames{get(hBasis, 'Value')};
                if strcmp(o.CoefVar, o.BasisVar)
                    errordlg('Coefficient and basis variables must differ.', 'Load', 'modal');
                    return
                end
            else
                o.Format = 'dense';
                o.ImageVar = numNames{get(hImg, 'Value')};
                o.DimOrder = lower(strrep(get(hOrder, 'String'), ' ', ''));
            end
            o.CropFactor = str2double(get(hCrop, 'String'));
            if ~isfinite(o.CropFactor) || o.CropFactor < 1
                errordlg('Crop factor must be >= 1.', 'Load', 'modal');
                return
            end
            o.Normalize = logical(get(hNorm, 'Value'));
        end
        o.AxisCodes = codes{get(hPlane, 'Value')};
        [pm, ~] = GraspData.codesToPermFlip(o.AxisCodes);
        pl = {'sagittal', 'coronal', 'axial'};
        o.Plane = pl{find(pm == 3)};
        v = [str2double(get(hVr, 'String')) str2double(get(hVc, 'String')) str2double(get(hVs, 'String'))];
        if any(~isfinite(v)) || any(v <= 0)
            errordlg('Voxel sizes must be positive numbers.', 'Load', 'modal');
            return
        end
        o.Voxel = v;
        o.dt = str2double(get(hDt, 'String'));
        o.t0 = str2double(get(hT0, 'String'));
        if ~isfinite(o.t0), o.t0 = 0; end
        out = o;
        delete(d);
    end

    function onCancel(~, ~)
        out = [];
        delete(d);
    end
end

function k = pickIdx(list, name, fallback)
k = find(strcmpi(list, name), 1);
if isempty(k), k = fallback; end
k = max(1, min(numel(list), k));
end

function s = ternary(c, a, b)
if c, s = a; else, s = b; end
end
