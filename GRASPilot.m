classdef GRASPilot < handle
%GRASPILOT 3-plane viewer and ROI analysis platform for dynamic GRASP DCE-MRI.
%
%   GRASPilot                       empty viewer (File > Open .mat)
%   GRASPilot('recon.mat')          open a reconstruction (load dialog)
%   GRASPilot(vk_all, basis_all)    subspace variables already in the workspace
%   GRASPilot(data1)                4D array from ReadData.m ([nx ny t nz])
%   GRASPilot(D)                    GRASPilotData object
%   app = GRASPilot(...)            handle for scripting
%
%   Name/value options: 'AxisCodes' (direction of increasing rows, cols,
%   slices; default 'LPS' = scanner x,y,z order of the GRASP reconstructions,
%   'PLS' = data1 slices already radiological), 'Voxel' [row col slice] mm,
%   'dt' s/frame, 't0' s, 'DimOrder' (dense input,
%   default 'xytz'), 'CropFactor', 'Normalize', 'CoefVar', 'BasisVar',
%   'ImageVar', 'BinSize' (average n consecutive frames), 'NoDialog'
%   (true skips the load dialog).
%
%   Mouse:  left = active tool, right-drag = window/level, shift- or
%           middle-drag = pan, wheel = slice, Ctrl+wheel = zoom.
%   Keys:   Left/Right frame -/+1 (Shift: 10), Up/Down slice, Home/End,
%           Space play, N/C/B/E tools, [ ] brush size, D toggle 3D,
%           P preview, S subtraction view, G grow from crosshair,
%           M parameter map overlay,
%           Ctrl+Z undo, Ctrl+Y or Ctrl+Shift+Z redo.
%
%   No toolboxes are required. Nothing is written except files you export.

    properties
        Data = [];
        Fig = [];
        H = struct();
        V = [];
        Tic = struct();
        Pos = [1 1 1];
        Frame = 1;
        CLim = [0 1];
        Tool = 'nav';
        BrushR = 3;
        Use3D = false;
        Exclusive = true;
        Constrain = false;
        Preview = false;
        ThrSrc = 2;             % index into ThrKeys (default: subtraction)
        SubFrames = [1 1];      % subtraction map = frame SubFrames(1) - frame SubFrames(2)
        ViewMode = 'frame';     % 'frame' or 'sub' (views show the subtraction map)
        CLimStore = struct('frame', [], 'sub', []);
        GrowR = 40;             % mm, max distance for grow-from-crosshair
        ParamKey = 'PE';        % parameter shown as a colour map
        ParamOn = false;        % parameter map overlay on the images
        ParamCmap = 'jet';
        ParamCLims = struct();  % colour range per parameter
        ParamAlpha = 0.65;
        ParamHideBelow = true;
        ParamHeadOnly = true;
        RefROIName = '';        % reference ROI for /SSS normalization ('' = first ROI named SSS*)
        IncludeRefStats = false;
        PW = struct();          % parameter window handles
        MW = struct();          % motion window handles
        MotionStop = false;     % stop request for the motion estimate
        CB = struct();          % colour bar in the sagittal panel
        ParamCache = struct('key', '', 'V', []);
        RefCache = [];
        HeadCache = struct('ver', -1, 'M', []);
        ThrLo = -Inf;
        ThrHi = Inf;
        BaseFrames = [1 5];
        Alpha = 0.35;
        ShowOutline = true;
        ShowCross = true;
        ShowCursorTIC = true;
        TicY = 1;
        TicStat = 1;
        TicX = 1;
        ROIs = struct('Name', {}, 'Color', {}, 'Visible', {}, 'Dirty', {}, 'Stats', {}, 'Count', {});
        Masks = {};
        ActiveROI = 0;
        UndoStack = {};
        RedoStack = {};
        Pending = [];
        Drag = struct('type', '');
        ActiveView = 1;
        ViewRect = zeros(3, 4);
        Layout = struct('h', {}, 'panel', {}, 'x', {}, 'top', {}, 'w', {}, 'ht', {});
        PlayTimer = [];
        Modified = false;
        LastDir = '';
        CursorTC = [];
        ViewImg = cell(1, 3);
        FrameImg = cell(1, 3);
        OvlRGB = cell(1, 3);
        OvlA = cell(1, 3);
        PrevA = cell(1, 3);
        ParamRGB = cell(1, 3);
        ParamA = cell(1, 3);
        LoadOpts = struct();
        Busy = false;
        Quiet = false;          % true: dialogs auto-answer with defaults (scripting/tests)
    end

    properties (Constant)
        Version = '1.0';
        ViewAxis = [3 2 1];
        ViewNames = {'Axial', 'Coronal', 'Sagittal'};
        ViewLabels = {{'R', 'L', 'A', 'P'}, {'R', 'L', 'S', 'I'}, {'A', 'P', 'S', 'I'}};
        ThrKeys = {'frame', 'sub', 'frameSub', 'frameRE', 'peakEnh', 'peakRE'};
        ThrNames = {'Signal at current frame', 'Subtraction: post - pre frame', ...
            'Current frame - baseline', 'Current frame rel. enhancement (%)', ...
            'Peak - baseline (max after baseline)', 'Peak rel. enhancement (%)'};
        Palette = [1.00 0.85 0.10; 0.20 0.90 1.00; 1.00 0.35 0.95; 0.40 1.00 0.35; ...
            1.00 0.60 0.15; 0.55 0.55 1.00; 0.60 1.00 0.85; 1.00 0.55 0.55; 0.85 0.65 1.00];
        SSSColor = [1.00 0.22 0.22];
        ParamKeys = {'PE', 'PE_SSS', 'RE', 'WIS', 'WIS_SSS', 'WOS', 'WOS_SSS', 'AUC', 'AUC_SSS', 'TTP', ...
            'BAT', 'T90', 'MITR', 'REE', 'SER', 'WOP', 'PLAT', 'CT'};
        ParamLabels = {'Peak enhancement', 'Peak enhancement / SSS', 'Relative peak enhancement', ...
            'Wash-in slope', 'Wash-in slope / SSS', 'Washout slope', 'Washout slope / SSS', ...
            'Area under curve', 'Area under curve / SSS', 'Time to peak', ...
            'Onset (bolus arrival) time', 'Time to 90% of peak', 'Max intensity-time ratio', ...
            'Rel. enhancement, early frame', 'Signal enhancement ratio', 'Washout percentage', ...
            'Plateau enhancement', 'Curve type (1 pers, 2 plat, 3 wash)'};
        CurveTypeNames = {'Curve type I (persistent) fraction', 'Curve type II (plateau) fraction', ...
            'Curve type III (washout) fraction'};
        CmapNames = {'jet', 'turbo', 'hot', 'gray'};
        CmapLabels = {'Jet (rainbow)', 'Turbo (rainbow)', 'Hot', 'Gray'};
        Col = struct('bg', [0.11 0.11 0.13], 'panel', [0.16 0.16 0.18], 'fg', [0.90 0.90 0.92], ...
            'accent', [0.40 0.72 1.00], 'dim', [0.62 0.62 0.66], 'cross', [0.25 0.75 1.00], ...
            'sel', [0.30 0.55 0.85], 'prev', [0.20 1.00 0.30]);
    end

    % =====================================================================
    methods

        function obj = GRASPilot(varargin)
            obj.LoadOpts = GRASPilot.defaultLoadOpts();
            obj.loadPrefs();
            obj.buildUI();
            set(obj.Fig, 'Visible', 'on');
            if nargin > 0
                try
                    obj.openInput(varargin{:});
                catch err
                    obj.showError(err, 'Could not load data');
                end
            end
        end

        % ------------------------------------------------------------ input
        function openInput(obj, varargin)
            a = varargin;
            if isa(a{1}, 'GRASPilotData')
                obj.setData(a{1});
                return
            end
            if ischar(a{1})
                opts = obj.mergeOpts(a(2:end));
                obj.openFile(a{1}, opts);
                return
            end
            if ~isnumeric(a{1})
                error('GRASPilot:input', 'Unrecognized input.');
            end
            sub = numel(a) >= 2 && isnumeric(a{2});
            if sub, nv = a(3:end); else, nv = a(2:end); end
            opts = obj.mergeOpts(nv);
            if ~opts.NoDialog
                g = graspilot_load_dialog('', opts);
                if isempty(g), return; end
                opts = g;
            end
            wb = obj.wbOpen('Preparing data...');
            try
                prog = @(f, m) obj.wbUpdate(wb, f, m);
                if sub
                    D = GRASPilotData.fromSubspace(a{1}, a{2}, 'CropFactor', opts.CropFactor, ...
                        'Normalize', opts.Normalize, 'Progress', prog);
                    D.SourceVars = 'vk_all + basis_all (workspace)';
                else
                    D = GRASPilotData.fromDense(a{1}, opts.DimOrder, 'CropFactor', opts.CropFactor, ...
                        'Normalize', opts.Normalize, 'Progress', prog);
                    D.SourceVars = sprintf('4D array (%s, workspace)', opts.DimOrder);
                end
            catch err
                obj.wbClose(wb);
                rethrow(err);
            end
            obj.wbClose(wb);
            D.SourceFile = '(workspace)';
            obj.applyLoadGeometry(D, opts);
            obj.setData(D);
        end

        function openFile(obj, file, opts)
            if nargin < 3, opts = obj.LoadOpts; end
            if ~opts.NoDialog
                o = graspilot_load_dialog(file, opts);
                if isempty(o), return; end
                opts = o;
            end
            wb = obj.wbOpen('Loading...');
            try
                D = GRASPilotData.fromFile(file, 'Format', opts.Format, 'CoefVar', opts.CoefVar, ...
                    'BasisVar', opts.BasisVar, 'ImageVar', opts.ImageVar, 'DimOrder', opts.DimOrder, ...
                    'CropFactor', opts.CropFactor, 'Normalize', opts.Normalize, ...
                    'Progress', @(f, m) obj.wbUpdate(wb, f, m));
            catch err
                obj.wbClose(wb);
                rethrow(err);
            end
            obj.wbClose(wb);
            obj.applyLoadGeometry(D, opts);
            obj.LastDir = fileparts(file);
            obj.setData(D);
        end

        function applyLoadGeometry(obj, D, opts)
            if isfield(opts, 'AxisCodes') && ~isempty(opts.AxisCodes)
                D.setAxisCodes(opts.AxisCodes);
            else
                D.setOrientation(opts.Plane, opts.SwapInPlane, opts.Flip);
            end
            D.VoxNative = opts.Voxel(:)';
            D.dt = opts.dt;
            D.t0 = opts.t0;
            if isfield(opts, 'BinSize') && ~isempty(opts.BinSize) && opts.BinSize > 1
                D.setBinning(opts.BinSize);
            end
            obj.LoadOpts.AxisCodes = D.AxisCodes;
            obj.LoadOpts.Plane = D.Plane;
            obj.LoadOpts.Voxel = opts.Voxel(:)';
            obj.LoadOpts.dt = opts.dt;
            obj.LoadOpts.t0 = opts.t0;
            obj.savePrefs();
        end

        function setData(obj, D)
            keep = false;
            if ~isempty(obj.Data) && ~isempty(obj.Masks) && any(cellfun(@(m) any(m(:)), obj.Masks)) ...
                    && isequal(obj.Data.NativeSize, D.NativeSize)
                b = obj.ask('Keep the current ROIs for the new series (same matrix size)?', 'ROIs', ...
                    {'Keep', 'Discard'}, 'Keep');
                keep = strcmp(b, 'Keep');
            end
            if keep
                nat = cellfun(@(m) obj.Data.canon2native(m), obj.Masks, 'UniformOutput', false);
            end
            obj.stopPlay();
            obj.Data = D;
            if keep
                obj.Masks = cellfun(@(m) D.native2canon(m), nat, 'UniformOutput', false);
                for r = 1:numel(obj.ROIs), obj.ROIs(r).Dirty = true; end
            else
                obj.ROIs = obj.ROIs([]);
                obj.Masks = {};
                obj.ActiveROI = 0;
            end
            obj.UndoStack = {}; obj.RedoStack = {};
            obj.CursorTC = [];
            n = D.canonSize();
            obj.Pos = max(1, round((n + 1) / 2));
            obj.Frame = max(1, min(D.NT, round(0.75 * D.NT)));
            obj.BaseFrames = [1 max(1, min(5, D.NT))];
            set(obj.H.base1, 'String', num2str(obj.BaseFrames(1)));
            set(obj.H.base2, 'String', num2str(obj.BaseFrames(2)));
            obj.setSubFrames(D.NT, 1, false);
            set(obj.H.binEdit, 'String', num2str(D.BinSize));
            obj.ViewMode = 'frame';
            set(obj.H.btnSub, 'Value', 0);
            obj.CLimStore = struct('frame', [], 'sub', []);
            obj.ParamOn = false;
            obj.ParamCLims = struct();
            obj.ParamCache = struct('key', '', 'V', []);
            obj.RefCache = [];
            obj.updateColorbar();
            obj.refreshMotionWindow();
            obj.setupTimeSlider();
            obj.resetViewsGeometry();
            obj.autoContrast(false);
            obj.updateDataInfo();
            obj.refreshROIList();
            obj.refreshAll();
            obj.updateThrSliderRange();
            obj.Modified = keep;
            [~, nm] = fileparts(D.SourceFile);
            set(obj.Fig, 'Name', ['GRASPilot - ' nm]);
            obj.setStatus(sprintf('Loaded %s. %s', nm, D.describe()));
            obj.checkOrientation();
        end

        function ok = checkOrientation(obj)
            % Heads are longer front-to-back than left-to-right: warn when the
            % canonical L-R extent exceeds the A-P extent (view rotated 90 deg).
            ok = true;
            V = obj.Data.canonVolume(obj.Frame);
            step = max(1, floor(numel(V) / 200000));
            sv = sort(double(V(1:step:end)));
            m = V > 0.15 * sv(max(1, round(0.99 * numel(sv))));
            if nnz(m) < 100, return; end
            vox = obj.Data.canonVox();
            ext = zeros(1, 2);
            for a = 1:2
                prof = squeeze(any(any(m, 3), 3 - a));
                idx = find(prof);
                ext(a) = (idx(end) - idx(1) + 1) * vox(a);
            end
            if ext(1) > 1.08 * ext(2)
                ok = false;
                obj.setStatus(sprintf(['Orientation check: the head is wider left-right (%.0f mm) than ' ...
                    'front-back (%.0f mm), so the views are probably rotated 90 degrees. ' ...
                    'Use View > Geometry (axis codes, Rotate 90).'], ext(1), ext(2)));
            end
        end

        % ------------------------------------------------------------ frame / position
        function setFrame(obj, f, force)
            if isempty(obj.Data), return; end
            if nargin < 3, force = false; end
            if isnan(f), f = obj.Frame; end
            f = max(1, min(obj.Data.NT, round(f)));
            if f == obj.Frame && ~force
                obj.updateTimeWidgets();
                return
            end
            obj.Frame = f;
            obj.updateTimeWidgets();
            for v = 1:3
                obj.refreshView(v, 'image');
            end
            obj.updateTicFrameLine();
            if obj.thrUsesFrame(), obj.updateThrInfo(); end
        end

        function setPos(obj, c)
            n = obj.Data.canonSize();
            c = max(1, min(n, round(c(:)')));
            old = obj.Pos;
            obj.Pos = c;
            for v = 1:3
                a = obj.ViewAxis(v);
                if old(a) ~= c(a)
                    obj.refreshView(v, 'all');
                else
                    obj.refreshView(v, 'cross');
                end
            end
        end

        function moveSlice(obj, v, d)
            if isempty(obj.Data), return; end
            c = obj.Pos;
            a = obj.ViewAxis(v);
            c(a) = c(a) + d;
            obj.setPos(c);
        end

        % ------------------------------------------------------------ tools
        function setTool(obj, name)
            obj.Tool = name;
            fns = fieldnames(obj.H.tool);
            dflt = get(0, 'defaultUicontrolBackgroundColor');
            for i = 1:numel(fns)
                on = strcmp(fns{i}, name);
                bgc = dflt;
                if on, bgc = obj.Col.sel; end
                set(obj.H.tool.(fns{i}), 'Value', double(on), 'BackgroundColor', bgc);
            end
            obj.hideCursors();
            hints = struct('nav', 'Navigate: click to move the crosshair; the cursor voxel curve is shown in the TIC plot.', ...
                'circle', 'Circle: press at the centre and drag to the radius (sphere when 3D is on).', ...
                'brush', 'Brush: paint into the active ROI ([ and ] change the radius).', ...
                'erase', 'Eraser: remove voxels from the active ROI.');
            obj.setStatus(hints.(name));
        end

        function setBrushR(obj, r)
            if ~isfinite(r), r = obj.BrushR; end
            r = min(50, max(0.25, r));
            obj.BrushR = r;
            set(obj.H.rSlider, 'Value', min(25, max(0.5, r)));
            set(obj.H.rEdit, 'String', sprintf('%.2g', r));
        end

        function set3D(obj, on)
            obj.Use3D = logical(on);
            set(obj.H.chk3D, 'Value', double(obj.Use3D));
        end

        % ------------------------------------------------------------ ROIs
        function r = addROI(obj, name, color)
            r = 0;
            if ~obj.requireData(), return; end
            if nargin < 3 || isempty(color), color = obj.nextColor(); end
            r = numel(obj.ROIs) + 1;
            obj.ROIs(r) = struct('Name', name, 'Color', color(:)', 'Visible', true, 'Dirty', true, ...
                'Stats', [], 'Count', 0);
            obj.Masks{r} = false(obj.Data.canonSize());
            obj.ActiveROI = r;
            obj.refreshROIList();
            obj.setActiveROI(r);
            obj.setStatus(sprintf('Added ROI "%s". Draw with Circle/Brush; Eraser removes.', name));
        end

        function drawReady(obj)
            if strcmp(obj.Tool, 'nav'), obj.setTool('circle'); end
        end

        function addSSS(obj)
            base = 'SSS (control)';
            nm = base; k = 2;
            while any(strcmp({obj.ROIs.Name}, nm))
                nm = sprintf('%s %d', base, k); k = k + 1;
            end
            obj.addROI(nm, obj.SSSColor);
            obj.drawReady();
        end

        function addTumor(obj)
            k = 1;
            while any(strcmp({obj.ROIs.Name}, sprintf('Tumor %d', k))), k = k + 1; end
            obj.addROI(sprintf('Tumor %d', k));
            obj.drawReady();
        end

        function addCustom(obj)
            if ~obj.requireData(), return; end
            a = inputdlg('ROI name:', 'New ROI', 1, {sprintf('ROI %d', numel(obj.ROIs) + 1)});
            if isempty(a) || isempty(strtrim(a{1})), return; end
            obj.addROI(strtrim(a{1}));
            obj.drawReady();
        end

        function setActiveROI(obj, r)
            n = numel(obj.ROIs);
            if n == 0
                obj.ActiveROI = 0;
                return
            end
            r = max(1, min(n, round(r)));
            obj.ActiveROI = r;
            set(obj.H.roiPop, 'Value', r);
            set(obj.H.swatch, 'BackgroundColor', obj.ROIs(r).Color);
            obj.updateThrInfo();
        end

        function deleteROI(obj)
            r = obj.ActiveROI;
            if r == 0, return; end
            if nnz(obj.Masks{r}) > 0
                b = obj.ask(sprintf('Delete ROI "%s" (%d voxels)?', obj.ROIs(r).Name, nnz(obj.Masks{r})), ...
                    'Delete ROI', {'Delete', 'Cancel'}, 'Delete');
                if ~strcmp(b, 'Delete'), return; end
            end
            obj.ROIs(r) = [];
            obj.Masks(r) = [];
            obj.UndoStack = GRASPilot.reindexStack(obj.UndoStack, r);
            obj.RedoStack = GRASPilot.reindexStack(obj.RedoStack, r);
            obj.ActiveROI = min(r, numel(obj.ROIs));
            obj.Modified = true;
            obj.refreshROIList();
            obj.setActiveROI(obj.ActiveROI);
            obj.refreshOverlays();
            obj.refreshTIC();
        end

        function renameROI(obj)
            r = obj.ActiveROI;
            if r == 0, return; end
            a = inputdlg('ROI name:', 'Rename ROI', 1, {obj.ROIs(r).Name});
            if isempty(a) || isempty(strtrim(a{1})), return; end
            obj.ROIs(r).Name = strtrim(a{1});
            obj.refreshROIList();
            obj.refreshTIC();
        end

        function recolorROI(obj)
            r = obj.ActiveROI;
            if r == 0, return; end
            c = uisetcolor(obj.ROIs(r).Color, 'ROI color');
            if numel(c) ~= 3, return; end
            obj.ROIs(r).Color = c;
            obj.refreshROIList();
            obj.setActiveROI(r);
            obj.refreshOverlays();
            obj.refreshTIC();
        end

        function toggleROI(obj)
            r = obj.ActiveROI;
            if r == 0, return; end
            obj.ROIs(r).Visible = ~obj.ROIs(r).Visible;
            obj.refreshROIList();
            obj.refreshOverlays();
            obj.refreshTIC();
        end

        function clearROI(obj, scope)
            r = obj.ActiveROI;
            if r == 0 || isempty(obj.Data), return; end
            if strcmp(scope, 'slice')
                lin = obj.sliceLin(obj.ActiveView);
            else
                lin = find(obj.Masks{r});
            end
            obj.beginAction();
            obj.setVoxels(r, lin, false);
            obj.commitAction('Clear');
            obj.afterMaskChange();
        end

        % ------------------------------------------------------------ drawing primitives
        function paintShape(obj, v, X, Y, rmm, erase)
            r = obj.ActiveROI;
            if r == 0, return; end
            if obj.Use3D
                lin = obj.sphereIndices(obj.screen2canonCont(v, X, Y), rmm);
            else
                lin = obj.diskIndices(v, X, Y, rmm);
            end
            if isempty(lin), return; end
            if erase
                obj.setVoxels(r, lin, false);
            else
                if obj.Constrain
                    lin = lin(obj.inWindow(obj.thrValues(lin)));
                end
                obj.addVoxels(r, lin);
            end
        end

        function addVoxels(obj, r, lin)
            obj.setVoxels(r, lin, true);
            if obj.Exclusive
                for q = 1:numel(obj.Masks)
                    if q ~= r, obj.setVoxels(q, lin, false); end
                end
            end
        end

        function setVoxels(obj, r, lin, val)
            if isempty(lin), return; end
            lin = lin(obj.Masks{r}(lin) ~= val);
            if isempty(lin), return; end
            obj.Masks{r}(lin) = val;
            if ~isstruct(obj.Pending)
                obj.Pending = struct('roi', {}, 'idx', {}, 'val', {});
            end
            obj.Pending(end + 1) = struct('roi', r, 'idx', uint32(lin(:)), 'val', logical(val));
            obj.ROIs(r).Dirty = true;
        end

        function beginAction(obj)
            obj.Pending = struct('roi', {}, 'idx', {}, 'val', {});
        end

        function commitAction(obj, label)
            P = obj.Pending;
            obj.Pending = [];
            if ~isstruct(P) || isempty(P), return; end
            keys = 2 * [P.roi] + double([P.val]);
            uk = unique(keys);
            A = struct('roi', {}, 'idx', {}, 'val', {});
            for k = uk
                sel = P(keys == k);
                A(end + 1) = struct('roi', sel(1).roi, 'idx', vertcat(sel.idx), 'val', sel(1).val); %#ok<AGROW>
            end
            obj.UndoStack{end + 1} = struct('label', label, 'changes', A);
            if numel(obj.UndoStack) > 100, obj.UndoStack(1) = []; end
            obj.RedoStack = {};
            obj.Modified = true;
        end

        function undo(obj)
            if isempty(obj.UndoStack), obj.setStatus('Nothing to undo.'); return; end
            a = obj.UndoStack{end};
            obj.UndoStack(end) = [];
            obj.applyChanges(a.changes, true);
            obj.RedoStack{end + 1} = a;
            obj.afterMaskChange();
            obj.setStatus(['Undo: ' a.label]);
        end

        function redo(obj)
            if isempty(obj.RedoStack), obj.setStatus('Nothing to redo.'); return; end
            a = obj.RedoStack{end};
            obj.RedoStack(end) = [];
            obj.applyChanges(a.changes, false);
            obj.UndoStack{end + 1} = a;
            obj.afterMaskChange();
            obj.setStatus(['Redo: ' a.label]);
        end

        % ------------------------------------------------------------ threshold
        function setThrSource(obj, k)
            if ischar(k), k = find(strcmp(obj.ThrKeys, k), 1); end
            obj.ThrSrc = k;
            set(obj.H.thrSrc, 'Value', k);
            if obj.thrNeedsMaps() && (obj.Preview || obj.Constrain)
                if ~obj.ensureMaps(false), return; end
            end
            obj.updateThrSliderRange();
            obj.refreshPreviews();
            obj.updateThrInfo();
        end

        function setThreshold(obj, lo, hi)
            if nargin < 3, hi = obj.ThrHi; end
            if isnan(lo), lo = -Inf; end
            if isnan(hi), hi = Inf; end
            obj.ThrLo = lo; obj.ThrHi = hi;
            set(obj.H.thrLo, 'String', GRASPilot.num2s(lo));
            set(obj.H.thrHi, 'String', GRASPilot.num2s(hi));
            if isfinite(lo)
                mn = get(obj.H.thrSlider, 'Min'); mx = get(obj.H.thrSlider, 'Max');
                set(obj.H.thrSlider, 'Value', min(mx, max(mn, lo)));
            end
            obj.refreshPreviews();
            obj.updateThrInfo();
        end

        function setPreview(obj, on)
            on = logical(on);
            if on && obj.thrNeedsMaps() && ~obj.ensureMaps(false), on = false; end
            obj.Preview = on;
            set(obj.H.chkPrev, 'Value', double(on));
            obj.refreshPreviews();
        end

        function setConstrain(obj, on)
            on = logical(on);
            if on && obj.thrNeedsMaps() && ~obj.ensureMaps(false), on = false; end
            obj.Constrain = on;
            set(obj.H.chkCons, 'Value', double(on));
        end

        function n = applyThreshold(obj, scope)
            n = 0;
            r = obj.ActiveROI;
            if r == 0 || ~any(obj.Masks{r}(:))
                obj.inform('Draw a region in the active ROI first, then apply the threshold to keep only voxels inside the window.', 'Threshold');
                return
            end
            if ~isfinite(obj.ThrLo) && ~isfinite(obj.ThrHi)
                obj.inform('Set a Low and/or High threshold first (or press Otsu).', 'Threshold');
                return
            end
            if ~obj.ensureThresholdReady(), return; end
            if strcmp(scope, 'slice')
                lin = obj.sliceLin(obj.ActiveView);
                lin = lin(obj.Masks{r}(lin));
            else
                lin = find(obj.Masks{r});
            end
            if isempty(lin)
                obj.setStatus('The active ROI has no voxels on this slice.');
                return
            end
            keep = obj.inWindow(obj.thrValues(lin));
            obj.beginAction();
            obj.setVoxels(r, lin(~keep), false);
            obj.commitAction('Threshold');
            obj.afterMaskChange();
            n = nnz(keep);
            obj.setStatus(sprintf('Threshold kept %d of %d voxels of %s (%s).', n, numel(lin), ...
                obj.ROIs(r).Name, scope));
        end

        function th = suggestThreshold(obj)
            th = NaN;
            if ~obj.requireData() || ~obj.ensureThresholdReady(), return; end
            r = obj.ActiveROI;
            if r > 0 && any(obj.Masks{r}(:))
                vals = obj.thrValues(find(obj.Masks{r})); %#ok<FNDSB>
                src = sprintf('ROI %s', obj.ROIs(r).Name);
            else
                M = obj.thrSliceMap(obj.ActiveView);
                vals = M(:);
                src = 'current slice';
            end
            th = GRASPilot.otsu(double(vals));
            if ~isfinite(th), return; end
            obj.setThreshold(th, Inf);
            obj.setStatus(sprintf('Otsu threshold %.4g from %s (%s). Use Keep: 3D to apply.', th, src, ...
                obj.ThrNames{obj.ThrSrc}));
        end

        function setSubFrames(obj, post, pre, refresh)
            % subtraction map = frame post - frame pre (default last - first)
            if nargin < 4, refresh = true; end
            T = obj.Data.NT;
            if ~isfinite(post), post = T; end
            if ~isfinite(pre), pre = 1; end
            post = max(1, min(T, round(post)));
            pre = max(1, min(T, round(pre)));
            obj.SubFrames = [post pre];
            set(obj.H.subPost, 'String', num2str(post));
            set(obj.H.subPre, 'String', num2str(pre));
            if ~refresh, return; end
            if strcmp(obj.ViewMode, 'sub')
                obj.CLimStore.sub = [];
                obj.autoContrast(true);
            end
            obj.refreshAll();
            obj.updateThrSliderRange();
        end

        function A = subVolume(obj)
            A = obj.Data.canonSubtraction(obj.SubFrames(1), obj.SubFrames(2));
        end

        function setViewMode(obj, mode)
            if ~obj.requireData(), set(obj.H.btnSub, 'Value', 0); return; end
            if strcmp(mode, obj.ViewMode), set(obj.H.btnSub, 'Value', double(strcmp(mode, 'sub'))); return; end
            obj.CLimStore.(obj.ViewMode) = obj.CLim;
            obj.ViewMode = mode;
            set(obj.H.btnSub, 'Value', double(strcmp(mode, 'sub')));
            for v = 1:3
                obj.refreshView(v, 'all');
            end
            if isempty(obj.CLimStore.(mode))
                obj.autoContrast(true);
            else
                obj.setCLim(obj.CLimStore.(mode));
            end
            if strcmp(mode, 'sub')
                obj.setStatus(sprintf(['Showing subtraction map: frame %d - frame %d. Threshold it (Map = ' ...
                    'Subtraction, Otsu/Low), then Grow from crosshair or Keep: 3D.'], obj.SubFrames));
            else
                obj.setStatus('Showing the current frame.');
            end
        end

        function n = growFromCrosshair(obj)
            % 3D (6-connected) region of in-window voxels connected to the
            % crosshair, within GrowR mm, added to the active ROI.
            n = 0;
            if ~obj.requireData() || ~obj.ensureThresholdReady(), return; end
            if ~isfinite(obj.ThrLo) && ~isfinite(obj.ThrHi)
                obj.inform('Set a threshold window first (Low/High or Otsu).', 'Grow from crosshair');
                return
            end
            nn = obj.Data.canonSize();
            vox = obj.Data.canonVox();
            seed = obj.Pos;
            rv = ceil(obj.GrowR ./ vox);
            lo = max(1, seed - rv);
            hi = min(nn, seed + rv);
            [I, J, K] = ndgrid(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));
            lin = sub2ind(nn, I(:), J(:), K(:));
            M = reshape(obj.inWindow(obj.thrValues(lin)), hi - lo + 1);
            s0 = seed - lo + 1;
            if ~M(s0(1), s0(2), s0(3))
                obj.inform(['The crosshair voxel is outside the threshold window. Put the crosshair on ' ...
                    'enhancing tumor (Navigate tool) and try again.'], 'Grow from crosshair');
                return
            end
            R = GRASPilot.floodFill3D(M, s0);
            add = lin(R(:));
            if ~obj.ensureActiveROI(), return; end
            r = obj.ActiveROI;
            obj.beginAction();
            obj.addVoxels(r, add);
            obj.commitAction('Grow');
            obj.afterMaskChange();
            n = numel(add);
            obj.setStatus(sprintf('Grew %d voxels (%.3f mL) connected to the crosshair into %s (%s, max %g mm).', ...
                n, n * prod(vox) / 1000, obj.ROIs(r).Name, obj.ThrNames{obj.ThrSrc}, obj.GrowR));
        end

        function n = windowToROI(obj)
            % add every voxel inside the threshold window to the active ROI
            n = 0;
            if ~obj.requireData() || ~obj.ensureThresholdReady(), return; end
            if ~isfinite(obj.ThrLo) && ~isfinite(obj.ThrHi)
                obj.inform('Set a threshold window first (Low/High or Otsu).', 'All in window -> ROI');
                return
            end
            lin = find(obj.inWindow(obj.thrVolume()));
            if isempty(lin)
                obj.inform('No voxels are inside the threshold window.', 'All in window -> ROI');
                return
            end
            frac = numel(lin) / prod(obj.Data.canonSize());
            if frac > 0.05
                b = obj.ask(sprintf(['This adds %d voxels (%.1f%% of the volume), which usually includes ' ...
                    'vessels and other enhancing tissue. Continue?'], numel(lin), 100 * frac), ...
                    'All in window -> ROI', {'Add', 'Cancel'}, 'Add');
                if ~strcmp(b, 'Add'), return; end
            end
            if ~obj.ensureActiveROI(), return; end
            r = obj.ActiveROI;
            obj.beginAction();
            obj.addVoxels(r, lin);
            obj.commitAction('Window -> ROI');
            obj.afterMaskChange();
            n = numel(lin);
            obj.setStatus(sprintf('Added %d voxels inside the window to %s.', n, obj.ROIs(r).Name));
        end

        function applyBinning(obj, n, keepPartial)
            % average n consecutive frames (signal intensity); n = 1 restores the original series
            if ~obj.requireData(), return; end
            D = obj.Data;
            if nargin < 3 || isempty(keepPartial), keepPartial = D.BinKeepPartial; end
            if ~isfinite(n), n = D.BinSize; end
            n = max(1, min(D.NTRaw, round(n)));
            rawCur = D.BinIdx{obj.Frame}(1);
            rawBase = [D.BinIdx{obj.BaseFrames(1)}(1), D.BinIdx{obj.BaseFrames(2)}(end)];
            obj.stopPlay();
            D.setBinning(n, keepPartial);
            T = D.NT;
            obj.Frame = min(T, floor((rawCur - 1) / n) + 1);
            a = min(T, floor((rawBase(1) - 1) / n) + 1);
            b = max(a, min(T, floor(rawBase(2) / n)));
            obj.BaseFrames = [a b];
            set(obj.H.base1, 'String', num2str(a));
            set(obj.H.base2, 'String', num2str(b));
            obj.setSubFrames(T, 1, false);
            obj.CLimStore.sub = [];
            set(obj.H.binEdit, 'String', num2str(n));
            for r = 1:numel(obj.ROIs), obj.ROIs(r).Dirty = true; end
            obj.CursorTC = [];
            if obj.Preview && obj.thrNeedsMaps(), obj.ensureMaps(false); end
            obj.setupTimeSlider();
            obj.updateDataInfo();
            if strcmp(obj.ViewMode, 'sub'), obj.autoContrast(true); end
            obj.refreshAll();
            obj.updateThrSliderRange();
            dr = D.droppedFrames();
            note = '';
            if numel(dr) == 1
                note = sprintf(' Raw frame %d dropped (incomplete bin).', dr);
            elseif ~isempty(dr)
                note = sprintf(' Raw frames %d-%d dropped (incomplete bin).', dr(1), dr(end));
            end
            if D.hasTiming()
                ts = sprintf(', %.3g s/frame', D.effectiveDt());
            else
                ts = '';
            end
            obj.setStatus(sprintf('Temporal binning: %d raw frames -> %d frames of %d averaged%s.%s', ...
                D.NTRaw, T, n, ts, note));
        end

        function fillHoles(obj)
            r = obj.ActiveROI;
            if r == 0 || ~any(obj.Masks{r}(:)), return; end
            n = obj.Data.canonSize();
            a = obj.ViewAxis(obj.ActiveView);
            oa = setdiff(1:3, a);
            M = obj.Masks{r};
            add = zeros(0, 1);
            for q = 1:n(a)
                S = GRASPilot.sliceOfArray(M, a, q, n);
                if ~any(S(:)), continue; end
                rows = find(any(S, 2)); cols = find(any(S, 1));
                r1 = max(1, rows(1) - 1); r2 = min(size(S, 1), rows(end) + 1);
                c1 = max(1, cols(1) - 1); c2 = min(size(S, 2), cols(end) + 1);
                sub = S(r1:r2, c1:c2);
                nw = GRASPilot.fillHoles2D(sub) & ~sub;
                if any(nw(:))
                    [ii, jj] = find(nw);
                    s3 = zeros(numel(ii), 3);
                    s3(:, a) = q;
                    s3(:, oa(1)) = ii + r1 - 1;
                    s3(:, oa(2)) = jj + c1 - 1;
                    add = [add; sub2ind(n, s3(:, 1), s3(:, 2), s3(:, 3))]; %#ok<AGROW>
                end
            end
            obj.beginAction();
            obj.addVoxels(r, add);
            obj.commitAction('Fill holes');
            obj.afterMaskChange();
            obj.setStatus(sprintf('Filled %d voxels (2D, %s slices).', numel(add), lower(obj.ViewNames{obj.ActiveView})));
        end

        function ok = syncRefWindow(obj)
            % WashIn = 'sss': the wash-in period is that of the steepest rise
            % of the reference (SSS) ROI mean curve; keep it up to date
            ok = true;
            D = obj.Data;
            if isempty(D) || ~strcmpi(D.ParamOpts.WashIn, 'sss'), return; end
            w = [];
            r = obj.refIndex();
            if r > 0 && any(obj.Masks{r}(:))
                if obj.ROIs(r).Dirty || isempty(obj.ROIs(r).Stats), obj.computeStats(r); end
                w = graspilot_ref_window(obj.ROIs(r).Stats.mean, D.paramTimes(), ...
                    obj.BaseFrames(1):obj.BaseFrames(2), D.ParamOpts);
            end
            ok = ~isempty(w);
            if ~isequal(w, D.ParamOpts.WashInFrames), D.setParamOptions('WashInFrames', w); end
        end

        function ok = ensureMaps(obj, force)
            ok = true;
            if isempty(obj.Data), ok = false; return; end
            obj.syncRefWindow();
            bf = obj.BaseFrames(1):obj.BaseFrames(2);
            if ~force && obj.Data.hasMaps(bf), return; end
            wb = obj.wbOpen('Computing enhancement maps...');
            try
                obj.Data.computeMaps(bf, @(f, m) obj.wbUpdate(wb, f, m));
            catch err
                ok = false;
                obj.wbClose(wb);
                obj.showError(err, 'Could not compute enhancement maps');
                return
            end
            obj.wbClose(wb);
            obj.ParamCache = struct('key', '', 'V', []);
            obj.RefCache = [];
            obj.updateThrSliderRange();
            obj.refreshPreviews();
            if obj.ParamOn, obj.refreshParamDisplay(); end
            obj.refreshParamWindow();
        end

        % ------------------------------------------------------------ semi-quantitative parameters
        function ok = computeParams(obj, varargin)
            % (Re)compute the curve-parameter maps. Options as in
            % GRASPilotData.setParamOptions: 'WashIn', 'WashOut', 'AUCWindow'.
            ok = false;
            if ~obj.requireData(), return; end
            if ~isempty(varargin), obj.Data.setParamOptions(varargin{:}); end
            if ~obj.syncRefWindow()
                obj.inform(['Wash-in in the SSS wash-in period needs a reference ROI in the superior ' ...
                    'sagittal sinus with a rising curve: press Auto SSS (or draw an ROI whose name ' ...
                    'contains SSS). Wash-in slopes are empty until then.'], 'Parameters');
            end
            ok = obj.ensureMaps(false);
        end

        function V = paramVolume(obj, key)
            % canonical map of a parameter (key from ParamKeys); [] if not available
            V = [];
            if isempty(obj.Data) || ~obj.Data.hasMaps(obj.BaseFrames(1):obj.BaseFrames(2)), return; end
            [name, fld, norm] = GRASPilot.paramInfo(key);
            ck = sprintf('%s_%d', key, obj.Data.MapsVersion);
            ref = [];
            if norm
                [ref, r] = obj.refParams();
                if isempty(ref), return; end
                ck = sprintf('%s_%d', ck, r);
            end
            if strcmp(obj.ParamCache.key, ck)
                V = obj.ParamCache.V;
                return
            end
            V = obj.Data.canonMap(name);
            if norm
                rv = ref.(fld);
                if isfinite(rv) && rv ~= 0
                    V = V / single(rv);
                else
                    V = NaN(size(V), 'single');
                end
            end
            obj.ParamCache = struct('key', ck, 'V', V);
        end

        function r = refIndex(obj)
            % reference ROI: RefROIName if set, else the first ROI whose name contains SSS
            r = 0;
            if isempty(obj.ROIs), return; end
            names = {obj.ROIs.Name};
            if ~isempty(obj.RefROIName)
                r = find(strcmp(names, obj.RefROIName), 1);
                if isempty(r), r = 0; end
                return
            end
            k = find(~cellfun(@isempty, regexpi(names, 'sss', 'once')), 1);
            if ~isempty(k), r = k; end
        end

        function [ref, r] = refParams(obj)
            % curve parameters of the reference (SSS) ROI mean curve
            ref = []; r = obj.refIndex();
            if r == 0 || ~any(obj.Masks{r}(:)), return; end
            if ~isempty(obj.RefCache) && obj.RefCache.r == r && obj.RefCache.ver == obj.Data.MapsVersion
                ref = obj.RefCache.P;
                return
            end
            if obj.ROIs(r).Dirty || isempty(obj.ROIs(r).Stats), obj.computeStats(r); end
            ref = graspilot_curve_params(obj.ROIs(r).Stats.mean, obj.Data.paramTimes(), ...
                obj.BaseFrames(1):obj.BaseFrames(2), obj.Data.ParamOpts);
            obj.RefCache = struct('r', r, 'ver', obj.Data.MapsVersion, 'P', ref);
        end

        function u = paramUnit(obj, key)
            [~, ~, norm] = GRASPilot.paramInfo(key);
            if norm, u = 'ratio'; return; end
            ts = 's';
            if ~obj.Data.hasTiming(), ts = 'frame'; end
            switch key
                case {'PE', 'PLAT'}, u = 'a.u.';
                case {'RE', 'REE', 'WOP'}, u = '%';
                case {'WIS', 'WOS', 'MITR'}, u = ['a.u./' ts];
                case 'AUC', u = ['a.u.*' ts];
                case {'TTP', 'BAT', 'T90'}, u = ts;
                case 'SER', u = 'ratio';
                case 'CT', u = 'type';
            end
        end

        function [vals, cval] = paramValues(obj, key, lin, cp)
            % voxel values of a parameter at canonical indices lin, and the value
            % from an ROI mean curve (cp = graspilot_curve_params output)
            vals = []; cval = NaN;
            [name, fld, norm] = GRASPilot.paramInfo(key);
            A = obj.Data.canonMap(name);
            vals = double(A(lin(:)));
            if nargin > 3 && ~isempty(cp), cval = cp.(fld); end
            if norm
                ref = obj.refParams();
                if isempty(ref) || ~isfinite(ref.(fld)) || ref.(fld) == 0
                    vals = NaN(size(vals)); cval = NaN;
                    return
                end
                vals = vals / ref.(fld);
                cval = cval / ref.(fld);
            end
        end

        function M = headMask(obj)
            if obj.HeadCache.ver == obj.Data.MapsVersion && ~isempty(obj.HeadCache.M)
                M = obj.HeadCache.M;
                return
            end
            B = obj.Data.canonMap('Baseline');
            step = max(1, floor(numel(B) / 300000));
            M = B > GRASPilot.otsu(double(B(1:step:end)));
            obj.HeadCache = struct('ver', obj.Data.MapsVersion, 'M', M);
        end

        function setParamDisplay(obj, on)
            on = logical(on);
            if on
                if ~obj.ensureMaps(false), on = false; end
            end
            if on && isempty(obj.paramVolume(obj.ParamKey))
                obj.inform(['Normalized maps need a reference ROI in the superior sagittal sinus: press ' ...
                    'Auto SSS (or draw an ROI whose name contains SSS).'], 'Parameter map');
                on = false;
            end
            obj.ParamOn = on;
            obj.refreshParamDisplay();
            obj.refreshOverlays();
            if on
                obj.setStatus(sprintf('Showing %s [%s] as a colour map (%.4g to %.4g). M toggles it.', ...
                    obj.ParamLabels{strcmp(obj.ParamKeys, obj.ParamKey)}, obj.paramUnit(obj.ParamKey), ...
                    obj.paramCLim()));
            end
        end

        function setParamKey(obj, key)
            if isnumeric(key), key = obj.ParamKeys{key}; end
            obj.ParamKey = key;
            if obj.ParamOn && isempty(obj.paramVolume(key))
                obj.setParamDisplay(false);
            end
            obj.refreshParamDisplay();
        end

        function c = paramCLim(obj, key)
            if nargin < 2, key = obj.ParamKey; end
            if isfield(obj.ParamCLims, key)
                c = obj.ParamCLims.(key);
            else
                c = obj.autoParamCLim('head', key, false);
            end
        end

        function setParamCLim(obj, c, key)
            if nargin < 3, key = obj.ParamKey; end
            if numel(c) ~= 2 || any(~isfinite(c)), return; end
            if c(2) <= c(1), c(2) = c(1) + max(abs(c(1)) * 1e-3, 1e-6); end
            obj.ParamCLims.(key) = c;
            obj.refreshParamDisplay();
        end

        function c = autoParamCLim(obj, scope, key, apply)
            % 1st-99.5th percentile inside the head ('head': keeps small
            % enhancing lesions in range) or 2nd-98th inside the tumor ROIs ('rois')
            if nargin < 3 || isempty(key), key = obj.ParamKey; end
            if nargin < 4, apply = true; end
            c = [0 1];
            V = obj.paramVolume(key);
            if isempty(V), return; end
            if strcmp(scope, 'rois')
                m = false(size(V));
                rr = obj.refIndex();
                for r = 1:numel(obj.Masks)
                    if r ~= rr, m = m | obj.Masks{r}; end
                end
                if ~any(m(:)), m = obj.headMask(); end
            else
                m = obj.headMask();
            end
            v = double(V(m));
            v = sort(v(isfinite(v)));
            pr = [0.02 0.98];
            if ~strcmp(scope, 'rois'), pr = [0.01 0.995]; end
            if numel(v) >= 2
                c = [v(max(1, round(pr(1) * numel(v)))) v(max(1, round(pr(2) * numel(v))))];
                if c(2) <= c(1), c = [v(1) v(end)]; end
                if c(2) <= c(1), c = c(1) + [-0.5 0.5]; end
            end
            if apply
                obj.ParamCLims.(key) = c;
                obj.refreshParamDisplay();
            end
        end

        function refreshParamDisplay(obj)
            if isempty(obj.Data), return; end
            for v = 1:3
                obj.computeParamOverlay(v);
                obj.composeView(v);
            end
            obj.updateColorbar();
            obj.syncParamWindow();
        end

        function [hdr, rows] = paramStatsTable(obj)
            % descriptive statistics of every parameter in each (tumor) ROI
            hdr = {'ROI', 'Parameter', 'Unit', 'N', 'Mean', 'SD', 'Median', 'Q1', 'Q3', 'P10', 'P90', ...
                'Min', 'Max', 'ROI mean curve'};
            rows = cell(0, numel(hdr));
            if isempty(obj.Data) || ~obj.Data.hasMaps(obj.BaseFrames(1):obj.BaseFrames(2)), return; end
            rr = obj.refIndex();
            for r = 1:numel(obj.ROIs)
                if r == rr && ~obj.IncludeRefStats, continue; end
                lin = find(obj.Masks{r});
                if isempty(lin), continue; end
                if obj.ROIs(r).Dirty || isempty(obj.ROIs(r).Stats), obj.computeStats(r); end
                cp = graspilot_curve_params(obj.ROIs(r).Stats.mean, obj.Data.paramTimes(), ...
                    obj.BaseFrames(1):obj.BaseFrames(2), obj.Data.ParamOpts);
                for k = 1:numel(obj.ParamKeys)
                    key = obj.ParamKeys{k};
                    [vals, cval] = obj.paramValues(key, lin, cp);
                    d = GRASPilot.describeVals(vals);
                    rows(end + 1, :) = {obj.ROIs(r).Name, obj.ParamLabels{k}, obj.paramUnit(key), d.n, d.mean, ...
                        d.sd, d.median, d.q1, d.q3, d.p10, d.p90, d.min, d.max, cval}; %#ok<AGROW>
                end
                ct = obj.paramValues('CT', lin);
                nct = nnz(isfinite(ct));
                for c = 1:3                      % voxel fractions of the curve types
                    rows(end + 1, :) = {obj.ROIs(r).Name, obj.CurveTypeNames{c}, '%', nct, ...
                        100 * nnz(ct == c) / max(1, nct), NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, ...
                        100 * double(cp.CurveType == c)}; %#ok<AGROW>
                end
            end
        end

        function files = exportParams(obj, file)
            % parameter statistics and voxel values of the (tumor) ROIs -> .xlsx / .csv
            files = {};
            if ~obj.computeParams(), return; end
            [hdr, rows] = obj.paramStatsTable();
            if isempty(rows)
                obj.inform('There are no tumor ROI voxels to export.', 'Parameter export');
                return
            end
            rr = obj.refIndex();
            blocks = struct('id', {}, 'name', {}, 'lin', {}, 'values', {});
            for r = 1:numel(obj.ROIs)
                if r == rr && ~obj.IncludeRefStats, continue; end
                lin = find(obj.Masks{r});
                if isempty(lin), continue; end
                X = zeros(numel(lin), numel(obj.ParamKeys));
                for k = 1:numel(obj.ParamKeys)
                    X(:, k) = obj.paramValues(obj.ParamKeys{k}, lin);
                end
                blocks(end + 1) = struct('id', r, 'name', obj.ROIs(r).Name, 'lin', lin, 'values', X); %#ok<AGROW>
            end
            units = cellfun(@(k) obj.paramUnit(k), obj.ParamKeys, 'UniformOutput', false);
            files = graspilot_export_params(file, obj.Data, hdr, rows, blocks, obj.ParamKeys, units, ...
                obj.paramSettings());
            obj.setStatus(['Wrote ' strjoin(files, ', ')]);
        end

        function S = paramSettings(obj)
            D = obj.Data;
            o = D.ParamOpts;
            switch lower(o.WashIn)
                case 'max', wis = 'maximum up-slope between consecutive frames';
                case 'sss'
                    wis = 'maximum up-slope within the period of most rapid SSS signal increase: no reference ROI';
                    if numel(o.WashInFrames) == 2
                        tt = D.paramTimes();
                        wis = sprintf(['maximum up-slope within the period of most rapid SSS signal increase ' ...
                            '(frames %d-%d, %.4g-%.4g on the time axis; Lee, Berger et al., Neurosurgery 2023)'], ...
                            o.WashInFrames, tt(o.WashInFrames));
                    end
                otherwise, wis = 'dS/(TP - To): onset to peak (onset = last frame below 10% of dS)';
            end
            ke = D.analysisEnd(); tt = D.paramTimes();
            S = {'Baseline frames', sprintf('%d-%d', obj.BaseFrames);
                'Motion correction', D.motionDescription();
                'Analysis endpoint', sprintf('frames 1-%d of %d (last frame at %.4g %s); later frames ignored', ...
                ke, D.NT, tt(ke), GRASPilot.ternaryStr(D.hasTiming(), 's', 'frames'));
                'Early / late frame (after end of baseline)', sprintf('%g / %g (nearest frames; Inf = last analysed frame)', ...
                o.EarlyTime, o.LateTime);
                'Curve type', sprintf(['1 persistent, 2 plateau, 3 washout: change early -> late frame ' ...
                '> +%g%%, within, < -%g%% of the early signal'], o.CurveTol, o.CurveTol);
                'Other definitions', ['T90 = time to 90% of PE; MITR = PE/TTP; SER = (S_early - S0)/(S_late - S0); ' ...
                'washout % = 100 (Smax - S_last)/Smax; plateau = mean S - S0 of the last 3 analysed frames'];
                'Wash-in slope', wis;
                'Washout slope', GRASPilot.ternaryStr(strcmpi(o.WashOut, 'endpoints'), '(S_last - S_peak)/(t_last - t_peak)', ...
                'least-squares slope from the peak (at least the last 3 frames) to the end; negative = washout');
                'AUC window (s after end of baseline)', num2str(o.AUCWindow);
                'Time origin', 'end of the baseline frames (TTP and onset are measured from it)';
                'Time unit', GRASPilot.ternaryStr(D.hasTiming(), 's', 'frame (temporal resolution not set)')};
            [ref, r] = obj.refParams();
            if isempty(ref)
                S(end + 1, :) = {'Reference ROI', 'none (normalized values are NaN)'};
            else
                S(end + 1, :) = {'Reference ROI', sprintf('%s (%d voxels); ROI mean curve', obj.ROIs(r).Name, ...
                    nnz(obj.Masks{r}))};
                S(end + 1, :) = {'Reference PE', ref.PE};
                S(end + 1, :) = {'Reference wash-in slope', ref.WashIn};
                S(end + 1, :) = {'Reference washout slope', ref.WashOut};
                S(end + 1, :) = {'Reference AUC', ref.AUC};
            end
            S(end + 1, :) = {'Normalized (/SSS)', 'parameter divided by the same parameter of the reference ROI mean curve'};
            S(end + 1, :) = {'ROI mean curve column', 'parameter computed from the ROI mean time-intensity curve'};
        end

        % ------------------------------------------------------------ display
        function setCLim(obj, c)
            if any(~isfinite(c)), return; end
            if c(2) <= c(1), c(2) = c(1) + max(abs(c(1)) * 1e-3, 1e-6); end
            obj.CLim = c;
            for v = 1:3
                set(obj.V(v).Ax, 'CLim', c);
                obj.composeView(v);
            end
            set(obj.H.wMin, 'String', sprintf('%.4g', c(1)));
            set(obj.H.wMax, 'String', sprintf('%.4g', c(2)));
        end

        function autoContrast(obj, ~)
            if isempty(obj.Data), return; end
            if strcmp(obj.ViewMode, 'sub')
                Vv = obj.subVolume();
            else
                Vv = obj.Data.canonVolume(obj.Frame);
            end
            step = max(1, floor(numel(Vv) / 200000));
            s = sort(double(Vv(1:step:end)));
            s = s(isfinite(s));
            if isempty(s), return; end
            if strcmp(obj.ViewMode, 'sub')
                hi = s(max(1, round(0.998 * numel(s))));
            else
                hi = s(max(1, round(0.995 * numel(s))));
            end
            lo = 0;
            if s(1) < 0 && ~strcmp(obj.ViewMode, 'sub'), lo = s(max(1, round(0.005 * numel(s)))); end
            if hi <= lo, hi = max(s(end), lo + 1e-6); end
            obj.setCLim([lo hi]);
        end

        function resetZoom(obj)
            for v = 1:3
                sz = obj.V(v).Size;
                set(obj.V(v).Ax, 'XLim', [0.5 sz(2) + 0.5], 'YLim', [0.5 sz(1) + 0.5]);
            end
        end

        % ------------------------------------------------------------ geometry
        function geometryDialog(obj)
            if ~obj.requireData(), return; end
            D = obj.Data;
            W = 560; Hh = 420;
            ss = get(0, 'ScreenSize');
            bg = [0.94 0.94 0.94];
            d = figure('Name', 'Geometry and orientation', 'NumberTitle', 'off', 'MenuBar', 'none', ...
                'ToolBar', 'none', 'Resize', 'off', 'Units', 'pixels', 'Color', bg, ...
                'Position', [(ss(3) - W) / 2, (ss(4) - Hh) / 2, W, Hh], 'WindowStyle', 'modal');
            t = @(x, y, w, s) uicontrol(d, 'Style', 'text', 'Units', 'pixels', 'Position', [x y w 18], ...
                'String', s, 'HorizontalAlignment', 'left', 'BackgroundColor', bg);
            e = @(x, y, w, s) uicontrol(d, 'Style', 'edit', 'Units', 'pixels', 'Position', [x y w 22], ...
                'String', s, 'BackgroundColor', 'w');
            b = @(x, y, w, s, cb) uicontrol(d, 'Style', 'pushbutton', 'Units', 'pixels', ...
                'Position', [x y w 26], 'String', s, 'Callback', cb);
            letters = 'LRPASI';
            opts = {'L  (R -> L)', 'R  (L -> R)', 'P  (A -> P)', 'A  (P -> A)', 'S  (I -> S)', 'I  (S -> I)'};
            y = Hh - 30;
            t(12, y, W - 24, 'Direction in which the stored array index increases (axis codes):');
            y = y - 30;
            lab = {'Rows', 'Columns', 'Slices'};
            hC = GRASPilot.handleArray(3);
            for k = 1:3
                t(12 + (k - 1) * 180, y, 56, lab{k});
                hC(k) = uicontrol(d, 'Style', 'popupmenu', 'Units', 'pixels', 'Position', [68 + (k - 1) * 180, y + 2, 110, 22], ...
                    'String', opts, 'Value', find(letters == D.AxisCodes(k)), 'BackgroundColor', 'w');
            end
            y = y - 36;
            t(12, y + 4, 60, 'Presets');
            b(72, y, 150, 'LPS  scanner order', @(s, ev) setCodes('LPS', true));
            b(226, y, 110, 'PLS  as stored', @(s, ev) setCodes('PLS', true));
            b(340, y, 100, 'ILP  coronal', @(s, ev) setCodes('ILP', true));
            b(444, y, 100, 'IPL  sagittal', @(s, ev) setCodes('IPL', true));
            y = y - 32;
            t(12, y + 4, 60, 'Adjust');
            b(72, y, 150, 'Rotate axial 90 CW', @(s, ev) setCodes(GRASPilotData.rotateCodesCW(getCodes()), true));
            b(226, y, 110, 'Flip L-R', @(s, ev) setCodes(GRASPilotData.flipCodes(getCodes(), 'LR'), true));
            b(340, y, 100, 'Flip A-P', @(s, ev) setCodes(GRASPilotData.flipCodes(getCodes(), 'AP'), true));
            b(444, y, 100, 'Flip S-I', @(s, ev) setCodes(GRASPilotData.flipCodes(getCodes(), 'SI'), true));
            y = y - 40;
            t(12, y, 170, 'Pixel spacing rows (mm)');
            hVr = e(185, y + 2, 60, num2str(D.VoxNative(1)));
            t(255, y, 40, 'cols');
            hVc = e(295, y + 2, 60, num2str(D.VoxNative(2)));
            t(365, y, 50, 'slices');
            hVs = e(415, y + 2, 60, num2str(D.VoxNative(3)));
            y = y - 32;
            t(12, y, 170, 'Time per frame (s)');
            hDt = e(185, y + 2, 60, GRASPilot.num2s(D.dt));
            t(255, y, 110, 'Frame 1 at (s)');
            hT0 = e(365, y + 2, 60, num2str(D.t0));
            uicontrol(d, 'Style', 'text', 'Units', 'pixels', 'Position', [12 56 W - 24 84], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', [1 1 0.92], 'String', { ...
                'Target: axial view anterior up, patient right on the left (radiological);', ...
                'coronal and sagittal views superior up. Presets and Adjust apply at once.', ...
                'Left-right cannot be judged from symmetric anatomy: confirm it once with a', ...
                'lateralized finding on the clinical images. Remembered for the next files.'});
            b(W - 300, 14, 90, 'Apply', @(s, ev) apply(false));
            uicontrol(d, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [W - 200 14 90 30], ...
                'String', 'OK', 'FontWeight', 'bold', 'Callback', @(s, ev) apply(true));
            b(W - 100, 14, 88, 'Close', @(s, ev) delete(d));

            function c = getCodes()
                c = letters([get(hC(1), 'Value') get(hC(2), 'Value') get(hC(3), 'Value')]);
            end
            function setCodes(c, doApply)
                for kk = 1:3, set(hC(kk), 'Value', find(letters == c(kk))); end
                if doApply, apply(false); end
            end
            function apply(closeAfter)
                c = getCodes();
                try
                    GRASPilotData.codesToPermFlip(c);
                catch err
                    errordlg(err.message, 'Geometry', 'modal');
                    return
                end
                vox = [str2double(get(hVr, 'String')) str2double(get(hVc, 'String')) str2double(get(hVs, 'String'))];
                if any(~isfinite(vox)) || any(vox <= 0)
                    errordlg('Voxel sizes must be positive numbers.', 'Geometry', 'modal');
                    return
                end
                t0v = str2double(get(hT0, 'String'));
                if ~isfinite(t0v), t0v = 0; end
                obj.applyGeometry(c, vox, str2double(get(hDt, 'String')), t0v);
                if closeAfter, delete(d); end
            end
        end

        function applyGeometry(obj, varargin)
            % applyGeometry(codes, vox, dt, t0)  or
            % applyGeometry(plane, swap, flips, vox, dt, t0)
            D = obj.Data;
            if numel(varargin) == 6
                [pm, fl] = GRASPilotData.presetOrientation(varargin{1}, varargin{2}, varargin{3});
                codes = GRASPilotData.permFlipToCodes(pm, fl);
                varargin = varargin(4:6);
            else
                codes = upper(varargin{1});
                varargin = varargin(2:4);
            end
            [vox, dt, t0] = deal(varargin{:});
            if ~isfinite(dt) || dt <= 0, dt = NaN; end
            if ~strcmp(D.AxisCodes, codes)
                nat = cellfun(@(m) D.canon2native(m), obj.Masks, 'UniformOutput', false);
                D.setAxisCodes(codes);
                obj.Masks = cellfun(@(m) D.native2canon(m), nat, 'UniformOutput', false);
                obj.UndoStack = {}; obj.RedoStack = {};
                obj.Pos = max(1, round((D.canonSize() + 1) / 2));
            end
            D.VoxNative = vox(:)';
            D.dt = dt;
            D.t0 = t0;
            obj.LoadOpts.AxisCodes = D.AxisCodes;
            obj.LoadOpts.Plane = D.Plane;
            obj.LoadOpts.Voxel = D.VoxNative;
            obj.LoadOpts.dt = dt;
            obj.LoadOpts.t0 = t0;
            obj.savePrefs();
            obj.resetViewsGeometry();
            obj.updateDataInfo();
            obj.updateTimeWidgets();
            obj.refreshROIList();
            obj.refreshAll();
            if obj.checkOrientation()
                obj.setStatus(sprintf('Orientation %s (rows, cols, slices point %s, %s, %s).', D.AxisCodes, ...
                    D.AxisCodes(1), D.AxisCodes(2), D.AxisCodes(3)));
            end
        end

        function setTemporalResolution(obj, dt)
            if ~obj.requireData(), return; end
            if ~isfinite(dt) || dt <= 0, dt = NaN; end
            obj.Data.dt = dt;
            obj.LoadOpts.dt = dt;
            obj.savePrefs();
            obj.updateDataInfo();
            obj.updateTimeWidgets();
            obj.refreshTIC();
        end

        function r = autoSSS(obj)
            % detect the posterior superior sagittal sinus from above the torcula
            % toward the vertex and put it in the ROI "SSS (auto)" (replaced on each run)
            r = 0;
            if ~obj.requireData(), return; end
            wb = obj.wbOpen('Detecting the superior sagittal sinus...');
            try
                [lin, info] = graspilot_detect_sss(obj.Data, 'BaseFrames', obj.BaseFrames(1):obj.BaseFrames(2), ...
                    'Progress', @(f, m) obj.wbUpdate(wb, f, m));
            catch err
                obj.wbClose(wb);
                obj.showError(err, 'SSS detection failed');
                return
            end
            obj.wbClose(wb);
            if isempty(lin)
                obj.inform(info.Message, 'Auto SSS');
                return
            end
            k = find(strcmp({obj.ROIs.Name}, 'SSS (auto)'), 1);
            if isempty(k)
                k = obj.addROI('SSS (auto)', obj.SSSColor);
            else
                obj.setActiveROI(k);
                obj.ROIs(k).Visible = true;
            end
            obj.beginAction();
            obj.setVoxels(k, find(obj.Masks{k}), false);
            obj.addVoxels(k, lin);
            obj.commitAction('Auto SSS');
            obj.setPos(round(info.Centroid));
            obj.afterMaskChange();
            if obj.ShowCursorTIC, obj.updateCursorTC(false); end
            obj.setStatus([info.Message ' Check it on the three views; edit or undo as needed.']);
            r = k;
        end

        % ------------------------------------------------------------ export / import
        function exportImageNifti(obj, file, frames)
            if nargin < 3 || isempty(frames), frames = 1:obj.Data.NT; end
            wb = obj.wbOpen('Writing NIfTI...');
            try
                graspilot_write_dce_nifti(obj.Data, file, frames, @(f, m) obj.wbUpdate(wb, f, m));
            catch err
                obj.wbClose(wb);
                obj.showError(err, 'NIfTI export failed');
                return
            end
            obj.wbClose(wb);
            obj.setStatus(sprintf('Wrote %s (%d frame(s)).', file, numel(frames)));
        end

        function exportMaskNifti(obj, file, separate)
            if nargin < 3, separate = false; end
            D = obj.Data;
            n = D.canonSize();
            nR = numel(obj.ROIs);
            if nR > 255, cls = 'uint16'; else, cls = 'uint8'; end
            L = zeros(n, cls);
            cnt = zeros(n, 'uint8');
            for r = 1:nR
                L(obj.Masks{r}) = r;
                cnt = cnt + uint8(obj.Masks{r});
            end
            hdr = struct('Dim', n, 'Vox', D.canonVox(), 'Affine', D.affine(), 'Datatype', cls, ...
                'Descrip', 'GRASPilot ROI labels', 'CalRange', [0 nR]);
            graspilot_nifti_write(file, L, hdr);
            base = regexprep(file, '(\.nii\.gz|\.nii)$', '', 'ignorecase');
            obj.writeLabelTables(base);
            if separate
                for r = 1:nR
                    if ~any(obj.Masks{r}(:)), continue; end
                    nm = regexprep(obj.ROIs(r).Name, '[^A-Za-z0-9_\-]+', '_');
                    h1 = hdr; h1.Datatype = 'uint8'; h1.CalRange = [0 1];
                    graspilot_nifti_write(sprintf('%s_roi%02d_%s.nii.gz', base, r, nm), uint8(obj.Masks{r}), h1);
                end
            end
            msg = sprintf('Wrote %s (+ label tables).', file);
            if any(cnt(:) > 1)
                msg = [msg ' Overlapping voxels were assigned to the higher ROI number.'];
            end
            obj.setStatus(msg);
        end

        function files = exportTimeseries(obj, file)
            files = {};
            names = {obj.ROIs.Name};
            wb = obj.wbOpen('Exporting time series...');
            try
                files = graspilot_export_timeseries(obj.Data, obj.Masks, names, file, ...
                    'BaseFrames', obj.BaseFrames(1):obj.BaseFrames(2), 'Progress', @(f, m) obj.wbUpdate(wb, f, m));
            catch err
                obj.wbClose(wb);
                obj.showError(err, 'Time-series export failed');
                return
            end
            obj.wbClose(wb);
            obj.setStatus(['Wrote ' strjoin(files, ', ')]);
        end

        function exportMaps(obj, folder, base)
            if ~obj.ensureMaps(false), return; end
            D = obj.Data;
            hdr = struct('Dim', D.canonSize(), 'Vox', D.canonVox(), 'Affine', D.affine(), ...
                'Datatype', 'float32', 'Descrip', 'GRASPilot map');
            if D.BinSize > 1, base = sprintf('%s_%dframes', base, D.NT); end
            nw = 0;
            for nm = {'Baseline', 'Max'}
                graspilot_nifti_write(fullfile(folder, sprintf('%s_%s.nii.gz', base, nm{1})), single(D.canonMap(nm{1})), hdr);
                nw = nw + 1;
            end
            for k = 1:numel(obj.ParamKeys)
                V = obj.paramVolume(obj.ParamKeys{k});
                if isempty(V), continue; end
                graspilot_nifti_write(fullfile(folder, sprintf('%s_%s.nii.gz', base, obj.ParamKeys{k})), single(V), hdr);
                nw = nw + 1;
            end
            obj.setStatus(sprintf('Wrote %d maps to %s (baseline frames %d-%d).', nw, folder, ...
                obj.BaseFrames(1), obj.BaseFrames(2)));
        end

        function saveROIs(obj, file)
            D = obj.Data;
            S = struct();
            S.Format = 'GRASPilot ROI session';
            S.Version = obj.Version;
            S.Saved = datestr(now, 'yyyy-mm-dd HH:MM:SS');
            S.SourceFile = D.SourceFile;
            S.NativeSize = D.NativeSize;
            S.AxisCodes = D.AxisCodes;
            S.Plane = D.Plane;
            S.SwapInPlane = D.SwapInPlane;
            S.UserFlip = D.UserFlip;
            S.VoxNative = D.VoxNative;
            S.dt = D.dt;
            S.t0 = D.t0;
            S.BaseFrames = obj.BaseFrames;
            S.BinSize = D.BinSize;
            S.BinKeepPartial = D.BinKeepPartial;
            S.SubFrames = obj.SubFrames;
            S.ParamOpts = D.ParamOpts;
            S.RefROIName = obj.RefROIName;
            S.MotionCorrected = D.MotionCorrected;
            if D.MotionCorrected
                S.Motion = D.MotionApplied;
            elseif ~isempty(D.Motion)
                S.Motion = D.Motion;
            end
            S.ROINames = {obj.ROIs.Name};
            S.ROIColors = reshape([obj.ROIs.Color], 3, [])';
            S.ROIVisible = [obj.ROIs.Visible];
            S.ROIIndexNative = cell(1, numel(obj.ROIs));
            for r = 1:numel(obj.ROIs)
                S.ROIIndexNative{r} = uint32(find(D.canon2native(obj.Masks{r})));
            end
            save(file, '-struct', 'S', '-v7');
            obj.Modified = false;
            obj.setStatus(['Saved ROIs to ' file]);
        end

        function loadROIs(obj, file)
            D = obj.Data;
            S = load(file);
            if ~isfield(S, 'ROIIndexNative')
                error('GRASPilot:session', 'Not a GRASPilot ROI file.');
            end
            if ~isequal(S.NativeSize(:)', D.NativeSize)
                error('GRASPilot:session', 'ROI file is for a %s matrix; data is %s.', ...
                    num2str(S.NativeSize(:)'), num2str(D.NativeSize));
            end
            mode = 'Replace';
            if ~isempty(obj.ROIs)
                mode = obj.ask('Replace the current ROIs or append?', 'Load ROIs', ...
                    {'Replace', 'Append', 'Cancel'}, 'Replace');
                if isempty(mode) || strcmp(mode, 'Cancel'), return; end
            end
            if isfield(S, 'AxisCodes')
                sCodes = S.AxisCodes;
            else
                [pm, fl] = GRASPilotData.presetOrientation(S.Plane, logical(S.SwapInPlane), logical(S.UserFlip(:)'));
                sCodes = GRASPilotData.permFlipToCodes(pm, fl);
            end
            geoDiff = ~strcmp(sCodes, D.AxisCodes) || any(abs(S.VoxNative(:)' - D.VoxNative) > 1e-6) ...
                || ~isequaln(S.dt, D.dt);
            sBin = 1;
            if isfield(S, 'BinSize'), sBin = S.BinSize; end
            binDiff = sBin ~= D.BinSize;
            applied = false;
            if geoDiff || binDiff
                b = obj.ask(['Apply the geometry and timing saved with the ROIs (orientation, voxel size, ' ...
                    'temporal resolution, frame binning)? The ROI voxels are correct either way.'], ...
                    'Load ROIs', {'Apply', 'Keep current'}, 'Apply');
                applied = strcmp(b, 'Apply');
                if applied && geoDiff
                    obj.applyGeometry(sCodes, S.VoxNative, S.dt, S.t0);
                end
                if applied && binDiff
                    kp = false;
                    if isfield(S, 'BinKeepPartial'), kp = logical(S.BinKeepPartial); end
                    obj.applyBinning(sBin, kp);
                end
            end
            if isfield(S, 'MotionCorrected') && isfield(S, 'Motion') && ~isempty(S.Motion) && ...
                    size(S.Motion.Params, 1) == D.NTRaw
                if S.MotionCorrected && ~D.MotionCorrected
                    % the ROIs were analysed on motion-corrected frames: redo the same correction
                    obj.applyMotion(S.Motion, true);
                elseif ~S.MotionCorrected && D.MotionCorrected
                    obj.revertMotion();
                elseif ~D.MotionCorrected
                    D.setMotion(S.Motion);
                end
            end
            if strcmp(mode, 'Replace')
                obj.ROIs = obj.ROIs([]);
                obj.Masks = {};
            end
            names = S.ROINames;
            if ischar(names), names = {names}; end
            for r = 1:numel(S.ROIIndexNative)
                Mn = false(D.NativeSize);
                Mn(double(S.ROIIndexNative{r})) = true;
                k = obj.addROI(names{r}, S.ROIColors(r, :));
                obj.Masks{k} = D.native2canon(Mn);
                obj.ROIs(k).Visible = logical(S.ROIVisible(r));
            end
            if isfield(S, 'BaseFrames') && (applied || ~binDiff)
                bfr = max(1, min(D.NT, S.BaseFrames));
                obj.BaseFrames = [bfr(1) max(bfr)];
                set(obj.H.base1, 'String', num2str(obj.BaseFrames(1)));
                set(obj.H.base2, 'String', num2str(obj.BaseFrames(2)));
            end
            if isfield(S, 'SubFrames') && (applied || ~binDiff)
                obj.setSubFrames(S.SubFrames(1), S.SubFrames(2), false);
            end
            if isfield(S, 'ParamOpts'), D.setParamOptions(S.ParamOpts); end
            if isfield(S, 'RefROIName'), obj.RefROIName = S.RefROIName; end
            obj.UndoStack = {}; obj.RedoStack = {};
            obj.afterMaskChange();
            obj.Modified = false;
            obj.setStatus(sprintf('Loaded %d ROI(s) from %s', numel(S.ROIIndexNative), file));
        end

        function importMaskNifti(obj, file)
            D = obj.Data;
            [img, h] = graspilot_nifti_read(file);
            img = img(:, :, :, 1);
            [img, msg] = graspilot_reorient_to_ref(img, h.Affine, D.affine(), D.canonSize());
            img = round(double(img));
            labels = unique(img(img ~= 0));
            if isempty(labels)
                obj.inform('The mask is empty.', 'Import mask');
                return
            end
            [names, cols] = GRASPilot.readLabelTable(regexprep(file, '(\.nii\.gz|\.nii)$', '', 'ignorecase'));
            mode = 'Replace';
            if ~isempty(obj.ROIs)
                mode = obj.ask('Replace the current ROIs or append?', 'Import mask', ...
                    {'Replace', 'Append', 'Cancel'}, 'Replace');
                if isempty(mode) || strcmp(mode, 'Cancel'), return; end
            end
            if strcmp(mode, 'Replace')
                obj.ROIs = obj.ROIs([]);
                obj.Masks = {};
            end
            for i = 1:numel(labels)
                L = labels(i);
                nm = sprintf('Label %d', L); c = [];
                if L <= numel(names) && ~isempty(names{L}), nm = names{L}; c = cols(L, :); end
                k = obj.addROI(nm, c);
                obj.Masks{k} = (img == L);
            end
            obj.UndoStack = {}; obj.RedoStack = {};
            obj.afterMaskChange();
            s = sprintf('Imported %d label(s) from %s.', numel(labels), file);
            if ~isempty(msg), s = [s ' ' msg]; end
            obj.setStatus(s);
        end

        function close(obj, force)
            if nargin < 2, force = false; end
            if ~force && ~obj.Quiet && obj.Modified && ~isempty(obj.Masks) && any(cellfun(@(m) any(m(:)), obj.Masks))
                b = obj.ask('Save ROIs before closing?', 'GRASPilot', {'Save', 'Discard', 'Cancel'}, 'Save');
                if isempty(b) || strcmp(b, 'Cancel'), return; end
                if strcmp(b, 'Save') && ~obj.onSaveROIs(), return; end
            end
            obj.stopPlay();
            try
                if ~isempty(obj.PlayTimer) && isvalid(obj.PlayTimer), delete(obj.PlayTimer); end
            catch
            end
            if isfield(obj.PW, 'Fig') && ishghandle(obj.PW.Fig), delete(obj.PW.Fig); end
            if isfield(obj.MW, 'Fig') && ishghandle(obj.MW.Fig), delete(obj.MW.Fig); end
            if ~isempty(obj.Fig) && ishghandle(obj.Fig), delete(obj.Fig); end
        end
    end

    % =====================================================================
    methods (Hidden)

        % ------------------------------------------------------------ UI build
        function buildUI(obj)
            C = obj.Col;
            ss = get(0, 'ScreenSize');
            W = min(1720, max(1200, ss(3) - 80));
            Hh = min(1020, max(770, ss(4) - 110));
            pos = [max(1, round((ss(3) - W) / 2)), max(30, round((ss(4) - Hh) / 2)), W, Hh];
            obj.Fig = figure('Name', 'GRASPilot', 'NumberTitle', 'off', 'MenuBar', 'none', ...
                'ToolBar', 'none', 'Color', C.bg, 'Units', 'pixels', 'Position', pos, 'Visible', 'off', ...
                'CloseRequestFcn', @(s, e) obj.close(), ...
                'WindowButtonDownFcn', @(s, e) obj.onMouseDown(), ...
                'WindowButtonMotionFcn', @(s, e) obj.onMouseMove(), ...
                'WindowButtonUpFcn', @(s, e) obj.onMouseUp(), ...
                'WindowScrollWheelFcn', @(s, e) obj.onScroll(e), ...
                'KeyPressFcn', @(s, e) obj.onKey(e));
            try, set(obj.Fig, 'DockControls', 'off'); catch, end
            obj.buildMenus();

            obj.H.left = uipanel('Parent', obj.Fig, 'Units', 'pixels', 'BorderType', 'none', ...
                'BackgroundColor', C.panel);
            obj.V = repmat(struct('Panel', [], 'Ax', [], 'Img', [], 'CrossH', [], ...
                'CrossV', [], 'Cursor', [], 'Title', [], 'Lbl', [], 'Size', [2 2]), 1, 3);
            for v = 1:3
                P = uipanel('Parent', obj.Fig, 'Units', 'pixels', 'BorderType', 'none', 'BackgroundColor', [0 0 0]);
                ax = axes('Parent', P, 'Units', 'normalized', 'Position', [0.005 0.005 0.99 0.905], ...
                    'Color', [0 0 0], 'XTick', [], 'YTick', [], 'YDir', 'reverse', 'Visible', 'off', ...
                    'CLim', [0 1], 'NextPlot', 'add');
                colormap(ax, gray(256));
                obj.V(v).Panel = P;
                obj.V(v).Ax = ax;
                obj.V(v).Img = image('Parent', ax, 'CData', zeros(2, 2, 3), 'HitTest', 'off');
                obj.V(v).CrossH = line('Parent', ax, 'XData', [NaN NaN], 'YData', [NaN NaN], 'Color', C.cross, ...
                    'LineWidth', 0.5, 'HitTest', 'off');
                obj.V(v).CrossV = line('Parent', ax, 'XData', [NaN NaN], 'YData', [NaN NaN], 'Color', C.cross, ...
                    'LineWidth', 0.5, 'HitTest', 'off');
                obj.V(v).Cursor = line('Parent', ax, 'XData', NaN, 'YData', NaN, 'Color', [1 1 0.4], ...
                    'LineWidth', 1, 'HitTest', 'off', 'Visible', 'off');
                obj.V(v).Title = uicontrol(P, 'Style', 'text', 'Units', 'pixels', 'Position', [4 4 200 18], ...
                    'String', obj.ViewNames{v}, 'BackgroundColor', [0 0 0], 'ForegroundColor', [0.95 0.95 0.95], ...
                    'HorizontalAlignment', 'left', 'FontSize', 9, 'FontWeight', 'bold');
                lb = GRASPilot.handleArray(4);
                for k = 1:4
                    lb(k) = uicontrol(P, 'Style', 'text', 'Units', 'pixels', 'Position', [1 1 14 16], ...
                        'String', obj.ViewLabels{v}{k}, 'BackgroundColor', [0 0 0], ...
                        'ForegroundColor', [0.95 0.75 0.30], 'FontSize', 9, 'FontWeight', 'bold');
                end
                obj.V(v).Lbl = lb;
                if v == 3
                    obj.CB.Ax = axes('Parent', P, 'Units', 'pixels', 'Position', [10 10 12 100], 'Visible', 'off', ...
                        'YDir', 'normal', 'XTick', [], 'YTick', []);
                    obj.CB.Img = image('Parent', obj.CB.Ax, 'CData', zeros(2, 1, 3), 'XData', [1 1], ...
                        'YData', [1 128], 'Visible', 'off', 'HitTest', 'off');
                    obj.CB.Hi = uicontrol(P, 'Style', 'text', 'Units', 'pixels', 'Position', [1 1 60 14], ...
                        'BackgroundColor', [0 0 0], 'ForegroundColor', [1 1 1], 'FontSize', 8, ...
                        'HorizontalAlignment', 'right', 'Visible', 'off');
                    obj.CB.Lo = uicontrol(P, 'Style', 'text', 'Units', 'pixels', 'Position', [1 1 60 14], ...
                        'BackgroundColor', [0 0 0], 'ForegroundColor', [1 1 1], 'FontSize', 8, ...
                        'HorizontalAlignment', 'right', 'Visible', 'off');
                    obj.CB.Title = uicontrol(P, 'Style', 'text', 'Units', 'pixels', 'Position', [1 1 120 14], ...
                        'BackgroundColor', [0 0 0], 'ForegroundColor', [1 1 1], 'FontSize', 8, ...
                        'HorizontalAlignment', 'right', 'Visible', 'off');
                end
            end
            obj.H.time = uipanel('Parent', obj.Fig, 'Units', 'pixels', 'BorderType', 'none', 'BackgroundColor', C.panel);
            obj.H.ticPanel = uipanel('Parent', obj.Fig, 'Units', 'pixels', 'BorderType', 'none', 'BackgroundColor', C.panel);
            obj.H.roiPanel = uipanel('Parent', obj.Fig, 'Units', 'pixels', 'BorderType', 'none', 'BackgroundColor', C.panel);
            obj.H.status = uicontrol(obj.Fig, 'Style', 'text', 'Units', 'pixels', 'Position', [1 1 400 20], ...
                'String', 'Open a GRASP reconstruction (.mat) to begin.', 'HorizontalAlignment', 'left', ...
                'BackgroundColor', C.bg, 'ForegroundColor', C.dim, 'FontSize', 9);
            obj.buildLeft();
            obj.buildTime();
            obj.buildTicPanel();
            obj.buildROIPanel();
            obj.setTool('nav');
            obj.setStatus('Open a GRASP reconstruction (.mat) to begin: File > Open, or the Open button.');
            try
                set(obj.Fig, 'SizeChangedFcn', @(s, e) obj.doLayout());
            catch
                set(obj.Fig, 'ResizeFcn', @(s, e) obj.doLayout());
            end
            obj.doLayout();
        end

        function buildMenus(obj)
            m = uimenu(obj.Fig, 'Label', 'File');
            uimenu(m, 'Label', 'Open .mat...', 'Callback', @(s, e) obj.onOpen());
            uimenu(m, 'Label', 'Save ROIs (session)...', 'Separator', 'on', 'Callback', @(s, e) obj.onSaveROIs());
            uimenu(m, 'Label', 'Load ROIs (session)...', 'Callback', @(s, e) obj.onLoadROIs());
            uimenu(m, 'Label', 'Import mask NIfTI...', 'Callback', @(s, e) obj.onImportMask());
            uimenu(m, 'Label', 'Export image -> NIfTI...', 'Separator', 'on', 'Callback', @(s, e) obj.onExportImage());
            uimenu(m, 'Label', 'Export mask -> NIfTI...', 'Callback', @(s, e) obj.onExportMask());
            uimenu(m, 'Label', 'Export ROI voxel time series -> spreadsheet...', 'Callback', @(s, e) obj.onExportTS());
            uimenu(m, 'Label', 'Export enhancement maps -> NIfTI...', 'Callback', @(s, e) obj.onExportMaps());
            uimenu(m, 'Label', 'Save TIC plot as PNG...', 'Callback', @(s, e) obj.onSavePNG());
            uimenu(m, 'Label', 'Close', 'Separator', 'on', 'Callback', @(s, e) obj.close());
            m = uimenu(obj.Fig, 'Label', 'Edit');
            uimenu(m, 'Label', 'Undo   (Ctrl+Z)', 'Callback', @(s, e) obj.undo());
            uimenu(m, 'Label', 'Redo   (Ctrl+Y)', 'Callback', @(s, e) obj.redo());
            uimenu(m, 'Label', 'Clear active ROI on slice', 'Separator', 'on', 'Callback', @(s, e) obj.clearROI('slice'));
            uimenu(m, 'Label', 'Clear active ROI (all slices)', 'Callback', @(s, e) obj.clearROI('all'));
            uimenu(m, 'Label', 'Fill holes in active ROI', 'Callback', @(s, e) obj.fillHoles());
            uimenu(m, 'Label', 'Detect superior sagittal sinus (auto ROI)', 'Separator', 'on', ...
                'Callback', @(s, e) obj.autoSSS());
            uimenu(m, 'Label', 'Grow from crosshair   (G)', 'Callback', @(s, e) obj.growFromCrosshair());
            uimenu(m, 'Label', 'All voxels in threshold window -> active ROI', 'Callback', @(s, e) obj.windowToROI());
            m = uimenu(obj.Fig, 'Label', 'View');
            uimenu(m, 'Label', 'Geometry and orientation...', 'Callback', @(s, e) obj.geometryDialog());
            uimenu(m, 'Label', 'Show subtraction map   (S)', 'Callback', @(s, e) obj.toggleSubView());
            uimenu(m, 'Label', 'Temporal binning...', 'Callback', @(s, e) obj.onBinningDialog());
            uimenu(m, 'Label', 'Auto contrast', 'Separator', 'on', 'Callback', @(s, e) obj.autoContrast(true));
            uimenu(m, 'Label', 'Reset zoom', 'Callback', @(s, e) obj.resetZoom());
            m = uimenu(obj.Fig, 'Label', 'Analysis');
            uimenu(m, 'Label', 'Semi-quantitative parameters (wash-in, washout, AUC, PE)...', ...
                'Callback', @(s, e) obj.openParamWindow());
            uimenu(m, 'Label', 'Show / hide parameter map   (M)', 'Callback', @(s, e) obj.setParamDisplay(~obj.ParamOn));
            uimenu(m, 'Label', 'Export parameter statistics and voxel values...', 'Callback', @(s, e) obj.onExportParams());
            uimenu(m, 'Label', 'Motion check and correction...', 'Separator', 'on', ...
                'Callback', @(s, e) obj.openMotionWindow());
            m = uimenu(obj.Fig, 'Label', 'Help');
            uimenu(m, 'Label', 'Mouse and keyboard...', 'Callback', @(s, e) obj.showHelp());
            uimenu(m, 'Label', 'About', 'Callback', @(s, e) obj.inform(sprintf(['GRASPilot %s\n' ...
                'Frames are synthesized from the subspace coefficients exactly as in ReadData.m ' ...
                '(TempPCASub adjoint, CropImg, global-max normalization).\nNo toolboxes required; ' ...
                'runs offline.'], obj.Version), 'About'));
        end

        function h = ctl(obj, parent, x, top, w, ht, style, str, varargin)
            pbg = get(parent, 'BackgroundColor');
            args = {'Parent', parent, 'Units', 'pixels', 'Position', [max(x, 1) 1 max(w, 10) ht], ...
                'Style', style, 'String', str, 'FontSize', 9};
            switch style
                case 'text'
                    args = [args, {'BackgroundColor', pbg, 'ForegroundColor', obj.Col.fg, 'HorizontalAlignment', 'left'}];
                case 'checkbox'
                    args = [args, {'BackgroundColor', pbg, 'ForegroundColor', obj.Col.fg}];
                case {'edit', 'popupmenu'}
                    args = [args, {'BackgroundColor', [1 1 1], 'ForegroundColor', [0 0 0]}];
                case 'slider'
                    args = [args, {'BackgroundColor', [0.85 0.85 0.85]}];
            end
            h = uicontrol(args{:}, varargin{:});
            if any(strcmp(style, {'pushbutton', 'togglebutton', 'checkbox'}))
                set(h, 'KeyPressFcn', @(s, e) obj.onKey(e));
            end
            obj.reg(h, parent, x, top, w, ht);
        end

        function reg(obj, h, parent, x, top, w, ht)
            obj.Layout(end + 1) = struct('h', h, 'panel', parent, 'x', x, 'top', top, 'w', w, 'ht', ht);
        end

        function buildLeft(obj)
            P = obj.H.left; C = obj.Col;
            hdr = @(top, s) obj.ctl(P, 10, top, -10, 18, 'text', s, 'FontWeight', 'bold', 'ForegroundColor', C.accent);
            hdr(8, 'DATA');
            obj.ctl(P, 10, 28, 135, 26, 'pushbutton', 'Open .mat...', 'Callback', @(s, e) obj.onOpen());
            obj.ctl(P, 150, 28, 140, 26, 'pushbutton', 'Geometry...', 'Callback', @(s, e) obj.geometryDialog());
            obj.H.info = obj.ctl(P, 10, 58, -10, 44, 'text', 'No data loaded.', 'FontSize', 8, 'ForegroundColor', C.dim);
            obj.ctl(P, 10, 108, 86, 18, 'text', 'Average every');
            obj.H.binEdit = obj.ctl(P, 98, 106, 36, 22, 'edit', '1', ...
                'TooltipString', 'Combine consecutive frames by averaging signal intensity (1 = original)', ...
                'Callback', @(s, e) obj.applyBinning(str2double(get(s, 'String'))));
            obj.ctl(P, 138, 108, 48, 18, 'text', 'frames');
            obj.ctl(P, 188, 106, 102, 22, 'pushbutton', 'Apply binning', ...
                'TooltipString', 'Lower the temporal resolution: each frame = mean of n consecutive frames', ...
                'Callback', @(s, e) obj.applyBinning(str2double(get(obj.H.binEdit, 'String'))));

            hdr(136, 'SEGMENTATION');
            names = {'nav', 'circle', 'brush', 'erase'};
            labels = {'Navigate', 'Circle', 'Brush', 'Eraser'};
            xs = [10 81 152 223];
            for i = 1:4
                nm = names{i};
                obj.H.tool.(nm) = obj.ctl(P, xs(i), 156, 66, 26, 'togglebutton', labels{i}, ...
                    'Callback', @(s, e) obj.setTool(nm));
            end
            obj.ctl(P, 10, 191, 80, 18, 'text', 'Radius (mm)');
            obj.H.rSlider = obj.ctl(P, 92, 191, 140, 18, 'slider', '', 'Min', 0.5, 'Max', 25, ...
                'Value', obj.BrushR, 'SliderStep', [0.5 2] / 24.5, 'Callback', @(s, e) obj.setBrushR(get(s, 'Value')));
            obj.H.rEdit = obj.ctl(P, 238, 189, 52, 22, 'edit', num2str(obj.BrushR), ...
                'Callback', @(s, e) obj.setBrushR(str2double(get(s, 'String'))));
            obj.H.chk3D = obj.ctl(P, 10, 216, 135, 20, 'checkbox', '3D (sphere)', 'Value', obj.Use3D, ...
                'TooltipString', 'Circle, brush and eraser act in 3D (key D)', 'Callback', @(s, e) obj.set3D(get(s, 'Value')));
            obj.H.chkExcl = obj.ctl(P, 150, 216, 140, 20, 'checkbox', 'Exclusive ROIs', 'Value', obj.Exclusive, ...
                'TooltipString', 'Painting removes the voxels from other ROIs (no overlap)', ...
                'Callback', @(s, e) obj.setFlag('Exclusive', get(s, 'Value')));
            obj.ctl(P, 10, 240, 66, 24, 'pushbutton', 'Undo', 'Callback', @(s, e) obj.undo());
            obj.ctl(P, 81, 240, 66, 24, 'pushbutton', 'Redo', 'Callback', @(s, e) obj.redo());
            obj.ctl(P, 152, 240, 66, 24, 'pushbutton', 'Clear slice', ...
                'TooltipString', 'Clear the active ROI on the current slice of the last-used view', ...
                'Callback', @(s, e) obj.clearROI('slice'));
            obj.ctl(P, 223, 240, 66, 24, 'pushbutton', 'Clear ROI', 'Callback', @(s, e) obj.clearROI('all'));

            hdr(272, 'THRESHOLD / SUBTRACTION');
            obj.ctl(P, 10, 293, 40, 18, 'text', 'Map');
            obj.H.thrSrc = obj.ctl(P, 50, 291, -10, 22, 'popupmenu', obj.ThrNames, 'Value', obj.ThrSrc, ...
                'Callback', @(s, e) obj.setThrSource(get(s, 'Value')));
            obj.ctl(P, 10, 321, 70, 18, 'text', 'Subtraction');
            obj.ctl(P, 82, 321, 28, 18, 'text', 'post');
            obj.H.subPost = obj.ctl(P, 110, 319, 38, 22, 'edit', '1', ...
                'TooltipString', 'Post-contrast frame (default: last)', ...
                'Callback', @(s, e) obj.setSubFrames(str2double(get(s, 'String')), obj.SubFrames(2)));
            obj.ctl(P, 151, 321, 32, 18, 'text', '- pre');
            obj.H.subPre = obj.ctl(P, 184, 319, 38, 22, 'edit', '1', ...
                'TooltipString', 'Pre-contrast frame (default: first)', ...
                'Callback', @(s, e) obj.setSubFrames(obj.SubFrames(1), str2double(get(s, 'String'))));
            obj.H.btnSub = obj.ctl(P, 226, 319, 64, 22, 'togglebutton', 'Show', ...
                'TooltipString', 'Display the subtraction map instead of the current frame (key S)', ...
                'Callback', @(s, e) obj.toggleSubView());
            obj.ctl(P, 10, 349, 100, 18, 'text', 'Baseline frames');
            obj.H.base1 = obj.ctl(P, 112, 347, 40, 22, 'edit', '1', 'Callback', @(s, e) obj.onBaseline());
            obj.ctl(P, 155, 349, 20, 18, 'text', 'to');
            obj.H.base2 = obj.ctl(P, 176, 347, 40, 22, 'edit', '5', 'Callback', @(s, e) obj.onBaseline());
            obj.ctl(P, 222, 347, 68, 22, 'pushbutton', 'Maps', ...
                'TooltipString', 'Recompute baseline / peak enhancement maps', 'Callback', @(s, e) obj.ensureMaps(true));
            obj.ctl(P, 10, 377, 32, 18, 'text', 'Low');
            obj.H.thrLo = obj.ctl(P, 42, 375, 74, 22, 'edit', '-Inf', 'Callback', @(s, e) obj.onThrEdit());
            obj.ctl(P, 122, 377, 34, 18, 'text', 'High');
            obj.H.thrHi = obj.ctl(P, 156, 375, 74, 22, 'edit', 'Inf', 'Callback', @(s, e) obj.onThrEdit());
            obj.ctl(P, 236, 375, 54, 22, 'pushbutton', 'Otsu', ...
                'TooltipString', 'Suggest Low from the active ROI (Otsu)', 'Callback', @(s, e) obj.suggestThreshold());
            obj.H.thrSlider = obj.ctl(P, 10, 403, -10, 18, 'slider', '', 'Min', 0, 'Max', 1, 'Value', 0, ...
                'Callback', @(s, e) obj.onThrSlider());
            try
                addlistener(obj.H.thrSlider, 'ContinuousValueChange', @(s, e) obj.onThrSlider());
            catch
            end
            obj.H.chkPrev = obj.ctl(P, 10, 425, 90, 20, 'checkbox', 'Preview', 'Value', 0, ...
                'TooltipString', 'Show voxels inside the threshold window (green) - key P', ...
                'Callback', @(s, e) obj.setPreview(get(s, 'Value')));
            obj.H.chkCons = obj.ctl(P, 105, 425, 185, 20, 'checkbox', 'Constrain brush/circle', 'Value', 0, ...
                'TooltipString', 'Brush and circle add only voxels inside the threshold window', ...
                'Callback', @(s, e) obj.setConstrain(get(s, 'Value')));
            obj.ctl(P, 10, 449, 92, 24, 'pushbutton', 'Keep: slice', ...
                'TooltipString', 'Remove voxels outside the window from the active ROI on this slice', ...
                'Callback', @(s, e) obj.applyThreshold('slice'));
            obj.ctl(P, 106, 449, 92, 24, 'pushbutton', 'Keep: 3D', ...
                'TooltipString', 'Remove voxels outside the window from the whole active ROI', ...
                'Callback', @(s, e) obj.applyThreshold('3d'));
            obj.ctl(P, 202, 449, 88, 24, 'pushbutton', 'Fill holes', ...
                'TooltipString', 'Fill enclosed holes slice by slice (plane of the last-used view)', ...
                'Callback', @(s, e) obj.fillHoles());
            obj.ctl(P, 10, 477, 130, 24, 'pushbutton', 'Grow from crosshair', ...
                'TooltipString', ['Add the 3D region of in-window voxels connected to the crosshair ' ...
                'voxel to the active ROI (key G)'], 'Callback', @(s, e) obj.growFromCrosshair());
            obj.ctl(P, 142, 480, 24, 18, 'text', 'max');
            obj.H.growR = obj.ctl(P, 166, 478, 32, 22, 'edit', num2str(obj.GrowR), ...
                'TooltipString', 'Maximum distance from the crosshair (mm)', 'Callback', @(s, e) obj.onGrowR());
            obj.ctl(P, 199, 480, 22, 18, 'text', 'mm');
            obj.ctl(P, 222, 477, 68, 24, 'pushbutton', 'All -> ROI', ...
                'TooltipString', 'Add every voxel inside the threshold window (whole volume) to the active ROI', ...
                'Callback', @(s, e) obj.windowToROI());
            obj.H.thrInfo = obj.ctl(P, 10, 505, -10, 18, 'text', '', 'FontSize', 8, 'ForegroundColor', C.dim);

            hdr(530, 'DISPLAY');
            obj.ctl(P, 10, 552, 50, 18, 'text', 'Window');
            obj.H.wMin = obj.ctl(P, 62, 550, 60, 22, 'edit', '0', 'Callback', @(s, e) obj.onWindowEdit());
            obj.H.wMax = obj.ctl(P, 126, 550, 60, 22, 'edit', '1', 'Callback', @(s, e) obj.onWindowEdit());
            obj.ctl(P, 190, 550, 48, 22, 'pushbutton', 'Auto', 'Callback', @(s, e) obj.autoContrast(true));
            obj.ctl(P, 242, 550, 48, 22, 'pushbutton', 'Fit', 'TooltipString', 'Reset zoom', ...
                'Callback', @(s, e) obj.resetZoom());
            obj.ctl(P, 10, 578, 100, 18, 'text', 'ROI opacity');
            obj.H.alpha = obj.ctl(P, 112, 578, -10, 18, 'slider', '', 'Min', 0, 'Max', 1, 'Value', obj.Alpha, ...
                'Callback', @(s, e) obj.setFlag('Alpha', get(s, 'Value')));
            obj.H.chkOutline = obj.ctl(P, 10, 600, 80, 20, 'checkbox', 'Outline', 'Value', obj.ShowOutline, ...
                'Callback', @(s, e) obj.setFlag('ShowOutline', get(s, 'Value')));
            obj.H.chkCross = obj.ctl(P, 95, 600, 90, 20, 'checkbox', 'Crosshair', 'Value', obj.ShowCross, ...
                'Callback', @(s, e) obj.setFlag('ShowCross', get(s, 'Value')));
            obj.H.chkCursor = obj.ctl(P, 190, 600, 100, 20, 'checkbox', 'Cursor TIC', 'Value', obj.ShowCursorTIC, ...
                'TooltipString', 'Plot the time course of the crosshair voxel', ...
                'Callback', @(s, e) obj.setFlag('ShowCursorTIC', get(s, 'Value')));

            hdr(628, 'EXPORT');
            obj.ctl(P, 10, 648, 137, 24, 'pushbutton', 'Image -> NIfTI', 'Callback', @(s, e) obj.onExportImage());
            obj.ctl(P, 152, 648, 138, 24, 'pushbutton', 'Mask -> NIfTI', 'Callback', @(s, e) obj.onExportMask());
            obj.ctl(P, 10, 676, -10, 24, 'pushbutton', 'ROI voxel time series -> spreadsheet', ...
                'Callback', @(s, e) obj.onExportTS());
            obj.ctl(P, 10, 704, 90, 24, 'pushbutton', 'Save ROIs', 'Callback', @(s, e) obj.onSaveROIs());
            obj.ctl(P, 104, 704, 90, 24, 'pushbutton', 'Load ROIs', 'Callback', @(s, e) obj.onLoadROIs());
            obj.ctl(P, 198, 704, 92, 24, 'pushbutton', 'Import mask', 'Callback', @(s, e) obj.onImportMask());
        end

        function buildTime(obj)
            P = obj.H.time;
            obj.ctl(P, 6, 7, 30, 22, 'pushbutton', '|<', 'TooltipString', 'First frame (Home)', ...
                'Callback', @(s, e) obj.setFrame(1));
            obj.ctl(P, 38, 7, 30, 22, 'pushbutton', '<', 'TooltipString', 'Previous frame (Left)', ...
                'Callback', @(s, e) obj.setFrame(obj.Frame - 1));
            obj.H.slider = obj.ctl(P, 72, 9, -552, 18, 'slider', '', 'Min', 0, 'Max', 1, 'Value', 0, ...
                'Callback', @(s, e) obj.onSlider());
            try
                addlistener(obj.H.slider, 'ContinuousValueChange', @(s, e) obj.onSlider());
            catch
            end
            obj.ctl(P, -546, 7, 30, 22, 'pushbutton', '>', 'TooltipString', 'Next frame (Right)', ...
                'Callback', @(s, e) obj.setFrame(obj.Frame + 1));
            obj.ctl(P, -514, 7, 30, 22, 'pushbutton', '>|', 'TooltipString', 'Last frame (End)', ...
                'Callback', @(s, e) obj.setFrame(Inf));
            obj.ctl(P, -480, 9, 48, 18, 'text', 'Frame', 'HorizontalAlignment', 'right');
            obj.H.frameEdit = obj.ctl(P, -428, 7, 46, 22, 'edit', '1', ...
                'Callback', @(s, e) obj.setFrame(str2double(get(s, 'String')), true));
            obj.H.frameTotal = obj.ctl(P, -378, 9, 44, 18, 'text', '/ 0');
            obj.ctl(P, -332, 9, 48, 18, 'text', 's/frame', 'HorizontalAlignment', 'right');
            obj.H.dtEdit = obj.ctl(P, -282, 7, 46, 22, 'edit', '', ...
                'TooltipString', ['Temporal resolution of the reconstructed frames in seconds ' ...
                '(spokes per frame x TR). Blank = unknown: curves are plotted per frame.'], ...
                'Callback', @(s, e) obj.setTemporalResolution(str2double(get(s, 'String'))));
            obj.H.timeText = obj.ctl(P, -232, 9, 106, 18, 'text', '', 'ForegroundColor', obj.Col.accent, ...
                'FontWeight', 'bold');
            obj.H.btnPlay = obj.ctl(P, -122, 7, 54, 22, 'togglebutton', 'Play', 'Callback', @(s, e) obj.togglePlay());
            obj.ctl(P, -66, 9, 24, 18, 'text', 'fps');
            obj.H.fps = obj.ctl(P, -42, 7, 36, 22, 'edit', '8', 'Callback', @(s, e) obj.onFps());
        end

        function buildTicPanel(obj)
            P = obj.H.ticPanel; C = obj.Col;
            obj.ctl(P, 10, 7, 32, 18, 'text', 'TIC', 'FontWeight', 'bold', 'ForegroundColor', C.accent);
            obj.ctl(P, 44, 7, 12, 18, 'text', 'Y');
            obj.H.ticY = obj.ctl(P, 58, 5, 150, 22, 'popupmenu', ...
                {'Signal (a.u.)', 'Relative enhancement (%)', 'S / S0', 'S - S0'}, ...
                'Callback', @(s, e) obj.setFlag('TicY', get(s, 'Value')));
            obj.ctl(P, 214, 7, 28, 18, 'text', 'Stat');
            obj.H.ticStat = obj.ctl(P, 244, 5, 112, 22, 'popupmenu', {'Mean +/- SD', 'Median [IQR]', 'Mean only'}, ...
                'Callback', @(s, e) obj.setFlag('TicStat', get(s, 'Value')));
            obj.ctl(P, 362, 7, 12, 18, 'text', 'X');
            obj.H.ticX = obj.ctl(P, 376, 5, 78, 22, 'popupmenu', {'Time (s)', 'Frame'}, ...
                'Callback', @(s, e) obj.setFlag('TicX', get(s, 'Value')));
            obj.ctl(P, -176, 5, 98, 22, 'pushbutton', 'Parameters...', ...
                'TooltipString', 'Wash-in/washout slopes, AUC, peak enhancement: maps, ROI statistics, export', ...
                'Callback', @(s, e) obj.openParamWindow());
            obj.ctl(P, -74, 5, 66, 22, 'pushbutton', 'PNG...', 'TooltipString', 'Save the plot as PNG', ...
                'Callback', @(s, e) obj.onSavePNG());
            ax = axes('Parent', P, 'Units', 'pixels', 'Position', [60 36 300 120], 'Color', [0.07 0.07 0.08], ...
                'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'Box', 'on', 'FontSize', 8, 'NextPlot', 'add');
            try
                set(ax, 'XGrid', 'on', 'YGrid', 'on', 'GridColor', [0.6 0.6 0.6], 'GridAlpha', 0.25);
            catch
                grid(ax, 'on');
            end
            obj.Tic.Ax = ax;
            obj.Tic.Items = [];
            obj.Tic.FrameLine = [];
            obj.Tic.CursorLine = [];
            obj.Tic.S0Cursor = NaN;
            xlabel(ax, 'Time (s)', 'Color', [0.85 0.85 0.85]);
            ylabel(ax, 'Signal (a.u.)', 'Color', [0.85 0.85 0.85]);
        end

        function buildROIPanel(obj)
            P = obj.H.roiPanel; C = obj.Col;
            obj.ctl(P, 8, 7, 40, 18, 'text', 'ROIs', 'FontWeight', 'bold', 'ForegroundColor', C.accent);
            obj.ctl(P, 50, 5, 76, 22, 'pushbutton', 'Auto SSS', ...
                'TooltipString', ['Detect the posterior superior sagittal sinus from just above the torcula ' ...
                'up toward the vertex (excludes transverse and straight sinuses) and place the ROI "SSS (auto)"'], ...
                'Callback', @(s, e) obj.autoSSS());
            obj.ctl(P, 130, 5, 64, 22, 'pushbutton', '+ SSS', ...
                'TooltipString', 'Add an empty superior sagittal sinus control ROI to draw by hand', ...
                'Callback', @(s, e) obj.addSSS());
            obj.ctl(P, 198, 5, 70, 22, 'pushbutton', '+ Tumor', 'Callback', @(s, e) obj.addTumor());
            obj.ctl(P, 272, 5, 80, 22, 'pushbutton', '+ Custom...', 'Callback', @(s, e) obj.addCustom());
            obj.ctl(P, 8, 35, 44, 18, 'text', 'Active');
            obj.H.swatch = obj.ctl(P, 54, 35, 18, 18, 'text', '', 'BackgroundColor', [0.3 0.3 0.3]);
            obj.H.roiPop = obj.ctl(P, 78, 33, -8, 22, 'popupmenu', {'(none)'}, ...
                'Callback', @(s, e) obj.setActiveROI(get(s, 'Value')));
            obj.ctl(P, 8, 61, 78, 22, 'pushbutton', 'Rename', 'Callback', @(s, e) obj.renameROI());
            obj.ctl(P, 90, 61, 70, 22, 'pushbutton', 'Color', 'Callback', @(s, e) obj.recolorROI());
            obj.ctl(P, 164, 61, 92, 22, 'pushbutton', 'Show/Hide', 'Callback', @(s, e) obj.toggleROI());
            obj.ctl(P, 260, 61, 92, 22, 'pushbutton', 'Delete', 'Callback', @(s, e) obj.deleteROI());
            obj.H.roiTable = uitable('Parent', P, 'Units', 'pixels', 'Position', [8 8 340 100], ...
                'ColumnName', {'#', 'Name', 'Show', 'Voxels', 'mL'}, ...
                'ColumnFormat', {'numeric', 'char', 'logical', 'numeric', 'numeric'}, ...
                'ColumnEditable', [false true true false false], 'ColumnWidth', {28, 140, 42, 64, 56}, ...
                'RowName', [], 'Data', cell(0, 5), 'FontSize', 9, ...
                'CellEditCallback', @(s, e) obj.onTableEdit(e), ...
                'CellSelectionCallback', @(s, e) obj.onTableSelect(e));
            obj.reg(obj.H.roiTable, P, 8, 89, -8, -8);
        end

        function doLayout(obj)
            if ~isfield(obj.H, 'roiPanel') || ~ishghandle(obj.H.roiPanel), return; end
            fp = get(obj.Fig, 'Position');
            W = fp(3); Hh = fp(4);
            LW = 300; SB = 22; roiW = 360;
            set(obj.H.left, 'Position', [1 SB + 1 LW max(1, Hh - SB)]);
            rx = LW + 2; rw = max(300, W - LW - 1);
            bottomH = max(200, round(0.34 * (Hh - SB)));
            timeH = 36;
            vy = SB + bottomH + timeH;
            vh = max(60, Hh - vy);
            vw = floor(rw / 3);
            for v = 1:3
                x = rx + (v - 1) * vw;
                w = vw;
                if v == 3, w = rw - 2 * vw; end
                p = [x, vy + 1, max(20, w - 2), vh - 1];
                set(obj.V(v).Panel, 'Position', p);
                obj.ViewRect(v, :) = p;
                pw = p(3); ph = p(4);
                set(obj.V(v).Title, 'Position', [6 ph - 21 max(20, pw - 12) 18]);
                L = obj.V(v).Lbl;
                mid = 0.005 * ph + 0.905 * ph / 2;
                set(L(1), 'Position', [3, mid - 8, 14, 16]);
                set(L(2), 'Position', [max(3, pw - 17), mid - 8, 14, 16]);
                set(L(3), 'Position', [max(3, pw / 2 - 7), 0.91 * ph - 17, 14, 16]);
                set(L(4), 'Position', [max(3, pw / 2 - 7), 3, 14, 16]);
                if v == 3 && isfield(obj.CB, 'Ax') && ishghandle(obj.CB.Ax)
                    h0 = max(40, round(0.45 * ph));
                    y0 = round(0.25 * ph);
                    set(obj.CB.Ax, 'Position', [max(1, pw - 34), y0, 12, h0]);
                    set(obj.CB.Hi, 'Position', [max(1, pw - 98), y0 + h0 - 12, 60, 14]);
                    set(obj.CB.Lo, 'Position', [max(1, pw - 98), y0 - 2, 60, 14]);
                    set(obj.CB.Title, 'Position', [max(1, pw - 158), y0 + h0 + 6, 150, 14]);
                end
            end
            set(obj.H.time, 'Position', [rx, SB + bottomH + 1, rw, timeH]);
            set(obj.H.ticPanel, 'Position', [rx, SB + 1, max(60, rw - roiW - 2), bottomH]);
            set(obj.H.roiPanel, 'Position', [rx + rw - roiW, SB + 1, roiW, bottomH]);
            set(obj.H.status, 'Position', [6 1 max(50, W - 12) 20]);
            for i = 1:numel(obj.Layout)
                L = obj.Layout(i);
                if ~ishghandle(L.h) || ~ishghandle(L.panel), continue; end
                pp = get(L.panel, 'Position');
                pw = pp(3); ph = pp(4);
                x = L.x;
                if x < 0, x = pw + x; end
                w = L.w;
                if w <= 0, w = pw - x + w; end
                ht = L.ht;
                if ht <= 0, ht = ph - L.top + ht; end
                set(L.h, 'Position', [max(1, x), ph - L.top - ht, max(w, 4), max(ht, 4)]);
            end
            tp = get(obj.H.ticPanel, 'Position');
            set(obj.Tic.Ax, 'Position', [60, 36, max(50, tp(3) - 60 - 16), max(40, tp(4) - 36 - 40)]);
        end

        % ------------------------------------------------------------ view geometry and rendering
        function resetViewsGeometry(obj)
            n = obj.Data.canonSize();
            for v = 1:3
                switch obj.ViewAxis(v)
                    case 3, sz = [n(2) n(1)];
                    case 2, sz = [n(3) n(1)];
                    case 1, sz = [n(3) n(2)];
                end
                obj.V(v).Size = sz;
                [sx, sy] = obj.viewSpacing(v);
                ax = obj.V(v).Ax;
                set(obj.V(v).Img, 'CData', zeros([sz 3]), 'XData', [1 sz(2)], 'YData', [1 sz(1)]);
                obj.ViewImg{v} = zeros(sz, 'single');
                obj.FrameImg{v} = zeros(sz, 'single');
                obj.OvlRGB{v} = []; obj.OvlA{v} = []; obj.PrevA{v} = [];
                set(ax, 'XLim', [0.5 sz(2) + 0.5], 'YLim', [0.5 sz(1) + 0.5], 'DataAspectRatio', [sy sx 1], ...
                    'CLim', obj.CLim);
            end
        end

        function [sx, sy] = viewSpacing(obj, v)
            vox = obj.Data.canonVox();
            switch obj.ViewAxis(v)
                case 3, sx = vox(1); sy = vox(2);
                case 2, sx = vox(1); sy = vox(3);
                case 1, sx = vox(2); sy = vox(3);
            end
        end

        function img = toScreen(obj, v, S)
            img = S.';
            if obj.ViewAxis(v) ~= 3, img = flipud(img); end
        end

        function img = screenSliceOf(obj, v, A)
            a = obj.ViewAxis(v);
            img = obj.toScreen(v, GRASPilot.sliceOfArray(A, a, obj.Pos(a), obj.Data.canonSize()));
        end

        function c = screen2canonCont(obj, v, X, Y)
            n = obj.Data.canonSize();
            switch obj.ViewAxis(v)
                case 3, c = [X, Y, obj.Pos(3)];
                case 2, c = [X, obj.Pos(2), n(3) - Y + 1];
                case 1, c = [obj.Pos(1), X, n(3) - Y + 1];
            end
        end

        function p = canon2screen(obj, v, c)
            n = obj.Data.canonSize();
            switch obj.ViewAxis(v)
                case 3, p = [c(1), c(2)];
                case 2, p = [c(1), n(3) - c(3) + 1];
                case 1, p = [c(2), n(3) - c(3) + 1];
            end
        end

        function lin = screen2lin(obj, v, X, Y)
            n = obj.Data.canonSize();
            X = X(:); Y = Y(:);
            o = ones(numel(X), 1);
            switch obj.ViewAxis(v)
                case 3, lin = sub2ind(n, X, Y, obj.Pos(3) * o);
                case 2, lin = sub2ind(n, X, obj.Pos(2) * o, n(3) - Y + 1);
                case 1, lin = sub2ind(n, obj.Pos(1) * o, X, n(3) - Y + 1);
            end
        end

        function lin = sliceLin(obj, v)
            n = obj.Data.canonSize();
            a = obj.ViewAxis(v);
            q = obj.Pos(a);
            switch a
                case 3
                    [I, J] = ndgrid(1:n(1), 1:n(2));
                    lin = sub2ind(n, I(:), J(:), q * ones(numel(I), 1));
                case 2
                    [I, K] = ndgrid(1:n(1), 1:n(3));
                    lin = sub2ind(n, I(:), q * ones(numel(I), 1), K(:));
                case 1
                    [J, K] = ndgrid(1:n(2), 1:n(3));
                    lin = sub2ind(n, q * ones(numel(J), 1), J(:), K(:));
            end
        end

        function lin = diskIndices(obj, v, X0, Y0, rmm)
            [sx, sy] = obj.viewSpacing(v);
            sz = obj.V(v).Size;
            h = sz(1); w = sz(2);
            xs = max(1, floor(X0 - rmm / sx)):min(w, ceil(X0 + rmm / sx));
            ys = max(1, floor(Y0 - rmm / sy)):min(h, ceil(Y0 + rmm / sy));
            [XX, YY] = meshgrid(xs, ys);
            in = ((XX - X0) * sx) .^ 2 + ((YY - Y0) * sy) .^ 2 <= rmm ^ 2;
            XX = XX(in); YY = YY(in);
            if isempty(XX)
                xr = round(X0); yr = round(Y0);
                if xr < 1 || xr > w || yr < 1 || yr > h
                    lin = [];
                    return
                end
                XX = xr; YY = yr;
            end
            lin = obj.screen2lin(v, XX, YY);
        end

        function lin = sphereIndices(obj, c, rmm)
            n = obj.Data.canonSize();
            vox = obj.Data.canonVox();
            rng = cell(1, 3);
            for a = 1:3
                rng{a} = max(1, floor(c(a) - rmm / vox(a))):min(n(a), ceil(c(a) + rmm / vox(a)));
            end
            if any(cellfun(@isempty, rng))
                lin = [];
                return
            end
            [I, J, K] = ndgrid(rng{1}, rng{2}, rng{3});
            in = ((I - c(1)) * vox(1)) .^ 2 + ((J - c(2)) * vox(2)) .^ 2 + ((K - c(3)) * vox(3)) .^ 2 <= rmm ^ 2;
            if ~any(in(:))
                ci = max(1, min(n, round(c)));
                lin = sub2ind(n, ci(1), ci(2), ci(3));
                return
            end
            lin = sub2ind(n, I(in), J(in), K(in));
        end

        function refreshAll(obj)
            if isempty(obj.Data), return; end
            for v = 1:3
                obj.refreshView(v, 'all');
            end
            obj.updateTimeWidgets();
            obj.refreshTIC();
            obj.updateThrInfo();
        end

        function refreshView(obj, v, what)
            if isempty(obj.Data), return; end
            a = obj.ViewAxis(v);
            comp = false;
            if any(strcmp(what, {'all', 'image'}))
                obj.FrameImg{v} = obj.toScreen(v, obj.Data.canonSlice(a, obj.Pos(a), obj.Frame));
                if strcmp(obj.ViewMode, 'sub')
                    obj.ViewImg{v} = obj.screenSliceOf(v, obj.subVolume());
                else
                    obj.ViewImg{v} = obj.FrameImg{v};
                end
                obj.computePreview(v);
                obj.computeParamOverlay(v);
                comp = true;
            end
            if any(strcmp(what, {'all', 'overlay'}))
                obj.computeOverlay(v);
                comp = true;
            end
            if comp, obj.composeView(v); end
            if any(strcmp(what, {'all', 'cross'}))
                obj.refreshCross(v);
                n = obj.Data.canonSize();
                tt = sprintf('%s   %d / %d', obj.ViewNames{v}, obj.Pos(a), n(a));
                if strcmp(obj.ViewMode, 'sub')
                    tt = sprintf('%s   [SUB %d - %d]', tt, obj.SubFrames(1), obj.SubFrames(2));
                end
                set(obj.V(v).Title, 'String', tt);
            end
        end

        function refreshOverlays(obj)
            for v = 1:3
                obj.refreshOverlay(v);
            end
        end

        function refreshOverlay(obj, v)
            obj.computeOverlay(v);
            obj.composeView(v);
        end

        function computeOverlay(obj, v)
            sz = obj.V(v).Size;
            rgb = zeros(sz(1), sz(2), 3);
            A = zeros(sz(1), sz(2));
            for r = 1:numel(obj.ROIs)
                if ~obj.ROIs(r).Visible, continue; end
                m = obj.screenSliceOf(v, obj.Masks{r});
                if ~any(m(:)), continue; end
                col = obj.ROIs(r).Color;
                for k = 1:3
                    ch = rgb(:, :, k);
                    ch(m) = col(k);
                    rgb(:, :, k) = ch;
                end
                fa = obj.Alpha;
                if obj.ParamOn, fa = 0; end      % outlines only over a parameter map
                if obj.ShowOutline || obj.ParamOn
                    e = GRASPilot.maskEdge(m);
                    A(m) = fa;
                    A(e) = 1;
                else
                    A(m) = fa;
                end
            end
            obj.OvlRGB{v} = rgb;
            obj.OvlA{v} = A;
        end

        function refreshPreviews(obj)
            for v = 1:3
                obj.refreshPreview(v);
            end
        end

        function refreshPreview(obj, v)
            obj.computePreview(v);
            obj.composeView(v);
        end

        function computePreview(obj, v)
            obj.PrevA{v} = [];
            if ~obj.Preview || isempty(obj.Data), return; end
            M = obj.thrSliceMap(v);
            if isempty(M), return; end
            obj.PrevA{v} = 0.38 * double(obj.inWindow(M));
        end

        function computeParamOverlay(obj, v)
            obj.ParamA{v} = []; obj.ParamRGB{v} = [];
            if ~obj.ParamOn || isempty(obj.Data), return; end
            V = obj.paramVolume(obj.ParamKey);
            if isempty(V), return; end
            M = double(obj.screenSliceOf(v, V));
            c = obj.paramCLim();
            ok = isfinite(M);
            if obj.ParamHideBelow, ok = ok & M >= c(1); end
            if obj.ParamHeadOnly, ok = ok & obj.screenSliceOf(v, obj.headMask()); end
            u = (M - c(1)) / (c(2) - c(1));
            u(~isfinite(u)) = 0;
            u = min(1, max(0, u));
            cm = GRASPilot.colormapByName(obj.ParamCmap, 256);
            idx = 1 + round(u * 255);
            obj.ParamRGB{v} = reshape(cm(idx(:), :), [size(M) 3]);
            obj.ParamA{v} = obj.ParamAlpha * double(ok);
        end

        function updateColorbar(obj)
            if ~isfield(obj.CB, 'Ax') || ~ishghandle(obj.CB.Ax), return; end
            if ~obj.ParamOn || isempty(obj.Data)
                set([obj.CB.Ax obj.CB.Img], 'Visible', 'off');
                set([obj.CB.Hi obj.CB.Lo obj.CB.Title], 'Visible', 'off');
                return
            end
            cm = GRASPilot.colormapByName(obj.ParamCmap, 128);
            set(obj.CB.Img, 'CData', reshape(cm, [128 1 3]), 'Visible', 'on');
            set(obj.CB.Ax, 'Visible', 'off', 'XLim', [0.5 1.5], 'YLim', [0.5 128.5]);
            c = obj.paramCLim();
            k = find(strcmp(obj.ParamKeys, obj.ParamKey));
            set(obj.CB.Hi, 'String', sprintf('%.3g', c(2)), 'Visible', 'on');
            set(obj.CB.Lo, 'String', sprintf('%.3g', c(1)), 'Visible', 'on');
            set(obj.CB.Title, 'String', sprintf('%s [%s]', obj.ParamKeys{k}, obj.paramUnit(obj.ParamKey)), ...
                'Visible', 'on');
        end

        function composeView(obj, v)
            % grayscale (window/level) + threshold preview + ROI overlay -> one RGB image
            S = obj.ViewImg{v};
            if isempty(S) || isempty(obj.V) || ~ishghandle(obj.V(v).Img), return; end
            c = obj.CLim;
            g = (double(S) - c(1)) / (c(2) - c(1));
            g(~isfinite(g)) = 0;
            g = min(1, max(0, g));
            rgb = cat(3, g, g, g);
            pa = obj.ParamA{v};
            if obj.ParamOn && ~isempty(pa) && isequal(size(pa), size(g)) && any(pa(:))
                PR = obj.ParamRGB{v};
                for k = 1:3
                    rgb(:, :, k) = rgb(:, :, k) .* (1 - pa) + PR(:, :, k) .* pa;
                end
            end
            a = obj.PrevA{v};
            if ~isempty(a) && isequal(size(a), size(g)) && any(a(:))
                pc = obj.Col.prev;
                for k = 1:3
                    rgb(:, :, k) = rgb(:, :, k) .* (1 - a) + pc(k) * a;
                end
            end
            A = obj.OvlA{v};
            if ~isempty(A) && isequal(size(A), size(g)) && any(A(:))
                O = obj.OvlRGB{v};
                for k = 1:3
                    rgb(:, :, k) = rgb(:, :, k) .* (1 - A) + O(:, :, k) .* A;
                end
            end
            set(obj.V(v).Img, 'CData', rgb);
        end

        function refreshCross(obj, v)
            sz = obj.V(v).Size;
            p = obj.canon2screen(v, obj.Pos);
            vis = 'off';
            if obj.ShowCross, vis = 'on'; end
            set(obj.V(v).CrossH, 'XData', [0.5 sz(2) + 0.5], 'YData', [p(2) p(2)], 'Visible', vis);
            set(obj.V(v).CrossV, 'XData', [p(1) p(1)], 'YData', [0.5 sz(1) + 0.5], 'Visible', vis);
        end

        function hideCursors(obj)
            for v = 1:3
                if ishghandle(obj.V(v).Cursor), set(obj.V(v).Cursor, 'Visible', 'off'); end
            end
        end

        function drawCircleLine(obj, v, X, Y, rmm, col)
            [sx, sy] = obj.viewSpacing(v);
            th = linspace(0, 2 * pi, 61);
            set(obj.V(v).Cursor, 'XData', X + rmm / sx * cos(th), 'YData', Y + rmm / sy * sin(th), ...
                'Color', col, 'Visible', 'on');
        end

        % ------------------------------------------------------------ threshold maps
        function k = thrKey(obj)
            k = obj.ThrKeys{obj.ThrSrc};
        end

        function tf = thrNeedsMaps(obj)
            tf = any(strcmp(obj.thrKey(), {'frameSub', 'frameRE', 'peakEnh', 'peakRE'}));
        end

        function tf = thrUsesFrame(obj)
            tf = any(strcmp(obj.thrKey(), {'frame', 'frameSub', 'frameRE'}));
        end

        function M = thrSliceMap(obj, v)
            % threshold map on the current slice of view v (screen orientation)
            M = [];
            if isempty(obj.Data) || isempty(obj.FrameImg{v}), return; end
            F = double(obj.FrameImg{v});
            key = obj.thrKey();
            switch key
                case 'frame'
                    M = F;
                    return
                case 'sub'
                    M = double(obj.screenSliceOf(v, obj.subVolume()));
                    return
            end
            bf = obj.BaseFrames(1):obj.BaseFrames(2);
            if ~obj.Data.hasMaps(bf), return; end
            switch key
                case 'frameSub'
                    M = F - double(obj.screenSliceOf(v, obj.Data.canonMap('Baseline')));
                case 'frameRE'
                    B = double(obj.screenSliceOf(v, obj.Data.canonMap('Baseline')));
                    M = 100 * (F - B) ./ B;
                    M(~(B > 0)) = NaN;
                case 'peakEnh'
                    M = double(obj.screenSliceOf(v, obj.Data.canonMap('PeakEnh')));
                case 'peakRE'
                    M = double(obj.screenSliceOf(v, obj.Data.canonMap('RelEnh')));
            end
        end

        function vals = thrValues(obj, lin)
            lin = lin(:);
            switch obj.thrKey()
                case 'frame'
                    vals = double(obj.Data.valuesAtFrame(lin, obj.Frame));
                case 'sub'
                    S = obj.subVolume();
                    vals = double(S(lin));
                case 'frameSub'
                    B = obj.Data.canonMap('Baseline');
                    vals = double(obj.Data.valuesAtFrame(lin, obj.Frame)) - double(B(lin));
                case 'frameRE'
                    B = double(obj.Data.canonMap('Baseline'));
                    b = B(lin);
                    vals = 100 * (double(obj.Data.valuesAtFrame(lin, obj.Frame)) - b) ./ b;
                    vals(~(b > 0)) = NaN;
                case 'peakEnh'
                    M = obj.Data.canonMap('PeakEnh');
                    vals = double(M(lin));
                case 'peakRE'
                    M = obj.Data.canonMap('RelEnh');
                    vals = double(M(lin));
            end
        end

        function V = thrVolume(obj)
            % threshold map over the whole volume (canonical)
            switch obj.thrKey()
                case 'frame'
                    V = obj.Data.canonVolume(obj.Frame);
                case 'sub'
                    V = obj.subVolume();
                case 'frameSub'
                    V = obj.Data.canonVolume(obj.Frame) - obj.Data.canonMap('Baseline');
                case 'frameRE'
                    B = obj.Data.canonMap('Baseline');
                    V = 100 * (obj.Data.canonVolume(obj.Frame) - B) ./ B;
                    V(~(B > 0)) = NaN;
                case 'peakEnh'
                    V = obj.Data.canonMap('PeakEnh');
                case 'peakRE'
                    V = obj.Data.canonMap('RelEnh');
            end
        end

        function tf = inWindow(obj, vals)
            tf = vals >= obj.ThrLo & vals <= obj.ThrHi;
        end

        function ok = ensureThresholdReady(obj)
            ok = true;
            if obj.thrNeedsMaps(), ok = obj.ensureMaps(false); end
        end

        function updateThrSliderRange(obj)
            if isempty(obj.Data), return; end
            vals = [];
            for v = 1:3
                M = obj.thrSliceMap(v);
                if ~isempty(M), vals = [vals; M(isfinite(M))]; end %#ok<AGROW>
            end
            if isempty(vals)
                lo = 0; hi = 1;
            else
                s = sort(vals);
                lo = s(max(1, round(0.001 * numel(s))));
                hi = s(max(1, round(0.999 * numel(s))));
                if hi <= lo, hi = lo + 1; end
            end
            val = obj.ThrLo;
            if ~isfinite(val), val = lo; end
            GRASPilot.setSlider(obj.H.thrSlider, lo, hi, val);
        end

        function updateThrInfo(obj)
            if ~isfield(obj.H, 'thrInfo'), return; end
            r = obj.ActiveROI;
            s = '';
            if ~isempty(obj.Data) && r > 0 && r <= numel(obj.Masks) && (isfinite(obj.ThrLo) || isfinite(obj.ThrHi))
                lin = find(obj.Masks{r});
                if ~isempty(lin) && (~obj.thrNeedsMaps() || obj.Data.hasMaps(obj.BaseFrames(1):obj.BaseFrames(2)))
                    k = nnz(obj.inWindow(obj.thrValues(lin)));
                    s = sprintf('Window keeps %d / %d voxels of %s (%.0f%%)', k, numel(lin), ...
                        obj.ROIs(r).Name, 100 * k / numel(lin));
                end
            end
            set(obj.H.thrInfo, 'String', s);
        end

        % ------------------------------------------------------------ TIC
        function x = ticX(obj)
            if obj.TicX == 1 && obj.Data.hasTiming()
                x = obj.Data.frameTimes();
            else
                x = 1:obj.Data.NT;
            end
        end

        function computeStats(obj, r)
            lin = find(obj.Masks{r});
            st = struct('n', numel(lin), 'mean', [], 'sd', [], 'median', [], 'q1', [], 'q3', []);
            if st.n > 0
                TC = double(obj.Data.timecourses(lin));
                st.mean = mean(TC, 1);
                if st.n > 1, st.sd = std(TC, 0, 1); else, st.sd = zeros(1, size(TC, 2)); end
                Ts = sort(TC, 1);
                st.median = GRASPilot.quant(Ts, 0.5);
                st.q1 = GRASPilot.quant(Ts, 0.25);
                st.q3 = GRASPilot.quant(Ts, 0.75);
            end
            obj.ROIs(r).Stats = st;
            obj.ROIs(r).Dirty = false;
            obj.ROIs(r).Count = st.n;
        end

        function [y, lo, hi] = ticTransform(obj, y, lo, hi)
            bf = obj.BaseFrames(1):obj.BaseFrames(2);
            bf = bf(bf <= numel(y));
            s0 = mean(y(bf));
            switch obj.TicY
                case 2, f = @(z) 100 * (z - s0) / s0;
                case 3, f = @(z) z / s0;
                case 4, f = @(z) z - s0;
                otherwise, f = @(z) z;
            end
            y = f(y);
            if ~isempty(lo), lo = f(lo); hi = f(hi); end
        end

        function refreshTIC(obj)
            ax = obj.Tic.Ax;
            it = obj.Tic.Items;
            delete(it(ishghandle(it)));
            obj.Tic.Items = [];
            obj.Tic.FrameLine = [];
            obj.Tic.CursorLine = [];
            try, legend(ax, 'off'); catch, end
            if isempty(obj.Data), return; end
            x = obj.ticX();
            items = GRASPilot.handleArray(0);
            hl = GRASPilot.handleArray(0);
            names = {};
            ymin = Inf; ymax = -Inf;
            for r = 1:numel(obj.ROIs)
                if ~obj.ROIs(r).Visible, continue; end
                if obj.ROIs(r).Dirty || isempty(obj.ROIs(r).Stats), obj.computeStats(r); end
                st = obj.ROIs(r).Stats;
                if st.n == 0, continue; end
                col = obj.ROIs(r).Color;
                switch obj.TicStat
                    case 1, y = st.mean; lo = st.mean - st.sd; hi = st.mean + st.sd;
                    case 2, y = st.median; lo = st.q1; hi = st.q3;
                    otherwise, y = st.mean; lo = []; hi = [];
                end
                [y, lo, hi] = obj.ticTransform(y, lo, hi);
                if ~isempty(lo) && numel(x) > 1
                    p = patch('Parent', ax, 'XData', [x fliplr(x)], 'YData', [lo fliplr(hi)], ...
                        'FaceColor', col, 'FaceAlpha', 0.15, 'EdgeColor', 'none', 'HitTest', 'off');
                    items(end + 1) = p; %#ok<AGROW>
                    ymin = min([ymin lo]); ymax = max([ymax hi]);
                end
                l = line('Parent', ax, 'XData', x, 'YData', y, 'Color', col, 'LineWidth', 1.8, ...
                    'Marker', 'none', 'HitTest', 'off');
                items(end + 1) = l; %#ok<AGROW>
                hl(end + 1) = l; %#ok<AGROW>
                names{end + 1} = sprintf('%s (n=%d)', obj.ROIs(r).Name, st.n); %#ok<AGROW>
                ymin = min([ymin y]); ymax = max([ymax y]);
            end
            if obj.ShowCursorTIC && ~isempty(obj.CursorTC)
                y = obj.ticTransform(obj.CursorTC, [], []);
                l = line('Parent', ax, 'XData', x, 'YData', y, 'Color', [0.75 0.75 0.75], 'LineWidth', 1, ...
                    'LineStyle', '--', 'HitTest', 'off');
                items(end + 1) = l;
                hl(end + 1) = l;
                names{end + 1} = sprintf('Cursor voxel [%d %d %d]', obj.Pos);
                obj.Tic.CursorLine = l;
                ymin = min([ymin y]); ymax = max([ymax y]);
            end
            if ~isfinite(ymin) || ~isfinite(ymax)
                ymin = 0; ymax = 1;
            end
            pad = 0.06 * max(ymax - ymin, eps);
            yl = [ymin - pad, ymax + pad];
            if yl(2) <= yl(1), yl = yl(1) + [-1 1]; end
            if numel(x) > 1, xl = [x(1) x(end)]; else, xl = x(1) + [-0.5 0.5]; end
            set(ax, 'XLim', xl, 'YLim', yl);
            xf = x(obj.Frame);
            fl = line('Parent', ax, 'XData', [xf xf], 'YData', yl, 'Color', [1 1 1] * 0.85, ...
                'LineStyle', ':', 'LineWidth', 1, 'HitTest', 'off');
            items(end + 1) = fl;
            ke = obj.Data.analysisEnd();
            if ke < obj.Data.NT                 % analysis endpoint
                items(end + 1) = line('Parent', ax, 'XData', [x(ke) x(ke)], 'YData', yl, ...
                    'Color', [1 0.6 0.2], 'LineStyle', '--', 'LineWidth', 1, 'HitTest', 'off');
            end
            obj.Tic.FrameLine = fl;
            obj.Tic.Items = items;
            if obj.TicX == 1 && obj.Data.hasTiming(), xs = 'Time (s)'; else, xs = 'Frame'; end
            ylabs = {'Signal (a.u.)', 'Relative enhancement (%)', 'S / S0', 'S - S0'};
            xlabel(ax, xs, 'Color', [0.85 0.85 0.85]);
            ylabel(ax, ylabs{obj.TicY}, 'Color', [0.85 0.85 0.85]);
            if ~isempty(hl)
                try
                    legend(ax, hl, names, 'Interpreter', 'none', 'TextColor', [0.92 0.92 0.92], ...
                        'Color', [0.14 0.14 0.16], 'EdgeColor', [0.4 0.4 0.4], 'Location', 'northwest', ...
                        'FontSize', 8, 'AutoUpdate', 'off');
                catch
                    try, legend(ax, hl, names, 'Location', 'northwest'); catch, end
                end
            end
        end

        function updateTicFrameLine(obj)
            fl = obj.Tic.FrameLine;
            if isempty(fl) || ~ishghandle(fl), return; end
            x = obj.ticX();
            xf = x(obj.Frame);
            set(fl, 'XData', [xf xf]);
        end

        function updateCursorTC(obj, light)
            if ~obj.ShowCursorTIC, return; end
            n = obj.Data.canonSize();
            lin = sub2ind(n, obj.Pos(1), obj.Pos(2), obj.Pos(3));
            obj.CursorTC = double(obj.Data.timecourses(lin));
            cl = obj.Tic.CursorLine;
            if light && ~isempty(cl) && ishghandle(cl)
                set(cl, 'YData', obj.ticTransform(obj.CursorTC, [], []));
            else
                obj.refreshTIC();
            end
        end

        function ticClick(obj)
            cp = get(obj.Tic.Ax, 'CurrentPoint');
            [~, f] = min(abs(obj.ticX() - cp(1, 1)));
            obj.setFrame(f);
        end

        % ------------------------------------------------------------ mouse / keyboard
        function [v, X, Y] = pointerView(obj)
            v = 0; X = NaN; Y = NaN;
            cp = get(obj.Fig, 'CurrentPoint');
            for k = 1:3
                r = obj.ViewRect(k, :);
                if cp(1) >= r(1) && cp(1) <= r(1) + r(3) && cp(2) >= r(2) && cp(2) <= r(2) + r(4)
                    v = k;
                    [X, Y] = obj.axesPoint(k);
                    return
                end
            end
        end

        function [X, Y] = axesPoint(obj, v)
            cp = get(obj.V(v).Ax, 'CurrentPoint');
            X = cp(1, 1); Y = cp(1, 2);
        end

        function tf = insideImage(obj, v, X, Y)
            sz = obj.V(v).Size;
            tf = X >= 0.5 && X <= sz(2) + 0.5 && Y >= 0.5 && Y <= sz(1) + 0.5;
        end

        function tf = pointerInTIC(obj)
            tf = false;
            try
                r = getpixelposition(obj.Tic.Ax, true);
                cp = get(obj.Fig, 'CurrentPoint');
                tf = cp(1) >= r(1) && cp(1) <= r(1) + r(3) && cp(2) >= r(2) && cp(2) <= r(2) + r(4);
            catch
            end
        end

        function onMouseDown(obj)
            if isempty(obj.Data) || obj.Busy, return; end
            sel = get(obj.Fig, 'SelectionType');
            [v, X, Y] = obj.pointerView();
            if v == 0
                if strcmp(sel, 'normal') && obj.pointerInTIC()
                    obj.Drag = struct('type', 'tic');
                    obj.ticClick();
                end
                return
            end
            obj.ActiveView = v;
            switch sel
                case 'alt'
                    obj.Drag = struct('type', 'wl', 'v', v, 'p0', get(obj.Fig, 'CurrentPoint'), 'clim0', obj.CLim);
                case 'extend'
                    obj.Drag = struct('type', 'pan', 'v', v, 'p0', [X Y]);
                case 'normal'
                    if ~obj.insideImage(v, X, Y), return; end
                    switch obj.Tool
                        case 'nav'
                            obj.Drag = struct('type', 'nav', 'v', v);
                            obj.navTo(v, X, Y, true);
                        case {'brush', 'erase'}
                            if ~obj.ensureActiveROI(), return; end
                            er = strcmp(obj.Tool, 'erase');
                            if ~er && obj.Constrain && ~obj.ensureThresholdReady(), return; end
                            obj.beginAction();
                            obj.Drag = struct('type', 'paint', 'v', v, 'last', [X Y], 'erase', er);
                            obj.paintShape(v, X, Y, obj.BrushR, er);
                            obj.refreshOverlays();
                        case 'circle'
                            if ~obj.ensureActiveROI(), return; end
                            if obj.Constrain && ~obj.ensureThresholdReady(), return; end
                            obj.Drag = struct('type', 'circle', 'v', v, 'c', [X Y], 'r', 0);
                            obj.drawCircleLine(v, X, Y, 0, obj.ROIs(obj.ActiveROI).Color);
                    end
            end
        end

        function onMouseMove(obj)
            if isempty(obj.Data), return; end
            d = obj.Drag;
            if ~isstruct(d) || ~isfield(d, 'type') || isempty(d.type)
                obj.hover();
                return
            end
            switch d.type
                case 'nav'
                    [X, Y] = obj.axesPoint(d.v);
                    obj.navTo(d.v, X, Y, true);
                case 'paint'
                    [X, Y] = obj.axesPoint(d.v);
                    [sx, sy] = obj.viewSpacing(d.v);
                    p0 = d.last;
                    dist = hypot((X - p0(1)) * sx, (Y - p0(2)) * sy);
                    step = max(0.4 * obj.BrushR, 0.5 * min(sx, sy));
                    ns = max(1, ceil(dist / step));
                    for s = 1:ns
                        p = p0 + ([X Y] - p0) * s / ns;
                        obj.paintShape(d.v, p(1), p(2), obj.BrushR, d.erase);
                    end
                    obj.Drag.last = [X Y];
                    obj.refreshOverlays();
                    obj.drawCircleLine(d.v, X, Y, obj.BrushR, obj.cursorColor());
                case 'circle'
                    [X, Y] = obj.axesPoint(d.v);
                    [sx, sy] = obj.viewSpacing(d.v);
                    r = hypot((X - d.c(1)) * sx, (Y - d.c(2)) * sy);
                    obj.Drag.r = r;
                    obj.drawCircleLine(d.v, d.c(1), d.c(2), r, obj.ROIs(obj.ActiveROI).Color);
                    s = '';
                    if obj.Use3D, s = ' (sphere)'; end
                    obj.setStatus(sprintf('Circle radius %.1f mm%s', r, s));
                case 'wl'
                    cp = get(obj.Fig, 'CurrentPoint');
                    dd = cp - d.p0;
                    w0 = max(diff(d.clim0), eps);
                    c0 = mean(d.clim0);
                    w = w0 * exp(dd(1) / 150);
                    c = c0 - dd(2) / 250 * w0;
                    obj.setCLim([c - w / 2, c + w / 2]);
                case 'pan'
                    ax = obj.V(d.v).Ax;
                    cp = get(ax, 'CurrentPoint');
                    delta = d.p0 - cp(1, 1:2);
                    set(ax, 'XLim', get(ax, 'XLim') + delta(1), 'YLim', get(ax, 'YLim') + delta(2));
                case 'tic'
                    obj.ticClick();
            end
        end

        function onMouseUp(obj)
            d = obj.Drag;
            obj.Drag = struct('type', '');
            if ~isstruct(d) || ~isfield(d, 'type') || isempty(d.type), return; end
            switch d.type
                case 'paint'
                    lbl = 'Brush';
                    if d.erase, lbl = 'Erase'; end
                    obj.commitAction(lbl);
                    obj.afterMaskChange();
                case 'circle'
                    obj.hideCursors();
                    [sx, sy] = obj.viewSpacing(d.v);
                    r = max(d.r, 0.5 * min(sx, sy));
                    obj.beginAction();
                    obj.paintShape(d.v, d.c(1), d.c(2), r, false);
                    obj.commitAction('Circle');
                    obj.afterMaskChange();
                case 'nav'
                    if obj.ShowCursorTIC, obj.refreshTIC(); end
            end
        end

        function navTo(obj, v, X, Y, light)
            n = obj.Data.canonSize();
            c = max(1, min(n, round(obj.screen2canonCont(v, X, Y))));
            obj.setPos(c);
            obj.updateCursorTC(light);
            obj.statusVoxel(v, c);
        end

        function hover(obj)
            [v, X, Y] = obj.pointerView();
            if v == 0 || ~obj.insideImage(v, X, Y)
                obj.hideCursors();
                return
            end
            for k = 1:3
                if k ~= v, set(obj.V(k).Cursor, 'Visible', 'off'); end
            end
            if any(strcmp(obj.Tool, {'brush', 'erase'}))
                obj.drawCircleLine(v, X, Y, obj.BrushR, obj.cursorColor());
            else
                set(obj.V(v).Cursor, 'Visible', 'off');
            end
            n = obj.Data.canonSize();
            obj.statusVoxel(v, max(1, min(n, round(obj.screen2canonCont(v, X, Y)))));
        end

        function statusVoxel(obj, v, c)
            sz = obj.V(v).Size;
            p = round(obj.canon2screen(v, c));
            val = NaN;
            if p(2) >= 1 && p(2) <= sz(1) && p(1) >= 1 && p(1) <= sz(2) && ~isempty(obj.ViewImg{v})
                val = obj.ViewImg{v}(p(2), p(1));
            end
            n = obj.Data.canonSize();
            lin = sub2ind(n, c(1), c(2), c(3));
            inr = {};
            for r = 1:numel(obj.Masks)
                if obj.Masks{r}(lin), inr{end + 1} = obj.ROIs(r).Name; end %#ok<AGROW>
            end
            rs = '';
            if ~isempty(inr), rs = ['   in: ' strjoin(inr, ', ')]; end
            tl = '';
            if obj.Data.hasTiming()
                ft = obj.Data.frameTimes();
                tl = sprintf('  t = %.1f s', ft(obj.Frame));
            end
            lab = 'S';
            if strcmp(obj.ViewMode, 'sub'), lab = sprintf('Sub(%d-%d)', obj.SubFrames); end
            ps = '';
            if obj.ParamOn
                PV = obj.paramVolume(obj.ParamKey);
                if ~isempty(PV), ps = sprintf('   %s = %.4g %s', obj.ParamKey, PV(lin), obj.paramUnit(obj.ParamKey)); end
            end
            obj.setStatus(sprintf('%s  [i j k] = [%d %d %d]   %s = %.4g   frame %d%s%s%s', obj.ViewNames{v}, ...
                c(1), c(2), c(3), lab, val, obj.Frame, tl, ps, rs));
        end

        function c = cursorColor(obj)
            c = [1 1 0.4];
            if strcmp(obj.Tool, 'erase')
                c = [1 0.45 0.45];
            elseif obj.ActiveROI > 0
                c = obj.ROIs(obj.ActiveROI).Color;
            end
        end

        function onScroll(obj, evt)
            if isempty(obj.Data) || obj.Busy, return; end
            [v, X, Y] = obj.pointerView();
            if v == 0, return; end
            obj.ActiveView = v;
            cnt = evt.VerticalScrollCount;
            mods = get(obj.Fig, 'CurrentModifier');
            if any(strcmp(mods, 'control')) || any(strcmp(mods, 'command'))
                obj.zoomView(v, X, Y, 1.2 ^ cnt);
            else
                obj.moveSlice(v, -cnt);
            end
        end

        function zoomView(obj, v, X, Y, f)
            ax = obj.V(v).Ax;
            sz = obj.V(v).Size;
            xl = X + (get(ax, 'XLim') - X) * f;
            yl = Y + (get(ax, 'YLim') - Y) * f;
            if diff(xl) >= sz(2) && diff(yl) >= sz(1)
                xl = [0.5 sz(2) + 0.5]; yl = [0.5 sz(1) + 0.5];
            end
            if diff(xl) < 4 || diff(yl) < 4, return; end
            set(ax, 'XLim', xl, 'YLim', yl);
        end

        function onKey(obj, evt)
            if isempty(obj.Data) || obj.Busy, return; end
            k = evt.Key;
            m = evt.Modifier;
            ctrl = any(strcmp(m, 'control')) || any(strcmp(m, 'command'));
            shift = any(strcmp(m, 'shift'));
            step = 1;
            if shift, step = 10; end
            switch k
                case 'rightarrow', obj.setFrame(obj.Frame + step);
                case 'leftarrow', obj.setFrame(obj.Frame - step);
                case 'uparrow', obj.moveSlice(obj.ActiveView, step);
                case 'downarrow', obj.moveSlice(obj.ActiveView, -step);
                case 'pageup', obj.setFrame(obj.Frame + 10);
                case 'pagedown', obj.setFrame(obj.Frame - 10);
                case 'home', obj.setFrame(1);
                case 'end', obj.setFrame(obj.Data.NT);
                case 'space', obj.togglePlay();
                case 'z'
                    if ctrl && shift, obj.redo(); elseif ctrl, obj.undo(); end
                case 'y'
                    if ctrl, obj.redo(); end
                case 'n', if ~ctrl, obj.setTool('nav'); end
                case 'c', if ~ctrl, obj.setTool('circle'); end
                case 'b', if ~ctrl, obj.setTool('brush'); end
                case 'e', if ~ctrl, obj.setTool('erase'); end
                case 'd', if ~ctrl, obj.set3D(~obj.Use3D); end
                case 'p', if ~ctrl, obj.setPreview(~obj.Preview); end
                case 's', if ~ctrl, obj.toggleSubView(); end
                case 'g', if ~ctrl, obj.growFromCrosshair(); end
                case 'm', if ~ctrl, obj.setParamDisplay(~obj.ParamOn); end
                case 'leftbracket', obj.setBrushR(obj.BrushR / 1.25);
                case 'rightbracket', obj.setBrushR(obj.BrushR * 1.25);
            end
        end

        % ------------------------------------------------------------ control callbacks
        function openParamWindow(obj)
            if ~obj.requireData(), return; end
            if isfield(obj.PW, 'Fig') && ishghandle(obj.PW.Fig)
                figure(obj.PW.Fig);
                obj.refreshParamWindow();
                return
            end
            obj.Layout = obj.Layout(arrayfun(@(L) ishghandle(L.h), obj.Layout));
            C = obj.Col;
            W = 1220; Hh = 860;
            ss = get(0, 'ScreenSize');
            f = figure('Name', 'Semi-quantitative DCE parameters', 'NumberTitle', 'off', 'MenuBar', 'none', ...
                'ToolBar', 'none', 'Units', 'pixels', 'Color', C.bg, ...
                'Position', [max(1, round((ss(3) - W) / 2)), max(40, round((ss(4) - Hh) / 2)), W, Hh], ...
                'KeyPressFcn', @(s, e) obj.onKey(e));
            P = uipanel('Parent', f, 'Units', 'pixels', 'Position', [1 1 W Hh], 'BorderType', 'none', ...
                'BackgroundColor', C.panel);
            obj.PW = struct('Fig', f, 'Panel', P);
            hdr = @(x, s) obj.ctl(P, x, 10, 340, 18, 'text', s, 'FontWeight', 'bold', 'ForegroundColor', C.accent);
            x1 = 12; x2 = 400; x3 = 788;
            % curve parameters
            hdr(x1, 'CURVE PARAMETERS');
            obj.ctl(P, x1, 36, 104, 18, 'text', 'Baseline frames');
            obj.PW.base1 = obj.ctl(P, x1 + 108, 34, 40, 22, 'edit', '', 'Callback', @(s, e) obj.onParamBaseline());
            obj.ctl(P, x1 + 152, 36, 20, 18, 'text', 'to');
            obj.PW.base2 = obj.ctl(P, x1 + 174, 34, 40, 22, 'edit', '', 'Callback', @(s, e) obj.onParamBaseline());
            obj.ctl(P, x1, 64, 96, 18, 'text', 'Wash-in slope');
            obj.PW.washIn = obj.ctl(P, x1 + 100, 62, 240, 22, 'popupmenu', ...
                {'dS / (TP - To), onset to peak', 'Maximum up-slope', 'Max slope in SSS wash-in period'}, ...
                'TooltipString', ['SSS wash-in period: frames around the steepest rise of the SSS reference curve ' ...
                '(Lee, Berger et al. 2023); the tumor slope is its maximum slope in that period'], ...
                'Callback', @(s, e) obj.onParamOptions());
            obj.ctl(P, x1, 92, 96, 18, 'text', 'Washout slope');
            obj.PW.washOut = obj.ctl(P, x1 + 100, 90, 240, 22, 'popupmenu', ...
                {'Regression, peak to end', 'Peak to last frame (2 points)'}, 'Callback', @(s, e) obj.onParamOptions());
            obj.ctl(P, x1, 120, 100, 18, 'text', 'AUC window (s)');
            obj.PW.aucW = obj.ctl(P, x1 + 100, 118, 50, 22, 'edit', '', ...
                'TooltipString', 'Seconds after the end of the baseline; blank = whole curve; 60 or 90 = IAUC60/IAUC90', ...
                'Callback', @(s, e) obj.onParamOptions());
            obj.ctl(P, x1 + 156, 120, 220, 18, 'text', 'blank = whole curve; 60 = IAUC60', 'ForegroundColor', C.dim);
            obj.ctl(P, x1, 148, 100, 18, 'text', 'Analyze up to');
            obj.PW.endV = obj.ctl(P, x1 + 100, 146, 50, 22, 'edit', '', ...
                'TooltipString', 'Analysis endpoint: later frames are ignored for all parameters; blank = all frames', ...
                'Callback', @(s, e) obj.onParamOptions());
            obj.PW.endU = obj.ctl(P, x1 + 156, 146, 80, 22, 'popupmenu', {'s', 'frames'}, ...
                'Callback', @(s, e) obj.onParamOptions());
            obj.ctl(P, x1 + 242, 148, 130, 18, 'text', 'blank = all frames', 'ForegroundColor', C.dim);
            obj.ctl(P, x1, 176, 100, 18, 'text', 'Early / late (s)');
            obj.PW.early = obj.ctl(P, x1 + 100, 174, 50, 22, 'edit', '', ...
                'TooltipString', 'Early frame, time after the end of the baseline (SER, early enhancement, curve type)', ...
                'Callback', @(s, e) obj.onParamOptions());
            obj.PW.late = obj.ctl(P, x1 + 156, 174, 50, 22, 'edit', '', ...
                'TooltipString', 'Late frame, time after the end of the baseline; blank = last analysed frame', ...
                'Callback', @(s, e) obj.onParamOptions());
            obj.ctl(P, x1 + 212, 176, 160, 18, 'text', 'late blank = last frame', 'ForegroundColor', C.dim);
            obj.ctl(P, x1, 204, 96, 18, 'text', 'SSS reference');
            obj.PW.ref = obj.ctl(P, x1 + 100, 202, 240, 22, 'popupmenu', {'(auto)'}, ...
                'Callback', @(s, e) obj.onParamRef());
            obj.ctl(P, x1, 232, 170, 28, 'pushbutton', 'Compute parameters', 'FontWeight', 'bold', ...
                'Callback', @(s, e) obj.onParamCompute());
            obj.PW.info = obj.ctl(P, x1, 264, 376, 58, 'text', '', 'FontSize', 8, 'ForegroundColor', C.dim);
            % parameter map
            hdr(x2, 'PARAMETER MAP');
            obj.ctl(P, x2, 36, 70, 18, 'text', 'Parameter');
            obj.PW.key = obj.ctl(P, x2 + 74, 34, 292, 22, 'popupmenu', obj.ParamLabels, ...
                'Callback', @(s, e) obj.setParamKey(get(s, 'Value')));
            obj.PW.show = obj.ctl(P, x2, 62, 160, 24, 'togglebutton', 'Show on images (M)', ...
                'Callback', @(s, e) obj.setParamDisplay(get(s, 'Value')));
            obj.PW.cmap = obj.ctl(P, x2 + 170, 63, 196, 22, 'popupmenu', obj.CmapLabels, ...
                'Callback', @(s, e) obj.onParamCmap(get(s, 'Value')));
            obj.ctl(P, x2, 94, 30, 18, 'text', 'Min');
            obj.PW.cmin = obj.ctl(P, x2 + 32, 92, 70, 22, 'edit', '', 'Callback', @(s, e) obj.onParamCLimEdit());
            obj.ctl(P, x2 + 106, 94, 32, 18, 'text', 'Max');
            obj.PW.cmax = obj.ctl(P, x2 + 140, 92, 70, 22, 'edit', '', 'Callback', @(s, e) obj.onParamCLimEdit());
            obj.ctl(P, x2 + 214, 92, 74, 22, 'pushbutton', 'Auto head', ...
                'TooltipString', '1st-99.5th percentile inside the head', 'Callback', @(s, e) obj.autoParamCLim('head'));
            obj.ctl(P, x2 + 292, 92, 74, 22, 'pushbutton', 'Auto ROIs', ...
                'TooltipString', '2nd-98th percentile inside the tumor ROIs', 'Callback', @(s, e) obj.autoParamCLim('rois'));
            ax = axes('Parent', P, 'Units', 'pixels', 'Position', [x2 100 366 14], 'Visible', 'off');
            obj.PW.cbImg = image('Parent', ax, 'CData', zeros(1, 2, 3), 'HitTest', 'off');
            obj.PW.cbAx = ax;
            obj.reg(ax, P, x2, 122, 366, 14);
            obj.PW.cbLo = obj.ctl(P, x2, 138, 110, 16, 'text', '', 'FontSize', 8);
            obj.PW.cbHi = obj.ctl(P, x2 + 216, 138, 150, 16, 'text', '', 'FontSize', 8, 'HorizontalAlignment', 'right');
            obj.ctl(P, x2, 162, 60, 18, 'text', 'Opacity');
            obj.PW.alpha = obj.ctl(P, x2 + 64, 162, 302, 18, 'slider', '', 'Min', 0, 'Max', 1, 'Value', obj.ParamAlpha, ...
                'Callback', @(s, e) obj.onParamFlag('ParamAlpha', get(s, 'Value')));
            obj.PW.hide = obj.ctl(P, x2, 188, 150, 20, 'checkbox', 'Hide below min', 'Value', obj.ParamHideBelow, ...
                'Callback', @(s, e) obj.onParamFlag('ParamHideBelow', get(s, 'Value')));
            obj.PW.head = obj.ctl(P, x2 + 160, 188, 180, 20, 'checkbox', 'Inside head only', 'Value', obj.ParamHeadOnly, ...
                'Callback', @(s, e) obj.onParamFlag('ParamHeadOnly', get(s, 'Value')));
            % statistics and export
            hdr(x3, 'TUMOR ROI STATISTICS');
            obj.PW.inclRef = obj.ctl(P, x3, 36, 300, 20, 'checkbox', 'Include the reference (SSS) ROI', ...
                'Value', obj.IncludeRefStats, 'Callback', @(s, e) obj.onParamFlag('IncludeRefStats', get(s, 'Value')));
            obj.ctl(P, x3, 62, 80, 24, 'pushbutton', 'Refresh', 'Callback', @(s, e) obj.refreshParamWindow());
            obj.ctl(P, x3 + 86, 62, 334, 24, 'pushbutton', 'Export statistics + voxel values...', ...
                'Callback', @(s, e) obj.onExportParams());
            obj.ctl(P, x3, 92, 420, 24, 'pushbutton', 'Export parameter maps (NIfTI)...', ...
                'Callback', @(s, e) obj.onExportMaps());
            obj.ctl(P, x3, 124, 420, 112, 'text', { ...
                'S0 = mean of the baseline frames; times from the end of the baseline', ...
                'PE = Smax - S0;  relative PE = 100 PE / S0', ...
                'Wash-in = PE / (TP - To), or the maximum up-slope', ...
                'Washout = slope from the peak to the end (negative = washout)', ...
                'AUC = area of S - S0 after the baseline (window)', ...
                '/SSS = divided by the value of the SSS ROI mean curve', ...
                'ROI mean curve = parameter of the ROI mean TIC'}, 'FontSize', 8, 'ForegroundColor', C.dim);
            obj.PW.table = uitable('Parent', P, 'Units', 'pixels', 'Position', [12 12 W - 24 400], ...
                'ColumnName', {'ROI', 'Parameter', 'Unit', 'N', 'Mean', 'SD', 'Median', 'Q1', 'Q3', 'P10', ...
                'P90', 'Min', 'Max', 'ROI curve'}, 'RowName', [], 'Data', cell(0, 14), 'FontSize', 9, ...
                'ColumnWidth', {110, 195, 55, 50, 78, 78, 78, 78, 78, 78, 78, 78, 78, 78});
            obj.reg(obj.PW.table, P, 12, 328, -12, -12);
            set(f, 'SizeChangedFcn', @(s, e) obj.layoutParamWindow());
            obj.layoutParamWindow();
            obj.refreshParamWindow();
        end

        function layoutParamWindow(obj)
            if ~isfield(obj.PW, 'Fig') || ~ishghandle(obj.PW.Fig), return; end
            fp = get(obj.PW.Fig, 'Position');
            set(obj.PW.Panel, 'Position', [1 1 fp(3) fp(4)]);
            obj.layoutRecords(obj.PW.Panel);
        end

        function layoutRecords(obj, panel)
            for i = 1:numel(obj.Layout)
                L = obj.Layout(i);
                if ~ishghandle(L.h) || ~ishghandle(L.panel), continue; end
                if nargin > 1 && ~isequal(L.panel, panel), continue; end
                pp = get(L.panel, 'Position');
                pw = pp(3); ph = pp(4);
                x = L.x;
                if x < 0, x = pw + x; end
                w = L.w;
                if w <= 0, w = pw - x + w; end
                ht = L.ht;
                if ht <= 0, ht = ph - L.top + ht; end
                set(L.h, 'Position', [max(1, x), ph - L.top - ht, max(w, 4), max(ht, 4)]);
            end
        end

        function refreshParamWindow(obj)
            if ~isfield(obj.PW, 'Fig') || ~ishghandle(obj.PW.Fig) || isempty(obj.Data), return; end
            D = obj.Data;
            set(obj.PW.base1, 'String', num2str(obj.BaseFrames(1)));
            set(obj.PW.base2, 'String', num2str(obj.BaseFrames(2)));
            o = D.ParamOpts;
            set(obj.PW.washIn, 'Value', max(1, find(strcmpi({'onset', 'max', 'sss'}, o.WashIn))));
            if isfinite(o.EndFrame)
                set(obj.PW.endV, 'String', num2str(o.EndFrame)); set(obj.PW.endU, 'Value', 2);
            elseif isfinite(o.EndTime)
                set(obj.PW.endV, 'String', num2str(o.EndTime)); set(obj.PW.endU, 'Value', 1);
            else
                set(obj.PW.endV, 'String', '');
            end
            set(obj.PW.early, 'String', num2str(o.EarlyTime));
            if isfinite(o.LateTime), set(obj.PW.late, 'String', num2str(o.LateTime));
            else, set(obj.PW.late, 'String', ''); end
            set(obj.PW.washOut, 'Value', 1 + strcmpi(o.WashOut, 'endpoints'));
            if isfinite(o.AUCWindow), set(obj.PW.aucW, 'String', num2str(o.AUCWindow));
            else, set(obj.PW.aucW, 'String', ''); end
            names = {obj.ROIs.Name};
            r = obj.refIndex();
            items = [{'(auto: first ROI named SSS)'}, names];
            val = 1;
            if ~isempty(obj.RefROIName) && r > 0, val = r + 1; end
            set(obj.PW.ref, 'String', items, 'Value', val);
            ready = D.hasMaps(obj.BaseFrames(1):obj.BaseFrames(2));
            if ~ready
                msgI = {'Parameters not computed for the current settings:', 'press Compute parameters.'};
            else
                ts = 's';
                if ~D.hasTiming(), ts = 'frames (enter s/frame for seconds)'; end
                [ref, rr] = obj.refParams();
                if isempty(ref)
                    msgI = {sprintf('Computed (time unit: %s).', ts), ...
                        'No SSS reference ROI: /SSS values unavailable (press Auto SSS).'};
                else
                    msgI = {sprintf('Computed (time unit: %s). Reference: %s (%d voxels)', ts, ...
                        obj.ROIs(rr).Name, nnz(obj.Masks{rr})), sprintf('PE %.3g, wash-in %.3g, washout %.3g, AUC %.3g', ...
                        ref.PE, ref.WashIn, ref.WashOut, ref.AUC)};
                end
            end
            ke = D.analysisEnd(); tp = D.paramTimes();
            if ke < D.NT
                msgI{end + 1} = sprintf('Analysis limited to frames 1-%d of %d (last at %.4g %s).', ke, D.NT, ...
                    tp(ke), GRASPilot.ternaryStr(D.hasTiming(), 's', 'frames'));
            end
            if strcmpi(o.WashIn, 'sss') && numel(o.WashInFrames) == 2
                msgI{end + 1} = sprintf('SSS wash-in period: frames %d-%d (%.4g-%.4g).', o.WashInFrames, ...
                    tp(o.WashInFrames));
            end
            set(obj.PW.info, 'String', msgI);
            [~, rows] = obj.paramStatsTable();
            dsp = rows;
            for i = 1:size(rows, 1)
                for j = 4:size(rows, 2)
                    if j == 4
                        dsp{i, j} = sprintf('%d', rows{i, j});
                    else
                        dsp{i, j} = sprintf('%.4g', rows{i, j});
                        if numel(dsp{i, j}) > 9, dsp{i, j} = sprintf('%.3g', rows{i, j}); end
                        if isnan(rows{i, j}), dsp{i, j} = ''; end
                    end
                end
            end
            set(obj.PW.table, 'Data', dsp);
            obj.syncParamWindow();
        end

        function syncParamWindow(obj)
            if ~isfield(obj.PW, 'Fig') || ~ishghandle(obj.PW.Fig) || isempty(obj.Data), return; end
            k = find(strcmp(obj.ParamKeys, obj.ParamKey));
            set(obj.PW.key, 'Value', k);
            set(obj.PW.show, 'Value', double(obj.ParamOn));
            set(obj.PW.cmap, 'Value', find(strcmp(obj.CmapNames, obj.ParamCmap)));
            c = obj.paramCLim();
            set(obj.PW.cmin, 'String', sprintf('%.4g', c(1)));
            set(obj.PW.cmax, 'String', sprintf('%.4g', c(2)));
            cm = GRASPilot.colormapByName(obj.ParamCmap, 128);
            set(obj.PW.cbImg, 'CData', reshape(cm, [1 128 3]), 'XData', [1 128], 'YData', [1 1]);
            set(obj.PW.cbAx, 'XLim', [0.5 128.5], 'YLim', [0.5 1.5]);
            u = obj.paramUnit(obj.ParamKey);
            set(obj.PW.cbLo, 'String', sprintf('%.4g %s', c(1), u));
            set(obj.PW.cbHi, 'String', sprintf('%.4g %s', c(2), u));
        end

        function onParamBaseline(obj)
            set(obj.H.base1, 'String', get(obj.PW.base1, 'String'));
            set(obj.H.base2, 'String', get(obj.PW.base2, 'String'));
            obj.onBaseline();
            obj.refreshParamWindow();
        end

        function onParamOptions(obj)
            wi = {'onset', 'max', 'sss'}; wo = {'regression', 'endpoints'};
            ev = str2double(get(obj.PW.endV, 'String'));
            ef = Inf; et = Inf;
            if isfinite(ev) && ev > 0
                if get(obj.PW.endU, 'Value') == 2 || ~obj.Data.hasTiming()
                    ef = ev;                      % no timing: the endpoint is in frames
                else
                    et = ev;
                end
            end
            obj.Data.setParamOptions('EndFrame', ef, 'EndTime', et, ...
                'EarlyTime', str2double(get(obj.PW.early, 'String')), ...
                'LateTime', str2double(get(obj.PW.late, 'String')));
            aw = str2double(get(obj.PW.aucW, 'String'));
            if ~isfinite(aw) || aw <= 0, aw = Inf; end
            obj.Data.setParamOptions('WashIn', wi{get(obj.PW.washIn, 'Value')}, ...
                'WashOut', wo{get(obj.PW.washOut, 'Value')}, 'AUCWindow', aw);
            obj.syncRefWindow();
            obj.refreshTIC();
            obj.refreshParamWindow();
            if obj.ParamOn, obj.refreshParamDisplay(); end
        end

        function onParamRef(obj)
            v = get(obj.PW.ref, 'Value');
            if v == 1
                obj.RefROIName = '';
            else
                obj.RefROIName = obj.ROIs(v - 1).Name;
            end
            obj.RefCache = [];
            obj.ParamCache = struct('key', '', 'V', []);
            obj.refreshParamDisplay();
            obj.refreshParamWindow();
        end

        function onParamCompute(obj)
            obj.onParamOptions();
            if obj.computeParams()
                obj.refreshParamDisplay();
                obj.refreshParamWindow();
            end
        end

        function onParamCmap(obj, k)
            obj.ParamCmap = obj.CmapNames{k};
            obj.refreshParamDisplay();
        end

        function onParamCLimEdit(obj)
            c = [str2double(get(obj.PW.cmin, 'String')) str2double(get(obj.PW.cmax, 'String'))];
            if all(isfinite(c)), obj.setParamCLim(c); else, obj.syncParamWindow(); end
        end

        function onParamFlag(obj, name, val)
            switch name
                case 'ParamAlpha', obj.ParamAlpha = val;
                otherwise, obj.(name) = logical(val);
            end
            if strcmp(name, 'IncludeRefStats')
                obj.refreshParamWindow();
            else
                obj.refreshParamDisplay();
            end
        end

        function onExportParams(obj)
            if ~obj.requireData(), return; end
            suf = '_dce_parameters.xlsx';
            file = obj.putFile({'*.xlsx', 'Excel workbook (*.xlsx)'; '*.csv', 'CSV files (*.csv)'}, ...
                'Export parameter statistics and voxel values', suf);
            if isempty(file), return; end
            [~, ~, ext] = fileparts(file);
            if isempty(ext), file = [file '.xlsx']; end
            try
                files = obj.exportParams(file);
            catch err
                obj.showError(err, 'Parameter export failed');
                return
            end
            if ~isempty(files)
                obj.inform(sprintf('Wrote:\n%s', strjoin(files, char(10))), 'Parameter export');
            end
        end

        % ------------------------------------------------------------ motion
        function M = estimateMotion(obj, varargin)
            % Rigid frame-by-frame motion estimate of the raw frames
            % (graspilot_motion_estimate); options are passed on. Warns when
            % excess motion is found. Reference: mean of the baseline frames.
            M = [];
            if ~obj.requireData(), return; end
            D = obj.Data;
            rf = [D.BinIdx{obj.BaseFrames(1):obj.BaseFrames(2)}];
            opts = [{'RefFrames', rf}, varargin];
            if isfield(obj.MW, 'Fig') && ishghandle(obj.MW.Fig)
                a = str2double(get(obj.MW.warnAbs, 'String'));
                r = str2double(get(obj.MW.warnRel, 'String'));
                if isfinite(a) && a > 0, opts = [opts, {'WarnAbsMM', a}]; end
                if isfinite(r) && r > 0, opts = [opts, {'WarnRelMM', r}]; end
                if get(obj.MW.refPop, 'Value') == 2
                    k = round(str2double(get(obj.MW.refFrame, 'String')));
                    if isfinite(k) && k >= 1 && k <= D.NTRaw, opts = [opts, {'Reference', k}]; end
                end
            end
            obj.MotionStop = false;
            wb = obj.wbOpen('Estimating head motion...');
            obj.setMotionButtons(true);
            try
                M = graspilot_motion_estimate(D, opts{:}, 'Progress', @(f, m) obj.motionProgress(wb, f, m), ...
                    'Cancel', @() obj.MotionStop);
            catch err
                obj.wbClose(wb);
                obj.setMotionButtons(false);
                obj.showError(err, 'Motion estimate failed');
                return
            end
            obj.wbClose(wb);
            obj.setMotionButtons(false);
            D.setMotion(M);
            obj.updateDataInfo();
            obj.refreshMotionWindow();
            pre = '';
            if D.MotionCorrected, pre = 'Residual motion after correction. '; end
            obj.setStatus([pre M.Message]);
            if M.Excess
                obj.inform(sprintf(['%s%s\n\nCheck the frames listed (Motion window). Motion during the ' ...
                    'acquisition also degrades the reconstruction itself; Apply correction realigns the frames ' ...
                    'but cannot remove blurring or ghosting.'], pre, M.Message), 'Motion warning');
            end
        end

        function motionProgress(obj, wb, f, m)
            obj.wbUpdate(wb, f, m);
            if isfield(obj.MW, 'Fig') && ishghandle(obj.MW.Fig)
                set(obj.MW.info, 'String', m);
            end
            drawnow;
        end

        function ok = applyMotion(obj, M, quiet)
            % realign the frames with motion estimate M (default: the last
            % estimate); ROIs stay in place (frames move to the reference)
            ok = false;
            if ~obj.requireData(), return; end
            D = obj.Data;
            if nargin < 2 || isempty(M), M = D.Motion; end
            if nargin < 3, quiet = false; end
            if isempty(M)
                obj.inform('Estimate the motion first.', 'Motion correction');
                return
            end
            if D.MotionCorrected
                obj.inform('The frames are already motion-corrected (Revert first).', 'Motion correction');
                return
            end
            if ~quiet && ~obj.Quiet
                mb = D.motionCorrectionBytes() / 2^20;
                b = obj.ask(sprintf(['Realign all %d frames to the reference (%s)?\n\nThe corrected frames are ' ...
                    'kept in memory as a full 4D series: about %.0f MB (the original data are kept for Revert). ' ...
                    'Images, curves, maps and exports then use the corrected frames.'], D.NTRaw, M.Reference, mb), ...
                    'Motion correction', {'Apply', 'Cancel'}, 'Apply');
                if ~strcmp(b, 'Apply'), return; end
            end
            wb = obj.wbOpen('Motion correction...');
            try
                D.applyMotionCorrection(M, @(f, m) obj.wbUpdate(wb, f, m));
            catch err
                obj.wbClose(wb);
                obj.showError(err, 'Motion correction failed');
                return
            end
            obj.wbClose(wb);
            obj.afterFramesChanged();
            obj.setStatus(sprintf('Motion correction applied: %d frames realigned to the %s.', D.NTRaw, M.Reference));
            ok = true;
        end

        function revertMotion(obj)
            if ~obj.requireData() || ~obj.Data.MotionCorrected, return; end
            M = obj.Data.MotionApplied;
            obj.Data.revertMotionCorrection();
            obj.Data.setMotion(M);
            obj.afterFramesChanged();
            obj.setStatus('Motion correction removed: original frames.');
        end

        function afterFramesChanged(obj)
            % the frame source changed (motion correction): refresh everything
            obj.stopPlay();
            for r = 1:numel(obj.ROIs), obj.ROIs(r).Dirty = true; end
            obj.CursorTC = [];
            obj.ParamCache = struct('key', '', 'V', []);
            obj.RefCache = [];
            obj.HeadCache = struct('ver', -1, 'M', []);
            if obj.Preview && obj.thrNeedsMaps(), obj.ensureMaps(false); end
            if obj.ParamOn, obj.ensureMaps(false); end
            obj.updateDataInfo();
            obj.refreshAll();
            obj.updateThrSliderRange();
            obj.refreshParamDisplay();
            obj.refreshParamWindow();
            obj.refreshMotionWindow();
        end

        function files = exportMotion(obj, file)
            % motion parameters and displacements per raw frame (.csv)
            files = {};
            D = obj.Data;
            M = D.Motion;
            what = 'estimate';
            if isempty(M), M = D.MotionApplied; end
            if D.MotionCorrected
                if isequal(M, D.MotionApplied)
                    what = 'estimate applied for the motion correction';
                else
                    what = 'residual motion of the corrected frames';
                end
            end
            if isempty(M), error('GRASPilot:motion', 'No motion estimate.'); end
            fid = fopen(file, 'w');
            if fid < 0, error('GRASPilot:motion', 'Cannot write %s.', file); end
            fprintf(fid, ['RawFrame,Time_s,TransL_mm,TransP_mm,TransS_mm,RotLR_deg,RotAP_deg,RotIS_deg,' ...
                'MeanDispRef_mm,MaxDispRef_mm,MeanDispPrev_mm,MaxDispPrev_mm,Excess\n']);
            tt = D.t0 + (0:D.NTRaw - 1) * D.dt;
            for f = 1:D.NTRaw
                v = [M.Params(f, :) M.AbsMean(f) M.AbsMax(f) M.RelMean(f) M.RelMax(f)];
                c = arrayfun(@(x) GRASPilot.ternaryStr(isfinite(x), sprintf('%.4f', x), ''), v, 'UniformOutput', false);
                ts = '';
                if isfinite(tt(f)), ts = sprintf('%.3f', tt(f)); end
                fprintf(fid, '%d,%s,%s,%d\n', f, ts, strjoin(c, ','), any(M.ExcessFrames == f));
            end
            fprintf(fid, '# Motion %s. %s\n', what, M.Message);
            fprintf(fid, '# Translations: head motion along L (R->L), P (A->P), S (I->S); rotations about the L-R, A-P, I-S axes through the head centre [%.1f %.1f %.1f] mm (canonical voxel order).\n', M.Centre);
            fprintf(fid, '# Motion correction: %s\n', D.motionDescription());
            fclose(fid);
            files = {file};
        end

        function openMotionWindow(obj)
            if ~obj.requireData(), return; end
            if isfield(obj.MW, 'Fig') && ishghandle(obj.MW.Fig)
                figure(obj.MW.Fig);
                obj.refreshMotionWindow();
                return
            end
            obj.Layout = obj.Layout(arrayfun(@(L) ishghandle(L.h), obj.Layout));
            C = obj.Col;
            W = 980; Hh = 720;
            ss = get(0, 'ScreenSize');
            f = figure('Name', 'Motion check and correction', 'NumberTitle', 'off', 'MenuBar', 'none', ...
                'ToolBar', 'none', 'Units', 'pixels', 'Color', C.bg, ...
                'Position', [max(1, round((ss(3) - W) / 2)), max(40, round((ss(4) - Hh) / 2)), W, Hh]);
            P = uipanel('Parent', f, 'Units', 'pixels', 'Position', [1 1 W Hh], 'BorderType', 'none', ...
                'BackgroundColor', C.panel);
            obj.MW = struct('Fig', f, 'Panel', P);
            vmm = max(1, prod(obj.Data.canonVox()) ^ (1 / 3));
            obj.ctl(P, 12, 10, 300, 18, 'text', 'RIGID MOTION, FRAME BY FRAME', 'FontWeight', 'bold', ...
                'ForegroundColor', C.accent);
            obj.ctl(P, 12, 36, 70, 18, 'text', 'Reference');
            obj.MW.refPop = obj.ctl(P, 86, 34, 210, 22, 'popupmenu', {'Mean of the baseline frames', 'Raw frame:'});
            obj.MW.refFrame = obj.ctl(P, 300, 34, 46, 22, 'edit', '1');
            obj.ctl(P, 360, 36, 176, 18, 'text', 'Warn if mean displacement >');
            obj.MW.warnAbs = obj.ctl(P, 538, 34, 46, 22, 'edit', sprintf('%.2g', vmm), ...
                'TooltipString', 'mm from the reference (default: one voxel)');
            obj.ctl(P, 588, 36, 140, 18, 'text', 'mm from reference, or');
            obj.MW.warnRel = obj.ctl(P, 730, 34, 46, 22, 'edit', sprintf('%.2g', vmm / 2), ...
                'TooltipString', 'mm between consecutive frames (default: half a voxel)');
            obj.ctl(P, 780, 36, 160, 18, 'text', 'mm between frames');
            obj.MW.btnEst = obj.ctl(P, 12, 64, 150, 26, 'pushbutton', 'Estimate motion', 'FontWeight', 'bold', ...
                'Callback', @(s, e) obj.estimateMotion());
            obj.MW.btnStop = obj.ctl(P, 168, 64, 60, 26, 'pushbutton', 'Stop', 'Enable', 'off', ...
                'Callback', @(s, e) obj.stopMotion());
            obj.MW.btnApply = obj.ctl(P, 236, 64, 150, 26, 'pushbutton', 'Apply correction', ...
                'Callback', @(s, e) obj.applyMotion());
            obj.MW.btnRevert = obj.ctl(P, 392, 64, 150, 26, 'pushbutton', 'Revert to original', ...
                'Callback', @(s, e) obj.revertMotion());
            obj.MW.btnExport = obj.ctl(P, 548, 64, 150, 26, 'pushbutton', 'Export (.csv)...', ...
                'Callback', @(s, e) obj.onExportMotion());
            obj.MW.info = obj.ctl(P, 12, 98, -12, 66, 'text', '', 'FontSize', 9);
            ttl = {'Translation (mm)', 'Rotation (deg)', 'Mean head displacement (mm)'};
            obj.MW.Ax = GRASPilot.handleArray(3);
            for i = 1:3
                ax = axes('Parent', P, 'Units', 'pixels', 'Position', [60 100 400 120], 'Color', [0.07 0.07 0.08], ...
                    'XColor', [0.8 0.8 0.8], 'YColor', [0.8 0.8 0.8], 'Box', 'on', 'FontSize', 8, 'NextPlot', 'add');
                title(ax, ttl{i}, 'Color', [0.85 0.85 0.85], 'FontSize', 9);
                obj.MW.Ax(i) = ax;
            end
            set(f, 'SizeChangedFcn', @(s, e) obj.layoutMotionWindow());
            obj.layoutMotionWindow();
            obj.refreshMotionWindow();
        end

        function layoutMotionWindow(obj)
            if ~isfield(obj.MW, 'Fig') || ~ishghandle(obj.MW.Fig), return; end
            fp = get(obj.MW.Fig, 'Position');
            set(obj.MW.Panel, 'Position', [1 1 fp(3) fp(4)]);
            obj.layoutRecords(obj.MW.Panel);
            top = 190; gap = 46; bottom = 40;
            h = max(60, (fp(4) - top - bottom - 2 * gap) / 3);
            for i = 1:3
                y = fp(4) - top - i * h - (i - 1) * gap;
                set(obj.MW.Ax(i), 'Position', [60 y max(100, fp(3) - 80) h]);
            end
        end

        function setMotionButtons(obj, busy)
            if ~isfield(obj.MW, 'Fig') || ~ishghandle(obj.MW.Fig), return; end
            on = {'on', 'off'};
            set(obj.MW.btnStop, 'Enable', on{2 - busy});
            set([obj.MW.btnEst obj.MW.btnApply obj.MW.btnRevert obj.MW.btnExport], 'Enable', on{1 + busy});
            if ~busy
                D = obj.Data;
                if D.MotionCorrected, set(obj.MW.btnApply, 'Enable', 'off'); end
                if ~D.MotionCorrected, set(obj.MW.btnRevert, 'Enable', 'off'); end
                if isempty(D.Motion) && ~D.MotionCorrected
                    set([obj.MW.btnApply obj.MW.btnExport], 'Enable', 'off');
                end
            end
        end

        function stopMotion(obj)
            obj.MotionStop = true;
        end

        function refreshMotionWindow(obj)
            if ~isfield(obj.MW, 'Fig') || ~ishghandle(obj.MW.Fig) || isempty(obj.Data), return; end
            D = obj.Data;
            obj.setMotionButtons(false);
            M = D.Motion;
            hdr = '';
            if D.MotionCorrected
                hdr = 'FRAMES ARE MOTION-CORRECTED. ';
                if isempty(M) || isequal(M, D.MotionApplied)
                    M = D.MotionApplied;
                    hdr = [hdr 'Applied estimate: '];
                else
                    hdr = [hdr 'Residual motion of the corrected frames: '];
                end
            end
            for i = 1:3, cla(obj.MW.Ax(i)); end
            if isempty(M)
                set(obj.MW.info, 'String', {['Estimate motion registers every raw frame rigidly to the reference ' ...
                    '(normalized mutual information, robust to contrast enhancement).'], ...
                    'About 1-3 s per frame. Apply correction then realigns the frames (Revert restores them).'}, ...
                    'ForegroundColor', obj.Col.fg);
                return
            end
            col = [0.95 0.40 0.35; 0.40 0.85 0.45; 0.40 0.65 1.00];
            x = 1:D.NTRaw;
            xl = 'Raw frame';
            if D.hasTiming(), x = D.t0 + (0:D.NTRaw - 1) * D.dt; xl = 'Time (s)'; end
            lbl = {{'L (R->L)', 'P (A->P)', 'S (I->S)'}, {'about L-R', 'about A-P', 'about I-S'}};
            for i = 1:2
                ax = obj.MW.Ax(i);
                hl = zeros(1, 3);
                for j = 1:3
                    hl(j) = plot(ax, x, M.Params(:, 3 * (i - 1) + j), '-', 'Color', col(j, :), 'LineWidth', 1.2);
                end
                obj.motionLegend(ax, hl, lbl{i});
            end
            ax = obj.MW.Ax(3);
            h1 = plot(ax, x, M.AbsMean, '-', 'Color', [1 0.75 0.2], 'LineWidth', 1.5);
            h2 = plot(ax, x, M.RelMean, '-', 'Color', [0.7 0.7 1], 'LineWidth', 1.2);
            plot(ax, [x(1) x(end)], M.WarnAbsMM * [1 1], '--', 'Color', [1 0.75 0.2]);
            plot(ax, [x(1) x(end)], M.WarnRelMM * [1 1], ':', 'Color', [0.7 0.7 1]);
            ex = M.ExcessFrames;
            if ~isempty(ex)
                plot(ax, x(ex), M.AbsMean(ex), 'o', 'Color', [1 0.3 0.3], 'MarkerFaceColor', [1 0.3 0.3], 'MarkerSize', 4);
            end
            obj.motionLegend(ax, [h1 h2], {'from reference', 'from previous frame'});
            xlabel(ax, xl, 'Color', [0.85 0.85 0.85]);
            for i = 1:3
                if numel(x) > 1, set(obj.MW.Ax(i), 'XLim', [x(1) x(end)]); end
            end
            fc = obj.Col.fg;
            if M.Excess, fc = [1 0.55 0.45]; end
            fp = get(obj.MW.Fig, 'Position');
            set(obj.MW.info, 'String', GRASPilot.wrapText([hdr M.Message], floor((fp(3) - 30) / 6.6)), ...
                'ForegroundColor', fc);
        end

        function motionLegend(~, ax, h, lbl)
            try
                legend(ax, h, lbl, 'TextColor', [0.9 0.9 0.9], 'Color', [0.1 0.1 0.12], 'Location', 'northwest');
            catch
                legend(ax, h, lbl, 'Location', 'northwest');
            end
        end

        function onExportMotion(obj)
            if ~obj.requireData(), return; end
            file = obj.putFile({'*.csv', 'CSV file (*.csv)'}, 'Export motion parameters', '_motion.csv');
            if isempty(file), return; end
            try
                obj.exportMotion(file);
            catch err
                obj.showError(err, 'Motion export failed');
                return
            end
            obj.setStatus(['Wrote ' file]);
        end

        function toggleSubView(obj)
            if strcmp(obj.ViewMode, 'sub')
                obj.setViewMode('frame');
            else
                obj.setViewMode('sub');
            end
        end

        function onGrowR(obj)
            r = str2double(get(obj.H.growR, 'String'));
            if isfinite(r) && r > 0, obj.GrowR = r; end
            set(obj.H.growR, 'String', num2str(obj.GrowR));
        end

        function onBinningDialog(obj)
            if ~obj.requireData(), return; end
            D = obj.Data;
            a = inputdlg({sprintf(['Average every n consecutive frames (1 = original, %d raw frames):'], D.NTRaw), ...
                'Keep an incomplete last bin? (0 = drop, 1 = keep)'}, 'Temporal binning', 1, ...
                {num2str(D.BinSize), num2str(double(D.BinKeepPartial))});
            if isempty(a), return; end
            obj.applyBinning(str2double(a{1}), str2double(a{2}) == 1);
        end

        function onSlider(obj)
            if isempty(obj.Data), return; end
            obj.setFrame(get(obj.H.slider, 'Value'));
        end

        function onFps(obj)
            if ~isempty(obj.PlayTimer)
                try
                    if isvalid(obj.PlayTimer) && strcmp(get(obj.PlayTimer, 'Running'), 'on')
                        obj.stopPlay(); obj.startPlay();
                    end
                catch
                end
            end
        end

        function setFlag(obj, name, val)
            switch name
                case {'Exclusive', 'ShowOutline', 'ShowCross', 'ShowCursorTIC'}
                    obj.(name) = logical(val);
                otherwise
                    obj.(name) = val;
            end
            switch name
                case {'Alpha', 'ShowOutline'}
                    obj.refreshOverlays();
                case 'ShowCross'
                    set(obj.H.chkCross, 'Value', double(obj.ShowCross));
                    if ~isempty(obj.Data), for v = 1:3, obj.refreshCross(v); end, end
                case 'ShowCursorTIC'
                    if ~isempty(obj.Data)
                        if obj.ShowCursorTIC, obj.updateCursorTC(false); else, obj.refreshTIC(); end
                    end
                case {'TicY', 'TicStat', 'TicX'}
                    hn = struct('TicY', 'ticY', 'TicStat', 'ticStat', 'TicX', 'ticX');
                    set(obj.H.(hn.(name)), 'Value', val);
                    obj.refreshTIC();
            end
        end

        function onBaseline(obj)
            if isempty(obj.Data), return; end
            T = obj.Data.NT;
            a = round(str2double(get(obj.H.base1, 'String')));
            b = round(str2double(get(obj.H.base2, 'String')));
            if ~isfinite(a), a = obj.BaseFrames(1); end
            if ~isfinite(b), b = obj.BaseFrames(2); end
            a = max(1, min(T, a)); b = max(a, min(T, b));
            set(obj.H.base1, 'String', num2str(a));
            set(obj.H.base2, 'String', num2str(b));
            if isequal([a b], obj.BaseFrames), return; end
            obj.BaseFrames = [a b];
            if obj.thrNeedsMaps() && (obj.Preview || obj.Constrain)
                obj.ensureMaps(false);
            end
            obj.refreshPreviews();
            obj.refreshTIC();
            obj.updateThrInfo();
        end

        function onThrEdit(obj)
            lo = str2double(get(obj.H.thrLo, 'String'));
            hi = str2double(get(obj.H.thrHi, 'String'));
            obj.setThreshold(lo, hi);
        end

        function onThrSlider(obj)
            obj.setThreshold(get(obj.H.thrSlider, 'Value'), obj.ThrHi);
        end

        function onWindowEdit(obj)
            c = [str2double(get(obj.H.wMin, 'String')) str2double(get(obj.H.wMax, 'String'))];
            if any(~isfinite(c))
                obj.setCLim(obj.CLim);
            else
                obj.setCLim(c);
            end
        end

        function onTableEdit(obj, e)
            if isempty(e.Indices), return; end
            r = e.Indices(1); c = e.Indices(2);
            if r > numel(obj.ROIs), return; end
            switch c
                case 2
                    if ischar(e.NewData) && ~isempty(strtrim(e.NewData))
                        obj.ROIs(r).Name = strtrim(e.NewData);
                    end
                case 3
                    obj.ROIs(r).Visible = logical(e.NewData);
            end
            obj.refreshROIList();
            obj.refreshOverlays();
            obj.refreshTIC();
        end

        function onTableSelect(obj, e)
            if isempty(e.Indices), return; end
            r = e.Indices(1, 1);
            if r >= 1 && r <= numel(obj.ROIs) && r ~= obj.ActiveROI
                obj.setActiveROI(r);
            end
        end

        function onOpen(obj)
            [f, p] = uigetfile({'*.mat', 'MATLAB data (*.mat)'}, 'Open GRASP reconstruction', obj.startDir());
            if isequal(f, 0), return; end
            try
                obj.openFile(fullfile(p, f), obj.LoadOpts);
            catch err
                obj.showError(err, 'Could not load data');
            end
        end

        function onExportImage(obj)
            if ~obj.requireData(), return; end
            c = obj.ask('Export which frames?', 'Image -> NIfTI', ...
                {'All frames (4D)', 'Current frame (3D)', 'Cancel'}, 'All frames (4D)');
            if isempty(c) || strcmp(c, 'Cancel'), return; end
            if strncmp(c, 'All', 3)
                frames = 1:obj.Data.NT; suf = '_dce.nii.gz';
            else
                frames = obj.Frame; suf = sprintf('_frame%03d.nii.gz', obj.Frame);
            end
            file = obj.putFile({'*.nii.gz', 'Compressed NIfTI (*.nii.gz)'; '*.nii', 'NIfTI (*.nii)'}, ...
                'Save image as', suf);
            if isempty(file), return; end
            obj.exportImageNifti(GRASPilot.fixNii(file), frames);
        end

        function onExportMask(obj)
            if ~obj.requireROIs(), return; end
            file = obj.putFile({'*.nii.gz', 'Compressed NIfTI (*.nii.gz)'; '*.nii', 'NIfTI (*.nii)'}, ...
                'Save mask as', '_mask.nii.gz');
            if isempty(file), return; end
            b = obj.ask('Also save one binary mask per ROI?', 'Mask -> NIfTI', ...
                {'Label map only', 'Label map + per-ROI masks'}, 'Label map only');
            if isempty(b), return; end
            try
                obj.exportMaskNifti(GRASPilot.fixNii(file), strcmp(b, 'Label map + per-ROI masks'));
            catch err
                obj.showError(err, 'Mask export failed');
            end
        end

        function onExportTS(obj)
            if ~obj.requireROIs(), return; end
            suf = '_roi_timeseries.xlsx';
            if obj.Data.BinSize > 1, suf = sprintf('_roi_timeseries_%dframes.xlsx', obj.Data.NT); end
            file = obj.putFile({'*.xlsx', 'Excel workbook (*.xlsx)'; '*.csv', 'CSV files (*.csv)'}, ...
                'Save ROI voxel time series', suf);
            if isempty(file), return; end
            [~, ~, ext] = fileparts(file);
            if isempty(ext), file = [file '.xlsx']; end
            files = obj.exportTimeseries(file);
            if ~isempty(files)
                obj.inform(sprintf('Wrote:\n%s', strjoin(files, char(10))), 'Time-series export');
            end
        end

        function onExportMaps(obj)
            if ~obj.requireData(), return; end
            folder = uigetdir(obj.startDir(), 'Folder for enhancement maps');
            if isequal(folder, 0), return; end
            [~, b] = fileparts(obj.Data.SourceFile);
            if isempty(b) || strcmp(b, '(workspace)'), b = 'grasp'; end
            try
                obj.exportMaps(folder, b);
            catch err
                obj.showError(err, 'Map export failed');
            end
        end

        function onSavePNG(obj)
            if ~obj.requireData(), return; end
            file = obj.putFile({'*.png', 'PNG image (*.png)'}, 'Save TIC plot', '_tic.png');
            if isempty(file), return; end
            try
                if exist('exportgraphics', 'file') == 2
                    exportgraphics(obj.Tic.Ax, file, 'Resolution', 200, 'BackgroundColor', 'current');
                else
                    f2 = figure('Visible', 'off', 'Color', [0.07 0.07 0.08], 'Position', [100 100 900 450]);
                    lg = findobj(obj.Fig, 'Type', 'legend');
                    copyobj([lg(:); obj.Tic.Ax], f2);
                    set(findobj(f2, 'Type', 'axes'), 'Units', 'normalized', 'Position', [0.08 0.12 0.88 0.8]);
                    set(f2, 'InvertHardcopy', 'off');
                    print(f2, file, '-dpng', '-r200');
                    delete(f2);
                end
                obj.setStatus(['Saved ' file]);
            catch err
                obj.showError(err, 'Could not save the plot');
            end
        end

        function ok = onSaveROIs(obj)
            ok = false;
            if ~obj.requireROIs(), return; end
            file = obj.putFile({'*.mat', 'ROI session (*.mat)'}, 'Save ROIs', '_rois.mat');
            if isempty(file), return; end
            try
                obj.saveROIs(file);
                ok = true;
            catch err
                obj.showError(err, 'Could not save ROIs');
            end
        end

        function onLoadROIs(obj)
            if ~obj.requireData(), return; end
            [f, p] = uigetfile({'*.mat', 'ROI session (*.mat)'}, 'Load ROIs', obj.startDir());
            if isequal(f, 0), return; end
            try
                obj.loadROIs(fullfile(p, f));
            catch err
                obj.showError(err, 'Could not load ROIs');
            end
        end

        function onImportMask(obj)
            if ~obj.requireData(), return; end
            [f, p] = uigetfile({'*.nii;*.nii.gz', 'NIfTI (*.nii, *.nii.gz)'}, 'Import mask', obj.startDir());
            if isequal(f, 0), return; end
            try
                obj.importMaskNifti(fullfile(p, f));
            catch err
                obj.showError(err, 'Could not import mask');
            end
        end

        function showHelp(obj)
            s = {'MOUSE (over any view)'
                '  Left button      active tool (Navigate / Circle / Brush / Eraser)'
                '  Right-drag       window/level (left-right = width, up-down = level)'
                '  Shift-drag       pan (middle button also works)'
                '  Wheel            change slice;  Ctrl+wheel = zoom'
                '  Click TIC plot   jump to that frame'
                ''
                'KEYS'
                '  Left / Right     previous / next frame (Shift = 10)'
                '  Home / End       first / last frame;  PgUp / PgDn = 10 frames'
                '  Up / Down        slice in the last-used view'
                '  Space            play / pause'
                '  N C B E          Navigate, Circle, Brush, Eraser'
                '  [  ]             brush radius down / up'
                '  D                toggle 3D (sphere) drawing'
                '  P                toggle threshold preview'
                '  S                show subtraction map (post - pre frame)'
                '  G                grow from crosshair (threshold window)'
                '  M                parameter colour map on / off'
                '  Ctrl+Z / Ctrl+Y  undo / redo (Ctrl+Shift+Z also redoes)'
                ''
                'ENHANCING TUMOR FROM THE SUBTRACTION MAP (last - first frame)'
                '  1. Map = Subtraction (post = last, pre = first), Show to view it'
                '  2. + Tumor, set Low (type a value, drag the slider, or draw a rough'
                '     sphere and press Otsu); tick Preview to see the green window'
                '  3. Navigate to the enhancing rim, then Grow from crosshair (connected'
                '     voxels within max mm), or draw a sphere and press Keep: 3D'
                '  All -> ROI adds every in-window voxel (includes vessels).'
                ''
                'ORIENTATION'
                '  View > Geometry: axis codes (direction of increasing rows, cols, slices),'
                '  presets (LPS = scanner order of the GRASP recon), Rotate 90 / Flip buttons'
                ''
                'SSS CONTROL ROI'
                '  Auto SSS: posterior superior sagittal sinus from just above the torcula'
                '  up toward the vertex (up to 50 mm above its thickest part)'
                '  (transverse and straight sinuses excluded); check and edit as needed'
                ''
                'MOTION (Analysis > Motion check and correction)'
                '  Estimate: rigid registration of every raw frame to the baseline mean;'
                '  warns when the mean head displacement exceeds the thresholds;'
                '  Apply correction realigns the frames (Revert restores them)'
                ''
                'SEMI-QUANTITATIVE PARAMETERS (Analysis menu or Parameters... button)'
                '  Analyze up to: limit the analysis to the first N s or frames'
                '  (orange dashed line on the curves)'
                '  Peak enhancement, wash-in and washout slopes, AUC, TTP, absolute and'
                '  normalized to the SSS ROI mean curve; colour maps (M), ROI statistics,'
                '  export of statistics and voxel values'
                ''
                'TEMPORAL BINNING'
                '  Average every n frames + Apply binning: each frame becomes the mean'
                '  signal of n consecutive frames; curves, maps and exports follow.'};
            if obj.Quiet
                fprintf('%s\n', s{:});
            else
                msgbox(s, 'Mouse and keyboard');
            end
        end

        % ------------------------------------------------------------ playback
        function togglePlay(obj)
            running = false;
            try
                running = ~isempty(obj.PlayTimer) && isvalid(obj.PlayTimer) && strcmp(get(obj.PlayTimer, 'Running'), 'on');
            catch
            end
            if running
                obj.stopPlay();
            else
                obj.startPlay();
            end
        end

        function startPlay(obj)
            if isempty(obj.Data) || obj.Data.NT < 2
                set(obj.H.btnPlay, 'Value', 0);
                return
            end
            fps = str2double(get(obj.H.fps, 'String'));
            if ~isfinite(fps) || fps <= 0, fps = 8; end
            fps = min(fps, 30);
            try
                if isempty(obj.PlayTimer) || ~isvalid(obj.PlayTimer)
                    obj.PlayTimer = timer('ExecutionMode', 'fixedSpacing', 'BusyMode', 'drop', ...
                        'TimerFcn', @(s, e) obj.playStep());
                end
                set(obj.PlayTimer, 'Period', max(0.001, round(1000 / fps) / 1000));
                start(obj.PlayTimer);
                set(obj.H.btnPlay, 'Value', 1, 'String', 'Pause');
            catch
                set(obj.H.btnPlay, 'Value', 0, 'String', 'Play');
                obj.setStatus('Playback is not available in this environment.');
            end
        end

        function stopPlay(obj)
            try
                if ~isempty(obj.PlayTimer) && isvalid(obj.PlayTimer), stop(obj.PlayTimer); end
            catch
            end
            if isfield(obj.H, 'btnPlay') && ishghandle(obj.H.btnPlay)
                set(obj.H.btnPlay, 'Value', 0, 'String', 'Play');
            end
        end

        function playStep(obj)
            if isempty(obj.Fig) || ~ishghandle(obj.Fig) || isempty(obj.Data)
                obj.stopPlay();
                return
            end
            if obj.Busy || (isstruct(obj.Drag) && isfield(obj.Drag, 'type') && ~isempty(obj.Drag.type)), return; end
            f = obj.Frame + 1;
            if f > obj.Data.NT, f = 1; end
            obj.setFrame(f);
            drawnow;
        end

        % ------------------------------------------------------------ misc helpers
        function setupTimeSlider(obj)
            T = obj.Data.NT;
            h = obj.H.slider;
            if T > 1
                GRASPilot.setSlider(h, 1, T, obj.Frame);
                set(h, 'SliderStep', [1 / (T - 1), min(1, 10 / (T - 1))], 'Enable', 'on');
            else
                GRASPilot.setSlider(h, 0, 1, 0);
                set(h, 'Enable', 'off');
            end
            set(obj.H.frameTotal, 'String', sprintf('/ %d', T));
        end

        function updateTimeWidgets(obj)
            if isempty(obj.Data), return; end
            if obj.Data.NT > 1 && abs(get(obj.H.slider, 'Value') - obj.Frame) > 1e-9
                set(obj.H.slider, 'Value', obj.Frame);
            end
            set(obj.H.frameEdit, 'String', num2str(obj.Frame));
            if obj.Data.hasTiming()
                ft = obj.Data.frameTimes();
                set(obj.H.timeText, 'String', sprintf('t = %.2f s', ft(obj.Frame)));
                set(obj.H.dtEdit, 'String', sprintf('%.4g', obj.Data.dt));
            else
                set(obj.H.timeText, 'String', '(time unknown)');
                set(obj.H.dtEdit, 'String', '');
            end
        end

        function updateDataInfo(obj)
            D = obj.Data;
            [~, nm, ext] = fileparts(D.SourceFile);
            v = D.VoxNative;
            n = D.NativeSize;
            if D.hasTiming(), dts = sprintf('%.3g s/frame', D.effectiveDt()); else, dts = 's/frame unknown'; end
            if strcmp(D.Mode, 'subspace'), md = sprintf('subspace K=%d', D.NK); else, md = 'dense'; end
            if D.MotionCorrected
                md = 'motion-corrected';
            elseif ~isempty(D.Motion) && D.Motion.Excess
                md = [md ', MOTION WARNING'];
            end
            if D.BinSize > 1
                ts = sprintf('%d frames (%dx binned from %d)', D.NT, D.BinSize, D.NTRaw);
            else
                ts = sprintf('%d frames', D.NT);
            end
            s = sprintf('%s%s\n%dx%dx%d %s (%s), %.3gx%.3gx%.3g mm, %s\n%s, %s, %.0f MB', nm, ext, n(1), n(2), n(3), ...
                D.Plane, D.AxisCodes, v(1), v(2), v(3), md, ts, dts, D.memoryBytes() / 2^20);
            set(obj.H.info, 'String', s);
        end

        function refreshROIList(obj)
            n = numel(obj.ROIs);
            if n == 0
                set(obj.H.roiPop, 'String', {'(none)'}, 'Value', 1);
                set(obj.H.swatch, 'BackgroundColor', [0.3 0.3 0.3]);
                set(obj.H.roiTable, 'Data', cell(0, 5));
                obj.ActiveROI = 0;
                return
            end
            vmL = 1;
            if ~isempty(obj.Data), vmL = prod(obj.Data.canonVox()) / 1000; end
            str = cell(1, n);
            data = cell(n, 5);
            bgc = zeros(n, 3);
            for r = 1:n
                cnt = nnz(obj.Masks{r});
                obj.ROIs(r).Count = cnt;
                str{r} = sprintf('%d: %s', r, obj.ROIs(r).Name);
                data(r, :) = {r, obj.ROIs(r).Name, logical(obj.ROIs(r).Visible), cnt, round(cnt * vmL * 1e4) / 1e4};
                bgc(r, :) = 0.45 * obj.ROIs(r).Color + 0.55;
            end
            ar = max(1, min(n, obj.ActiveROI));
            obj.ActiveROI = ar;
            set(obj.H.roiPop, 'String', str, 'Value', ar);
            set(obj.H.swatch, 'BackgroundColor', obj.ROIs(ar).Color);
            set(obj.H.roiTable, 'Data', data);
            if ~exist('OCTAVE_VERSION', 'builtin')
                try
                    set(obj.H.roiTable, 'BackgroundColor', bgc, 'RowStriping', 'on', 'ForegroundColor', [0 0 0]);
                catch
                end
            end
        end

        function afterMaskChange(obj)
            obj.Modified = true;
            obj.RefCache = [];
            obj.ParamCache = struct('key', '', 'V', []);
            obj.syncRefWindow();
            obj.refreshROIList();
            if obj.ParamOn
                [~, ~, nrm] = GRASPilot.paramInfo(obj.ParamKey);
                if nrm, obj.refreshParamDisplay(); end
            end
            obj.refreshOverlays();
            obj.refreshTIC();
            obj.updateThrInfo();
            obj.refreshParamWindow();
        end

        function applyChanges(obj, ch, inverse)
            for c = numel(ch):-1:1
                e = ch(c);
                if e.roi > numel(obj.Masks), continue; end
                val = e.val;
                if inverse, val = ~val; end
                obj.Masks{e.roi}(double(e.idx)) = val;
                obj.ROIs(e.roi).Dirty = true;
            end
        end

        function ok = ensureActiveROI(obj)
            ok = true;
            if isempty(obj.ROIs)
                obj.addROI('ROI 1');
            end
            if obj.ActiveROI == 0, obj.setActiveROI(1); end
            if ~obj.ROIs(obj.ActiveROI).Visible
                obj.ROIs(obj.ActiveROI).Visible = true;
                obj.refreshROIList();
            end
        end

        function c = nextColor(obj)
            P = obj.Palette;
            if isempty(obj.ROIs)
                c = P(1, :);
                return
            end
            used = reshape([obj.ROIs.Color], 3, [])';
            for i = 1:size(P, 1)
                if ~any(all(abs(used - repmat(P(i, :), size(used, 1), 1)) < 1e-3, 2))
                    c = P(i, :);
                    return
                end
            end
            c = P(mod(numel(obj.ROIs), size(P, 1)) + 1, :);
        end

        function writeLabelTables(obj, base)
            vmL = prod(obj.Data.canonVox()) / 1000;
            fid = fopen([base '_labels.txt'], 'w');
            fprintf(fid, '# GRASPilot label table\n# label\tname\tR\tG\tB\tvoxels\tvolume_mL\n');
            for r = 1:numel(obj.ROIs)
                c = round(255 * obj.ROIs(r).Color);
                k = nnz(obj.Masks{r});
                fprintf(fid, '%d\t%s\t%d\t%d\t%d\t%d\t%.4f\n', r, obj.ROIs(r).Name, c, k, k * vmL);
            end
            fclose(fid);
            fid = fopen([base '_itksnap_labels.txt'], 'w');
            fprintf(fid, ['################################################\n' ...
                '# ITK-SnAP Label Description File\n' ...
                '#    IDX   -R-  -G-  -B-  -A--  VIS MSH  LABEL\n' ...
                '################################################\n']);
            fprintf(fid, '%5d %5d %4d %4d %8.2f %2d %2d    "%s"\n', 0, 0, 0, 0, 0, 0, 0, 'Clear Label');
            for r = 1:numel(obj.ROIs)
                c = round(255 * obj.ROIs(r).Color);
                fprintf(fid, '%5d %5d %4d %4d %8.2f %2d %2d    "%s"\n', r, c(1), c(2), c(3), 1, 1, 1, ...
                    strrep(obj.ROIs(r).Name, '"', ''''));
            end
            fclose(fid);
        end

        function file = putFile(obj, filt, ttl, suffix)
            file = '';
            [~, b] = fileparts(obj.Data.SourceFile);
            if isempty(b) || strcmp(b, '(workspace)'), b = 'grasp'; end
            [f, p] = uiputfile(filt, ttl, fullfile(obj.startDir(), [b suffix]));
            if isequal(f, 0), return; end
            file = fullfile(p, f);
            obj.LastDir = p;
            obj.savePrefs();
        end

        function d = startDir(obj)
            d = obj.LastDir;
            if isempty(d) || exist(d, 'dir') ~= 7, d = pwd; end
        end

        function ok = requireData(obj)
            ok = ~isempty(obj.Data);
            if ~ok, obj.inform('Open a reconstruction first.', 'GRASPilot'); end
        end

        function ok = requireROIs(obj)
            ok = obj.requireData();
            if ok
                ok = ~isempty(obj.Masks) && any(cellfun(@(m) any(m(:)), obj.Masks));
                if ~ok, obj.inform('There are no ROI voxels yet.', 'GRASPilot'); end
            end
        end

        function b = ask(obj, q, ttl, opts, def)
            if obj.Quiet
                b = def;
                return
            end
            if numel(opts) == 2
                b = questdlg(q, ttl, opts{1}, opts{2}, def);
            else
                b = questdlg(q, ttl, opts{1}, opts{2}, opts{3}, def);
            end
            if isempty(b), b = ''; end
        end

        function inform(obj, msg, ttl)
            if obj.Quiet
                fprintf('[%s] %s\n', ttl, msg);
            else
                msgbox(msg, ttl, 'modal');
            end
        end

        function showError(obj, err, ttl)
            fprintf(2, '%s: %s\n', ttl, err.message);
            for k = 1:min(4, numel(err.stack))
                fprintf(2, '    at %s (line %d)\n', err.stack(k).name, err.stack(k).line);
            end
            if obj.Quiet, return; end
            errordlg(sprintf('%s:\n\n%s', ttl, err.message), 'GRASPilot', 'modal');
        end

        function setStatus(obj, s)
            if isfield(obj.H, 'status') && ishghandle(obj.H.status)
                set(obj.H.status, 'String', s);
            end
        end

        function wb = wbOpen(obj, msg)
            obj.Busy = true;
            wb = [];
            if obj.Quiet, return; end
            try
                wb = waitbar(0, msg, 'Name', 'GRASPilot');
                set(findall(wb, 'Type', 'text'), 'Interpreter', 'none');
            catch
                wb = [];
            end
        end

        function wbUpdate(~, wb, f, msg)
            if ~isempty(wb) && ishghandle(wb)
                waitbar(max(0, min(1, f)), wb, msg);
            end
        end

        function wbClose(obj, wb)
            obj.Busy = false;
            if ~isempty(wb) && ishghandle(wb), delete(wb); end
        end

        function o = mergeOpts(obj, nv)
            o = obj.LoadOpts;
            fn = fieldnames(o);
            if mod(numel(nv), 2) ~= 0
                error('GRASPilot:opts', 'Options must be name/value pairs.');
            end
            given = {};
            for i = 1:2:numel(nv)
                m = find(strcmpi(fn, nv{i}), 1);
                if isempty(m)
                    error('GRASPilot:opts', 'Unknown option "%s".', nv{i});
                end
                o.(fn{m}) = nv{i + 1};
                given{end + 1} = fn{m}; %#ok<AGROW>
            end
            % an explicit Plane/SwapInPlane/Flip (old style) overrides the axis codes
            if ~any(strcmp(given, 'AxisCodes')) && any(ismember({'Plane', 'SwapInPlane', 'Flip'}, given))
                o.AxisCodes = '';
            end
        end

        function loadPrefs(obj)
            try
                if ispref('GRASPilot', 'Geometry')
                    g = getpref('GRASPilot', 'Geometry');
                    f = {'AxisCodes', 'Voxel', 'dt', 't0'};
                    for i = 1:numel(f)
                        if isfield(g, f{i}), obj.LoadOpts.(f{i}) = g.(f{i}); end
                    end
                end
                if ispref('GRASPilot', 'LastDir')
                    obj.LastDir = getpref('GRASPilot', 'LastDir');
                end
                if ~ispref('GRASPilot', 'Geometry') && ispref('GraspDCEViewer', 'Geometry')
                    % settings saved under the tool's former name
                    g = getpref('GraspDCEViewer', 'Geometry');
                    f = {'AxisCodes', 'Voxel', 'dt', 't0'};
                    for i = 1:numel(f)
                        if isfield(g, f{i}), obj.LoadOpts.(f{i}) = g.(f{i}); end
                    end
                end
            catch
            end
        end

        function savePrefs(obj)
            if obj.Quiet, return; end
            try
                g = struct('AxisCodes', obj.LoadOpts.AxisCodes, 'Voxel', obj.LoadOpts.Voxel, ...
                    'dt', obj.LoadOpts.dt, 't0', obj.LoadOpts.t0);
                setpref('GRASPilot', 'Geometry', g);
                setpref('GRASPilot', 'LastDir', obj.LastDir);
            catch
            end
        end
    end

    % =====================================================================
    methods (Static)
        function o = defaultLoadOpts()
            o = struct('Format', 'auto', 'CoefVar', 'vk_all', 'BasisVar', 'basis_all', 'ImageVar', '', ...
                'DimOrder', 'xytz', 'CropFactor', [], 'Normalize', [], 'Plane', 'axial', ...
                'SwapInPlane', false, 'Flip', [false false false], 'Voxel', [1 1 1], 'dt', NaN, ...
                't0', 0, 'NoDialog', false, 'BinSize', 1, 'AxisCodes', 'LPS');
        end

        function L = wrapText(s, w)
            % word-wrap a string into a cell array of lines of <= w characters
            L = {};
            words = regexp(s, '\s+', 'split');
            cur = '';
            for i = 1:numel(words)
                if isempty(cur)
                    cur = words{i};
                elseif numel(cur) + 1 + numel(words{i}) <= w
                    cur = [cur ' ' words{i}]; %#ok<AGROW>
                else
                    L{end + 1} = cur; %#ok<AGROW>
                    cur = words{i};
                end
            end
            if ~isempty(cur), L{end + 1} = cur; end
        end

        function s = ternaryStr(c, a, b)
            if c, s = a; else, s = b; end
        end

        function [name, fld, norm] = paramInfo(key)
            % map name in GRASPilotData, field of graspilot_curve_params, normalized flag
            norm = numel(key) > 4 && strcmp(key(end - 3:end), '_SSS');
            if norm, key = key(1:end - 4); end
            switch key
                case 'PE',  name = 'PeakEnh'; fld = 'PE';
                case 'RE',  name = 'RelEnh';  fld = 'RE';
                case 'WIS', name = 'WashIn';  fld = 'WashIn';
                case 'WOS', name = 'WashOut'; fld = 'WashOut';
                case 'AUC', name = 'AUC';     fld = 'AUC';
                case 'TTP', name = 'TTP';     fld = 'TTP';
                case 'BAT', name = 'Onset';   fld = 'Onset';
                case 'T90', name = 'T90';     fld = 'T90';
                case 'MITR', name = 'MITR';   fld = 'MITR';
                case 'REE', name = 'REearly'; fld = 'REearly';
                case 'SER', name = 'SER';     fld = 'SER';
                case 'WOP', name = 'WashoutPct'; fld = 'WashoutPct';
                case 'PLAT', name = 'Plateau'; fld = 'Plateau';
                case 'CT',  name = 'CurveType'; fld = 'CurveType';
                otherwise
                    error('GRASPilot:param', 'Unknown parameter %s.', key);
            end
        end

        function d = describeVals(v)
            v = sort(double(v(isfinite(v))));
            n = numel(v);
            d = struct('n', n, 'mean', NaN, 'sd', NaN, 'median', NaN, 'q1', NaN, 'q3', NaN, ...
                'p10', NaN, 'p90', NaN, 'min', NaN, 'max', NaN);
            if n == 0, return; end
            v = v(:);
            d.mean = mean(v);
            if n > 1, d.sd = std(v); else, d.sd = 0; end
            q = @(p) GRASPilot.quant(v, p);
            d.median = q(0.5); d.q1 = q(0.25); d.q3 = q(0.75); d.p10 = q(0.1); d.p90 = q(0.9);
            d.min = v(1); d.max = v(end);
        end

        function cm = colormapByName(name, n)
            x = linspace(0, 1, n)';
            switch lower(name)
                case 'turbo'      % polynomial approximation of Google's Turbo
                    r = 0.13572138 + x .* (4.61539260 + x .* (-42.66032258 + x .* (132.13108234 + x .* (-152.94239396 + x * 59.28637943))));
                    g = 0.09140261 + x .* (2.19418839 + x .* (4.84296658 + x .* (-14.18503333 + x .* (4.27729857 + x * 2.82956604))));
                    b = 0.10667330 + x .* (12.64194608 + x .* (-60.58204836 + x .* (110.36276771 + x .* (-89.90310912 + x * 27.34824973))));
                case 'hot'
                    r = 3 * x; g = 3 * x - 1; b = 3 * x - 2;
                case 'gray'
                    r = x; g = x; b = x;
                otherwise         % jet
                    r = min(4 * x - 1.5, -4 * x + 4.5);
                    g = min(4 * x - 0.5, -4 * x + 3.5);
                    b = min(4 * x + 0.5, -4 * x + 2.5);
            end
            cm = min(1, max(0, [r g b]));
        end

        function h = handleArray(n)
            if exist('gobjects', 'builtin') == 5 || exist('gobjects', 'file') == 2
                h = gobjects(1, n);
            else
                h = zeros(1, n);
            end
        end

        function S = sliceOfArray(A, a, q, n)
            switch a
                case 1, S = reshape(A(q, :, :), n(2), n(3));
                case 2, S = reshape(A(:, q, :), n(1), n(3));
                case 3, S = A(:, :, q);
            end
        end

        function e = maskEdge(m)
            [h, w] = size(m);
            up = [false(1, w); m(1:end - 1, :)];
            dn = [m(2:end, :); false(1, w)];
            lf = [false(h, 1), m(:, 1:end - 1)];
            rt = [m(:, 2:end), false(h, 1)];
            e = m & ~(up & dn & lf & rt);
        end

        function R = floodFill3D(M, s0)
            % 6-connected region of true voxels in M containing s0; grows
            % inside a box that expands by one voxel per iteration.
            sz = [size(M, 1) size(M, 2) size(M, 3)];
            R = false(sz);
            R(s0(1), s0(2), s0(3)) = true;
            b0 = s0(:)'; b1 = s0(:)';
            cnt = 1;
            while true
                a0 = max(1, b0 - 1); a1 = min(sz, b1 + 1);
                i1 = a0(1):a1(1); i2 = a0(2):a1(2); i3 = a0(3):a1(3);
                Rs = R(i1, i2, i3);
                Rn = Rs;
                Rn(2:end, :, :) = Rn(2:end, :, :) | Rs(1:end - 1, :, :);
                Rn(1:end - 1, :, :) = Rn(1:end - 1, :, :) | Rs(2:end, :, :);
                Rn(:, 2:end, :) = Rn(:, 2:end, :) | Rs(:, 1:end - 1, :);
                Rn(:, 1:end - 1, :) = Rn(:, 1:end - 1, :) | Rs(:, 2:end, :);
                Rn(:, :, 2:end) = Rn(:, :, 2:end) | Rs(:, :, 1:end - 1);
                Rn(:, :, 1:end - 1) = Rn(:, :, 1:end - 1) | Rs(:, :, 2:end);
                Rn = Rn & M(i1, i2, i3);
                c = nnz(Rn);
                if c == cnt, break; end
                cnt = c;
                R(i1, i2, i3) = Rn;
                b0 = a0; b1 = a1;
            end
        end

        function F = fillHoles2D(M)
            bg = ~M;
            [h, w] = size(M);
            R = false(h, w);
            R([1 h], :) = bg([1 h], :);
            R(:, [1 w]) = R(:, [1 w]) | bg(:, [1 w]);
            while true
                Rn = R;
                Rn(2:end, :) = Rn(2:end, :) | R(1:end - 1, :);
                Rn(1:end - 1, :) = Rn(1:end - 1, :) | R(2:end, :);
                Rn(:, 2:end) = Rn(:, 2:end) | R(:, 1:end - 1);
                Rn(:, 1:end - 1) = Rn(:, 1:end - 1) | R(:, 2:end);
                Rn = Rn & bg;
                if isequal(Rn, R), break; end
                R = Rn;
            end
            F = ~R;
        end

        function th = otsu(vals)
            vals = vals(isfinite(vals));
            th = NaN;
            if numel(vals) < 2, return; end
            lo = min(vals); hi = max(vals);
            if hi <= lo, th = lo; return; end
            nb = 256;
            idx = min(nb, floor((vals - lo) / (hi - lo) * nb) + 1);
            counts = accumarray(idx(:), 1, [nb 1]);
            p = counts / sum(counts);
            centers = lo + ((1:nb)' - 0.5) * (hi - lo) / nb;
            om = cumsum(p);
            mu = cumsum(p .* centers);
            muT = mu(end);
            sb = (muT * om - mu) .^ 2 ./ (om .* (1 - om));
            sb(~isfinite(sb)) = 0;
            [~, k] = max(sb);
            th = lo + k * (hi - lo) / nb;
        end

        function q = quant(Ts, p)
            n = size(Ts, 1);
            pos = 1 + (n - 1) * p;
            a = floor(pos); b = ceil(pos); w = pos - a;
            q = Ts(a, :) * (1 - w) + Ts(b, :) * w;
        end

        function S = reindexStack(S, r)
            out = {};
            for i = 1:numel(S)
                a = S{i};
                ch = a.changes;
                if ~isempty(ch)
                    ch = ch([ch.roi] ~= r);
                    for j = 1:numel(ch)
                        if ch(j).roi > r, ch(j).roi = ch(j).roi - 1; end
                    end
                end
                if ~isempty(ch)
                    a.changes = ch;
                    out{end + 1} = a; %#ok<AGROW>
                end
            end
            S = out;
        end

        function setSlider(h, lo, hi, val)
            cur = get(h, 'Value');
            set(h, 'Min', min([lo cur]), 'Max', max([hi cur]));
            set(h, 'Value', min(max(val, lo), hi));
            set(h, 'Min', lo, 'Max', hi);
        end

        function s = num2s(x)
            if isinf(x)
                if x > 0, s = 'Inf'; else, s = '-Inf'; end
            else
                s = sprintf('%.4g', x);
            end
        end

        function f = fixNii(f)
            if isempty(regexpi(f, '\.nii(\.gz)?$', 'once'))
                f = [f '.nii.gz'];
            end
        end

        function [names, cols] = readLabelTable(base)
            names = {}; cols = zeros(0, 3);
            f = [base '_labels.txt'];
            if exist(f, 'file') ~= 2, return; end
            fid = fopen(f, 'r');
            while true
                l = fgetl(fid);
                if ~ischar(l), break; end
                if isempty(l) || l(1) == '#', continue; end
                p = regexp(l, '\t', 'split');
                if numel(p) < 5, continue; end
                k = str2double(p{1});
                if ~isfinite(k) || k < 1, continue; end
                names{k} = p{2}; %#ok<AGROW>
                cols(k, :) = [str2double(p{3}) str2double(p{4}) str2double(p{5})] / 255; %#ok<AGROW>
            end
            fclose(fid);
        end
    end
end
