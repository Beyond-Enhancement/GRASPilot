classdef GraspDCEViewer < handle
%GRASPDCEVIEWER 3-plane viewer and ROI analysis platform for dynamic GRASP DCE-MRI.
%
%   GraspDCEViewer                       empty viewer (File > Open .mat)
%   GraspDCEViewer('recon.mat')          open a reconstruction (load dialog)
%   GraspDCEViewer(vk_all, basis_all)    subspace variables already in the workspace
%   GraspDCEViewer(data1)                4D array from ReadData.m ([nx ny t nz])
%   GraspDCEViewer(D)                    GraspData object
%   app = GraspDCEViewer(...)            handle for scripting
%
%   Name/value options: 'Plane' ('axial'|'coronal'|'sagittal'), 'Voxel'
%   [row col slice] mm, 'dt' s/frame, 't0' s, 'DimOrder' (dense input,
%   default 'xytz'), 'CropFactor', 'Normalize', 'CoefVar', 'BasisVar',
%   'ImageVar', 'NoDialog' (true skips the load dialog).
%
%   Mouse:  left = active tool, right-drag = window/level, shift- or
%           middle-drag = pan, wheel = slice, Ctrl+wheel = zoom.
%   Keys:   Left/Right frame -/+1 (Shift: 10), Up/Down slice, Home/End,
%           Space play, N/C/B/E tools, [ ] brush size, D toggle 3D,
%           P preview, Ctrl+Z undo, Ctrl+Y or Ctrl+Shift+Z redo.
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
        ThrSrc = 1;
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
        OvlRGB = cell(1, 3);
        OvlA = cell(1, 3);
        PrevA = cell(1, 3);
        LoadOpts = struct();
        Busy = false;
        Quiet = false;          % true: dialogs auto-answer with defaults (scripting/tests)
    end

    properties (Constant)
        Version = '1.0';
        ViewAxis = [3 2 1];
        ViewNames = {'Axial', 'Coronal', 'Sagittal'};
        ViewLabels = {{'R', 'L', 'A', 'P'}, {'R', 'L', 'S', 'I'}, {'A', 'P', 'S', 'I'}};
        ThrNames = {'Signal at current frame', 'Current frame - baseline', ...
            'Current frame rel. enhancement (%)', 'Peak - baseline (max over time)', ...
            'Peak rel. enhancement (%)'};
        Palette = [1.00 0.85 0.10; 0.20 0.90 1.00; 1.00 0.35 0.95; 0.40 1.00 0.35; ...
            1.00 0.60 0.15; 0.55 0.55 1.00; 0.60 1.00 0.85; 1.00 0.55 0.55; 0.85 0.65 1.00];
        SSSColor = [1.00 0.22 0.22];
        Col = struct('bg', [0.11 0.11 0.13], 'panel', [0.16 0.16 0.18], 'fg', [0.90 0.90 0.92], ...
            'accent', [0.40 0.72 1.00], 'dim', [0.62 0.62 0.66], 'cross', [0.25 0.75 1.00], ...
            'sel', [0.30 0.55 0.85], 'prev', [0.20 1.00 0.30]);
    end

    % =====================================================================
    methods

        function obj = GraspDCEViewer(varargin)
            obj.LoadOpts = GraspDCEViewer.defaultLoadOpts();
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
            if isa(a{1}, 'GraspData')
                obj.setData(a{1});
                return
            end
            if ischar(a{1})
                opts = obj.mergeOpts(a(2:end));
                obj.openFile(a{1}, opts);
                return
            end
            if ~isnumeric(a{1})
                error('GraspDCEViewer:input', 'Unrecognized input.');
            end
            sub = numel(a) >= 2 && isnumeric(a{2});
            if sub, nv = a(3:end); else, nv = a(2:end); end
            opts = obj.mergeOpts(nv);
            if ~opts.NoDialog
                g = gdce_load_dialog('', opts);
                if isempty(g), return; end
                opts = g;
            end
            wb = obj.wbOpen('Preparing data...');
            try
                prog = @(f, m) obj.wbUpdate(wb, f, m);
                if sub
                    D = GraspData.fromSubspace(a{1}, a{2}, 'CropFactor', opts.CropFactor, ...
                        'Normalize', opts.Normalize, 'Progress', prog);
                    D.SourceVars = 'vk_all + basis_all (workspace)';
                else
                    D = GraspData.fromDense(a{1}, opts.DimOrder, 'CropFactor', opts.CropFactor, ...
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
                o = gdce_load_dialog(file, opts);
                if isempty(o), return; end
                opts = o;
            end
            wb = obj.wbOpen('Loading...');
            try
                D = GraspData.fromFile(file, 'Format', opts.Format, 'CoefVar', opts.CoefVar, ...
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
            D.setOrientation(opts.Plane, opts.SwapInPlane, opts.Flip);
            D.VoxNative = opts.Voxel(:)';
            D.dt = opts.dt;
            D.t0 = opts.t0;
            obj.LoadOpts.Plane = opts.Plane;
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
            obj.setupTimeSlider();
            obj.resetViewsGeometry();
            obj.autoContrast(false);
            obj.updateDataInfo();
            obj.refreshROIList();
            obj.refreshAll();
            obj.updateThrSliderRange();
            obj.Modified = keep;
            [~, nm] = fileparts(D.SourceFile);
            set(obj.Fig, 'Name', ['GRASP DCE Viewer - ' nm]);
            obj.setStatus(sprintf('Loaded %s. %s', nm, D.describe()));
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
            if obj.ThrSrc <= 3, obj.updateThrInfo(); end
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
            obj.UndoStack = GraspDCEViewer.reindexStack(obj.UndoStack, r);
            obj.RedoStack = GraspDCEViewer.reindexStack(obj.RedoStack, r);
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
            obj.ThrSrc = k;
            set(obj.H.thrSrc, 'Value', k);
            if k >= 2 && (obj.Preview || obj.Constrain)
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
            set(obj.H.thrLo, 'String', GraspDCEViewer.num2s(lo));
            set(obj.H.thrHi, 'String', GraspDCEViewer.num2s(hi));
            if isfinite(lo)
                mn = get(obj.H.thrSlider, 'Min'); mx = get(obj.H.thrSlider, 'Max');
                set(obj.H.thrSlider, 'Value', min(mx, max(mn, lo)));
            end
            obj.refreshPreviews();
            obj.updateThrInfo();
        end

        function setPreview(obj, on)
            on = logical(on);
            if on && obj.ThrSrc >= 2 && ~obj.ensureMaps(false), on = false; end
            obj.Preview = on;
            set(obj.H.chkPrev, 'Value', double(on));
            obj.refreshPreviews();
        end

        function setConstrain(obj, on)
            on = logical(on);
            if on && obj.ThrSrc >= 2 && ~obj.ensureMaps(false), on = false; end
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
            th = GraspDCEViewer.otsu(double(vals));
            if ~isfinite(th), return; end
            obj.setThreshold(th, Inf);
            obj.setStatus(sprintf('Otsu threshold %.4g from %s (%s). Use Keep: 3D to apply.', th, src, ...
                obj.ThrNames{obj.ThrSrc}));
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
                S = GraspDCEViewer.sliceOfArray(M, a, q, n);
                if ~any(S(:)), continue; end
                rows = find(any(S, 2)); cols = find(any(S, 1));
                r1 = max(1, rows(1) - 1); r2 = min(size(S, 1), rows(end) + 1);
                c1 = max(1, cols(1) - 1); c2 = min(size(S, 2), cols(end) + 1);
                sub = S(r1:r2, c1:c2);
                nw = GraspDCEViewer.fillHoles2D(sub) & ~sub;
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

        function ok = ensureMaps(obj, force)
            ok = true;
            if isempty(obj.Data), ok = false; return; end
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
            obj.updateThrSliderRange();
            obj.refreshPreviews();
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
            Vv = obj.Data.canonVolume(obj.Frame);
            step = max(1, floor(numel(Vv) / 200000));
            s = sort(double(Vv(1:step:end)));
            hi = s(max(1, round(0.995 * numel(s))));
            lo = 0;
            if s(1) < 0, lo = s(max(1, round(0.005 * numel(s)))); end
            if hi <= lo, hi = s(end); end
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
            W = 470; Hh = 330;
            ss = get(0, 'ScreenSize');
            d = figure('Name', 'Geometry and orientation', 'NumberTitle', 'off', 'MenuBar', 'none', ...
                'ToolBar', 'none', 'Resize', 'off', 'Units', 'pixels', 'Color', [0.94 0.94 0.94], ...
                'Position', [(ss(3) - W) / 2, (ss(4) - Hh) / 2, W, Hh], 'WindowStyle', 'modal');
            bg = [0.94 0.94 0.94];
            t = @(x, y, w, s) uicontrol(d, 'Style', 'text', 'Units', 'pixels', 'Position', [x y w 18], ...
                'String', s, 'HorizontalAlignment', 'left', 'BackgroundColor', bg);
            e = @(x, y, w, s) uicontrol(d, 'Style', 'edit', 'Units', 'pixels', 'Position', [x y w 22], ...
                'String', s, 'BackgroundColor', 'w');
            planes = {'axial', 'coronal', 'sagittal'};
            t(12, Hh - 34, 160, 'Native (slice) plane');
            hP = uicontrol(d, 'Style', 'popupmenu', 'Units', 'pixels', 'Position', [175 Hh - 32 120 22], ...
                'String', {'Axial', 'Coronal', 'Sagittal'}, 'Value', find(strcmp(planes, D.Plane)), ...
                'BackgroundColor', 'w');
            hS = uicontrol(d, 'Style', 'checkbox', 'Units', 'pixels', 'Position', [305 Hh - 32 160 20], ...
                'String', 'Swap rows/cols', 'Value', D.SwapInPlane, 'BackgroundColor', bg);
            hF = GraspDCEViewer.handleArray(3);
            lab = {'Flip L-R', 'Flip A-P', 'Flip S-I'};
            for i = 1:3
                hF(i) = uicontrol(d, 'Style', 'checkbox', 'Units', 'pixels', ...
                    'Position', [12 + (i - 1) * 110, Hh - 62, 105, 20], 'String', lab{i}, ...
                    'Value', D.UserFlip(i), 'BackgroundColor', bg);
            end
            t(12, Hh - 98, 160, 'Pixel spacing rows (mm)');
            hVr = e(175, Hh - 96, 60, num2str(D.VoxNative(1)));
            t(245, Hh - 98, 40, 'cols');
            hVc = e(285, Hh - 96, 60, num2str(D.VoxNative(2)));
            t(12, Hh - 126, 160, 'Slice spacing (mm)');
            hVs = e(175, Hh - 124, 60, num2str(D.VoxNative(3)));
            t(12, Hh - 162, 160, 'Temporal resolution (s)');
            hDt = e(175, Hh - 160, 60, num2str(D.dt));
            t(245, Hh - 162, 100, 'Frame 1 at (s)');
            hT0 = e(345, Hh - 160, 60, num2str(D.t0));
            uicontrol(d, 'Style', 'text', 'Units', 'pixels', 'Position', [12 60 W - 24 70], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', [1 1 0.92], 'String', ...
                ['Check against anatomy: the axial view should show anterior at the top and the ' ...
                'patient''s right on the left (radiological), coronal and sagittal superior at the top. ' ...
                'Use the flips to correct. ROIs follow the change. Exports use this geometry.']);
            uicontrol(d, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [W - 300 14 90 30], ...
                'String', 'Apply', 'Callback', @(s, ev) apply(false));
            uicontrol(d, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [W - 200 14 90 30], ...
                'String', 'OK', 'FontWeight', 'bold', 'Callback', @(s, ev) apply(true));
            uicontrol(d, 'Style', 'pushbutton', 'Units', 'pixels', 'Position', [W - 100 14 88 30], ...
                'String', 'Close', 'Callback', @(s, ev) delete(d));

            function apply(closeAfter)
                vox = [str2double(get(hVr, 'String')) str2double(get(hVc, 'String')) str2double(get(hVs, 'String'))];
                if any(~isfinite(vox)) || any(vox <= 0)
                    errordlg('Voxel sizes must be positive numbers.', 'Geometry', 'modal');
                    return
                end
                t0v = str2double(get(hT0, 'String'));
                if ~isfinite(t0v), t0v = 0; end
                fl = logical([get(hF(1), 'Value') get(hF(2), 'Value') get(hF(3), 'Value')]);
                obj.applyGeometry(planes{get(hP, 'Value')}, logical(get(hS, 'Value')), fl, vox, ...
                    str2double(get(hDt, 'String')), t0v);
                if closeAfter, delete(d); end
            end
        end

        function applyGeometry(obj, plane, swap, flips, vox, dt, t0)
            D = obj.Data;
            orientChanged = ~strcmp(D.Plane, lower(plane)) || D.SwapInPlane ~= logical(swap) || ...
                any(D.UserFlip ~= logical(flips(:)'));
            if orientChanged
                nat = cellfun(@(m) D.canon2native(m), obj.Masks, 'UniformOutput', false);
                D.setOrientation(plane, swap, flips);
                obj.Masks = cellfun(@(m) D.native2canon(m), nat, 'UniformOutput', false);
                obj.UndoStack = {}; obj.RedoStack = {};
                obj.Pos = max(1, round((D.canonSize() + 1) / 2));
            end
            D.VoxNative = vox(:)';
            D.dt = dt;
            D.t0 = t0;
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
        end

        % ------------------------------------------------------------ export / import
        function exportImageNifti(obj, file, frames)
            if nargin < 3 || isempty(frames), frames = 1:obj.Data.NT; end
            wb = obj.wbOpen('Writing NIfTI...');
            try
                gdce_write_dce_nifti(obj.Data, file, frames, @(f, m) obj.wbUpdate(wb, f, m));
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
                'Descrip', 'GRASP DCE Viewer ROI labels', 'CalRange', [0 nR]);
            gdce_nifti_write(file, L, hdr);
            base = regexprep(file, '(\.nii\.gz|\.nii)$', '', 'ignorecase');
            obj.writeLabelTables(base);
            if separate
                for r = 1:nR
                    if ~any(obj.Masks{r}(:)), continue; end
                    nm = regexprep(obj.ROIs(r).Name, '[^A-Za-z0-9_\-]+', '_');
                    h1 = hdr; h1.Datatype = 'uint8'; h1.CalRange = [0 1];
                    gdce_nifti_write(sprintf('%s_roi%02d_%s.nii.gz', base, r, nm), uint8(obj.Masks{r}), h1);
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
                files = gdce_export_timeseries(obj.Data, obj.Masks, names, file, ...
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
            names = {'Baseline', 'Max', 'PeakEnh', 'RelEnh', 'TTP', 'AUC'};
            hdr = struct('Dim', D.canonSize(), 'Vox', D.canonVox(), 'Affine', D.affine(), ...
                'Datatype', 'float32', 'Descrip', 'GRASP DCE Viewer map');
            for i = 1:numel(names)
                A = D.canonMap(names{i});
                if strcmp(names{i}, 'TTP') && D.hasTiming()
                    A = D.t0 + (A - 1) * D.dt;
                end
                gdce_nifti_write(fullfile(folder, sprintf('%s_%s.nii.gz', base, names{i})), single(A), hdr);
            end
            obj.setStatus(sprintf('Wrote %d maps to %s (baseline frames %d-%d).', numel(names), folder, ...
                obj.BaseFrames(1), obj.BaseFrames(2)));
        end

        function saveROIs(obj, file)
            D = obj.Data;
            S = struct();
            S.Format = 'GraspDCEViewer ROI session';
            S.Version = obj.Version;
            S.Saved = datestr(now, 'yyyy-mm-dd HH:MM:SS');
            S.SourceFile = D.SourceFile;
            S.NativeSize = D.NativeSize;
            S.Plane = D.Plane;
            S.SwapInPlane = D.SwapInPlane;
            S.UserFlip = D.UserFlip;
            S.VoxNative = D.VoxNative;
            S.dt = D.dt;
            S.t0 = D.t0;
            S.BaseFrames = obj.BaseFrames;
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
                error('GraspDCEViewer:session', 'Not a GRASP DCE Viewer ROI file.');
            end
            if ~isequal(S.NativeSize(:)', D.NativeSize)
                error('GraspDCEViewer:session', 'ROI file is for a %s matrix; data is %s.', ...
                    num2str(S.NativeSize(:)'), num2str(D.NativeSize));
            end
            mode = 'Replace';
            if ~isempty(obj.ROIs)
                mode = obj.ask('Replace the current ROIs or append?', 'Load ROIs', ...
                    {'Replace', 'Append', 'Cancel'}, 'Replace');
                if isempty(mode) || strcmp(mode, 'Cancel'), return; end
            end
            geoDiff = ~strcmp(S.Plane, D.Plane) || logical(S.SwapInPlane) ~= D.SwapInPlane || ...
                any(logical(S.UserFlip(:)') ~= D.UserFlip) || any(abs(S.VoxNative(:)' - D.VoxNative) > 1e-6) ...
                || ~isequaln(S.dt, D.dt);
            if geoDiff
                b = obj.ask(['Apply the geometry saved with the ROIs (orientation, voxel size, timing)? ' ...
                    'The ROI voxels are correct either way.'], 'Load ROIs', {'Apply', 'Keep current'}, 'Apply');
                if strcmp(b, 'Apply')
                    obj.applyGeometry(S.Plane, logical(S.SwapInPlane), logical(S.UserFlip(:)'), ...
                        S.VoxNative, S.dt, S.t0);
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
            if isfield(S, 'BaseFrames')
                obj.BaseFrames = S.BaseFrames;
                set(obj.H.base1, 'String', num2str(S.BaseFrames(1)));
                set(obj.H.base2, 'String', num2str(S.BaseFrames(2)));
            end
            obj.UndoStack = {}; obj.RedoStack = {};
            obj.afterMaskChange();
            obj.Modified = false;
            obj.setStatus(sprintf('Loaded %d ROI(s) from %s', numel(S.ROIIndexNative), file));
        end

        function importMaskNifti(obj, file)
            D = obj.Data;
            [img, h] = gdce_nifti_read(file);
            img = img(:, :, :, 1);
            [img, msg] = gdce_reorient_to_ref(img, h.Affine, D.affine(), D.canonSize());
            img = round(double(img));
            labels = unique(img(img ~= 0));
            if isempty(labels)
                obj.inform('The mask is empty.', 'Import mask');
                return
            end
            [names, cols] = GraspDCEViewer.readLabelTable(regexprep(file, '(\.nii\.gz|\.nii)$', '', 'ignorecase'));
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
                b = obj.ask('Save ROIs before closing?', 'GRASP DCE Viewer', {'Save', 'Discard', 'Cancel'}, 'Save');
                if isempty(b) || strcmp(b, 'Cancel'), return; end
                if strcmp(b, 'Save') && ~obj.onSaveROIs(), return; end
            end
            obj.stopPlay();
            try
                if ~isempty(obj.PlayTimer) && isvalid(obj.PlayTimer), delete(obj.PlayTimer); end
            catch
            end
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
            Hh = min(1020, max(760, ss(4) - 110));
            pos = [max(1, round((ss(3) - W) / 2)), max(30, round((ss(4) - Hh) / 2)), W, Hh];
            obj.Fig = figure('Name', 'GRASP DCE Viewer', 'NumberTitle', 'off', 'MenuBar', 'none', ...
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
                lb = GraspDCEViewer.handleArray(4);
                for k = 1:4
                    lb(k) = uicontrol(P, 'Style', 'text', 'Units', 'pixels', 'Position', [1 1 14 16], ...
                        'String', obj.ViewLabels{v}{k}, 'BackgroundColor', [0 0 0], ...
                        'ForegroundColor', [0.95 0.75 0.30], 'FontSize', 9, 'FontWeight', 'bold');
                end
                obj.V(v).Lbl = lb;
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
            m = uimenu(obj.Fig, 'Label', 'View');
            uimenu(m, 'Label', 'Geometry and orientation...', 'Callback', @(s, e) obj.geometryDialog());
            uimenu(m, 'Label', 'Auto contrast', 'Callback', @(s, e) obj.autoContrast(true));
            uimenu(m, 'Label', 'Reset zoom', 'Callback', @(s, e) obj.resetZoom());
            m = uimenu(obj.Fig, 'Label', 'Help');
            uimenu(m, 'Label', 'Mouse and keyboard...', 'Callback', @(s, e) obj.showHelp());
            uimenu(m, 'Label', 'About', 'Callback', @(s, e) obj.inform(sprintf(['GRASP DCE Viewer %s\n' ...
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
            obj.H.info = obj.ctl(P, 10, 58, -10, 48, 'text', 'No data loaded.', 'FontSize', 8, 'ForegroundColor', C.dim);

            hdr(112, 'SEGMENTATION');
            names = {'nav', 'circle', 'brush', 'erase'};
            labels = {'Navigate', 'Circle', 'Brush', 'Eraser'};
            xs = [10 81 152 223];
            for i = 1:4
                nm = names{i};
                obj.H.tool.(nm) = obj.ctl(P, xs(i), 132, 66, 26, 'togglebutton', labels{i}, ...
                    'Callback', @(s, e) obj.setTool(nm));
            end
            obj.ctl(P, 10, 167, 80, 18, 'text', 'Radius (mm)');
            obj.H.rSlider = obj.ctl(P, 92, 167, 140, 18, 'slider', '', 'Min', 0.5, 'Max', 25, ...
                'Value', obj.BrushR, 'SliderStep', [0.5 2] / 24.5, 'Callback', @(s, e) obj.setBrushR(get(s, 'Value')));
            obj.H.rEdit = obj.ctl(P, 238, 165, 52, 22, 'edit', num2str(obj.BrushR), ...
                'Callback', @(s, e) obj.setBrushR(str2double(get(s, 'String'))));
            obj.H.chk3D = obj.ctl(P, 10, 192, 135, 20, 'checkbox', '3D (sphere)', 'Value', obj.Use3D, ...
                'TooltipString', 'Circle, brush and eraser act in 3D (key D)', 'Callback', @(s, e) obj.set3D(get(s, 'Value')));
            obj.H.chkExcl = obj.ctl(P, 150, 192, 140, 20, 'checkbox', 'Exclusive ROIs', 'Value', obj.Exclusive, ...
                'TooltipString', 'Painting removes the voxels from other ROIs (no overlap)', ...
                'Callback', @(s, e) obj.setFlag('Exclusive', get(s, 'Value')));
            obj.ctl(P, 10, 216, 66, 24, 'pushbutton', 'Undo', 'Callback', @(s, e) obj.undo());
            obj.ctl(P, 81, 216, 66, 24, 'pushbutton', 'Redo', 'Callback', @(s, e) obj.redo());
            obj.ctl(P, 152, 216, 66, 24, 'pushbutton', 'Clear slice', ...
                'TooltipString', 'Clear the active ROI on the current slice of the last-used view', ...
                'Callback', @(s, e) obj.clearROI('slice'));
            obj.ctl(P, 223, 216, 66, 24, 'pushbutton', 'Clear ROI', 'Callback', @(s, e) obj.clearROI('all'));

            hdr(250, 'THRESHOLD');
            obj.ctl(P, 10, 271, 40, 18, 'text', 'Map');
            obj.H.thrSrc = obj.ctl(P, 50, 269, -10, 22, 'popupmenu', obj.ThrNames, 'Value', obj.ThrSrc, ...
                'Callback', @(s, e) obj.setThrSource(get(s, 'Value')));
            obj.ctl(P, 10, 299, 100, 18, 'text', 'Baseline frames');
            obj.H.base1 = obj.ctl(P, 112, 297, 40, 22, 'edit', '1', 'Callback', @(s, e) obj.onBaseline());
            obj.ctl(P, 155, 299, 20, 18, 'text', 'to');
            obj.H.base2 = obj.ctl(P, 176, 297, 40, 22, 'edit', '5', 'Callback', @(s, e) obj.onBaseline());
            obj.ctl(P, 222, 297, 68, 22, 'pushbutton', 'Maps', ...
                'TooltipString', 'Recompute baseline / peak enhancement maps', 'Callback', @(s, e) obj.ensureMaps(true));
            obj.ctl(P, 10, 327, 32, 18, 'text', 'Low');
            obj.H.thrLo = obj.ctl(P, 42, 325, 74, 22, 'edit', '-Inf', 'Callback', @(s, e) obj.onThrEdit());
            obj.ctl(P, 122, 327, 34, 18, 'text', 'High');
            obj.H.thrHi = obj.ctl(P, 156, 325, 74, 22, 'edit', 'Inf', 'Callback', @(s, e) obj.onThrEdit());
            obj.ctl(P, 236, 325, 54, 22, 'pushbutton', 'Otsu', ...
                'TooltipString', 'Suggest Low from the active ROI (Otsu)', 'Callback', @(s, e) obj.suggestThreshold());
            obj.H.thrSlider = obj.ctl(P, 10, 353, -10, 18, 'slider', '', 'Min', 0, 'Max', 1, 'Value', 0, ...
                'Callback', @(s, e) obj.onThrSlider());
            try
                addlistener(obj.H.thrSlider, 'ContinuousValueChange', @(s, e) obj.onThrSlider());
            catch
            end
            obj.H.chkPrev = obj.ctl(P, 10, 377, 90, 20, 'checkbox', 'Preview', 'Value', 0, ...
                'TooltipString', 'Show voxels inside the threshold window (green) - key P', ...
                'Callback', @(s, e) obj.setPreview(get(s, 'Value')));
            obj.H.chkCons = obj.ctl(P, 105, 377, 185, 20, 'checkbox', 'Constrain brush/circle', 'Value', 0, ...
                'TooltipString', 'Brush and circle add only voxels inside the threshold window', ...
                'Callback', @(s, e) obj.setConstrain(get(s, 'Value')));
            obj.ctl(P, 10, 401, 92, 24, 'pushbutton', 'Keep: slice', ...
                'TooltipString', 'Remove voxels outside the window from the active ROI on this slice', ...
                'Callback', @(s, e) obj.applyThreshold('slice'));
            obj.ctl(P, 106, 401, 92, 24, 'pushbutton', 'Keep: 3D', ...
                'TooltipString', 'Remove voxels outside the window from the whole active ROI', ...
                'Callback', @(s, e) obj.applyThreshold('3d'));
            obj.ctl(P, 202, 401, 88, 24, 'pushbutton', 'Fill holes', ...
                'TooltipString', 'Fill enclosed holes slice by slice (plane of the last-used view)', ...
                'Callback', @(s, e) obj.fillHoles());
            obj.H.thrInfo = obj.ctl(P, 10, 429, -10, 18, 'text', '', 'FontSize', 8, 'ForegroundColor', C.dim);

            hdr(455, 'DISPLAY');
            obj.ctl(P, 10, 477, 50, 18, 'text', 'Window');
            obj.H.wMin = obj.ctl(P, 62, 475, 60, 22, 'edit', '0', 'Callback', @(s, e) obj.onWindowEdit());
            obj.H.wMax = obj.ctl(P, 126, 475, 60, 22, 'edit', '1', 'Callback', @(s, e) obj.onWindowEdit());
            obj.ctl(P, 190, 475, 48, 22, 'pushbutton', 'Auto', 'Callback', @(s, e) obj.autoContrast(true));
            obj.ctl(P, 242, 475, 48, 22, 'pushbutton', 'Fit', 'TooltipString', 'Reset zoom', ...
                'Callback', @(s, e) obj.resetZoom());
            obj.ctl(P, 10, 503, 100, 18, 'text', 'ROI opacity');
            obj.H.alpha = obj.ctl(P, 112, 503, -10, 18, 'slider', '', 'Min', 0, 'Max', 1, 'Value', obj.Alpha, ...
                'Callback', @(s, e) obj.setFlag('Alpha', get(s, 'Value')));
            obj.H.chkOutline = obj.ctl(P, 10, 527, 80, 20, 'checkbox', 'Outline', 'Value', obj.ShowOutline, ...
                'Callback', @(s, e) obj.setFlag('ShowOutline', get(s, 'Value')));
            obj.H.chkCross = obj.ctl(P, 95, 527, 90, 20, 'checkbox', 'Crosshair', 'Value', obj.ShowCross, ...
                'Callback', @(s, e) obj.setFlag('ShowCross', get(s, 'Value')));
            obj.H.chkCursor = obj.ctl(P, 190, 527, 100, 20, 'checkbox', 'Cursor TIC', 'Value', obj.ShowCursorTIC, ...
                'TooltipString', 'Plot the time course of the crosshair voxel', ...
                'Callback', @(s, e) obj.setFlag('ShowCursorTIC', get(s, 'Value')));

            hdr(557, 'EXPORT');
            obj.ctl(P, 10, 577, 137, 26, 'pushbutton', 'Image -> NIfTI', 'Callback', @(s, e) obj.onExportImage());
            obj.ctl(P, 152, 577, 138, 26, 'pushbutton', 'Mask -> NIfTI', 'Callback', @(s, e) obj.onExportMask());
            obj.ctl(P, 10, 607, -10, 26, 'pushbutton', 'ROI voxel time series -> spreadsheet', ...
                'Callback', @(s, e) obj.onExportTS());
            obj.ctl(P, 10, 637, 90, 26, 'pushbutton', 'Save ROIs', 'Callback', @(s, e) obj.onSaveROIs());
            obj.ctl(P, 104, 637, 90, 26, 'pushbutton', 'Load ROIs', 'Callback', @(s, e) obj.onLoadROIs());
            obj.ctl(P, 198, 637, 92, 26, 'pushbutton', 'Import mask', 'Callback', @(s, e) obj.onImportMask());
        end

        function buildTime(obj)
            P = obj.H.time;
            obj.ctl(P, 6, 7, 30, 22, 'pushbutton', '|<', 'TooltipString', 'First frame (Home)', ...
                'Callback', @(s, e) obj.setFrame(1));
            obj.ctl(P, 38, 7, 30, 22, 'pushbutton', '<', 'TooltipString', 'Previous frame (Left)', ...
                'Callback', @(s, e) obj.setFrame(obj.Frame - 1));
            obj.H.slider = obj.ctl(P, 72, 9, -466, 18, 'slider', '', 'Min', 0, 'Max', 1, 'Value', 0, ...
                'Callback', @(s, e) obj.onSlider());
            try
                addlistener(obj.H.slider, 'ContinuousValueChange', @(s, e) obj.onSlider());
            catch
            end
            obj.ctl(P, -460, 7, 30, 22, 'pushbutton', '>', 'TooltipString', 'Next frame (Right)', ...
                'Callback', @(s, e) obj.setFrame(obj.Frame + 1));
            obj.ctl(P, -428, 7, 30, 22, 'pushbutton', '>|', 'TooltipString', 'Last frame (End)', ...
                'Callback', @(s, e) obj.setFrame(Inf));
            obj.ctl(P, -394, 9, 48, 18, 'text', 'Frame', 'HorizontalAlignment', 'right');
            obj.H.frameEdit = obj.ctl(P, -342, 7, 46, 22, 'edit', '1', ...
                'Callback', @(s, e) obj.setFrame(str2double(get(s, 'String')), true));
            obj.H.frameTotal = obj.ctl(P, -292, 9, 44, 18, 'text', '/ 0');
            obj.H.timeText = obj.ctl(P, -246, 9, 120, 18, 'text', '', 'ForegroundColor', obj.Col.accent, ...
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
            obj.ctl(P, 50, 5, 92, 22, 'pushbutton', '+ SSS control', ...
                'TooltipString', 'Add a superior sagittal sinus control ROI', 'Callback', @(s, e) obj.addSSS());
            obj.ctl(P, 146, 5, 80, 22, 'pushbutton', '+ Tumor', 'Callback', @(s, e) obj.addTumor());
            obj.ctl(P, 230, 5, 122, 22, 'pushbutton', '+ Custom...', 'Callback', @(s, e) obj.addCustom());
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
            img = obj.toScreen(v, GraspDCEViewer.sliceOfArray(A, a, obj.Pos(a), obj.Data.canonSize()));
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
                obj.ViewImg{v} = obj.toScreen(v, obj.Data.canonSlice(a, obj.Pos(a), obj.Frame));
                obj.computePreview(v);
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
                set(obj.V(v).Title, 'String', sprintf('%s   %d / %d', obj.ViewNames{v}, obj.Pos(a), n(a)));
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
                if obj.ShowOutline
                    e = GraspDCEViewer.maskEdge(m);
                    A(m) = obj.Alpha;
                    A(e) = 1;
                else
                    A(m) = obj.Alpha;
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

        function composeView(obj, v)
            % grayscale (window/level) + threshold preview + ROI overlay -> one RGB image
            S = obj.ViewImg{v};
            if isempty(S) || isempty(obj.V) || ~ishghandle(obj.V(v).Img), return; end
            c = obj.CLim;
            g = (double(S) - c(1)) / (c(2) - c(1));
            g(~isfinite(g)) = 0;
            g = min(1, max(0, g));
            rgb = cat(3, g, g, g);
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
        function M = thrSliceMap(obj, v)
            M = [];
            if isempty(obj.Data) || isempty(obj.ViewImg{v}), return; end
            S = double(obj.ViewImg{v});
            if obj.ThrSrc == 1
                M = S;
                return
            end
            bf = obj.BaseFrames(1):obj.BaseFrames(2);
            if ~obj.Data.hasMaps(bf), return; end
            switch obj.ThrSrc
                case 2
                    M = S - double(obj.screenSliceOf(v, obj.Data.canonMap('Baseline')));
                case 3
                    B = double(obj.screenSliceOf(v, obj.Data.canonMap('Baseline')));
                    M = 100 * (S - B) ./ B;
                    M(~(B > 0)) = NaN;
                case 4
                    M = double(obj.screenSliceOf(v, obj.Data.canonMap('PeakEnh')));
                case 5
                    M = double(obj.screenSliceOf(v, obj.Data.canonMap('RelEnh')));
            end
        end

        function vals = thrValues(obj, lin)
            lin = lin(:);
            switch obj.ThrSrc
                case 1
                    vals = double(obj.Data.valuesAtFrame(lin, obj.Frame));
                case 2
                    B = obj.Data.canonMap('Baseline');
                    vals = double(obj.Data.valuesAtFrame(lin, obj.Frame)) - double(B(lin));
                case 3
                    B = double(obj.Data.canonMap('Baseline'));
                    b = B(lin);
                    vals = 100 * (double(obj.Data.valuesAtFrame(lin, obj.Frame)) - b) ./ b;
                    vals(~(b > 0)) = NaN;
                case 4
                    M = obj.Data.canonMap('PeakEnh');
                    vals = double(M(lin));
                case 5
                    M = obj.Data.canonMap('RelEnh');
                    vals = double(M(lin));
            end
        end

        function tf = inWindow(obj, vals)
            tf = vals >= obj.ThrLo & vals <= obj.ThrHi;
        end

        function ok = ensureThresholdReady(obj)
            ok = true;
            if obj.ThrSrc >= 2, ok = obj.ensureMaps(false); end
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
            GraspDCEViewer.setSlider(obj.H.thrSlider, lo, hi, val);
        end

        function updateThrInfo(obj)
            if ~isfield(obj.H, 'thrInfo'), return; end
            r = obj.ActiveROI;
            s = '';
            if ~isempty(obj.Data) && r > 0 && r <= numel(obj.Masks) && (isfinite(obj.ThrLo) || isfinite(obj.ThrHi))
                lin = find(obj.Masks{r});
                if ~isempty(lin) && (obj.ThrSrc == 1 || obj.Data.hasMaps(obj.BaseFrames(1):obj.BaseFrames(2)))
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
                st.median = GraspDCEViewer.quant(Ts, 0.5);
                st.q1 = GraspDCEViewer.quant(Ts, 0.25);
                st.q3 = GraspDCEViewer.quant(Ts, 0.75);
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
            items = GraspDCEViewer.handleArray(0);
            hl = GraspDCEViewer.handleArray(0);
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
            x = cp(1, 1);
            if obj.TicX == 1 && obj.Data.hasTiming()
                f = round((x - obj.Data.t0) / obj.Data.dt) + 1;
            else
                f = round(x);
            end
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
            obj.setStatus(sprintf('%s  [i j k] = [%d %d %d]   S = %.4g   frame %d%s%s', obj.ViewNames{v}, ...
                c(1), c(2), c(3), val, obj.Frame, tl, rs));
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
                case 'leftbracket', obj.setBrushR(obj.BrushR / 1.25);
                case 'rightbracket', obj.setBrushR(obj.BrushR * 1.25);
            end
        end

        % ------------------------------------------------------------ control callbacks
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
            if obj.ThrSrc >= 2 && (obj.Preview || obj.Constrain)
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
            obj.exportImageNifti(GraspDCEViewer.fixNii(file), frames);
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
                obj.exportMaskNifti(GraspDCEViewer.fixNii(file), strcmp(b, 'Label map + per-ROI masks'));
            catch err
                obj.showError(err, 'Mask export failed');
            end
        end

        function onExportTS(obj)
            if ~obj.requireROIs(), return; end
            file = obj.putFile({'*.xlsx', 'Excel workbook (*.xlsx)'; '*.csv', 'CSV files (*.csv)'}, ...
                'Save ROI voxel time series', '_roi_timeseries.xlsx');
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
                '  Ctrl+Z / Ctrl+Y  undo / redo (Ctrl+Shift+Z also redoes)'
                ''
                'ENHANCING TUMOR WORKFLOW'
                '  1. + Tumor, then draw a generous circle/sphere around the lesion'
                '  2. Map = rel. enhancement (current or peak), baseline = pre-contrast frames'
                '  3. Otsu (or type Low), check Preview, then Keep: 3D'
                '  Alternatively tick "Constrain brush/circle" and paint.'};
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
                GraspDCEViewer.setSlider(h, 1, T, obj.Frame);
                set(h, 'SliderStep', [1 / (T - 1), min(1, 10 / (T - 1))], 'Enable', 'on');
            else
                GraspDCEViewer.setSlider(h, 0, 1, 0);
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
            else
                set(obj.H.timeText, 'String', '(dt not set)');
            end
        end

        function updateDataInfo(obj)
            D = obj.Data;
            [~, nm, ext] = fileparts(D.SourceFile);
            v = D.VoxNative;
            if D.hasTiming(), dts = sprintf('%.3g s/frame', D.dt); else, dts = 'dt not set'; end
            s = sprintf('%s%s\n%s\n%s plane, %.3gx%.3gx%.3g mm, %s', nm, ext, D.describe(), D.Plane, ...
                v(1), v(2), v(3), dts);
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
            try
                set(obj.H.roiTable, 'BackgroundColor', bgc, 'RowStriping', 'on', 'ForegroundColor', [0 0 0]);
            catch
            end
        end

        function afterMaskChange(obj)
            obj.Modified = true;
            obj.refreshROIList();
            obj.refreshOverlays();
            obj.refreshTIC();
            obj.updateThrInfo();
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
            fprintf(fid, '# GRASP DCE Viewer label table\n# label\tname\tR\tG\tB\tvoxels\tvolume_mL\n');
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
            if ~ok, obj.inform('Open a reconstruction first.', 'GRASP DCE Viewer'); end
        end

        function ok = requireROIs(obj)
            ok = obj.requireData();
            if ok
                ok = ~isempty(obj.Masks) && any(cellfun(@(m) any(m(:)), obj.Masks));
                if ~ok, obj.inform('There are no ROI voxels yet.', 'GRASP DCE Viewer'); end
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
            errordlg(sprintf('%s:\n\n%s', ttl, err.message), 'GRASP DCE Viewer', 'modal');
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
                wb = waitbar(0, msg, 'Name', 'GRASP DCE Viewer');
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
                error('GraspDCEViewer:opts', 'Options must be name/value pairs.');
            end
            for i = 1:2:numel(nv)
                m = find(strcmpi(fn, nv{i}), 1);
                if isempty(m)
                    error('GraspDCEViewer:opts', 'Unknown option "%s".', nv{i});
                end
                o.(fn{m}) = nv{i + 1};
            end
        end

        function loadPrefs(obj)
            try
                if ispref('GraspDCEViewer', 'Geometry')
                    g = getpref('GraspDCEViewer', 'Geometry');
                    f = {'Plane', 'Voxel', 'dt', 't0'};
                    for i = 1:numel(f)
                        if isfield(g, f{i}), obj.LoadOpts.(f{i}) = g.(f{i}); end
                    end
                end
                if ispref('GraspDCEViewer', 'LastDir')
                    obj.LastDir = getpref('GraspDCEViewer', 'LastDir');
                end
            catch
            end
        end

        function savePrefs(obj)
            if obj.Quiet, return; end
            try
                g = struct('Plane', obj.LoadOpts.Plane, 'Voxel', obj.LoadOpts.Voxel, ...
                    'dt', obj.LoadOpts.dt, 't0', obj.LoadOpts.t0);
                setpref('GraspDCEViewer', 'Geometry', g);
                setpref('GraspDCEViewer', 'LastDir', obj.LastDir);
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
                't0', 0, 'NoDialog', false);
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
