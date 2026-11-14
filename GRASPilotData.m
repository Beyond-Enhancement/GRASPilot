classdef GRASPilotData < handle
%GRASPILOTDATA Data model for dynamic GRASP DCE-MRI image series.
%
%   D = GRASPilotData.fromFile(matFile, 'Name', Value, ...)
%   D = GRASPilotData.fromSubspace(vk_all, basis_all, 'Name', Value, ...)
%   D = GRASPilotData.fromDense(img, dimOrder, 'Name', Value, ...)
%
%   SUBSPACE MODE reproduces ReadData.m frame by frame:
%       data1(:,:,t,z) = abs( Bz * Phi_z' )        (TempPCASub adjoint)
%   with Bz    = reshape(vk_all(z,:,:,:), nx*ny, K)
%        Phi_z = basis_all(z,:,:) as a [T x K] matrix (single),
%   followed by the CropImg central crop (default 2 = remove 2x readout
%   oversampling) and division by the global maximum over all slices and
%   frames (taken before the crop, as in ReadData.m). Only the cropped
%   coefficients are kept and frames are synthesized on demand, so memory
%   is ~1/4 of vk_all and no 4D frame array is ever built.
%
%   DENSE MODE holds a magnitude series as single [nx ny nz T].
%
%   MOTION CORRECTION (applyMotionCorrection) realigns every raw frame to
%   the reference of a rigid motion estimate (graspilot_motion_estimate)
%   and keeps the realigned frames as a dense series; the original frame
%   source is kept for revertMotionCorrection.
%
%   TEMPORAL BINNING (setBinning) averages the magnitude of n consecutive
%   reconstructed frames into one frame. Every frame-based method (slices,
%   volumes, time courses, maps, frame times, exports) then works on the
%   binned series; NT is the binned frame count and NTRaw the original.
%
%   Index conventions
%     native    : [rows cols partitions] = data1(:,:,t,z) of ReadData.m
%     canonical : array used by the viewer and every export. Axis 1 runs
%                 patient Right->Left, axis 2 Anterior->Posterior, axis 3
%                 Inferior->Superior (LPS+ voxel order).
%
%   Orientation is given by 3-letter axis codes (setAxisCodes): the
%   direction in which native rows, cols and slices increase, e.g.
%     'LPS'  rows R->L, cols A->P, slices I->S (scanner x,y,z order; the
%            deep-learning GRASP reconstructions are stored this way)
%     'PLS'  rows A->P, cols R->L: data1(:,:,t,z) already radiological axial
%     'ILP'  coronal slices,  'IPL'  sagittal slices
%
%   Options (fromFile / fromSubspace / fromDense)
%     'Format'     'auto' | 'subspace' | 'dense'        (fromFile)
%     'CoefVar'    coefficient variable name   ('vk_all')
%     'BasisVar'   temporal basis variable name ('basis_all')
%     'ImageVar'   4D image variable name       (dense)
%     'DimOrder'   dimension order of dense array, e.g. 'xytz' (ReadData
%                  data1) or 'xyzt'
%     'CropFactor' central crop factor (default 2 subspace, 1 dense)
%     'Normalize'  divide by global max (default true subspace, false dense)
%     'Progress'   function handle progress(fraction, message)

    properties (SetAccess = private)
        Mode = '';                  % 'subspace' | 'dense'
        C = [];                     % subspace coefficients [nx ny K nz] single
        Phi = [];                   % temporal basis [T K nz] single
        D = [];                     % dense magnitude [nx ny nz T] single
        NativeSize = [0 0 0];
        NT = 0;                     % frames after temporal binning
        NTRaw = 0;                  % frames as reconstructed
        BinSize = 1;                % consecutive frames averaged per binned frame
        BinKeepPartial = false;     % keep an incomplete last bin
        BinIdx = {};                % raw frame indices of each binned frame
        BinW = [];                  % [NTRaw x NT] averaging matrix
        NK = 0;
        ScaleFactor = 1;            % normalization factor applied (1/max)
        CropFactor = 1;
        Normalized = false;
        Plane = 'axial';
        AxisCodes = 'PLS';          % direction of increasing native rows/cols/slices
        SwapInPlane = false;
        UserFlip = [false false false];
        Perm = [2 1 3];             % canonical axis a <- native dim Perm(a)
        Flip = [false false false]; % flip along canonical axis a
        Maps = [];                  % native-space enhancement maps
        CanonMapCache = struct();
        MapsVersion = 0;            % incremented whenever maps or orientation change
        Motion = [];                % last motion estimate (graspilot_motion_estimate)
        MotionCorrected = false;    % frames realigned by applyMotionCorrection
        MotionApplied = [];         % the estimate that was applied
        Orig = [];                  % original frame source while motion-corrected
        SrcAffine = [];             % native voxel (0-based) -> RAS mm of a NIfTI / DICOM source
        SourceType = 'mat';         % 'mat' | 'nifti' | 'dicom' | 'array'
    end

    properties
        VoxNative = [1 1 1];        % mm: [row spacing, col spacing, slice spacing]
        dt = NaN;                   % s per frame
        t0 = 0;                     % s, time of frame 1
        RawTimes = [];              % measured times of the raw frames (s from frame 1); [] = uniform dt
        SourceFile = '';
        SourceVars = '';
        ParamOpts = struct('WashIn', 'onset', 'WashOut', 'regression', 'AUCWindow', Inf, 'OnsetFrac', 0.1, ...
            'WashInFrames', [], 'EndFrame', Inf, 'EndTime', Inf, 'EarlyTime', 60, 'LateTime', Inf, ...
            'CurveTol', 10, 'RefSlopeFrac', 0.5);
    end

    % ------------------------------------------------------------------
    methods (Static)

        function obj = fromFile(matFile, varargin)
            p = GRASPilotData.parseOpts(varargin, GRASPilotData.defaultOpts());
            if exist(matFile, 'file') ~= 2
                error('GRASPilotData:file', 'File not found: %s', matFile);
            end
            info = whos('-file', matFile);
            names = {info.name};
            fmt = lower(p.Format);
            if strcmp(fmt, 'auto')
                if any(strcmp(names, p.CoefVar)) && any(strcmp(names, p.BasisVar))
                    fmt = 'subspace';
                elseif ~isempty(p.ImageVar) && any(strcmp(names, p.ImageVar))
                    fmt = 'dense';
                else
                    error('GRASPilotData:format', ['Could not find %s/%s in %s. Use the load ' ...
                        'dialog or pass ''Format'',''dense'',''ImageVar'',name.'], ...
                        p.CoefVar, p.BasisVar, matFile);
                end
            end
            GRASPilotData.report(p.Progress, 0, 'Reading .mat file (this can take a while for large files)...');
            switch fmt
                case 'subspace'
                    S = load(matFile, p.CoefVar, p.BasisVar);
                    if ~isfield(S, p.CoefVar) || ~isfield(S, p.BasisVar)
                        error('GRASPilotData:vars', 'Variables %s and %s are required.', p.CoefVar, p.BasisVar);
                    end
                    vk = S.(p.CoefVar); bs = S.(p.BasisVar);
                    S = []; %#ok<NASGU>
                    obj = GRASPilotData.fromSubspace(vk, bs, p);
                    obj.SourceVars = sprintf('%s + %s', p.CoefVar, p.BasisVar);
                case 'dense'
                    if isempty(p.ImageVar)
                        error('GRASPilotData:vars', 'ImageVar must be given for dense format.');
                    end
                    S = load(matFile, p.ImageVar);
                    img = S.(p.ImageVar);
                    S = []; %#ok<NASGU>
                    obj = GRASPilotData.fromDense(img, p.DimOrder, p);
                    obj.SourceVars = sprintf('%s (%s)', p.ImageVar, p.DimOrder);
                otherwise
                    error('GRASPilotData:format', 'Unknown format "%s".', fmt);
            end
            obj.SourceFile = matFile;
        end

        function obj = fromSubspace(vk_all, basis_all, varargin)
            p = GRASPilotData.parseOpts(varargin, GRASPilotData.defaultOpts());
            if isempty(p.CropFactor), p.CropFactor = 2; end
            if isempty(p.Normalize), p.Normalize = true; end
            if ndims(vk_all) > 4
                error('GRASPilotData:size', 'vk_all must be [slices nx ny K].');
            end
            nz = size(vk_all, 1); nxF = size(vk_all, 2); nyF = size(vk_all, 3); K = size(vk_all, 4);

            % temporal basis -> Phi [T K nz]; ReadData uses Phi_z = basis_all(z,:,:) as [T K]
            B = basis_all;
            if ndims(B) == 2 && size(B, 2) == K && ~(nz > 1 && size(B, 1) == nz && K == 1) %#ok<ISMAT>
                Phi = repmat(single(B), [1 1 nz]);        % one basis shared by all slices
            else
                if size(B, 1) ~= nz
                    error('GRASPilotData:basis', ['basis_all has %d entries along dim 1 but vk_all ' ...
                        'has %d slices.'], size(B, 1), nz);
                end
                if size(B, 3) == K
                    Phi = permute(single(B), [2 3 1]);    % [nz T K] -> [T K nz] (ReadData layout)
                elseif size(B, 2) == K
                    Phi = permute(single(B), [3 2 1]);    % [nz K T] -> [T K nz]
                else
                    error('GRASPilotData:basis', ['basis_all [%s] is not compatible with K = %d ' ...
                        'coefficients in vk_all.'], num2str(size(B)), K);
                end
            end
            T = size(Phi, 1);

            % central crop (CropImg.m)
            cf = p.CropFactor;
            nx = max(1, floor(nxF / cf)); ny = max(1, floor(nyF / cf));
            rr = floor((nxF - nx) / 2) + (1:nx);
            cc = floor((nyF - ny) / 2) + (1:ny);

            isCplx = ~isreal(vk_all);
            Cc = zeros(nx, ny, K, nz, 'single');
            if isCplx, Cc = complex(Cc, Cc); end
            P = nxF * nyF;
            blk = 65536;
            gmax = 0;
            for z = 1:nz
                Bz = reshape(single(vk_all(z, :, :, :)), P, K);
                if p.Normalize
                    PzH = Phi(:, :, z)';                  % [K T], conjugate transpose (TempPCASub)
                    for b0 = 1:blk:P
                        rows = b0:min(P, b0 + blk - 1);
                        m = max(max(abs(Bz(rows, :) * PzH)));
                        if m > gmax, gmax = m; end
                    end
                end
                Bz = reshape(Bz, nxF, nyF, K);
                Cc(:, :, :, z) = Bz(rr, cc, :);
                GRASPilotData.report(p.Progress, z / nz, sprintf('Preparing slice %d/%d', z, nz));
            end
            sf = 1;
            if p.Normalize && gmax > 0
                sf = 1 / double(gmax);
                for z = 1:nz
                    Cc(:, :, :, z) = Cc(:, :, :, z) * single(sf);
                end
            end

            obj = GRASPilotData();
            obj.Mode = 'subspace';
            obj.C = Cc;
            obj.Phi = Phi;
            obj.NativeSize = [nx ny nz];
            obj.NTRaw = T;
            obj.NK = K;
            obj.setBinning(1);
            obj.ScaleFactor = sf;
            obj.CropFactor = cf;
            obj.Normalized = logical(p.Normalize);
            obj.setOrientation('axial', false, [false false false]);
        end

        function obj = fromDense(img, dimOrder, varargin)
            p = GRASPilotData.parseOpts(varargin, GRASPilotData.defaultOpts());
            if isempty(p.CropFactor), p.CropFactor = 1; end
            if isempty(p.Normalize), p.Normalize = false; end
            if nargin < 2 || isempty(dimOrder), dimOrder = 'xytz'; end
            o = lower(strrep(dimOrder, ' ', ''));
            nd = ndims(img);
            if nd == 3 && numel(o) == 4
                o = o(o ~= 't');                           % 3D input, single frame
            end
            if numel(o) == 3
                o = [o 't'];
            end
            if numel(o) ~= 4 || ~all(sort(o) == sort('xyzt'))
                error('GRASPilotData:dimorder', 'DimOrder must be a permutation of xyzt (got "%s").', dimOrder);
            end
            pm = [find(o == 'x') find(o == 'y') find(o == 'z') find(o == 't')];
            if ~isa(img, 'single')
                if isreal(img), img = single(img); else, img = single(abs(img)); end
            elseif ~isreal(img)
                img = abs(img);
            end
            if ~isequal(pm, 1:4)
                img = permute(img, pm);
            end
            nxF = size(img, 1); nyF = size(img, 2);
            cf = p.CropFactor;
            if cf ~= 1
                nx = max(1, floor(nxF / cf)); ny = max(1, floor(nyF / cf));
                img = img(floor((nxF - nx) / 2) + (1:nx), floor((nyF - ny) / 2) + (1:ny), :, :);
            end
            sf = 1;
            if p.Normalize
                m = max(img(:));
                if m > 0
                    sf = 1 / double(m);
                    img = img * single(sf);
                end
            end
            obj = GRASPilotData();
            obj.Mode = 'dense';
            obj.D = img;
            obj.NativeSize = [size(img, 1) size(img, 2) size(img, 3)];
            obj.NTRaw = size(img, 4);
            obj.NK = 0;
            obj.setBinning(1);
            obj.ScaleFactor = sf;
            obj.CropFactor = cf;
            obj.Normalized = logical(p.Normalize);
            obj.setOrientation('axial', false, [false false false]);
        end

        function obj = fromVolume(V, A, varargin)
            % Dense series with scanner geometry.
            %   D = GRASPilotData.fromVolume(V, A, 'dt', s, 't0', s, 'SourceFile', name, 'SourceType', t)
            % V [nx ny nz T] (or 3D), A 4x4 voxel (0-based) -> RAS mm. The
            % anatomical orientation (axis codes) and voxel size are taken
            % from A; signal values are kept as they are (not normalized).
            q = struct('dt', NaN, 't0', 0, 'SourceFile', '', 'SourceType', 'array', 'Times', []);
            for i = 1:2:numel(varargin)
                fn = fieldnames(q);
                m = find(strcmpi(fn, varargin{i}), 1);
                if isempty(m), error('GRASPilotData:opt', 'Unknown option %s.', varargin{i}); end
                q.(fn{m}) = varargin{i + 1};
            end
            if ~isequal(size(A), [4 4]) || any(~isfinite(A(:))) || abs(det(A(1:3, 1:3))) < 1e-9
                error('GRASPilotData:affine', 'The image geometry (affine) is not valid.');
            end
            obj = GRASPilotData.fromDense(V, 'xyzt', 'CropFactor', 1, 'Normalize', false);
            obj.SrcAffine = double(A);
            obj.VoxNative = sqrt(sum(A(1:3, 1:3) .^ 2, 1));
            obj.setAxisCodes(GRASPilotData.codesFromAffine(A));
            if isfinite(q.dt) && q.dt > 0, obj.dt = q.dt; else, obj.dt = NaN; end
            obj.t0 = q.t0;
            if ~isempty(q.Times), obj.setFrameTimes(q.Times); end
            obj.SourceFile = q.SourceFile;
            obj.SourceType = q.SourceType;
            obj.SourceVars = q.SourceType;
        end

        function codes = codesFromAffine(A)
            % anatomical direction of increasing index for each array axis
            % (nearest of L/R, A/P, S/I), from a voxel -> RAS affine
            R = A(1:3, 1:3) ./ repmat(sqrt(sum(A(1:3, 1:3) .^ 2, 1)), 3, 1);
            W = abs(R);
            lett = ['RL'; 'AP'; 'SI'];
            codes = '   ';
            for k = 1:3                                  % strongest pairings first
                [~, i] = max(W(:));
                [w, d] = ind2sub([3 3], i);              % world axis w <-> array dim d
                codes(d) = lett(w, 1 + (R(w, d) < 0));
                W(w, :) = -1; W(:, d) = -1;
            end
        end

        function obj = fromNifti(files, varargin)
            % 4D NIfTI (.nii / .nii.gz), or several 3D files (one per frame,
            % in the given order; a cell array). Options: 'dt' (s/frame,
            % overrides the header), 'Progress'.
            q = struct('dt', [], 'Progress', []);
            for i = 1:2:numel(varargin)
                fn = fieldnames(q);
                m = find(strcmpi(fn, varargin{i}), 1);
                if isempty(m), error('GRASPilotData:opt', 'Unknown option %s.', varargin{i}); end
                q.(fn{m}) = varargin{i + 1};
            end
            if ischar(files), files = {files}; end
            V = [];
            for k = 1:numel(files)
                GRASPilotData.report(q.Progress, (k - 1) / numel(files), sprintf('Reading %s...', files{k}));
                [img, hdr] = graspilot_nifti_read(files{k});
                if ndims(img) > 4
                    img = reshape(img, size(img, 1), size(img, 2), size(img, 3), []);
                end
                if ~isreal(img), img = abs(img); end
                img = single(img);
                if k == 1
                    V = img; h1 = hdr;
                else
                    if ~isequal(size(img(:, :, :, 1)), size(V(:, :, :, 1))) || ...
                            max(abs(hdr.Affine(:) - h1.Affine(:))) > 1e-2
                        error('GRASPilotData:nifti', '%s does not have the same grid as %s.', files{k}, files{1});
                    end
                    V = cat(4, V, img);
                end
            end
            dt = NaN; src = 'header'; ftimes = [];
            if ~isempty(q.dt)
                dt = q.dt; src = 'given';
            else
                % sidecar of a GRASPilot / dcm2niix export, else the header
                js = regexprep(files{1}, '(\.nii\.gz|\.nii)$', '.json', 'ignorecase');
                if exist(js, 'file') == 2
                    t = fileread(js);
                    tk = regexp(t, '"TemporalResolution_s"\s*:\s*([-+0-9.eE]+)', 'tokens', 'once');
                    if ~isempty(tk), dt = str2double(tk{1}); src = 'sidecar'; end
                    tb = regexp(t, '"BinSize"\s*:\s*([0-9]+)', 'tokens', 'once');
                    tf = regexp(t, '"FrameTimes_s"\s*:\s*\[([^\]]*)\]', 'tokens', 'once');
                    if ~isempty(tf) && (isempty(tb) || str2double(tb{1}) == 1)
                        ftimes = str2double(regexp(tf{1}, ',', 'split'));
                    end
                end
                if ~(isfinite(dt) && dt > 0) && numel(files) == 1 && h1.dtSeconds >= 0.5 && h1.dtSeconds <= 600
                    dt = h1.dtSeconds;       % shorter values are usually the sequence TR, not the frame time
                end
            end
            if ~(isfinite(dt) && dt > 0), dt = NaN; src = 'not found'; end
            if numel(ftimes) ~= size(V, 4) || any(~isfinite(ftimes)) || any(diff(ftimes) <= 0), ftimes = []; end
            obj = GRASPilotData.fromVolume(V, h1.Affine, 'dt', dt, 'Times', ftimes, 'SourceFile', files{1}, ...
                'SourceType', 'nifti');
            obj.SourceVars = sprintf('NIfTI (s/frame from %s)', src);
            GRASPilotData.report(q.Progress, 1, 'Done');
        end

        function [obj, info] = fromDicom(src, varargin)
            % DICOM series (folder, file, or a scan from graspilot_dicom_read).
            % Options as graspilot_dicom_read ('Series', 'Progress', ...).
            [V, info] = graspilot_dicom_read(src, varargin{:});
            name = info.Folder;
            if ischar(src) && exist(src, 'file') == 2, name = src; end
            obj = GRASPilotData.fromVolume(V, info.Affine, 'dt', info.dt, 'Times', info.Times, ...
                'SourceFile', name, 'SourceType', 'dicom');
            obj.SourceVars = sprintf('DICOM series %g: %s', info.SeriesNumber, info.Series);
        end

        function p = defaultOpts()
            p = struct('Format', 'auto', 'CoefVar', 'vk_all', 'BasisVar', 'basis_all', ...
                'ImageVar', '', 'DimOrder', 'xytz', 'CropFactor', [], 'Normalize', [], ...
                'Progress', []);
        end

        function p = parseOpts(args, p)
            if numel(args) == 1 && isstruct(args{1})
                s = args{1}; f = fieldnames(s);
                fn = fieldnames(p);
                for i = 1:numel(f)
                    m = find(strcmpi(fn, f{i}), 1);
                    if ~isempty(m), p.(fn{m}) = s.(f{i}); end
                end
                return
            end
            if mod(numel(args), 2) ~= 0
                error('GRASPilotData:opts', 'Options must be name/value pairs.');
            end
            fn = fieldnames(p);
            for i = 1:2:numel(args)
                m = find(strcmpi(fn, args{i}), 1);
                if isempty(m)
                    error('GRASPilotData:opts', 'Unknown option "%s".', args{i});
                end
                p.(fn{m}) = args{i + 1};
            end
        end

        function s = ternary(c, a, b)
            if c, s = a; else, s = b; end
        end

        function report(fcn, frac, msg)
            if ~isempty(fcn)
                try
                    fcn(frac, msg);
                catch
                end
            end
        end

        function [perm, flp] = codesToPermFlip(codes)
            if ~ischar(codes) || numel(codes) ~= 3
                error('GRASPilotData:codes', 'Axis codes must be 3 letters, e.g. ''LPS''.');
            end
            codes = upper(codes);
            perm = zeros(1, 3); flp = false(1, 3);
            for d = 1:3
                switch codes(d)
                    case 'L', a = 1; f = false;
                    case 'R', a = 1; f = true;
                    case 'P', a = 2; f = false;
                    case 'A', a = 2; f = true;
                    case 'S', a = 3; f = false;
                    case 'I', a = 3; f = true;
                    otherwise
                        error('GRASPilotData:codes', 'Invalid axis code ''%s'' in %s (use L R A P S I).', codes(d), codes);
                end
                if perm(a) ~= 0
                    error('GRASPilotData:codes', 'Axis codes %s use the same anatomical axis twice.', codes);
                end
                perm(a) = d; flp(a) = f;
            end
        end

        function codes = permFlipToCodes(perm, flp)
            pos = 'LPS'; neg = 'RAI';
            codes = '   ';
            for a = 1:3
                if flp(a), codes(perm(a)) = neg(a); else, codes(perm(a)) = pos(a); end
            end
        end

        function codes = rotateCodesCW(codes)
            % rotate the axial display 90 degrees clockwise
            m = 'LPRA'; r = 'PRAL';
            for d = 1:3
                k = find(m == codes(d), 1);
                if ~isempty(k), codes(d) = r(k); end
            end
        end

        function codes = flipCodes(codes, which)
            % which: 'LR', 'AP' or 'SI'
            for d = 1:3
                if codes(d) == which(1)
                    codes(d) = which(2);
                elseif codes(d) == which(2)
                    codes(d) = which(1);
                end
            end
        end

        function [perm, flp] = presetOrientation(plane, swapInPlane, userFlip)
            % native dims: 1 = rows (top->bottom), 2 = cols (left->right), 3 = partitions
            switch lower(plane)
                case 'axial'      % rows A->P, cols R->L
                    perm = [2 1 3]; flp = [false false false];
                case 'coronal'    % rows S->I, cols R->L
                    perm = [2 3 1]; flp = [false false true];
                case 'sagittal'   % rows S->I, cols A->P
                    perm = [3 2 1]; flp = [false false true];
                otherwise
                    error('GRASPilotData:plane', 'Plane must be axial, coronal or sagittal.');
            end
            if swapInPlane
                a1 = find(perm == 1); a2 = find(perm == 2);
                perm([a1 a2]) = perm([a2 a1]);
            end
            flp = xor(flp, logical(userFlip(:)'));
        end
    end

    % ------------------------------------------------------------------
    methods

        function setOrientation(obj, plane, swapInPlane, userFlip)
            if nargin < 3 || isempty(swapInPlane), swapInPlane = false; end
            if nargin < 4 || isempty(userFlip), userFlip = [false false false]; end
            [obj.Perm, obj.Flip] = GRASPilotData.presetOrientation(plane, swapInPlane, userFlip);
            obj.Plane = lower(plane);
            obj.SwapInPlane = logical(swapInPlane);
            obj.UserFlip = logical(userFlip(:)');
            obj.AxisCodes = GRASPilotData.permFlipToCodes(obj.Perm, obj.Flip);
            obj.CanonMapCache = struct();
            obj.MapsVersion = obj.MapsVersion + 1;
        end

        function setAxisCodes(obj, codes)
            % codes: direction of increasing native rows, cols, slices,
            % one of L/R, A/P, S/I each, e.g. 'LPS' or 'PLS'
            [perm, flp] = GRASPilotData.codesToPermFlip(codes);
            planes = {'sagittal', 'coronal', 'axial'};
            plane = planes{find(perm == 3)};                       % anatomical axis of the slices
            for sw = 0:1
                for f = 0:7
                    uf = logical(bitget(f, 1:3));
                    [p2, f2] = GRASPilotData.presetOrientation(plane, sw, uf);
                    if isequal(p2, perm) && isequal(f2, flp)
                        obj.setOrientation(plane, sw, uf);
                        return
                    end
                end
            end
            error('GRASPilotData:codes', 'Could not represent axis codes %s.', codes);
        end

        function n = canonSize(obj)
            n = obj.NativeSize(obj.Perm);
        end

        function v = canonVox(obj)
            v = obj.VoxNative(obj.Perm);
        end

        function tf = hasTiming(obj)
            tf = isfinite(obj.dt) && obj.dt > 0;
        end

        function setBinning(obj, n, keepPartial)
            % Average n consecutive frames (magnitudes). An incomplete last
            % bin is dropped unless keepPartial is true.
            if nargin < 3 || isempty(keepPartial), keepPartial = false; end
            T0 = obj.NTRaw;
            n = max(1, min(T0, round(n)));
            nb = floor(T0 / n);
            if keepPartial && mod(T0, n) > 0, nb = nb + 1; end
            nb = max(1, nb);
            idx = cell(1, nb);
            W = zeros(T0, nb);
            for b = 1:nb
                f = (b - 1) * n + 1:min(T0, b * n);
                idx{b} = f;
                W(f, b) = 1 / numel(f);
            end
            obj.BinSize = n;
            obj.BinKeepPartial = logical(keepPartial);
            obj.BinIdx = idx;
            obj.BinW = W;
            obj.NT = nb;
            obj.Maps = [];
            obj.CanonMapCache = struct();
            obj.MapsVersion = obj.MapsVersion + 1;
        end

        function f = binFrames(obj, t)
            f = obj.BinIdx{t};
        end

        function f = droppedFrames(obj)
            % raw frames not used by any bin
            used = [obj.BinIdx{:}];
            f = setdiff(1:obj.NTRaw, used);
        end

        function d = effectiveDt(obj)
            d = obj.dt * obj.BinSize;
        end

        function t = frameTimes(obj)
            % time of each (binned) frame = mean time of its raw frames
            t = obj.rawTimes() * obj.BinW;
        end

        function t = rawTimes(obj)
            % times of the raw frames (s): measured times if known, else t0 + (f - 1) dt
            if ~obj.hasTiming()
                t = NaN(1, obj.NTRaw);
            elseif numel(obj.RawTimes) == obj.NTRaw
                t = obj.t0 + obj.RawTimes;
            else
                t = obj.t0 + (0:obj.NTRaw - 1) * obj.dt;
            end
        end

        function setFrameTimes(obj, t)
            % measured (possibly uneven) times of the raw frames, s; [] = uniform dt
            if isempty(t)
                obj.RawTimes = [];
                return
            end
            t = double(t(:)');
            if numel(t) ~= obj.NTRaw || any(~isfinite(t)) || any(diff(t) <= 0)
                error('GRASPilotData:times', 'Frame times must be %d increasing values.', obj.NTRaw);
            end
            obj.RawTimes = t - t(1);
            if obj.NTRaw > 1, obj.dt = median(diff(t)); end
        end

        function [tf, rng] = unevenTiming(obj)
            % true when the measured frame intervals differ by > 5% from their median
            tf = false; rng = [NaN NaN];
            if numel(obj.RawTimes) ~= obj.NTRaw || obj.NTRaw < 3, return; end
            d = diff(obj.RawTimes);
            rng = [min(d) max(d)];
            tf = max(abs(d - median(d))) > 0.05 * median(d);
        end

        function b = memoryBytes(obj)
            if strcmp(obj.Mode, 'subspace')
                b = numel(obj.C) * 4 * (1 + ~isreal(obj.C)) + numel(obj.Phi) * 4 * (1 + ~isreal(obj.Phi));
            else
                b = numel(obj.D) * 4;
            end
            if ~isempty(obj.Orig)
                O = obj.Orig;
                b = b + numel(O.C) * 4 * (1 + ~isreal(O.C)) + numel(O.Phi) * 4 * (1 + ~isreal(O.Phi)) + numel(O.D) * 4;
            end
        end

        % ---------- motion ----------
        function V = rawCanonVolume(obj, f)
            % canonical volume of raw (unbinned) frame f
            if strcmp(obj.Mode, 'dense')
                V = obj.native2canon(obj.D(:, :, :, f));
                return
            end
            nx = obj.NativeSize(1); ny = obj.NativeSize(2); nz = obj.NativeSize(3);
            K = obj.NK;
            V = zeros(nx, ny, nz, 'single');
            for z = 1:nz
                V(:, :, z) = reshape(abs(reshape(obj.C(:, :, :, z), nx * ny, K) * obj.Phi(f, :, z)'), nx, ny);
            end
            V = obj.native2canon(V);
        end

        function s = motionDescription(obj)
            % one line for exports / sidecars
            if obj.MotionCorrected
                M = obj.MotionApplied;
                s = sprintf(['applied: rigid realignment of every raw frame to the %s (normalized mutual ' ...
                    'information); before correction max mean head displacement %.2f mm from the reference, ' ...
                    '%.2f mm between frames'], M.Reference, max(M.AbsMean), max([0; M.RelMean(isfinite(M.RelMean))]));
            elseif ~isempty(obj.Motion)
                M = obj.Motion;
                s = sprintf(['not applied; estimated max mean head displacement %.2f mm from the %s, %.2f mm ' ...
                    'between frames%s'], max(M.AbsMean), M.Reference, max([0; M.RelMean(isfinite(M.RelMean))]), ...
                    GRASPilotData.ternary(M.Excess, ' (EXCESS MOTION WARNING)', ''));
            else
                s = 'not applied; motion not estimated';
            end
        end

        function setMotion(obj, M)
            % store a motion estimate (graspilot_motion_estimate)
            obj.Motion = M;
        end

        function b = motionCorrectionBytes(obj)
            % memory of the realigned (dense) series
            b = prod(obj.NativeSize) * obj.NTRaw * 4;
        end

        function applyMotionCorrection(obj, M, progress)
            % Realign every raw frame with the rigid transforms of motion
            % estimate M (default: obj.Motion): frame f is resampled at
            % y = Rot (x - Centre) + Centre + t (trilinear; outside the
            % volume the nearest edge value). Frames with sub-0.05-voxel
            % motion are copied. The result is kept as a dense series.
            if nargin < 2 || isempty(M), M = obj.Motion; end
            if nargin < 3, progress = []; end
            if isempty(M) || ~isfield(M, 'Params')
                error('GRASPilotData:motion', 'No motion estimate: run graspilot_motion_estimate first.');
            end
            if obj.MotionCorrected
                error('GRASPilotData:motion', 'Already motion-corrected: revert first.');
            end
            if size(M.Params, 1) ~= obj.NTRaw
                error('GRASPilotData:motion', 'The motion estimate has %d frames, the series %d.', ...
                    size(M.Params, 1), obj.NTRaw);
            end
            n = obj.canonSize();
            vox = obj.canonVox();
            [g1, g2, g3] = ndgrid(single(0:n(1) - 1) * vox(1), single(0:n(2) - 1) * vox(2), ...
                single(0:n(3) - 1) * vox(3));
            X = [g1(:) g2(:) g3(:)];
            clear g1 g2 g3
            Cc = single(M.Centre);
            out = zeros([obj.NativeSize obj.NTRaw], 'single');
            for f = 1:obj.NTRaw
                GRASPilotData.report(progress, (f - 1) / obj.NTRaw, ...
                    sprintf('Motion correction: frame %d of %d...', f, obj.NTRaw));
                V = obj.rawCanonVolume(f);
                pf = M.Params(f, :);
                if any(~isfinite(pf)), pf = zeros(1, 6); end
                Rm = single(graspilot_rigid_matrix(pf));
                Y = (X - repmat(Cc, size(X, 1), 1)) * Rm' + repmat(Cc + single(pf(1:3)), size(X, 1), 1);
                q = Y ./ repmat(single(vox), size(Y, 1), 1) + 1;
                dq = abs(q - X ./ repmat(single(vox), size(X, 1), 1) - 1);
                if max(dq(:)) > 0.05
                    for a = 1:3, q(:, a) = min(max(q(:, a), 1), n(a)); end
                    V = reshape(interpn(double(V), double(q(:, 1)), double(q(:, 2)), double(q(:, 3)), 'linear'), n);
                end
                out(:, :, :, f) = obj.canon2native(single(V));
            end
            obj.Orig = struct('Mode', obj.Mode, 'C', obj.C, 'Phi', obj.Phi, 'D', obj.D, 'NK', obj.NK);
            obj.Mode = 'dense';
            obj.D = out;
            obj.C = []; obj.Phi = [];
            obj.NK = 0;
            obj.MotionCorrected = true;
            obj.MotionApplied = M;
            obj.clearMaps();
            GRASPilotData.report(progress, 1, 'Motion correction: done');
        end

        function revertMotionCorrection(obj)
            % back to the original (uncorrected) frames
            if ~obj.MotionCorrected, return; end
            O = obj.Orig;
            obj.Mode = O.Mode; obj.C = O.C; obj.Phi = O.Phi; obj.D = O.D; obj.NK = O.NK;
            obj.Orig = [];
            obj.MotionCorrected = false;
            obj.MotionApplied = [];
            obj.clearMaps();
        end

        % ---------- orientation helpers ----------
        function A = native2canon(obj, A)
            A = permute(A, obj.Perm);
            for a = 1:3
                if obj.Flip(a), A = flip(A, a); end
            end
        end

        function A = canon2native(obj, A)
            for a = 1:3
                if obj.Flip(a), A = flip(A, a); end
            end
            A = ipermute(A, obj.Perm);
        end

        function nl = canon2nativeLin(obj, cl)
            cs = obj.canonSize();
            [i1, i2, i3] = ind2sub(cs, cl(:));
            sub = {i1, i2, i3};
            ns = cell(1, 3);
            for a = 1:3
                s = sub{a};
                if obj.Flip(a), s = cs(a) - s + 1; end
                ns{obj.Perm(a)} = s;
            end
            nl = sub2ind(obj.NativeSize, ns{1}, ns{2}, ns{3});
        end

        % ---------- canonical access ----------
        function S = canonSlice(obj, a, q, t)
            % 2D slice perpendicular to canonical axis a at index q, frame t.
            % Rows/cols follow the remaining canonical axes in ascending order.
            nd = obj.Perm(a);
            idx = q;
            if obj.Flip(a), idx = obj.NativeSize(nd) - q + 1; end
            S = obj.nativeSlice(nd, idx, t);
            oa = setdiff(1:3, a);
            if obj.Perm(oa(1)) > obj.Perm(oa(2)), S = S.'; end
            if obj.Flip(oa(1)), S = flipud(S); end
            if obj.Flip(oa(2)), S = fliplr(S); end
        end

        function V = canonVolume(obj, t)
            V = obj.native2canon(obj.nativeVolume(t));
        end

        function TC = timecourses(obj, canonLin)
            TC = obj.nativeTC(obj.canon2nativeLin(canonLin));
        end

        function v = valuesAtFrame(obj, canonLin, t)
            v = obj.nativeValuesAtFrame(obj.canon2nativeLin(canonLin), t);
        end

        % ---------- native access ----------
        function S = nativeSlice(obj, d, idx, t)
            nx = obj.NativeSize(1); ny = obj.NativeSize(2); nz = obj.NativeSize(3);
            f = obj.BinIdx{t};
            if strcmp(obj.Mode, 'dense')
                switch d
                    case 3, S = mean(obj.D(:, :, idx, f), 4);
                    case 1, S = reshape(mean(obj.D(idx, :, :, f), 4), ny, nz);
                    case 2, S = reshape(mean(obj.D(:, idx, :, f), 4), nx, nz);
                end
                return
            end
            K = obj.NK;
            switch d
                case 3
                    S = abs(reshape(obj.C(:, :, :, idx), nx * ny, K) * obj.Phi(f, :, idx)');
                    S = reshape(mean(S, 2), nx, ny);
                case {1, 2}
                    if d == 1
                        X = reshape(obj.C(idx, :, :, :), ny, K, nz); m = ny;
                    else
                        X = reshape(obj.C(:, idx, :, :), nx, K, nz); m = nx;
                    end
                    S = zeros(m, nz, 'single');
                    for j = 1:numel(f)
                        w = reshape(conj(obj.Phi(f(j), :, :)), 1, K, nz);
                        S = S + abs(reshape(sum(X .* w, 2), m, nz));
                    end
                    S = S / numel(f);
            end
        end

        function V = nativeVolume(obj, t)
            f = obj.BinIdx{t};
            if strcmp(obj.Mode, 'dense')
                V = mean(obj.D(:, :, :, f), 4);
                return
            end
            nx = obj.NativeSize(1); ny = obj.NativeSize(2); nz = obj.NativeSize(3);
            K = obj.NK;
            V = zeros(nx, ny, nz, 'single');
            for z = 1:nz
                S = abs(reshape(obj.C(:, :, :, z), nx * ny, K) * obj.Phi(f, :, z)');
                V(:, :, z) = reshape(mean(S, 2), nx, ny);
            end
        end

        function TC = nativeTC(obj, lin)
            lin = lin(:);
            if strcmp(obj.Mode, 'dense')
                D2 = reshape(obj.D, [], obj.NTRaw);
                TC = D2(lin, :);
            else
                [Cv, z] = obj.gatherCoefs(lin);
                TC = zeros(numel(lin), obj.NTRaw, 'single');
                [zu, ~, g] = unique(z);
                for ii = 1:numel(zu)
                    m = (g == ii);
                    TC(m, :) = abs(Cv(m, :) * obj.Phi(:, :, zu(ii))');
                end
            end
            if obj.BinSize > 1 || obj.NT ~= obj.NTRaw
                TC = single(double(TC) * obj.BinW);
            end
        end

        function v = nativeValuesAtFrame(obj, lin, t)
            lin = lin(:);
            f = obj.BinIdx{t};
            if strcmp(obj.Mode, 'dense')
                N = prod(obj.NativeSize);
                v = mean(obj.D(repmat(lin, 1, numel(f)) + repmat((f - 1) * N, numel(lin), 1)), 2);
                return
            end
            [Cv, z] = obj.gatherCoefs(lin);
            v = zeros(numel(lin), 1, 'single');
            [zu, ~, g] = unique(z);
            for ii = 1:numel(zu)
                m = (g == ii);
                v(m) = mean(abs(Cv(m, :) * obj.Phi(f, :, zu(ii))'), 2);
            end
        end

        % ---------- enhancement maps ----------
        function tf = hasMaps(obj, baseFrames)
            % maps exist for these baseline frames, parameter options and timing
            tf = ~isempty(obj.Maps) && isequal(obj.Maps.BaseFrames, baseFrames) && ...
                isequaln(obj.Maps.Key, obj.mapsKey());
        end

        function k = mapsKey(obj)
            k = {obj.ParamOpts, obj.dt, obj.t0, obj.BinSize, obj.NT, obj.RawTimes};
        end

        function t = paramTimes(obj)
            % time axis used for curve parameters: seconds if known, else frames
            if obj.hasTiming()
                t = obj.frameTimes();
            else
                t = 1:obj.NT;
            end
        end

        function setParamOptions(obj, varargin)
            % 'WashIn' ('onset'|'max'|'sss'), 'WashOut' ('regression'|'endpoints'),
            % 'AUCWindow' (s after the baseline, Inf = to the end), 'OnsetFrac',
            % 'EndFrame' / 'EndTime' (analysis endpoint in frames / s; Inf = all),
            % 'EarlyTime', 'LateTime' (s after the baseline), 'CurveTol' (%),
            % 'WashInFrames' ([a b], set from the reference curve for 'sss'),
            % 'RefSlopeFrac'. See graspilot_curve_params.
            p = GRASPilotData.parseOpts(varargin, obj.ParamOpts);
            for f = {'AUCWindow', 'EndFrame', 'EndTime', 'LateTime'}
                v = p.(f{1});
                if isempty(v) || ~(v > 0), p.(f{1}) = Inf; end
            end
            if isempty(p.EarlyTime) || ~(p.EarlyTime > 0), p.EarlyTime = 60; end
            if isfinite(p.EndFrame), p.EndFrame = max(2, round(p.EndFrame)); end
            p.WashInFrames = double(p.WashInFrames(:)');
            obj.ParamOpts = p;
        end

        function k = analysisEnd(obj)
            % last frame used for the curve parameters (EndFrame / EndTime)
            k = graspilot_end_frame(obj.paramTimes(), obj.ParamOpts);
        end

        function computeMaps(obj, baseFrames, progress)
            % Voxel-wise enhancement and semi-quantitative curve parameters
            % (graspilot_curve_params) on the (binned) series, slice by slice.
            if nargin < 3, progress = []; end
            T = obj.NT;
            bf = baseFrames(baseFrames >= 1 & baseFrames <= T);
            if isempty(bf), bf = 1; end
            nx = obj.NativeSize(1); ny = obj.NativeSize(2); nz = obj.NativeSize(3);
            tt = obj.paramTimes();
            M.BaseFrames = baseFrames;
            M.Key = obj.mapsKey();
            M.TimeUnit = 's';
            if ~obj.hasTiming(), M.TimeUnit = 'frame'; end
            names = {'Baseline', 'Max', 'TTP', 'Onset', 'WashIn', 'WashOut', 'AUC', 'T90', 'REearly', ...
                'SER', 'WashoutPct', 'Plateau', 'CurveType'};
            src = names; src{1} = 'S0'; src{2} = 'Smax';
            for i = 1:numel(names)
                M.(names{i}) = zeros(nx, ny, nz, 'single');
            end
            P = nx * ny;
            blk = 16384;
            for z = 1:nz
                if strcmp(obj.Mode, 'subspace')
                    S = abs(reshape(obj.C(:, :, :, z), P, obj.NK) * obj.Phi(:, :, z)');
                else
                    S = reshape(obj.D(:, :, z, :), P, obj.NTRaw);
                end
                if obj.BinSize > 1 || T ~= obj.NTRaw
                    S = single(double(S) * obj.BinW);
                end
                out = cell(1, numel(names));
                for i = 1:numel(names), out{i} = zeros(P, 1, 'single'); end
                for b0 = 1:blk:P
                    rows = b0:min(P, b0 + blk - 1);
                    R = graspilot_curve_params(S(rows, :), tt, bf, obj.ParamOpts);
                    for i = 1:numel(names)
                        out{i}(rows) = R.(src{i});
                    end
                end
                for i = 1:numel(names)
                    M.(names{i})(:, :, z) = reshape(out{i}, nx, ny);
                end
                GRASPilotData.report(progress, z / nz, sprintf('Enhancement maps: slice %d/%d', z, nz));
            end
            obj.Maps = M;
            obj.CanonMapCache = struct();
            obj.MapsVersion = obj.MapsVersion + 1;
        end

        function A = canonSubtraction(obj, post, pre)
            % Subtraction map frame(post) - frame(pre), canonical order (cached)
            key = sprintf('Sub_%d_%d', post, pre);
            if isfield(obj.CanonMapCache, key)
                A = obj.CanonMapCache.(key);
                return
            end
            fn = fieldnames(obj.CanonMapCache);
            for i = 1:numel(fn)
                if strncmp(fn{i}, 'Sub_', 4), obj.CanonMapCache = rmfield(obj.CanonMapCache, fn{i}); end
            end
            A = obj.native2canon(obj.nativeVolume(post) - obj.nativeVolume(pre));
            obj.CanonMapCache.(key) = A;
        end

        function clearMaps(obj)
            obj.Maps = [];
            obj.CanonMapCache = struct();
            obj.MapsVersion = obj.MapsVersion + 1;
        end

        function A = canonMap(obj, name)
            % name: Baseline | Max | PeakEnh | RelEnh | TTP | Onset | WashIn | WashOut | AUC
            if isempty(obj.Maps)
                error('GRASPilotData:maps', 'Enhancement maps have not been computed.');
            end
            if isfield(obj.CanonMapCache, name)
                A = obj.CanonMapCache.(name);
                return
            end
            M = obj.Maps;
            switch name
                case {'Baseline', 'Max', 'TTP', 'Onset', 'WashIn', 'WashOut', 'AUC', 'T90', 'REearly', ...
                        'SER', 'WashoutPct', 'Plateau', 'CurveType'}
                    A = M.(name);
                case 'MITR'
                    A = (M.Max - M.Baseline) ./ M.TTP;
                case 'PeakEnh'
                    A = M.Max - M.Baseline;
                case 'RelEnh'
                    A = 100 * (M.Max - M.Baseline) ./ M.Baseline;
                    A(~(M.Baseline > 0)) = NaN;
                otherwise
                    error('GRASPilotData:maps', 'Unknown map "%s".', name);
            end
            A = obj.native2canon(A);
            obj.CanonMapCache.(name) = A;
        end

        % ---------- geometry ----------
        function A = affine(obj)
            % voxel (0-based, canonical LPS order) -> RAS mm. Data read from
            % NIfTI or DICOM keep the scanner coordinates of the source
            % (exports overlay on the original images); otherwise the
            % volume is centred at 0.
            v = obj.canonVox(); n = obj.canonSize();
            if ~isempty(obj.SrcAffine)
                % canonical index c -> native index: native dim Perm(a) = c(a), reversed if Flip(a)
                M = zeros(4);
                M(4, 4) = 1;
                for a = 1:3
                    d = obj.Perm(a);
                    if obj.Flip(a)
                        M(d, a) = -1; M(d, 4) = n(a) - 1;
                    else
                        M(d, a) = 1;
                    end
                end
                A = obj.SrcAffine * M;
                return
            end
            A = [-v(1) 0 0  v(1) * (n(1) - 1) / 2;
                 0 -v(2) 0  v(2) * (n(2) - 1) / 2;
                 0 0  v(3) -v(3) * (n(3) - 1) / 2;
                 0 0 0 1];
        end

        function s = describe(obj)
            n = obj.NativeSize;
            if obj.BinSize > 1
                ts = sprintf('T=%d (bin %d of %d)', obj.NT, obj.BinSize, obj.NTRaw);
            else
                ts = sprintf('T=%d', obj.NT);
            end
            if strcmp(obj.Mode, 'subspace')
                s = sprintf('Subspace: %dx%dx%d, K=%d, %s, %.0f MB', n(1), n(2), n(3), obj.NK, ...
                    ts, obj.memoryBytes() / 2^20);
            else
                s = sprintf('Dense: %dx%dx%d, %s, %.0f MB', n(1), n(2), n(3), ts, ...
                    obj.memoryBytes() / 2^20);
            end
        end
    end

    methods (Access = private)
        function [Cv, z] = gatherCoefs(obj, lin)
            nx = obj.NativeSize(1); ny = obj.NativeSize(2); K = obj.NK;
            [x, y, z] = ind2sub(obj.NativeSize, lin);
            base = x + (y - 1) * nx + (z - 1) * nx * ny * K;
            idx = base + (0:K - 1) * nx * ny;              % [nvox K]
            Cv = obj.C(idx);
            if numel(lin) == 1, Cv = reshape(Cv, 1, K); end
            if K == 1, Cv = reshape(Cv, [], 1); end
        end
    end
end
