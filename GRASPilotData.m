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
        PKMaps = [];                % native-space pharmacokinetic model maps (computePK)
        AIFCache = [];              % last VCAI fit (pkAIF)
        CanonMapCache = struct();
        MapsVersion = 0;            % incremented whenever maps or orientation change
        Motion = [];                % last motion estimate (graspilot_motion_estimate)
        MotionCorrected = false;    % frames realigned by applyMotionCorrection
        MotionApplied = [];         % the estimate that was applied
        Orig = [];                  % original frame source while motion-corrected
        SrcAffine = [];             % native voxel (0-based) -> RAS mm of a NIfTI / DICOM source
        SourceType = 'mat';         % 'mat' | 'nifti' | 'dicom' | 'array'
        FastCopyFile = '';          % fast-open copy the data were read from ('' = the source itself)
    end

    properties
        VoxNative = [1 1 1];        % mm: [row spacing, col spacing, slice spacing]
        dt = NaN;                   % s per frame
        t0 = 0;                     % s, time of frame 1
        RawTimes = [];              % measured times of the raw frames (s from frame 1); [] = uniform dt
        SourceFile = '';
        SourceVars = '';
        ParamOpts = struct('WashIn', 'sss', 'WashOut', 'regression', 'AUCWindow', Inf, 'OnsetFrac', 0.1, ...
            'WashInFrames', [], 'EndFrame', Inf, 'EndTime', Inf, 'EarlyTime', 60, 'LateTime', Inf, ...
            'CurveTol', 10, 'RefSlopeFrac', 0.5);
        % pharmacokinetic models (setPKOptions / computePK): global T10 (ms) of tissue and
        % blood, TR (ms), flip angle (deg), relaxivity r1 (1/(mM s)), hematocrit, AIF source
        % ('ref' = reference (SSS) ROI | 'population' = Parker 2006), AIF shift (s), region
        % fitted ('head' | 'rois' | 'all'), LM refinement
        PKOpts = graspilot_pk_defaults();
    end

    properties (Constant)
        % key, model (graspilot_pk_fit), field, label, unit
        PKTable = {'KT_TM', 'tofts', 'Ktrans', 'Tofts Ktrans', '1/min';
            'VE_TM', 'tofts', 've', 'Tofts ve', 'fraction';
            'KEP_TM', 'tofts', 'kep', 'Tofts kep', '1/min';
            'R2_TM', 'tofts', 'R2', 'Tofts fit R2', 'R2';
            'KT_ETM', 'etofts', 'Ktrans', 'Extended Tofts Ktrans', '1/min';
            'VE_ETM', 'etofts', 've', 'Extended Tofts ve', 'fraction';
            'VP_ETM', 'etofts', 'vp', 'Extended Tofts vp', 'fraction';
            'KEP_ETM', 'etofts', 'kep', 'Extended Tofts kep', '1/min';
            'R2_ETM', 'etofts', 'R2', 'Extended Tofts fit R2', 'R2';
            'FP_2CXM', '2cxm', 'Fp', '2CXM plasma flow Fp', 'ml/100ml/min';
            'PS_2CXM', '2cxm', 'PS', '2CXM permeability-surface PS', 'ml/100ml/min';
            'VP_2CXM', '2cxm', 'vp', '2CXM vp', 'fraction';
            'VE_2CXM', '2cxm', 've', '2CXM ve', 'fraction';
            'KT_2CXM', '2cxm', 'Ktrans', '2CXM Ktrans = E Fp', '1/min';
            'R2_2CXM', '2cxm', 'R2', '2CXM fit R2', 'R2';
            'KEP_RATIO', 'ratio', 'kep', 'RATIO kep', '1/min';
            'KPE_RATIO', 'ratio', 'kpe', 'RATIO kpe = E/MTT', '1/min';
            'MTT_RATIO', 'ratio', 'MTT', 'RATIO vascular MTT', 's';
            'E_RATIO', 'ratio', 'E', 'RATIO extraction fraction E', 'fraction';
            'VEVP_RATIO', 'ratio', 'VeVp', 'RATIO ve/vp', 'ratio';
            'K1_RATIO', 'ratio', 'k1', 'RATIO one-compartment rate k', '1/min';
            'CLASS_RATIO', 'ratio', 'Class', 'RATIO class (1 2C, 2 fast, 3 slow, 4 indet.)', 'class';
            'LEAK_RATIO', 'ratio', 'Leak', 'RATIO leakage call (1 yes, 0 no)', 'call';
            'VP_RATIO', 'ratio', 'vp', 'RATIO apparent vp', 'fraction (apparent)';
            'KT_RATIO', 'ratio', 'Ktrans', 'RATIO apparent Ktrans', '1/min (apparent)';
            'FP_RATIO', 'ratio', 'Fp', 'RATIO apparent Fp', 'ml/100ml/min (apparent)';
            'VE_RATIO', 'ratio', 've', 'RATIO apparent ve', 'fraction (apparent)';
            'DELAY_RATIO', 'ratio', 'Delay', 'RATIO delay', 's';
            'R2_RATIO', 'ratio', 'R2', 'RATIO fit R2', 'R2'};
    end

    % ------------------------------------------------------------------
    methods (Static)

        function obj = fromFile(matFile, varargin)
            p = GRASPilotData.parseOpts(varargin, GRASPilotData.defaultOpts());
            if exist(matFile, 'file') ~= 2
                error('GRASPilotData:file', 'File not found: %s', matFile);
            end
            if any(strcmpi(p.Format, {'auto', 'subspace'}))
                fc = GRASPilotData.fastCopyFor(matFile, p);  % current fast-open copy: skip the large file
                if ~isempty(fc)
                    obj = GRASPilotData.fromFastCopy(fc, p.Progress);
                    obj.SourceFile = matFile;
                    return
                end
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
            if any(strcmp(names, 'graspilot_fast'))         % a fast-open copy itself
                obj = GRASPilotData.fromFastCopy(matFile, p.Progress);
                return
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

        function f = fastCopyName(matFile)
            % default name of the fast-open copy of a .mat reconstruction
            [d, b] = fileparts(matFile);
            f = fullfile(d, [b '_graspilot.mat']);
        end

        function fc = fastCopyFor(matFile, p)
            % the fast-open copy next to matFile if it is current: same source size
            % and date, same variables, crop factor and normalization; else ''
            fc = '';
            f = GRASPilotData.fastCopyName(matFile);
            if exist(f, 'file') ~= 2, return; end
            try
                S = load(f, 'graspilot_fast');
                g = S.graspilot_fast;
                src = dir(matFile);
                cf = p.CropFactor; if isempty(cf), cf = 2; end
                nm = p.Normalize; if isempty(nm), nm = true; end
                if g.SourceBytes == src.bytes && abs(g.SourceDate - src.datenum) < 1e-6 && ...
                        strcmp(g.CoefVar, p.CoefVar) && strcmp(g.BasisVar, p.BasisVar) && ...
                        g.CropFactor == cf && logical(g.Normalized) == logical(nm)
                    fc = f;
                end
            catch
            end
        end

        function obj = fromFastCopy(file, progress)
            % read a fast-open copy (GRASPilotData.saveFastCopy): cropped, normalized
            % single-precision coefficients and basis, no preparation needed
            if nargin < 2, progress = []; end
            GRASPilotData.report(progress, 0, 'Reading the fast-open copy...');
            S = load(file, 'graspilot_fast', 'C', 'Phi');
            g = S.graspilot_fast;
            obj = GRASPilotData();
            obj.Mode = 'subspace';
            obj.C = S.C;
            obj.Phi = S.Phi;
            n = size(S.C); n(end + 1:4) = 1;
            obj.NativeSize = n([1 2 4]);
            obj.NTRaw = size(S.Phi, 1);
            obj.NK = n(3);
            obj.setBinning(1);
            obj.ScaleFactor = g.ScaleFactor;
            obj.CropFactor = g.CropFactor;
            obj.Normalized = logical(g.Normalized);
            obj.setOrientation('axial', false, [false false false]);
            obj.SourceFile = g.SourceFile;
            obj.SourceVars = sprintf('%s + %s (fast-open copy)', g.CoefVar, g.BasisVar);
            obj.FastCopyFile = file;
            GRASPilotData.report(progress, 1, 'Fast-open copy read');
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
            % Global maximum of |B Phi'| over all slices and frames (ReadData.m
            % normalization), found exactly without evaluating every voxel: for a
            % voxel with coefficients b, |b Phi(t,:)'| <= ||b|| max_t ||Phi(t,:)||.
            % Pass 1 evaluates the voxels with the largest bounds of each slice;
            % pass 2 only those whose bound still exceeds the best value found.
            ub = cell(1, nz);
            for z = 1:nz
                Bz = reshape(single(vk_all(z, :, :, :)), P, K);
                if p.Normalize
                    PzH = Phi(:, :, z)';                  % [K T], conjugate transpose (TempPCASub)
                    cz = max(sqrt(sum(abs(double(Phi(:, :, z))) .^ 2, 2)));
                    ub{z} = sqrt(sum(abs(double(Bz)) .^ 2, 2)) * cz;
                    [~, o] = sort(ub{z}, 'descend');
                    top = o(1:min(P, 2048));
                    m = max(max(abs(Bz(top, :) * PzH)));
                    if m > gmax, gmax = m; end
                end
                Bz = reshape(Bz, nxF, nyF, K);
                Cc(:, :, :, z) = Bz(rr, cc, :);
                GRASPilotData.report(p.Progress, 0.7 * z / nz, sprintf('Preparing slice %d/%d', z, nz));
            end
            if p.Normalize
                for z = 1:nz
                    cand = find(ub{z} > gmax);
                    if isempty(cand), continue; end
                    Bz = reshape(single(vk_all(z, :, :, :)), P, K);
                    PzH = Phi(:, :, z)';
                    for b0 = 1:blk:numel(cand)
                        rows = cand(b0:min(numel(cand), b0 + blk - 1));
                        m = max(max(abs(Bz(rows, :) * PzH)));
                        if m > gmax, gmax = m; end
                    end
                    GRASPilotData.report(p.Progress, 0.7 + 0.3 * z / nz, sprintf('Normalization: slice %d/%d', z, nz));
                end
            end
            ub = []; %#ok<NASGU>
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

        function a = aifStruct(a)
            % AIF inputs as a struct (a numeric curve is the reference ROI signal)
            f = {'Ref', 'Art', 'Ring', 'Phase'};
            if ~isstruct(a), a = struct('Ref', a); end
            for i = 1:numel(f)
                if ~isfield(a, f{i}), a.(f{i}) = []; end
                a.(f{i}) = double(a.(f{i})(:)');
            end
            a = orderfields(a, f);
        end

        function tf = validCurve(a, ke)
            tf = ~isempty(a) && numel(a) >= ke && all(isfinite(a(1:ke)));
        end

        function q = quantile(v, p)
            v = sort(v(:));
            q = v(max(1, min(numel(v), round(p * numel(v)))));
        end

        function th = otsu(vals)
            % Otsu threshold of a vector
            vals = double(vals(isfinite(vals)));
            lo = min(vals); hi = max(vals);
            th = lo;
            if isempty(vals) || hi <= lo, return; end
            nb = 256;
            idx = min(nb, floor((vals - lo) / (hi - lo) * nb) + 1);
            cnt = accumarray(idx(:), 1, [nb 1]);
            pp = cnt / sum(cnt);
            ctr = lo + ((1:nb)' - 0.5) * (hi - lo) / nb;
            om = cumsum(pp);
            mu = cumsum(pp .* ctr);
            sb = (mu(end) * om - mu) .^ 2 ./ (om .* (1 - om));
            sb(~isfinite(sb)) = 0;
            [~, k] = max(sb);
            th = lo + k * (hi - lo) / nb;
        end

        function h = maskHash(M)
            % cheap identity of a logical mask (for the PK map key)
            if isempty(M), h = []; return; end
            lin = find(M);
            h = [numel(M), numel(lin), sum(mod(lin, 104729)), sum(lin .^ 2) / numel(M) ^ 2];
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
            obj.PKMaps = [];
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
                'SER', 'WashoutPct', 'Plateau', 'CurveType', 'WashInSSS', 'WashInMax', 'WashInOnset'};
            src = names; src{1} = 'S0'; src{2} = 'Smax';
            for i = 1:numel(names)
                M.(names{i}) = zeros(nx, ny, nz, 'single');
            end
            P = nx * ny;
            blk = 16384;
            for z = 1:nz
                S = obj.sliceTC(z);
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

        function [bf, info] = detectBaseline(obj)
            % Automatic baseline (pre-contrast) frames 1..be: bolus arrival on
            % the mean curve of the earliest-arriving tenth of the strongly enhancing voxels (arteries,
            % which enhance first). A constant-then-linear-rise model is fitted
            % to that curve up to its peak (bolus arrival time, e.g. Cheong et
            % al. 2003); the baseline ends before the knee and before any frame
            % already 5% of the rise above the baseline level (temporal
            % smoothing of the reconstruction). Falls back to frames
            % 1-min(5, T) without a clear bolus. info: Arrival (frame),
            % ArrivalTime, Curve, NVoxels, Message.
            T = obj.NT;
            bf = 1:max(1, min(5, T));
            info = struct('Arrival', NaN, 'ArrivalTime', NaN, 'Curve', [], 'NVoxels', 0, 'Message', '');
            if T < 4
                info.Message = 'Too few frames for automatic baseline detection: frames 1-5 (or all) used.';
                return
            end
            N = prod(obj.NativeSize);
            lin = (1:max(1, floor(N / 60000)):N)';
            TC = double(obj.nativeTC(lin));
            m = mean(TC, 2);
            head = m > GRASPilotData.otsu(m);
            E0 = max(TC, [], 2) - mean(TC(:, 1:2), 2);
            E0(~head) = -Inf;
            nh = nnz(head);
            if nh < 50
                info.Message = 'No head found for automatic baseline detection: frames 1-5 used.';
                return
            end
            % clearly enhancing voxels (vessels, enhancing tissue): at least 20% of the
            % 99.9th percentile of the enhancement in the head, at most the top 3%
            Eh = sort(E0(head), 'descend');
            ref = Eh(max(1, round(0.001 * nh)));
            nt = max(10, min([nh, round(0.03 * nh), nnz(E0 >= 0.2 * ref)]));
            [~, o] = sort(E0, 'descend');
            top = o(1:nt);
            % earliest-arriving tenth of them (arteries before veins and tissue):
            % first frame at 20% of each voxel's rise
            Tt = TC(top, :);
            b1 = mean(Tt(:, 1:2), 2);
            nrm = (Tt - repmat(b1, 1, T)) ./ repmat(max(Tt, [], 2) - b1, 1, T);
            [~, arr] = max(nrm >= 0.2, [], 2);
            early = top(arr <= max(min(arr), GRASPilotData.quantile(arr, 0.1)));
            if numel(early) < 10
                [~, oa] = sort(arr);
                early = top(oa(1:min(10, numel(oa))));
            end
            c = mean(TC(early, :), 1);
            nt = numel(early);
            info.Curve = c; info.NVoxels = nt;
            t = obj.paramTimes();
            [cmax, tp] = max(c);
            if tp < 3 || ~(cmax > c(1))
                info.Message = 'No bolus arrival found (peak in the first frames): frames 1-5 used.';
                return
            end
            sse = Inf(1, tp - 1);
            for k = 1:tp - 1
                X = [ones(tp, 1), max(0, t(1:tp)' - t(k))];
                beta = X \ c(1:tp)';
                if beta(2) > 0, sse(k) = sum((c(1:tp)' - X * beta) .^ 2); end
            end
            [~, k0] = min(sse);
            b0 = median(c(1:k0));
            sg = 0;
            if k0 >= 3, sg = std(c(1:k0 - 1)); end
            thr = b0 + max(0.03 * (cmax - b0), 3 * sg);
            be = k0;
            while be > 1 && c(be) > thr, be = be - 1; end
            bf = 1:be;
            info.Arrival = k0 + 1;
            info.ArrivalTime = t(min(T, k0 + 1));
            unit = 'frame'; if obj.hasTiming(), unit = 's'; end
            info.Message = sprintf(['Baseline frames 1-%d detected automatically: contrast arrives after frame %d ' ...
                '(%.4g %s) on the mean curve of the %d earliest of the most enhancing voxels.'], be, k0, t(k0), unit, nt);
        end

        function file = saveFastCopy(obj, file)
            % Save a fast-open copy of a .mat reconstruction: the cropped, normalized
            % coefficients and the basis in single precision, uncompressed (~1/8 of
            % a double-precision source with crop factor 2). Opening the source .mat
            % then reads it instead while it is current. Default name: <source>_graspilot.mat
            if ~strcmp(obj.SourceType, 'mat') || exist(obj.SourceFile, 'file') ~= 2
                error('GRASPilotData:fast', 'A fast-open copy needs data opened from a .mat file.');
            end
            if obj.MotionCorrected || ~strcmp(obj.Mode, 'subspace')
                error('GRASPilotData:fast', ['A fast-open copy holds the subspace coefficients: revert the ' ...
                    'motion correction first.']);
            end
            if nargin < 2 || isempty(file), file = GRASPilotData.fastCopyName(obj.SourceFile); end
            src = dir(obj.SourceFile);
            vars = regexp(obj.SourceVars, '^(\S+) \+ (\S+)', 'tokens', 'once');
            if isempty(vars), vars = {'vk_all', 'basis_all'}; end
            v = GRASPilot.versionInfo();
            w = whos('-file', obj.SourceFile);              % listing for the load dialog
            sw = struct('name', {w.name}, 'size', {w.size}, 'bytes', {w.bytes}, 'class', {w.class}, ...
                'complex', {w.complex});
            graspilot_fast = struct('Format', 'GRASPilot fast-open copy', 'Version', v.Text, 'SourceWhos', {sw}, ...
                'SourceFile', obj.SourceFile, 'SourceBytes', src.bytes, 'SourceDate', src.datenum, ...
                'CoefVar', vars{1}, 'BasisVar', vars{2}, 'CropFactor', obj.CropFactor, ...
                'Normalized', obj.Normalized, 'ScaleFactor', obj.ScaleFactor); %#ok<NASGU>
            C = obj.C; Phi = obj.Phi; %#ok<NASGU>
            if exist('OCTAVE_VERSION', 'builtin')
                save('-v6', file, 'graspilot_fast', 'C', 'Phi');
            else
                try
                    save(file, 'graspilot_fast', 'C', 'Phi', '-v7.3', '-nocompression');
                catch
                    save(file, 'graspilot_fast', 'C', 'Phi', '-v7.3');
                end
            end
        end

        function S = sliceTC(obj, z)
            % time courses of every voxel of native slice z: [nx*ny x NT] (binned series)
            P = obj.NativeSize(1) * obj.NativeSize(2);
            if strcmp(obj.Mode, 'subspace')
                S = abs(reshape(obj.C(:, :, :, z), P, obj.NK) * obj.Phi(:, :, z)');
            else
                S = reshape(obj.D(:, :, z, :), P, obj.NTRaw);
            end
            if obj.BinSize > 1 || obj.NT ~= obj.NTRaw
                S = single(double(S) * obj.BinW);
            end
        end

        % ---------- pharmacokinetic models ----------
        function setPKOptions(obj, varargin)
            % 'Models' (cell of 'tofts' | 'etofts' | '2cxm'), 'T10' (ms, global
            % tissue T10), 'T10Blood' (ms, AIF), 'TR' (ms), 'FlipAngle' (deg),
            % 'r1' (1/(mM s)), 'Hct', 'AIF' ('ref' | 'population'), 'AIFShift'
            % (s; positive = later), 'Region' ('head' | 'rois' | 'all'), 'Refine'
            p = GRASPilotData.parseOpts(varargin, obj.PKOpts);
            if ischar(p.Models), p.Models = {p.Models}; end
            p.Models = lower(p.Models(:)');
            allm = {'tofts', 'etofts', '2cxm', 'ratio'};
            bad = setdiff(p.Models, allm);
            if ~isempty(bad), error('GRASPilotData:pk', 'Unknown model "%s".', bad{1}); end
            p.Models = intersect(allm, p.Models, 'stable');
            pos = {'T10', 'T10Blood', 'TR', 'FlipAngle', 'r1'};
            for i = 1:numel(pos)
                v = p.(pos{i});
                if ~isnumeric(v) || ~isscalar(v) || ~(v > 0) || ~isfinite(v)
                    error('GRASPilotData:pk', '%s must be a positive number.', pos{i});
                end
            end
            if ~(p.FlipAngle < 180), error('GRASPilotData:pk', 'FlipAngle must be below 180 degrees.'); end
            if ~(p.Hct >= 0 && p.Hct < 1), error('GRASPilotData:pk', 'Hct must be in [0, 1).'); end
            if isempty(p.AIFShift) || ~isfinite(p.AIFShift), p.AIFShift = 0; end
            p.AIF = lower(p.AIF);
            if ~any(strcmp(p.AIF, {'ref', 'population', 'vcai', 'artery', 'auto'}))
                error('GRASPilotData:pk', 'AIF must be ''ref'', ''population'', ''vcai'', ''artery'' or ''auto''.');
            end
            p.VCAIAnchor = lower(p.VCAIAnchor);
            if ~any(strcmp(p.VCAIAnchor, {'phase', 'magnitude'}))
                error('GRASPilotData:pk', 'VCAIAnchor must be ''phase'' or ''magnitude''.');
            end
            p.VCAIRing = logical(p.VCAIRing);
            p.RingWidth = max(1, round(p.RingWidth));
            p.RatioDomain = lower(p.RatioDomain);
            if ~any(strcmp(p.RatioDomain, {'basis', 'frames'}))
                error('GRASPilotData:pk', 'RatioDomain must be ''basis'' or ''frames''.');
            end
            p.Region = lower(p.Region);
            if ~any(strcmp(p.Region, {'head', 'rois', 'all'}))
                error('GRASPilotData:pk', 'Region must be ''head'', ''rois'' or ''all''.');
            end
            p.Refine = logical(p.Refine);
            obj.PKOpts = p;
        end

        function k = pkKey(obj, baseFrames, aif, regionHash)
            % everything the PK maps depend on
            aif = GRASPilotData.aifStruct(aif);
            switch obj.PKOpts.AIF
                case 'population', aif = [];
                case 'ref', aif = aif.Ref;
                case {'artery', 'auto'}, aif = aif.Art;
            end
            k = {obj.PKOpts, baseFrames, obj.dt, obj.t0, obj.BinSize, obj.NT, obj.RawTimes, ...
                obj.analysisEnd(), aif, regionHash};
        end

        function tf = hasPK(obj, baseFrames, aif, regionHash)
            tf = ~isempty(obj.PKMaps) && isequaln(obj.PKMaps.Key, obj.pkKey(baseFrames, aif, regionHash));
        end

        function [Cp, t, Cb, info] = pkAIF(obj, baseFrames, aif)
            % plasma AIF (mM) on the analysed frames 1..analysisEnd and their
            % times t (s). aif: mean signal of the reference (SSS) ROI, or a
            % struct with fields Ref (SSS mean signal), Art (arterial ROI mean
            % signal), Ring (neighbouring tissue ring mean signal) and Phase
            % (SSS phase difference from baseline, rad), all over every frame.
            %   'ref'        Ref converted with the blood T10, Cp = Cb/(1 - Hct)
            %   'population' Parker 2006, injection at the last baseline frame
            %   'vcai'       venous-calibrated arterial input (graspilot_vcai):
            %                Art (and Ring) calibrated against the SSS phase
            %                (or magnitude) curve
            %   'artery', 'auto'  Art (arterial ROI, drawn or found by
            %                graspilot_detect_aif) converted with the blood T10
            % Shifted by AIFShift seconds (positive = later). info: VCAI fit.
            if ~obj.hasTiming()
                error('GRASPilotData:pk', ['Pharmacokinetic models need the frame times: enter the ' ...
                    'temporal resolution (s/frame).']);
            end
            aif = GRASPilotData.aifStruct(aif);
            o = obj.PKOpts;
            ke = obj.analysisEnd();
            tt = obj.paramTimes();
            t = tt(1:ke);
            bf = baseFrames(baseFrames >= 1 & baseFrames <= ke);
            if isempty(bf), bf = 1; end
            info = [];
            sec = @(x) x / 1000;
            switch o.AIF
                case 'population'
                    Cb = graspilot_pop_aif(t - t(max(bf)) - o.AIFShift);
                case 'ref'
                    Cb = obj.bloodConc(aif.Ref, bf, ke, 'reference ROI');
                case {'artery', 'auto'}
                    Cb = obj.bloodConc(aif.Art, bf, ke, 'arterial ROI');
                case 'vcai'
                    if ~GRASPilotData.validCurve(aif.Art, ke)
                        error('GRASPilotData:pk', ['The VCAI needs an arterial ROI (MCA or distal ICA; a ' ...
                            'name containing "artery", "AIF", "MCA" or "ICA").']);
                    end
                    ring = [];
                    if o.VCAIRing
                        if ~GRASPilotData.validCurve(aif.Ring, ke)
                            error('GRASPilotData:pk', 'The VCAI ring of neighbouring tissue is empty.');
                        end
                        ring = aif.Ring(1:ke);
                    end
                    anchor = 'phase';
                    if strcmp(o.VCAIAnchor, 'phase') && GRASPilotData.validCurve(aif.Phase, ke)
                        s = aif.Phase(1:ke);
                    else
                        anchor = 'magnitude';
                        s = obj.bloodConc(aif.Ref, bf, ke, 'reference (SSS) ROI');
                    end
                    vo = struct('T10Blood', sec(o.T10Blood), 'TR', sec(o.TR), 'FlipAngle', o.FlipAngle, 'r1', o.r1);
                    ck = {double(aif.Art(1:ke)), double(ring), double(s), t, bf, vo};
                    if ~isempty(obj.AIFCache) && isequaln(obj.AIFCache.Key, ck)
                        Cb = obj.AIFCache.Cb; info = obj.AIFCache.Info;
                    else
                        [Cb, info] = graspilot_vcai(aif.Art(1:ke), ring, s, t, bf, vo);
                        info.Anchor = anchor;
                        obj.AIFCache = struct('Key', {ck}, 'Cb', Cb, 'Info', info);
                    end
                    if ~info.Ok
                        error('GRASPilotData:pk', ['The VCAI fit failed: the arterial signal cannot be ' ...
                            'converted for any blood fraction (check TR, flip angle and blood T10).']);
                    end
            end
            if o.AIFShift ~= 0 && ~strcmp(o.AIF, 'population')
                tq = min(max(t - o.AIFShift, t(1)), t(end));
                Cb = interp1(t, Cb, tq, 'linear');
            end
            Cp = Cb / (1 - o.Hct);
        end

        function Cb = bloodConc(obj, a, bf, ke, what)
            % blood concentration of a mean vessel signal (blood T10)
            o = obj.PKOpts;
            if ~GRASPilotData.validCurve(a, ke)
                if strcmp(what, 'arterial ROI')
                    error('GRASPilotData:pk', ['The AIF needs an arterial ROI: press Auto AIF, or + Artery ' ...
                        'and draw a few voxels in the MCA or distal ICA.']);
                end
                error('GRASPilotData:pk', ['The AIF needs a reference ROI in the superior sagittal ' ...
                    'sinus (Auto SSS) or the population AIF.']);
            end
            a = double(a(:)');
            [Cb, valid] = graspilot_signal2conc(a(1:ke), mean(a(bf)), o.T10Blood / 1000, o.TR / 1000, ...
                o.FlipAngle, o.r1);
            if ~all(valid)
                error('GRASPilotData:pk', ['The %s signal rises %.3g-fold, more than the spoiled ' ...
                    'gradient echo equation allows for TR %g ms, flip angle %g deg and blood T10 %g ms: ' ...
                    'check these values.'], what, max(a(1:ke)) / mean(a(bf)), o.TR, o.FlipAngle, o.T10Blood);
            end
        end

        function tf = hasPhase(obj)
            % complex subspace data: phase time courses available (sinus phase for the VCAI)
            tf = strcmp(obj.Mode, 'subspace') && (~isreal(obj.C) || ~isreal(obj.Phi));
        end

        function ph = phaseCurve(obj, canonLin, baseFrames)
            % mean phase difference from the baseline (rad) of the voxels canonLin
            % [1 x NT]; [] without complex data. Complex values are averaged
            % over the raw frames of each binned frame.
            ph = [];
            if ~obj.hasPhase() || isempty(canonLin), return; end
            lin = obj.canon2nativeLin(canonLin(:));
            [Cv, z] = obj.gatherCoefs(lin);
            Z = complex(zeros(numel(lin), obj.NTRaw));
            [zu, ~, g] = unique(z);
            for ii = 1:numel(zu)
                m = (g == ii);
                Z(m, :) = double(Cv(m, :)) * double(obj.Phi(:, :, zu(ii)))';
            end
            Z = Z * obj.BinW;
            bf = baseFrames(baseFrames >= 1 & baseFrames <= obj.NT);
            if isempty(bf), bf = 1; end
            Z0 = mean(Z(:, bf), 2);
            ph = mean(angle(Z .* repmat(conj(Z0), 1, obj.NT)), 1);
        end

        function Q = pkBasis(obj, slices)
            % orthonormal real temporal basis [analysed frames x n] spanned by the
            % (binned) subspace basis of the native slices: RATIO fits the
            % coefficients on it. [] for dense data or RatioDomain 'frames'.
            Q = [];
            if ~strcmp(obj.Mode, 'subspace') || ~strcmp(obj.PKOpts.RatioDomain, 'basis'), return; end
            ke = obj.analysisEnd();
            zs = unique(slices(:)');
            X = zeros(ke, 0);
            for z = zs
                Pb = obj.BinW' * double(obj.Phi(:, :, z));
                Pb = Pb(1:ke, :);
                X = [X, real(Pb), imag(Pb)]; %#ok<AGROW>
            end
            X = X(:, any(X ~= 0, 1));
            if isempty(X), return; end
            [U, S] = svd(X, 'econ');
            sv = diag(S);
            Q = U(:, sv > 1e-6 * max(sv));
        end

        function C = pkConc(obj, S, baseFrames)
            % tissue concentration (mM) of signal curves S [N x >= analysisEnd]
            % on the analysed frames (global tissue T10, signal-ratio method)
            o = obj.PKOpts;
            ke = obj.analysisEnd();
            bf = baseFrames(baseFrames >= 1 & baseFrames <= ke);
            if isempty(bf), bf = 1; end
            S = double(S(:, 1:ke));
            C = graspilot_signal2conc(S, mean(S(:, bf), 2), o.T10 / 1000, o.TR / 1000, o.FlipAngle, o.r1);
        end

        function [R, Fit, C, Cp, t] = fitPKCurves(obj, S, baseFrames, aif, canonLin)
            % PK model fits of signal curves S [N x NT] (e.g. ROI mean curves).
            % canonLin: voxels of the curves (temporal basis for RATIO)
            if nargin < 5, canonLin = []; end
            [Cp, t] = obj.pkAIF(baseFrames, aif);
            C = obj.pkConc(S, baseFrames);
            o = obj.PKOpts;
            R = struct(); Fit = struct();
            m = setdiff(o.Models, {'ratio'}, 'stable');
            if ~isempty(m)
                [R, Fit] = graspilot_pk_fit(C, t, Cp, m, struct('Refine', o.Refine));
            end
            if any(strcmp(o.Models, 'ratio'))
                Q = [];
                if ~isempty(canonLin)
                    [~, ~, z] = ind2sub(obj.NativeSize, obj.canon2nativeLin(canonLin(:)));
                    Q = obj.pkBasis(z);
                end
                [R.ratio, Fit.ratio] = graspilot_ratio_fit(C, t, Cp, obj.ratioOpts(Q));
            end
        end

        function ro = ratioOpts(obj, Q)
            o = obj.PKOpts;
            ro = struct('Basis', Q, 'Refine', o.Refine, 'EThresh', o.RatioEThresh, 'KFast', o.RatioKFast, ...
                'KSlow', o.RatioKSlow, 'RhoMax', o.RatioRhoMax, 'DelayRange', o.RatioDelayRange);
        end

        function computePK(obj, baseFrames, aif, region, progress)
            % Voxel-wise Tofts / extended Tofts / 2CXM maps. region: canonical
            % logical mask of the voxels to fit ([] = all voxels)
            if nargin < 5, progress = []; end
            o = obj.PKOpts;
            if isempty(o.Models), error('GRASPilotData:pk', 'Select at least one model.'); end
            [Cp, t, ~, vinfo] = obj.pkAIF(baseFrames, aif);
            nx = obj.NativeSize(1); ny = obj.NativeSize(2); nz = obj.NativeSize(3);
            if isempty(region)
                nm = true(nx, ny, nz);
            else
                nm = obj.canon2native(logical(region));
            end
            keys = obj.PKTable(ismember(obj.PKTable(:, 2), o.Models), :);
            M = struct();
            M.Key = obj.pkKey(baseFrames, aif, GRASPilotData.maskHash(region));
            M.Keys = keys(:, 1)';
            M.Lin = uint32(find(nm));                 % fitted voxels (native), values per key in columns
            M.Values = NaN(numel(M.Lin), size(keys, 1), 'single');
            ntot = max(1, nnz(nm)); ndone = 0;
            blk = 8192;
            fo = struct('Refine', o.Refine);
            mConv = setdiff(o.Models, {'ratio'}, 'stable');
            doRatio = any(strcmp(o.Models, 'ratio'));
            nObs = NaN;
            tic0 = tic;
            for z = 1:nz
                rows = find(reshape(nm(:, :, z), [], 1));
                if isempty(rows), continue; end
                S = obj.sliceTC(z);
                S = S(rows, :);
                for b0 = 1:blk:numel(rows)
                    q = b0:min(numel(rows), b0 + blk - 1);
                    Cq = obj.pkConc(S(q, :), baseFrames);
                    R = struct();
                    if ~isempty(mConv), R = graspilot_pk_fit(Cq, t, Cp, mConv, fo); end
                    if doRatio
                        if b0 == 1, Qz = obj.pkBasis(z); end
                        R.ratio = graspilot_ratio_fit(Cq, t, Cp, obj.ratioOpts(Qz));
                        nObs = R.ratio.NObs(1);
                    end
                    for i = 1:size(keys, 1)
                        fld = keys{i, 2}; if strcmp(fld, '2cxm'), fld = 'cxm'; end
                        M.Values(ndone + (1:numel(q)), i) = R.(fld).(keys{i, 3});
                    end
                    ndone = ndone + numel(q);
                    GRASPilotData.report(progress, ndone / ntot, sprintf(['Pharmacokinetic models: slice %d/%d, ' ...
                        '%d of %d voxels'], z, nz, ndone, ntot));
                end
            end
            M.Info = struct('t', t, 'Cp', Cp, 'NVoxels', nnz(nm), 'Seconds', toc(tic0), ...
                'BaseFrames', baseFrames, 'VCAI', vinfo, 'RatioNObs', nObs);
            obj.PKMaps = M;
            fn = fieldnames(obj.CanonMapCache);
            for i = 1:numel(fn)
                if strncmp(fn{i}, 'PK_', 3), obj.CanonMapCache = rmfield(obj.CanonMapCache, fn{i}); end
            end
            obj.MapsVersion = obj.MapsVersion + 1;
        end

        function A = canonPK(obj, key)
            % canonical map of a PK parameter (key from PKTable); [] if not computed
            A = [];
            if isempty(obj.PKMaps), return; end
            k = find(strcmp(obj.PKMaps.Keys, key), 1);
            if isempty(k), return; end
            ck = ['PK_' key];
            if isfield(obj.CanonMapCache, ck)
                A = obj.CanonMapCache.(ck);
                return
            end
            A = NaN(obj.NativeSize, 'single');
            A(obj.PKMaps.Lin) = obj.PKMaps.Values(:, k);
            A = obj.native2canon(A);
            obj.CanonMapCache.(ck) = A;
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
            obj.PKMaps = [];
            obj.CanonMapCache = struct();
            obj.MapsVersion = obj.MapsVersion + 1;
        end

        function A = canonMap(obj, name)
            % name: Baseline | Max | PeakEnh | RelEnh | TTP | Onset | WashIn | WashOut | AUC
            % | WashInSSS | WashInMax | WashInOnset (each wash-in method) | ...
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
                        'SER', 'WashoutPct', 'Plateau', 'CurveType', 'WashInSSS', 'WashInMax', 'WashInOnset'}
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
