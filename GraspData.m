classdef GraspData < handle
%GRASPDATA Data model for dynamic GRASP DCE-MRI image series.
%
%   D = GraspData.fromFile(matFile, 'Name', Value, ...)
%   D = GraspData.fromSubspace(vk_all, basis_all, 'Name', Value, ...)
%   D = GraspData.fromDense(img, dimOrder, 'Name', Value, ...)
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
%   Index conventions
%     native    : [rows cols partitions] = data1(:,:,t,z) of ReadData.m
%     canonical : array used by the viewer and every export. Axis 1 runs
%                 patient Right->Left, axis 2 Anterior->Posterior, axis 3
%                 Inferior->Superior (LPS+ voxel order). setOrientation()
%                 defines the native->canonical permutation/flips.
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
        NT = 0;
        NK = 0;
        ScaleFactor = 1;            % normalization factor applied (1/max)
        CropFactor = 1;
        Normalized = false;
        Plane = 'axial';
        SwapInPlane = false;
        UserFlip = [false false false];
        Perm = [2 1 3];             % canonical axis a <- native dim Perm(a)
        Flip = [false false false]; % flip along canonical axis a
        Maps = [];                  % native-space enhancement maps
        CanonMapCache = struct();
    end

    properties
        VoxNative = [1 1 1];        % mm: [row spacing, col spacing, slice spacing]
        dt = NaN;                   % s per frame
        t0 = 0;                     % s, time of frame 1
        SourceFile = '';
        SourceVars = '';
    end

    % ------------------------------------------------------------------
    methods (Static)

        function obj = fromFile(matFile, varargin)
            p = GraspData.parseOpts(varargin, GraspData.defaultOpts());
            if exist(matFile, 'file') ~= 2
                error('GraspData:file', 'File not found: %s', matFile);
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
                    error('GraspData:format', ['Could not find %s/%s in %s. Use the load ' ...
                        'dialog or pass ''Format'',''dense'',''ImageVar'',name.'], ...
                        p.CoefVar, p.BasisVar, matFile);
                end
            end
            GraspData.report(p.Progress, 0, 'Reading .mat file (this can take a while for large files)...');
            switch fmt
                case 'subspace'
                    S = load(matFile, p.CoefVar, p.BasisVar);
                    if ~isfield(S, p.CoefVar) || ~isfield(S, p.BasisVar)
                        error('GraspData:vars', 'Variables %s and %s are required.', p.CoefVar, p.BasisVar);
                    end
                    vk = S.(p.CoefVar); bs = S.(p.BasisVar);
                    S = []; %#ok<NASGU>
                    obj = GraspData.fromSubspace(vk, bs, p);
                    obj.SourceVars = sprintf('%s + %s', p.CoefVar, p.BasisVar);
                case 'dense'
                    if isempty(p.ImageVar)
                        error('GraspData:vars', 'ImageVar must be given for dense format.');
                    end
                    S = load(matFile, p.ImageVar);
                    img = S.(p.ImageVar);
                    S = []; %#ok<NASGU>
                    obj = GraspData.fromDense(img, p.DimOrder, p);
                    obj.SourceVars = sprintf('%s (%s)', p.ImageVar, p.DimOrder);
                otherwise
                    error('GraspData:format', 'Unknown format "%s".', fmt);
            end
            obj.SourceFile = matFile;
        end

        function obj = fromSubspace(vk_all, basis_all, varargin)
            p = GraspData.parseOpts(varargin, GraspData.defaultOpts());
            if isempty(p.CropFactor), p.CropFactor = 2; end
            if isempty(p.Normalize), p.Normalize = true; end
            if ndims(vk_all) > 4
                error('GraspData:size', 'vk_all must be [slices nx ny K].');
            end
            nz = size(vk_all, 1); nxF = size(vk_all, 2); nyF = size(vk_all, 3); K = size(vk_all, 4);

            % temporal basis -> Phi [T K nz]; ReadData uses Phi_z = basis_all(z,:,:) as [T K]
            B = basis_all;
            if ndims(B) == 2 && size(B, 2) == K && ~(nz > 1 && size(B, 1) == nz && K == 1) %#ok<ISMAT>
                Phi = repmat(single(B), [1 1 nz]);        % one basis shared by all slices
            else
                if size(B, 1) ~= nz
                    error('GraspData:basis', ['basis_all has %d entries along dim 1 but vk_all ' ...
                        'has %d slices.'], size(B, 1), nz);
                end
                if size(B, 3) == K
                    Phi = permute(single(B), [2 3 1]);    % [nz T K] -> [T K nz] (ReadData layout)
                elseif size(B, 2) == K
                    Phi = permute(single(B), [3 2 1]);    % [nz K T] -> [T K nz]
                else
                    error('GraspData:basis', ['basis_all [%s] is not compatible with K = %d ' ...
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
                GraspData.report(p.Progress, z / nz, sprintf('Preparing slice %d/%d', z, nz));
            end
            sf = 1;
            if p.Normalize && gmax > 0
                sf = 1 / double(gmax);
                for z = 1:nz
                    Cc(:, :, :, z) = Cc(:, :, :, z) * single(sf);
                end
            end

            obj = GraspData();
            obj.Mode = 'subspace';
            obj.C = Cc;
            obj.Phi = Phi;
            obj.NativeSize = [nx ny nz];
            obj.NT = T;
            obj.NK = K;
            obj.ScaleFactor = sf;
            obj.CropFactor = cf;
            obj.Normalized = logical(p.Normalize);
            obj.setOrientation('axial', false, [false false false]);
        end

        function obj = fromDense(img, dimOrder, varargin)
            p = GraspData.parseOpts(varargin, GraspData.defaultOpts());
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
                error('GraspData:dimorder', 'DimOrder must be a permutation of xyzt (got "%s").', dimOrder);
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
            obj = GraspData();
            obj.Mode = 'dense';
            obj.D = img;
            obj.NativeSize = [size(img, 1) size(img, 2) size(img, 3)];
            obj.NT = size(img, 4);
            obj.NK = 0;
            obj.ScaleFactor = sf;
            obj.CropFactor = cf;
            obj.Normalized = logical(p.Normalize);
            obj.setOrientation('axial', false, [false false false]);
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
                error('GraspData:opts', 'Options must be name/value pairs.');
            end
            fn = fieldnames(p);
            for i = 1:2:numel(args)
                m = find(strcmpi(fn, args{i}), 1);
                if isempty(m)
                    error('GraspData:opts', 'Unknown option "%s".', args{i});
                end
                p.(fn{m}) = args{i + 1};
            end
        end

        function report(fcn, frac, msg)
            if ~isempty(fcn)
                try
                    fcn(frac, msg);
                catch
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
                    error('GraspData:plane', 'Plane must be axial, coronal or sagittal.');
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
            [obj.Perm, obj.Flip] = GraspData.presetOrientation(plane, swapInPlane, userFlip);
            obj.Plane = lower(plane);
            obj.SwapInPlane = logical(swapInPlane);
            obj.UserFlip = logical(userFlip(:)');
            obj.CanonMapCache = struct();
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

        function t = frameTimes(obj)
            if obj.hasTiming()
                t = obj.t0 + (0:obj.NT - 1) * obj.dt;
            else
                t = NaN(1, obj.NT);
            end
        end

        function b = memoryBytes(obj)
            if strcmp(obj.Mode, 'subspace')
                b = numel(obj.C) * 4 * (1 + ~isreal(obj.C)) + numel(obj.Phi) * 4 * (1 + ~isreal(obj.Phi));
            else
                b = numel(obj.D) * 4;
            end
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
            if strcmp(obj.Mode, 'dense')
                switch d
                    case 3, S = obj.D(:, :, idx, t);
                    case 1, S = reshape(obj.D(idx, :, :, t), ny, nz);
                    case 2, S = reshape(obj.D(:, idx, :, t), nx, nz);
                end
                return
            end
            K = obj.NK;
            w = conj(obj.Phi(t, :, :));                    % [1 K nz]
            switch d
                case 3
                    S = reshape(obj.C(:, :, :, idx), nx * ny, K) * reshape(w(1, :, idx), K, 1);
                    S = abs(reshape(S, nx, ny));
                case 1
                    X = reshape(obj.C(idx, :, :, :), ny, K, nz);
                    S = abs(reshape(sum(X .* reshape(w, 1, K, nz), 2), ny, nz));
                case 2
                    X = reshape(obj.C(:, idx, :, :), nx, K, nz);
                    S = abs(reshape(sum(X .* reshape(w, 1, K, nz), 2), nx, nz));
            end
        end

        function V = nativeVolume(obj, t)
            if strcmp(obj.Mode, 'dense')
                V = obj.D(:, :, :, t);
                return
            end
            nx = obj.NativeSize(1); ny = obj.NativeSize(2); nz = obj.NativeSize(3);
            K = obj.NK;
            V = zeros(nx, ny, nz, 'single');
            for z = 1:nz
                w = reshape(conj(obj.Phi(t, :, z)), K, 1);
                V(:, :, z) = abs(reshape(reshape(obj.C(:, :, :, z), nx * ny, K) * w, nx, ny));
            end
        end

        function TC = nativeTC(obj, lin)
            lin = lin(:);
            if strcmp(obj.Mode, 'dense')
                D2 = reshape(obj.D, [], obj.NT);
                TC = D2(lin, :);
                return
            end
            [Cv, z] = obj.gatherCoefs(lin);
            TC = zeros(numel(lin), obj.NT, 'single');
            [zu, ~, g] = unique(z);
            for ii = 1:numel(zu)
                m = (g == ii);
                TC(m, :) = abs(Cv(m, :) * obj.Phi(:, :, zu(ii))');
            end
        end

        function v = nativeValuesAtFrame(obj, lin, t)
            lin = lin(:);
            if strcmp(obj.Mode, 'dense')
                v = obj.D(lin + (t - 1) * prod(obj.NativeSize));
                return
            end
            [Cv, z] = obj.gatherCoefs(lin);
            v = zeros(numel(lin), 1, 'single');
            [zu, ~, g] = unique(z);
            for ii = 1:numel(zu)
                m = (g == ii);
                v(m) = abs(Cv(m, :) * reshape(conj(obj.Phi(t, :, zu(ii))), obj.NK, 1));
            end
        end

        % ---------- enhancement maps ----------
        function tf = hasMaps(obj, baseFrames)
            tf = ~isempty(obj.Maps) && isequal(obj.Maps.BaseFrames, baseFrames);
        end

        function computeMaps(obj, baseFrames, progress)
            if nargin < 3, progress = []; end
            T = obj.NT;
            bf = baseFrames(baseFrames >= 1 & baseFrames <= T);
            if isempty(bf), bf = 1; end
            nx = obj.NativeSize(1); ny = obj.NativeSize(2); nz = obj.NativeSize(3);
            M.BaseFrames = baseFrames;
            M.Baseline = zeros(nx, ny, nz, 'single');
            M.Max = zeros(nx, ny, nz, 'single');
            M.TTP = zeros(nx, ny, nz, 'single');
            M.AUC = zeros(nx, ny, nz, 'single');
            post = max(bf) + 1:T;
            for z = 1:nz
                if strcmp(obj.Mode, 'subspace')
                    S = abs(reshape(obj.C(:, :, :, z), nx * ny, obj.NK) * obj.Phi(:, :, z)');
                else
                    S = reshape(obj.D(:, :, z, :), nx * ny, T);
                end
                b = mean(S(:, bf), 2);
                [mx, imx] = max(S, [], 2);
                M.Baseline(:, :, z) = reshape(b, nx, ny);
                M.Max(:, :, z) = reshape(mx, nx, ny);
                M.TTP(:, :, z) = reshape(single(imx), nx, ny);
                if ~isempty(post)
                    M.AUC(:, :, z) = reshape(sum(S(:, post) - b, 2), nx, ny);
                end
                GraspData.report(progress, z / nz, sprintf('Enhancement maps: slice %d/%d', z, nz));
            end
            obj.Maps = M;
            obj.CanonMapCache = struct();
        end

        function clearMaps(obj)
            obj.Maps = [];
            obj.CanonMapCache = struct();
        end

        function A = canonMap(obj, name)
            % name: Baseline | Max | PeakEnh | RelEnh | TTP | AUC
            if isempty(obj.Maps)
                error('GraspData:maps', 'Enhancement maps have not been computed.');
            end
            if isfield(obj.CanonMapCache, name)
                A = obj.CanonMapCache.(name);
                return
            end
            M = obj.Maps;
            switch name
                case {'Baseline', 'Max', 'TTP', 'AUC'}
                    A = M.(name);
                case 'PeakEnh'
                    A = M.Max - M.Baseline;
                case 'RelEnh'
                    A = 100 * (M.Max - M.Baseline) ./ M.Baseline;
                    A(~(M.Baseline > 0)) = NaN;
                otherwise
                    error('GraspData:maps', 'Unknown map "%s".', name);
            end
            A = obj.native2canon(A);
            obj.CanonMapCache.(name) = A;
        end

        % ---------- geometry ----------
        function A = affine(obj)
            % voxel (0-based, canonical LPS order) -> RAS mm, volume centred at 0
            v = obj.canonVox(); n = obj.canonSize();
            A = [-v(1) 0 0  v(1) * (n(1) - 1) / 2;
                 0 -v(2) 0  v(2) * (n(2) - 1) / 2;
                 0 0  v(3) -v(3) * (n(3) - 1) / 2;
                 0 0 0 1];
        end

        function s = describe(obj)
            n = obj.NativeSize;
            if strcmp(obj.Mode, 'subspace')
                s = sprintf('Subspace: %dx%dx%d, K=%d, T=%d, %.0f MB', n(1), n(2), n(3), obj.NK, ...
                    obj.NT, obj.memoryBytes() / 2^20);
            else
                s = sprintf('Dense: %dx%dx%d, T=%d, %.0f MB', n(1), n(2), n(3), obj.NT, ...
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
