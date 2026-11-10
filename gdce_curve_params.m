function P = gdce_curve_params(S, t, bf, opts)
%GDCE_CURVE_PARAMS Semi-quantitative descriptors of DCE time-intensity curves.
%
%   P = gdce_curve_params(S, t, baseFrames)
%   P = gdce_curve_params(S, t, baseFrames, opts)
%
%   S           [N x T] signal curves (one per row: voxels or ROI means)
%   t           [1 x T] frame times (s), or 1:T when the timing is unknown
%   baseFrames  pre-contrast frames; the last one is taken as the end of the
%               baseline (contrast injection)
%   opts        struct, all fields optional:
%     WashIn     'onset' (default): dS / (TP - To), the line from the onset
%                to the peak (Khalifa et al. 2014); 'max': maximum up-slope
%                between consecutive frames (Li et al. 2024)
%     WashOut    'regression' (default): least-squares slope from the peak
%                to the last frame; 'endpoints': (S_last - S_peak) /
%                (t_last - t_peak). At least the last 3 frames are used, so
%                curves still rising at the end get a positive slope.
%     AUCWindow  seconds (or frames) after the end of the baseline over
%                which the AUC is integrated; Inf (default) = to the end.
%                60 or 90 give the initial AUC (IAUC60, IAUC90).
%     OnsetFrac  enhancement fraction defining the onset (default 0.1)
%
%   Output fields, each [N x 1] (NaN where undefined):
%     S0       baseline signal (mean over baseFrames)
%     Smax     peak signal after the baseline
%     PE       peak enhancement Smax - S0                      (a.u.)
%     RE       relative peak enhancement 100 (Smax - S0)/S0    (%)
%     TTP      time to peak from the end of the baseline       (s or frames)
%     Onset    onset time To from the end of the baseline: last frame before
%              the peak with enhancement < OnsetFrac x PE      (s or frames)
%     WashIn   wash-in slope                                   (a.u./s)
%     WashOut  washout slope; negative = washout               (a.u./s)
%     AUC      trapezoidal area of S - S0 from the end of the
%              baseline over AUCWindow                         (a.u. x s)
%     PeakFrame frame index of the peak

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'WashIn'), opts.WashIn = 'onset'; end
if ~isfield(opts, 'WashOut'), opts.WashOut = 'regression'; end
if ~isfield(opts, 'AUCWindow') || isempty(opts.AUCWindow), opts.AUCWindow = Inf; end
if ~isfield(opts, 'OnsetFrac'), opts.OnsetFrac = 0.1; end

S = double(S);
t = double(t(:)');
[N, T] = size(S);
bf = bf(bf >= 1 & bf <= T);
if isempty(bf), bf = 1; end
be = max(bf);
nanv = NaN(N, 1);
P = struct('S0', nanv, 'Smax', nanv, 'PE', nanv, 'RE', nanv, 'TTP', nanv, 'Onset', nanv, ...
    'WashIn', nanv, 'WashOut', nanv, 'AUC', nanv, 'PeakFrame', nanv);
S0 = mean(S(:, bf), 2);
P.S0 = S0;
if be >= T, return; end

% peak after the baseline
[Smax, k] = max(S(:, be + 1:T), [], 2);
ipk = k + be;
PE = Smax - S0;
P.Smax = Smax;
P.PE = PE;
RE = 100 * PE ./ S0;
RE(~(S0 > 0)) = NaN;
P.RE = RE;
tpk = t(ipk)';
P.TTP = tpk - t(be);
P.PeakFrame = ipk;

E = S - repmat(S0, 1, T);
J = 1:T;

% onset: last frame in [be, peak-1] with enhancement below OnsetFrac x PE
below = E < repmat(opts.OnsetFrac * PE, 1, T);
win = repmat(J, N, 1) >= be & repmat(J, N, 1) <= repmat(ipk - 1, 1, T);
io = max((below & win) .* repmat(J, N, 1), [], 2);
io(io < be) = be;
To = t(io)';
P.Onset = To - t(be);

% wash-in slope
switch lower(opts.WashIn)
    case 'max'
        dS = diff(S, 1, 2) ./ repmat(diff(t), N, 1);
        JJ = repmat(1:T - 1, N, 1);
        ok = JJ >= be & JJ <= repmat(ipk - 1, 1, T - 1);
        dS(~ok) = -Inf;
        wi = max(dS, [], 2);
        wi(~isfinite(wi)) = NaN;
        P.WashIn = wi;
    otherwise
        P.WashIn = PE ./ (tpk - To);
end

% washout slope from the peak (at least the last 3 frames) to the end
ks = min(ipk, max(be + 1, T - 2));
n = T - ks + 1;
lin = (1:N)' + (ks - 1) * N;
switch lower(opts.WashOut)
    case 'endpoints'
        P.WashOut = (S(:, T) - S(lin)) ./ (t(T) - t(ks)');
    otherwise
        CS = fliplr(cumsum(fliplr(S), 2));                 % suffix sums
        CTS = fliplr(cumsum(fliplr(S .* repmat(t, N, 1)), 2));
        ct = fliplr(cumsum(fliplr(t)));
        ct2 = fliplr(cumsum(fliplr(t .^ 2)));
        sy = CS(lin); sty = CTS(lin);
        st = ct(ks)'; st2 = ct2(ks)';
        den = n .* st2 - st .^ 2;
        wo = (n .* sty - st .* sy) ./ den;
        wo(n < 2 | den == 0) = NaN;
        P.WashOut = wo;
end
P.WashOut(n < 2) = NaN;

% area under the enhancement curve
ke = find(t <= t(be) + opts.AUCWindow, 1, 'last');
if isempty(ke) || ke <= be
    P.AUC = zeros(N, 1);
else
    Y = E(:, be:ke);
    dt = diff(t(be:ke));
    P.AUC = sum((Y(:, 1:end - 1) + Y(:, 2:end)) / 2 .* repmat(dt, N, 1), 2);
end
end
