function [w, info] = graspilot_ref_window(Sref, t, bf, opts)
%GRASPILOT_REF_WINDOW Period of most rapid signal increase of a reference (SSS) curve.
%
%   [w, info] = graspilot_ref_window(Sref, t, baseFrames, opts)
%
%   Sref   [1 x T] reference curve (mean of the superior sagittal sinus ROI)
%   w      [a b] frames: the wash-in period. The slopes between consecutive
%          frames a..b are those of the contiguous run around the steepest
%          rise of Sref (after the baseline, up to its peak) in which the
%          slope stays >= opts.RefSlopeFrac (default 0.5) of the steepest.
%          [] if the curve does not rise.
%   info   MaxSlope, MaxSlopeFrame (the slope from that frame to the next)
%
%   Used for WashIn = 'sss' in graspilot_curve_params: the wash-in slope of a
%   tissue curve is its maximum slope within the same period, following the
%   wash-in definition "period of maximally rapid increasing superior
%   sagittal sinus signal" of Lee, Berger et al., Neurosurgery 2023;92:497.
%   opts.EndFrame / opts.EndTime limit the frames as in graspilot_curve_params.

if nargin < 4 || isempty(opts), opts = struct(); end
if ~isfield(opts, 'RefSlopeFrac') || isempty(opts.RefSlopeFrac), opts.RefSlopeFrac = 0.5; end
w = []; info = struct('MaxSlope', NaN, 'MaxSlopeFrame', NaN);
Sref = double(Sref(:)');
t = double(t(:)');
Te = graspilot_end_frame(t, opts);
Sref = Sref(1:Te); t = t(1:Te);
T = numel(t);
bf = bf(bf >= 1 & bf <= T);
if isempty(bf), bf = 1; end
be = max(bf);
if be >= T, return; end
[~, k] = max(Sref(be + 1:T));
ipk = k + be;
d = diff(Sref) ./ diff(t);              % d(i): slope from frame i to i+1
idx = be:ipk - 1;
[m, j] = max(d(idx));
if ~(m > 0), return; end
i0 = idx(j);
a = i0; b = i0;
thr = opts.RefSlopeFrac * m;
while a > be && d(a - 1) >= thr, a = a - 1; end
while b < ipk - 1 && d(b + 1) >= thr, b = b + 1; end
w = [a, b + 1];
info.MaxSlope = m;
info.MaxSlopeFrame = i0;
end
