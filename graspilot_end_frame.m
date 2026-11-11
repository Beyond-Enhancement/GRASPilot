function Te = graspilot_end_frame(t, opts)
%GRASPILOT_END_FRAME Last frame analysed for the options EndFrame / EndTime.
%
%   Te = graspilot_end_frame(t, opts)   t: frame times; opts.EndFrame (frames)
%   and opts.EndTime (same axis as t); the earlier one applies, at least 2
%   frames are kept.
T = numel(t);
Te = T;
if isfield(opts, 'EndFrame') && ~isempty(opts.EndFrame) && isfinite(opts.EndFrame)
    Te = min(Te, floor(opts.EndFrame));
end
if isfield(opts, 'EndTime') && ~isempty(opts.EndTime) && isfinite(opts.EndTime)
    k = find(t <= opts.EndTime, 1, 'last');
    if isempty(k), k = 1; end
    Te = min(Te, k);
end
Te = max(min(2, T), Te);
end
