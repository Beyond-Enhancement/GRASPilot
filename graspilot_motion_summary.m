function s = graspilot_motion_summary(M)
%GRASPILOT_MOTION_SUMMARY One-paragraph summary of a motion estimate.
%
%   s = graspilot_motion_summary(M)    M from graspilot_motion_estimate
if isempty(M) || ~isfield(M, 'Params') || all(isnan(M.AbsMean))
    s = 'No motion estimate.';
    return
end
[am, fa] = max(M.AbsMean);
[rm, fr] = max(M.RelMean);
P = M.Params;
tr = max(abs(P(:, 1:3)), [], 1);
ro = max(abs(P(:, 4:6)), [], 1);
s = sprintf(['Reference: %s. Mean head displacement: max %.2f mm (frame %d) from the reference, ' ...
    'max %.2f mm (frame %d) between consecutive frames. Largest translations L %.2f, P %.2f, S %.2f mm; ' ...
    'rotations %.2f, %.2f, %.2f deg.'], M.Reference, am, fa, rm, fr, tr, ro);
if isnan(rm), s = strrep(s, sprintf('max NaN mm (frame %d) between consecutive frames', fr), 'n/a between frames'); end
if M.Excess
    ex = M.ExcessFrames;
    s = sprintf(['%s WARNING: excess motion in %d frame(s) (%s): mean displacement > %.2f mm from the ' ...
        'reference or > %.2f mm between consecutive frames.'], s, numel(ex), frameList(ex), M.WarnAbsMM, M.WarnRelMM);
else
    s = sprintf('%s No excess motion (thresholds %.2f mm from the reference, %.2f mm between frames).', ...
        s, M.WarnAbsMM, M.WarnRelMM);
end
if isfield(M, 'Complete') && ~M.Complete
    s = sprintf('%s Stopped early: %d frames registered.', s, numel(M.Frames));
end
end

function t = frameList(f)
% compact list of frame numbers: 3-7, 12
t = '';
i = 1;
while i <= numel(f)
    j = i;
    while j < numel(f) && f(j + 1) == f(j) + 1, j = j + 1; end
    if j > i, u = sprintf('%d-%d', f(i), f(j)); else, u = sprintf('%d', f(i)); end
    if isempty(t), t = u; else, t = [t ', ' u]; end %#ok<AGROW>
    i = j + 1;
end
end
