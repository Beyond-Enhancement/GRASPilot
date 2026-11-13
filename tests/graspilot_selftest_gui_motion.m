function nfail = graspilot_selftest_gui_motion(out)
%GRASPILOT_SELFTEST_GUI_MOTION Motion window, warning, correction and session on a moving phantom.
nfail = 0;
[D, tr] = graspilot_motion_phantom([1 1.5 -1 1.5 1 -2], 9:12);
D.SourceFile = 'motion_phantom.mat';
app = GRASPilot(D);
app.Quiet = true;
cleanup = onCleanup(@() app.close(true));
n = D.canonSize();
% ROI across the posterior SSS, skull and brain: its curves change when the head moves
M = false(n); M(46:54, 95:104, 36:44) = true;
app.addTumor();
app.beginAction(); app.setVoxels(1, find(M), true); app.commitAction('box'); app.afterMaskChange();
lin = find(M);
app.openMotionWindow();
okW0 = ishghandle(app.MW.Fig) && strcmp(get(app.MW.btnApply, 'Enable'), 'off');
Mo = app.estimateMotion();
info = get(app.H.info, 'String');
if iscell(info), info = strjoin(info, ' '); end
okE = ~isempty(Mo) && Mo.Excess && isequal(Mo.ExcessFrames, 9:12) && ~isempty(strfind(info, 'MOTION WARNING')) && ...
    ~isempty(strfind(get(app.MW.info, 'String'), 'WARNING')) && numel(get(app.MW.Ax(3), 'Children')) >= 4 && ...
    strcmp(get(app.MW.btnApply, 'Enable'), 'on');
nfail = nfail + check('motion window: estimate, warning, plots', okW0 && okE, numel(Mo.ExcessFrames));
if nargin > 0 && ~isempty(out)
    try, print(app.MW.Fig, '-dpng', '-r60', fullfile(out, 'motion_window.png')); catch, end
end
% ROI curve against the motion-free truth, before and after correction
app.computeStats(1);
still = double(tr.Still.timecourses(lin));
before = double(D.timecourses(lin));
app.applyMotion();
after = double(D.timecourses(lin));
d = before(:, 9:12) - still(:, 9:12); eb = sqrt(mean(d(:) .^ 2));
d = after(:, 9:12) - still(:, 9:12); ea = sqrt(mean(d(:) .^ 2));
app.computeStats(1);
okC = D.MotionCorrected && ea < 0.5 * eb && app.ROIs(1).Dirty == false && ...
    max(abs(app.ROIs(1).Stats.mean(:)' - mean(after, 1))) < 1e-5 && strcmp(get(app.MW.btnRevert, 'Enable'), 'on');
nfail = nfail + check('Apply correction: ROI voxel curves match the motion-free series', okC, ea / eb);
% exports carry the motion state
f1 = fullfile(out, 'motion.csv');
app.exportMotion(f1);
fid = fopen(f1); L = {}; while true, l = fgetl(fid); if ~ischar(l), break; end, L{end + 1} = l; end, fclose(fid); %#ok<AGROW>
fn = fullfile(out, 'mc_dce.nii.gz');
app.exportImageNifti(fn, 1:D.NT);
js = fileread(regexprep(fn, '\.nii\.gz$', '.json'));
okX = numel(L) == 1 + D.NTRaw + 3 && ~isempty(strfind(L{end}, 'applied')) && ~isempty(strfind(js, '"MotionCorrected":true'));
nfail = nfail + check('motion csv and NIfTI sidecar record the correction', okX, numel(L));
% session: reloading on uncorrected data re-applies the same correction
fs = fullfile(out, 'motion_session.mat');
app.saveROIs(fs);
D2 = graspilot_motion_phantom([1 1.5 -1 1.5 1 -2], 9:12);
app2 = GRASPilot(D2); app2.Quiet = true;
c2 = onCleanup(@() app2.close(true));
app2.loadROIs(fs);
okS = D2.MotionCorrected && isequal(D2.D, D.D) && isequal(find(app2.Masks{1}), lin);
nfail = nfail + check('ROI session re-applies the motion correction', okS, 0);
% residual check on the corrected frames, then revert
Mr = app.estimateMotion();
okR = ~Mr.Excess && max(Mr.AbsMean) < 0.2;
app.revertMotion();
okV = ~D.MotionCorrected && isequal(double(D.timecourses(lin)), before) && strcmp(get(app.MW.btnApply, 'Enable'), 'on');
nfail = nfail + check('residual motion after correction small; revert', okR && okV, max(Mr.AbsMean));
end

function f = check(name, cond, val)
if cond
    fprintf('  PASS  %-48s (%.3g)\n', name, val);
    f = 0;
else
    fprintf('  FAIL  %-48s (%.3g)\n', name, val);
    f = 1;
end
end
