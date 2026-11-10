function nfail = gdce_selftest_gui_sss(out)
%GDCE_SELFTEST_GUI_SSS Auto SSS button and orientation check on the venous phantom.
nfail = 0;
[D, tr] = gdce_sss_phantom();
app = GraspDCEViewer(D);
app.Quiet = true;
cleanup = onCleanup(@() app.close(true));
nfail = nfail + check('orientation check passes (LPS phantom)', app.checkOrientation(), 0);
app.addTumor();
k = app.autoSSS();
lin = find(app.Masks{k});
[ref, info] = gdce_detect_sss(D, 'BaseFrames', app.BaseFrames(1):app.BaseFrames(2));
nfail = nfail + check('Auto SSS button creates "SSS (auto)" ROI', k > 0 && strcmp(app.ROIs(k).Name, 'SSS (auto)') ...
    && isequal(sort(lin), sort(ref)) && all(tr.SSS(lin)), numel(lin));
nfail = nfail + check('crosshair moved to the SSS ROI', isequal(app.Pos, round(info.Centroid)), 0);
k2 = app.autoSSS();
nfail = nfail + check('second run replaces, does not duplicate', k2 == k && numel(app.ROIs) == 2 && ...
    nnz(app.Masks{k}) == numel(ref), numel(app.ROIs));
app.undo(); app.undo();
nfail = nfail + check('Auto SSS undoable', ~any(app.Masks{k}(:)), 0);
app.computeStats(k);
app.redo();
app.computeStats(k);
st = app.ROIs(k).Stats;
[~, pk] = max(st.mean);
nfail = nfail + check('SSS curve peaks in the venous phase (frame 6)', pk == 6, pk);
app.applyGeometry('PLS', [1 1 1], 4, 0);
nfail = nfail + check('rotated view flagged by orientation check', ~app.checkOrientation(), 0);
nfail = nfail + check('ROI follows axis-code change', isequal(D.canon2native(app.Masks{k}) , tr.SSS & ...
    D.canon2native(app.Masks{k})), 0);
app.applyGeometry('LPS', [1 1 1], 4, 0);
nfail = nfail + check('ROI back in place after restoring LPS', isequal(sort(find(app.Masks{k})), sort(ref)), 0);
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
