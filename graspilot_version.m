function v = graspilot_version()
%GRASPILOT_VERSION Version of GRASPilot and the exact revision of a git checkout.
%
%   v = graspilot_version()
%
%   v.Version   release number, e.g. '0.1.0': the latest GitHub release this
%               code is based on (tag v<Version>). Update it when tagging a
%               new release.
%   v.Revision  'git describe' of the checkout this file is in, e.g.
%               'v0.1.0' (the release itself), 'v0.1.0-12-gab9e173' (12
%               commits after it, commit ab9e173), with '-dirty' when files
%               were changed locally; '' when not run from a git checkout
%               (download, compiled app) or git is not available
%   v.URL       GitHub repository
%   v.Text      one line for About, exports and sessions, e.g.
%               '0.1.0 (v0.1.0-12-gab9e173)'
%   v.Details   human-readable description of the revision

persistent cache
if ~isempty(cache), v = cache; return; end
v = struct('Version', '0.1.0', 'Revision', '', ...
    'URL', 'https://github.com/Beyond-Enhancement/GRASPilot', 'Text', '', 'Details', '');
here = fileparts(mfilename('fullpath'));
if exist(fullfile(here, '.git'), 'dir') == 7 || exist(fullfile(here, '.git'), 'file') == 2
    try
        [st, out] = system(sprintf('git -C "%s" describe --tags --dirty --always', here));
        out = strtrim(out);
        if st == 0 && ~isempty(out) && isempty(regexp(out, '\s', 'once')), v.Revision = out; end
    catch
    end
end
v.Text = v.Version;
v.Details = sprintf(['Version %s; exact revision unknown (not run from a git checkout, or git not ' ...
    'available)'], v.Version);
r = v.Revision;
if ~isempty(r)
    v.Text = sprintf('%s (%s)', v.Version, r);
    dirty = ~isempty(regexp(r, '-dirty$', 'once'));
    r = regexprep(r, '-dirty$', '');
    tk = regexp(r, '^(.*)-(\d+)-g([0-9a-f]+)$', 'tokens', 'once');
    if ~isempty(tk)
        v.Details = sprintf('%s commit(s) after release %s (commit %s)', tk{2}, tk{1}, tk{3});
    elseif ~isempty(regexp(r, '^v?\d', 'once'))
        v.Details = sprintf('Release %s', r);
    else
        v.Details = sprintf('Commit %s (no release tag found in this checkout)', r);
    end
    if dirty, v.Details = [v.Details ', with local changes']; end
end
cache = v;
end
