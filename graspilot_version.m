function v = graspilot_version()
%GRASPILOT_VERSION Version of GRASPilot and the exact revision of a git checkout.
%
%   v = graspilot_version()      (same as GRASPilot.versionInfo)
%
%   v.Version   release number, e.g. '0.1.0': the latest GitHub release this
%               code is based on (tag v<Version>); set in GRASPilot.m
%               (versionInfo) when tagging a new release
%   v.Revision  'git describe' of the checkout, e.g. 'v0.1.0' (the release),
%               'v0.1.0-12-gab9e173' (12 commits after it), '-dirty' with
%               local changes; '' when not run from a git checkout
%   v.URL       GitHub repository
%   v.Text      one line for About, exports and sessions
%   v.Details   human-readable description of the revision
v = GRASPilot.versionInfo();
end
