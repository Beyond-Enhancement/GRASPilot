function app = graspilot_app(varargin)
%GRASPILOT_APP Launch GRASPilot (same arguments as the GRASPilot class).
%
%   graspilot_app
%   graspilot_app('recon.mat')
%
%   Also the entry point for a standalone build with MATLAB Compiler:
%       mcc -m graspilot_app.m -a GRASPilot.m -a GRASPilotData.m -a graspilot_*.m
a = GRASPilot(varargin{:});
if nargout > 0, app = a; end
if isdeployed
    waitfor(a.Fig);
end
end
