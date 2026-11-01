function app = gdce_app(varargin)
%GDCE_APP Launch the GRASP DCE Viewer (same arguments as GraspDCEViewer).
%
%   gdce_app
%   gdce_app('recon.mat')
%
%   Also the entry point for a standalone build with MATLAB Compiler:
%       mcc -m gdce_app.m -a GraspDCEViewer.m -a GraspData.m -a gdce_*.m
a = GraspDCEViewer(varargin{:});
if nargout > 0, app = a; end
if isdeployed
    waitfor(a.Fig);
end
end
