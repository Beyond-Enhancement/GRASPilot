function varargout = GRASPilot(varargin)
%GRASPILOT Launch GRASPilot, the GRASP DCE-MRI viewer and ROI analysis tool.
%
%   GRASPilot                          then File > Open .mat
%   GRASPilot('recon.mat')             opens the load dialog for this file
%   GRASPilot(vk_all, basis_all)       variables already in the workspace
%   app = GRASPilot(...)               returns the viewer object for scripting
%
%   Same arguments as the GraspDCEViewer class, which implements the tool.
app = GraspDCEViewer(varargin{:});
if nargout > 0, varargout{1} = app; end
end
