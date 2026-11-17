function o = graspilot_pk_defaults()
%GRASPILOT_PK_DEFAULTS Default pharmacokinetic model options (GRASPilotData.PKOpts).
%
%   Models     'tofts', 'etofts', '2cxm'
%   T10        1400 ms   global pre-contrast T1 of tissue (brain lesion / grey matter at 3 T)
%   T10Blood   1664 ms   pre-contrast T1 of blood at 3 T (Lu et al. 2004, Hct 0.42); AIF
%   TR         4 ms, FlipAngle 12 deg: typical GRASP VIBE; set them to the protocol
%   r1         4.5 /(mM s): set to the agent and field strength
%   Hct        0.45      large-vessel hematocrit
%   AIF        'ref'     mean signal of the reference (SSS) ROI; or 'population'
%   AIFShift   0 s
%   Region     'head'    voxels fitted: inside the head | 'rois' | 'all'
%   Refine     true      Levenberg-Marquardt after the grid search
o = struct('Models', {{'tofts', 'etofts', '2cxm'}}, 'T10', 1400, 'T10Blood', 1664, ...
    'TR', 4, 'FlipAngle', 12, 'r1', 4.5, 'Hct', 0.45, 'AIF', 'ref', 'AIFShift', 0, ...
    'Region', 'head', 'Refine', true);
end
