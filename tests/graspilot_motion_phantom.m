function [D, truth] = graspilot_motion_phantom(P, frames)
%GRASPILOT_MOTION_PHANTOM Venous phantom with known rigid head motion.
%
%   [D, truth] = graspilot_motion_phantom(P, frames)
%
%   The frames listed in "frames" of graspilot_sss_phantom are moved
%   rigidly with parameters P (1x6: mm, degrees; see
%   graspilot_motion_estimate) about the phantom head centre: the head point
%   x of the original appears at y = Rot (x - c) + c + t. truth.Params
%   [T x 6] holds the motion of every frame, truth.Centre the centre (mm),
%   truth.Still the motion-free series (GRASPilotData).
[S, ~] = graspilot_sss_phantom();
n = S.canonSize(); vox = S.canonVox(); T = S.NTRaw;
c = ([50 58 40] - 1) .* vox;
I = zeros([n T], 'single');
[g1, g2, g3] = ndgrid(0:n(1) - 1, 0:n(2) - 1, 0:n(3) - 1);
Z = [g1(:) g2(:) g3(:)] .* repmat(vox, numel(g1), 1);
clear g1 g2 g3
Pt = zeros(T, 6);
for f = 1:T
    V = double(S.rawCanonVolume(f));
    if any(frames == f)
        Pt(f, :) = P;
        Rm = graspilot_rigid_matrix(P);
        % moved(z) = original(x) with z = R (x - c) + c + t
        X = (Z - repmat(c + P(1:3), size(Z, 1), 1)) * Rm + repmat(c, size(Z, 1), 1);
        q = X ./ repmat(vox, size(X, 1), 1) + 1;
        V = reshape(interpn(V, q(:, 1), q(:, 2), q(:, 3), 'linear', 0), n);
    end
    I(:, :, :, f) = single(V);
end
D = GRASPilotData.fromDense(I, 'xyzt');
D.setAxisCodes('LPS');
D.VoxNative = vox;
D.dt = 4;
truth = struct('Params', Pt, 'Centre', c, 'Still', S);
end
