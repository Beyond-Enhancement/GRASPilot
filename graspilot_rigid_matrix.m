function R = graspilot_rigid_matrix(p)
%GRASPILOT_RIGID_MATRIX Rotation matrix of rigid motion parameters.
%
%   R = graspilot_rigid_matrix(p)    p = [tL tP tS rL rP rS] (mm, degrees)
%
%   R = Rz(rS) * Ry(rP) * Rx(rL): rotations about canonical axis 1 (L-R),
%   2 (A-P) and 3 (I-S), applied in that order. See graspilot_motion_estimate.
a = p(4) * pi / 180; b = p(5) * pi / 180; c = p(6) * pi / 180;
Rx = [1 0 0; 0 cos(a) -sin(a); 0 sin(a) cos(a)];
Ry = [cos(b) 0 sin(b); 0 1 0; -sin(b) 0 cos(b)];
Rz = [cos(c) -sin(c) 0; sin(c) cos(c) 0; 0 0 1];
R = Rz * Ry * Rx;
end
