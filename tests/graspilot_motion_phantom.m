function [D, truth] = graspilot_motion_phantom(P, frames, streaks, lesion)
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
%
%   graspilot_motion_phantom(P, frames, a) adds radial streak artifacts of
%   amplitude a (the tissue signal is ~0.35): in every axial slice, lines
%   through the centre of the field of view at golden-angle steps that
%   rotate from frame to frame, as in stack-of-stars radial GRASP. They
%   are most visible outside the head. The brain then also gets a texture
%   (smooth grey/white-matter-like pattern and ventricles) for the
%   registration to lock on, as real brains have.
%
%   graspilot_motion_phantom(P, frames, a, true) also adds a strongly, slowly
%   enhancing lesion in the left occipital lobe (one side only), whose
%   enhancement can bias the registration.
if nargin < 3, streaks = 0; end
if nargin < 4, lesion = false; end
[S, ~] = graspilot_sss_phantom();
n = S.canonSize(); vox = S.canonVox(); T = S.NTRaw;
c = ([50 58 40] - 1) .* vox;
I = zeros([n T], 'single');
[g1, g2, g3] = ndgrid(0:n(1) - 1, 0:n(2) - 1, 0:n(3) - 1);
Z = [g1(:) g2(:) g3(:)] .* repmat(vox, numel(g1), 1);
clear g1 g2 g3
Pt = zeros(T, 6);
tex = ones(n);
les = zeros(n);
lc = [0 0 0 0.2 0.5 0.7 0.85 0.95 1 1 1 1];
if streaks > 0 || lesion
    [X, Y, Zg] = ndgrid(1:n(1), 1:n(2), 1:n(3));
    qq = @(a, cc) ((X - cc(1)) / a(1)) .^ 2 + ((Y - cc(2)) / a(2)) .^ 2 + ((Zg - cc(3)) / a(3)) .^ 2;
    inner = qq([38 48 44], [50 58 40]) <= 1;
    randn('state', 21); %#ok<RAND>
    g = randn(n);
    k = exp(-(-6:6) .^ 2 / (2 * 2.5 ^ 2)); k = k / sum(k);
    g = convn(convn(convn(g, k(:), 'same'), k(:)', 'same'), reshape(k, 1, 1, []), 'same');
    g = g / std(g(inner));
    tex(inner) = 1 + 0.25 * tanh(g(inner));                       % grey/white-like pattern
    tex(inner & (qq([5 12 8], [45 52 42]) <= 1 | qq([5 12 8], [55 52 42]) <= 1)) = 0.3;   % ventricles
    if lesion, les = 0.5 * double(qq([10 10 10], [70 82 44]) <= 1); end
end
for f = 1:T
    V = double(S.rawCanonVolume(f)) .* tex + les * lc(min(f, numel(lc)));
    if any(frames == f)
        Pt(f, :) = P;
        Rm = graspilot_rigid_matrix(P);
        % moved(z) = original(x) with z = R (x - c) + c + t
        X = (Z - repmat(c + P(1:3), size(Z, 1), 1)) * Rm + repmat(c, size(Z, 1), 1);
        q = X ./ repmat(vox, size(X, 1), 1) + 1;
        V = reshape(interpn(V, q(:, 1), q(:, 2), q(:, 3), 'linear', 0), n);
    end
    if streaks > 0
        V = V + repmat(streakPattern(n, f, streaks), [1 1 n(3)]);
    end
    I(:, :, :, f) = single(V);
end
D = GRASPilotData.fromDense(I, 'xyzt');
D.setAxisCodes('LPS');
D.VoxNative = vox;
D.dt = 4;
Still = S;
if streaks > 0 || lesion
    Is = zeros([n T], 'single');
    for f = 1:T, Is(:, :, :, f) = single(double(S.rawCanonVolume(f)) .* tex + les * lc(min(f, numel(lc)))); end
    Still = GRASPilotData.fromDense(Is, 'xyzt'); Still.setAxisCodes('LPS'); Still.VoxNative = vox; Still.dt = 4;
end
truth = struct('Params', Pt, 'Centre', c, 'Still', Still);
end

function A = streakPattern(n, f, a)
% 6 radial streaks per frame, angles advancing by the golden angle
[x, y] = ndgrid((1:n(1)) - (n(1) + 1) / 2, (1:n(2)) - (n(2) + 1) / 2);
A = zeros(n(1), n(2));
ga = 111.246;
for m = 1:6
    th = mod((6 * (f - 1) + m) * ga, 180);
    d = abs(-sind(th) * x + cosd(th) * y);            % distance to the line through the centre
    A = A + a * (1 - 2 * mod(m, 2)) * exp(-d .^ 2 / (2 * 0.8 ^ 2));
end
end
