function [D, truth] = graspilot_sss_phantom(zmin, mode)
%GRASPILOT_SSS_PHANTOM Synthetic DCE head with dural sinuses for testing graspilot_detect_sss.
%
%   [D, truth] = graspilot_sss_phantom          full coverage
%   [D, truth] = graspilot_sss_phantom(zmin)    slab starting at slice zmin
%
%   1 mm isotropic, LPS voxel order (axis codes 'LPS'), 12 frames. Contains
%   the superior sagittal sinus (thickest posteriorly), torcula, transverse
%   sinuses (descending laterally), straight sinus (joining the torcula from
%   anterior-superior, ~2 mm from the SSS just above the torcula), a scalp
%   vein as bright as the sinuses and a brighter posterior artery near the
%   midline. truth holds the masks and the torcula slice.
%
%   [D, truth] = graspilot_sss_phantom(zmin, 'realistic') adds what real
%   post-contrast data have around the sinus: bridging cortical veins
%   draining into the SSS (one low posterior), an enhancing falx, tentorium
%   and dural shell, and a ring-enhancing tumor beside the posterior
%   midline; 2 mm slices; stored rotated in-plane and flipped S-I (native
%   axis codes 'PLI', already set). truth is in canonical (LPS) order.

if nargin < 1 || isempty(zmin), zmin = 1; end
real = nargin > 1 && strcmpi(mode, 'realistic');
n = [100 120 90];
[X, Y, Z] = ndgrid(1:n(1), 1:n(2), 1:n(3));
c = [50 58 40];
q = @(a) ((X - c(1)) / a(1)) .^ 2 + ((Y - c(2)) / a(2)) .^ 2 + ((Z - c(3)) / a(3)) .^ 2;
head = q([44 54 50]) <= 1;
inner = q([38 48 44]) <= 1;
skull = head & ~inner & q([41 51 47]) <= 1;
scalp = head & ~inner & ~skull;

% superior sagittal sinus: midline arc 5 mm inside the inner table
th = (-35:0.5:150)';
P = [50 * ones(size(th)), 58 + 43 * cosd(th), 40 + 39 * sind(th)];
r = 4.5 - 2 * max(0, th) / 150;
tor = P(1, :);
SSS = tube(n, P, r) | tube(n, tor, 6);
% transverse sinuses
ph = [(-60:0.5:-5)'; (5:0.5:60)'];
PT = [50 + 30 * sind(ph), 58 + 35.2 * cosd(ph), tor(3) - 6 * abs(sind(ph))];
TS = tube(n, PT, 3.5 * ones(size(ph)));
% straight sinus
s = (0:0.01:1)';
PS = repmat(tor, numel(s), 1) + s * ([50 70 42] - tor);
SS = tube(n, PS, 2.5 * ones(size(s)));
% scalp vein (outside the skull) and a posterior artery near the midline
zz = (25:0.5:60)';
PV = [50 * ones(size(zz)), 58 + 52.5 * sqrt(1 - ((zz - 40) / 50) .^ 2) - 1, zz];
SV = tube(n, PV, 1.5 * ones(size(zz))) & scalp;
s2 = (0:0.01:1)';
PA = repmat([44 72 28], numel(s2), 1) + s2 * ([47 96 46] - [44 72 28]);
AR = tube(n, PA, 1.2 * ones(size(s2)));

ven = [0 0 0 0.3 0.8 1.0 0.85 0.7 0.6 0.55 0.5 0.5];
dur = 0.45 * [0 0 0 0.2 0.5 0.7 0.8 0.85 0.9 0.92 0.94 0.95];     % dura, falx, tentorium, tumor (slow)
CV = false(n); DU = false(n); TU = false(n);
if real
    % bridging cortical veins on the brain surface, from the SSS laterally
    ai = [38 48 44] * 0.93;
    for tj = [20 45 65 85 105 125]
        P0 = [50, 58 + 43 * cosd(tj), 40 + 39 * sind(tj)];
        for sd = [-1 1]
            if tj == 20 && sd > 0, continue; end            % one low posterior vein
            u0 = (P0 - c) ./ ai; u0 = u0 / norm(u0);
            u1 = [sd * 0.9, 0.3 * cosd(tj), 0.3 * sind(tj)]; u1 = u1 / norm(u1);
            sv = (0:0.01:1)';
            U = (1 - sv) * u0 + sv * u1;
            U = U ./ repmat(sqrt(sum(U .^ 2, 2)), 1, 3);
            CV = CV | tube(n, repmat(c, numel(sv), 1) + U .* repmat(ai, numel(sv), 1), 1.4 * ones(size(sv)));
        end
    end
    CV = CV & inner;
    % falx (midline sheet below the SSS), tentorium (sloping sheet from the torcula), dural shell
    falx = abs(X - 50) <= 0.8 & inner & Z > c(3) - 5 & ~(q([24 32 26]) <= 1);
    tent = abs(Z - (tor(3) + 0.45 * (tor(2) - Y))) <= 0.8 & inner & Y > 55 & abs(X - 50) <= 34;
    shell = inner & ~(q([37 47 43]) <= 1);
    DU = falx | tent | shell;
    % ring-enhancing tumor 14 mm left of the midline, ~22 mm from the posterior inner table
    rt = sqrt((X - 64) .^ 2 + (Y - 84) .^ 2 + (Z - 36) .^ 2);
    TU = rt <= 9 & rt > 5;
end
art = [0 0 0 1.3 0.9 0.5 0.35 0.3 0.28 0.27 0.26 0.25];
tis = 0.04 * [0 0 0 0.2 0.5 0.7 0.8 0.85 0.9 0.92 0.94 0.95];
base = 0.35 * inner + 0.03 * skull + 0.5 * scalp;
veins = SSS | TS | SS | SV | CV;
base(veins | AR) = 0.3;
T = numel(ven);
I = zeros([n T], 'single');
randn('state', 11); %#ok<RAND>
for t = 1:T
    f = base + tis(t) * (inner | scalp) + ven(t) * veins + art(t) * AR + dur(t) * (DU & ~veins) + ...
        1.4 * dur(t) * TU;
    I(:, :, :, t) = single(f + 0.005 * randn(n));
end
k = zmin:n(3);
if real, k = k(1:2:end); end                          % 2 mm slices
I = I(:, :, k, :);
if real
    D = GRASPilotData.fromDense(flip(permute(I, [2 1 3 4]), 3), 'xyzt');   % rows A->P, cols R->L, slices S->I
    D.setAxisCodes('PLI');
    D.VoxNative = [1 1 2];
else
    D = GRASPilotData.fromDense(I, 'xyzt');
    D.setAxisCodes('LPS');
    D.VoxNative = [1 1 1];
end
D.dt = 4;
truth = struct('SSS', SSS(:, :, k), 'TS', TS(:, :, k) & ~SSS(:, :, k), 'SS', SS(:, :, k) & ~SSS(:, :, k), ...
    'SV', SV(:, :, k), 'AR', AR(:, :, k), 'TorculaZ', tor(3) - zmin + 1, 'CV', CV(:, :, k) & ~SSS(:, :, k), ...
    'Tumor', TU(:, :, k));
if real, truth.TorculaZ = find(k <= tor(3), 1, 'last'); end
end

function M = tube(n, P, r)
M = false(n);
for i = 1:size(P, 1)
    lo = max(1, floor(P(i, :) - r(i)));
    hi = min(n, ceil(P(i, :) + r(i)));
    if any(hi < lo), continue; end
    [x, y, z] = ndgrid(lo(1):hi(1), lo(2):hi(2), lo(3):hi(3));
    in = (x - P(i, 1)) .^ 2 + (y - P(i, 2)) .^ 2 + (z - P(i, 3)) .^ 2 <= r(i) ^ 2;
    M(sub2ind(n, x(in), y(in), z(in))) = true;
end
end
