function [D, tr] = graspilot_pk_phantom(noise)
%GRASPILOT_PK_PHANTOM Synthetic DCE series with known Tofts / extended Tofts / 2CXM parameters.
%
%   [D, tr] = graspilot_pk_phantom          noiseless
%   [D, tr] = graspilot_pk_phantom(sigma)   Gaussian noise (fraction of the tissue S0)
%
%   Dense LPS series, 32 x 32 x 6 voxels, 60 frames of 3 s, injection after
%   frame 5. Signals follow the spoiled gradient echo equation (TR 4 ms,
%   flip angle 12 deg, r1 4.5 /mM/s) with a tissue T10 of 1400 ms and a blood
%   T10 of 1664 ms (the GRASPilotData PK defaults). tr holds the masks
%   (canonical order) and the true parameters:
%     tr.SSS     vessel: blood concentration = population AIF (Parker 2006)
%     tr.Tofts   Ktrans 0.25 /min, ve 0.25 (vp 0)
%     tr.ETofts  Ktrans 0.10 /min, ve 0.20, vp 0.05
%     tr.CXM     Fp 50 ml/100ml/min, PS 12 ml/100ml/min, vp 0.08, ve 0.25
%   The forward model convolves on a 0.05 s grid (independent of the
%   fitting code). Hematocrit 0.45.

if nargin < 1, noise = 0; end
n = [32 32 6];
T = 60; dt = 3; t = (0:T - 1) * dt;
tinj = t(5) + 1e-6;                  % contrast arrives just after frame 5
TR = 4e-3; FA = 12; r1 = 4.5; T10 = 1.4; T10b = 1.664; Hct = 0.45;
spgr = @(M0, R1) M0 * sind(FA) * (1 - exp(-TR * R1)) ./ (1 - cosd(FA) * exp(-TR * R1));
Cb = graspilot_pop_aif(t - tinj);
Cp = Cb / (1 - Hct);
% fine-grid convolution of the piecewise-linear Cp with exp(-k t) (k in 1/min)
tf = 0:0.05:t(end);
Cpf = interp1(t, Cp, tf);
d = 0.05 / 60;
ex = @(k) interp1(tf, filter(d * [0.5 0.5], [1 -exp(-k * d)], Cpf), t);

V = zeros([n T]);
mk = @(i, j) local_box(n, i, j);
tr.Head = false(n); tr.Head(3:30, 3:30, :) = true;
tr.SSS = mk(5:8, 5:8);
tr.Tofts = mk(12:17, 5:10);
tr.ETofts = mk(12:17, 14:19);
tr.CXM = mk(22:27, 14:19);
tr.T10 = T10; tr.Hct = Hct; tr.t = t; tr.Cp = Cp;
% background brain tissue: no enhancement
Ctis = zeros(n(1) * n(2) * n(3), T);
% Tofts
tr.P.Tofts = struct('Ktrans', 0.25, 've', 0.25, 'kep', 1);
Ctis(tr.Tofts(:), :) = repmat(0.25 * ex(1), nnz(tr.Tofts), 1);
% extended Tofts
tr.P.ETofts = struct('Ktrans', 0.10, 've', 0.20, 'kep', 0.5, 'vp', 0.05);
Ctis(tr.ETofts(:), :) = repmat(0.10 * ex(0.5) + 0.05 * Cp, nnz(tr.ETofts), 1);
% 2CXM (rates per minute)
Fp = 0.5; PS = 0.12; vp = 0.08; ve = 0.25;
tr.P.CXM = struct('Fp', 100 * Fp, 'PS', 100 * PS, 'vp', vp, 've', ve, 'Ktrans', Fp * PS / (Fp + PS));
kpo = Fp / vp; kpe = PS / vp; kep = PS / ve; a = kpo + kpe; b = kep;
s = sqrt((a + b) ^ 2 - 4 * kpo * b);
Kp = (a + b + s) / 2; Km = (a + b - s) / 2; A = (Kp - (kep + kpe)) / (Kp - Km);
Ctis(tr.CXM(:), :) = repmat(Fp * (A * ex(Kp) + (1 - A) * ex(Km)), nnz(tr.CXM), 1);
S = spgr(1, 1 / T10 + r1 * Ctis);
S(tr.SSS(:), :) = repmat(spgr(1.3, 1 / T10b + r1 * Cb), nnz(tr.SSS), 1);
S(~tr.Head(:), :) = 0.02;
S = S / max(S(:));
if noise > 0
    s0 = median(S(tr.Head(:), 1));
    randn('state', 7); %#ok<RAND>
    S = abs(S + noise * s0 * randn(size(S)));
end
V = reshape(single(S), [n T]);
D = GRASPilotData.fromDense(V, 'xyzt');
D.setAxisCodes('LPS');
D.VoxNative = [1 1 1];
D.dt = dt; D.t0 = 0;
tr.BaseFrames = 1:5;
end

function M = local_box(n, i, j)
M = false(n);
M(i, j, 2:5) = true;
end
