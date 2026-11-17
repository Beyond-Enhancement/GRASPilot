function [D, tr] = graspilot_ratio_phantom(noise)
%GRASPILOT_RATIO_PHANTOM Complex subspace DCE phantom for RATIO and the VCAI.
%
%   [D, tr] = graspilot_ratio_phantom          noiseless
%   [D, tr] = graspilot_ratio_phantom(sigma)   complex Gaussian noise on the
%                                              coefficients (fraction of the tissue S0)
%
%   32 x 32 x 8 voxels (LPS), 120 frames of 2 s, injection after frame 5,
%   stored as subspace coefficients on a shared complex temporal basis
%   (vk_all, basis_all), so magnitude, phase and the basis are available.
%   Spoiled gradient echo signals (TR 4 ms, flip angle 12 deg, r1 4.5 /mM/s)
%   with tissue T10 1400 ms and blood T10 1664 ms (GRASPilotData defaults).
%     tr.Artery  arterial voxels: volume blood fraction 0.3 in tissue
%                (tr.WSignal = signal-weighted blood fraction)
%     tr.SSS     sinus: blood concentration = arterial curve dispersed with
%                exp(-t/tau)/tau, tau = 5 s, and delayed 6 s; phase = 0.5 rad/mM
%     tr.Lesion2C  2CXM, Fp 0.5, PS 0.12, vp 0.08, ve 0.25 (ml/ml/min),
%                delay 3 s: kep 0.387, kpe 1.21 /min, MTT 9.6 s, E 0.194
%     tr.LesionSlow  Tofts-like (one-compartment slow): k 0.5 /min
%     tr.LesionFast  vascular only, transit k 12 /min (one-compartment fast)
%   Every voxel has its own constant phase offset.

if nargin < 1, noise = 0; end
n = [32 32 8];
T = 120; dt = 2; t = (0:T - 1) * dt;
TR = 4e-3; FA = 12; r1 = 4.5; T10 = 1.4; T10b = 1.664; Hct = 0.45;
spgr = @(M0, R1) M0 * sind(FA) * (1 - exp(-TR * R1)) ./ (1 - cosd(FA) * exp(-TR * R1));
ca = graspilot_pop_aif(t - t(5) - 1e-6);                 % arterial blood (mM)
cp = ca / (1 - Hct);
% fine-grid convolution of a delayed piecewise-linear input with exp(-k t)
tf = 0:0.02:t(end); d = 0.02 / 60;
ex = @(u, k, del) interp1(tf, filter(d * [0.5 0.5], [1 -exp(-k * d)], ...
    interp1(t, u, max(tf - del, 0), 'linear', 0)), t);
tm = (t - t(1)) / 60;
% sinus: dispersed, delayed arterial blood
tau = 5; dl = 6;
cv = graspilot_expconv(ca, t, 1 / tau) / tau;
cv = interp1(t, cv, max(t - dl, 0), 'linear', 0);
% tissue (brain background): small vascular component and slow leakage
Cgm = 0.03 * cp + 0.005 * ex(cp, 0.3, 0);
% lesions
Fp = 0.5; PS = 0.12; vp = 0.08; ve = 0.25; del2 = 3;
kpo = Fp / vp; kpe = PS / vp; kep = PS / ve; a = kpo + kpe; b = kep;
s = sqrt((a + b) ^ 2 - 4 * kpo * b);
Kp = (a + b + s) / 2; Km = (a + b - s) / 2; A = (Kp - (kep + kpe)) / (Kp - Km);
C2 = Fp * (A * ex(cp, Kp, del2) + (1 - A) * ex(cp, Km, del2));
E = PS / (PS + Fp);
tr.P2C = struct('kep', kep * (1 - E), 'kpe', kpo * E, 'MTT', 60 / kpo, 'E', E, 'Delay', del2, ...
    'vp', vp, 'Ktrans', E * Fp);
Cslow = 0.2 * 0.5 * ex(cp, 0.5, 2);
Cfast = 0.06 * 12 * ex(cp, 12, 2);
tr.kSlow = 0.5; tr.kFast = 12;
% signals (rows = curve types)
St = spgr(1, 1 / T10 + r1 * Cgm);
Sb = spgr(1, 1 / T10b + r1 * ca);
fb = 0.3;
Sart = fb * Sb + (1 - fb) * St;
tr.WSignal = fb * Sb(1) / Sart(1);
types = [St; Sart; spgr(1.3, 1 / T10b + r1 * cv) .* exp(1i * 0.5 * cv); spgr(1, 1 / T10 + r1 * C2); ...
    spgr(1, 1 / T10 + r1 * Cslow); spgr(1, 1 / T10 + r1 * Cfast); 0.02 * ones(1, T)];
[U, S] = svd(types.', 'econ');
K = nnz(diag(S) > 1e-9 * S(1));
U = U(:, 1:K);
lab = 7 * ones(n);                                       % outside the head
head = false(n); head(3:30, 3:30, :) = true;
lab(head) = 1;
tr.Artery = false(n); tr.Artery(15:16, 6:7, 2:7) = true;
tr.SSS = false(n); tr.SSS(5:8, 24:27, 2:7) = true;
tr.Lesion2C = false(n); tr.Lesion2C(18:23, 14:19, 3:6) = true;
tr.LesionSlow = false(n); tr.LesionSlow(8:12, 12:16, 3:6) = true;
tr.LesionFast = false(n); tr.LesionFast(24:28, 6:10, 3:6) = true;
lab(tr.Artery) = 2; lab(tr.SSS) = 3; lab(tr.Lesion2C) = 4; lab(tr.LesionSlow) = 5; lab(tr.LesionFast) = 6;
tr.Head = head;
rand('state', 11); randn('state', 12); %#ok<RAND>
ph = exp(1i * 2 * pi * rand(n));
coef = (U' * types.').';                                 % [ntypes x K], z_col = U * coef.'
V = coef(lab(:), :) .* repmat(ph(:), 1, K);              % [nvox x K]
if noise > 0
    V = V + noise * St(1) / sqrt(2) * complex(randn(size(V)), randn(size(V)));
end
vk = permute(reshape(single(V), [n K]), [3 1 2 4]);      % [nz nx ny K]
D = GRASPilotData.fromSubspace(vk, single(conj(U)), 'CropFactor', 1, 'Normalize', false);
D.setAxisCodes('LPS');
D.VoxNative = [1 1 1];
D.dt = dt; D.t0 = 0;
tr.BaseFrames = 1:5;
tr.Tau = tau; tr.Delta = dl; tr.K = K; tr.t = t; tr.ca = ca;
end
