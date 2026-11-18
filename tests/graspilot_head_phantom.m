function [D, tr] = graspilot_head_phantom()
%GRASPILOT_HEAD_PHANTOM Head-shaped complex subspace DCE phantom (README screenshots, Auto AIF).
%
%   [D, tr] = graspilot_head_phantom()
%
%   128 x 128 x 30 voxels of 0.94 x 0.94 x 2 mm (LPS), 100 frames of 2 s,
%   complex subspace data (normalized like ReadData.m). Head with scalp and
%   brain; a ring-enhancing tumour (two-compartment exchange, PS varying
%   around the ring) with a necrotic core (tr.Ring, tr.Core), a slowly
%   enhancing lesion (tr.Slow), the superior sagittal sinus (tr.SSS:
%   dispersed, delayed venous blood with phase) and an MCA-like artery
%   with 30% blood volume fraction (tr.Art). Spoiled gradient echo signals
%   with the GRASPilotData PK defaults.
n = [128 128 30]; vox = [0.94 0.94 2];
T = 100; dt = 2; t = (0:T - 1) * dt;
TR = 4e-3; FA = 12; r1 = 4.5; T10 = 1.4; T10b = 1.664; Hct = 0.45;
spgr = @(M0, R1) M0 * sind(FA) * (1 - exp(-TR * R1)) ./ (1 - cosd(FA) * exp(-TR * R1));
ca = graspilot_pop_aif(t - t(5) - 1e-6); cp = ca / (1 - Hct);
tf = 0:0.02:t(end); d = 0.02 / 60;
ex = @(u, k, del) interp1(tf, filter(d * [0.5 0.5], [1 -exp(-k * d)], interp1(t, u, max(tf - del, 0), 'linear', 0)), t);
cxm = @(Fp, PS, vp, ve, del) cxmc(Fp, PS, vp, ve, del, cp, ex);
[X, Y, Z] = ndgrid((1:n(1)) * vox(1), (1:n(2)) * vox(2), (1:n(3)) * vox(3));
c = [60 62 30];
q = @(a) ((X - c(1)) / a(1)) .^ 2 + ((Y - c(2)) / a(2)) .^ 2 + ((Z - c(3)) / a(3)) .^ 2;
head = q([50 56 40]) <= 1; brain = q([45 51 36]) <= 1; scalp = head & ~brain;
% tumor: ring around a necrotic core, PS varying around the ring
tc = [c(1) + 17, c(2) - 13, c(3) + 4];
rr = sqrt((X - tc(1)) .^ 2 + (Y - tc(2)) .^ 2 + (Z - tc(3)) .^ 2);
ring = rr <= 11 & rr > 6; core = rr <= 6;
ang = atan2(Y - tc(2), X - tc(1));
lev = min(5, max(1, ceil((ang + pi) / (2 * pi) * 5)));
% slow lesion (Tofts-like)
r2 = sqrt((X - (c(1) - 20)) .^ 2 + (Y - (c(2) + 14)) .^ 2 + (Z - (c(3) + 2)) .^ 2);
slow = r2 <= 7;
% SSS: posterior-superior midline strip; artery: small lateral tube (MCA-like)
sss = abs(X - c(1)) <= 3 & brain & ~(q([41 47 32]) <= 1) & Y > c(2) + 15 & Z > c(3) + 4;
art = abs(X - (c(1) - 22)) <= 1.5 & abs(Y - (c(2) - 12)) <= 1.5 & abs(Z - c(3)) <= 6;
% curves
tau = 5; dl = 6;
cv = graspilot_expconv(ca, t, 1 / tau) / tau; cv = interp1(t, cv, max(t - dl, 0), 'linear', 0);
Cgm = 0.03 * cp + 0.002 * ex(cp, 0.2, 0);
St = spgr(1, 1 / T10 + r1 * Cgm);
Sb = spgr(1, 1 / T10b + r1 * ca);
types = [0.02 * ones(1, T); spgr(0.6, 1 / 0.9 + 0 * t); St; spgr(1, 1 / 2.2 + 0 * t); ...
    0.3 * Sb + 0.7 * St; spgr(1.3, 1 / T10b + r1 * cv) .* exp(1i * 0.5 * cv); ...
    spgr(1, 1 / T10 + r1 * 0.15 * 0.4 * ex(cp, 0.4, 3))];
PS = [0.03 0.06 0.1 0.15 0.22];
for k = 1:5
    types(end + 1, :) = spgr(1, 1 / T10 + r1 * cxm(0.45, PS(k), 0.07, 0.3, 3)); %#ok<AGROW>
end
lab = ones(n); lab(scalp) = 2; lab(brain) = 3; lab(core) = 4; lab(art) = 5; lab(sss) = 6; lab(slow) = 7;
for k = 1:5, lab(ring & lev == k) = 7 + k; end
[U, S] = svd(types.', 'econ'); K = nnz(diag(S) > 1e-9 * S(1)); U = U(:, 1:K);
rand('state', 3); randn('state', 4); %#ok<RAND>
ph = exp(1i * (2 * pi * 0.15 * (X / 120) + 0.3 * rand(n)));
V = ((U' * types.').') ; V = V(lab(:), :) .* repmat(ph(:), 1, K);
V = V + 0.004 * St(1) * complex(randn(size(V)), randn(size(V)));
vk = permute(reshape(single(V), [n K]), [3 1 2 4]);
D = GRASPilotData.fromSubspace(vk, single(conj(U)), 'CropFactor', 1, 'Normalize', true);
D.setAxisCodes('LPS'); D.VoxNative = vox; D.dt = dt; D.t0 = 0;
D.SourceFile = 'phantom_recon.mat';
tr = struct('Ring', ring, 'Core', core, 'Slow', slow, 'SSS', sss, 'Art', art, 'TumorCenter', tc ./ [vox(1) vox(2) vox(3)]);
end

function C = cxmc(Fp, PS, vp, ve, del, cp, ex)
kpo = Fp / vp; kpe = PS / vp; kep = PS / ve; a = kpo + kpe; b = kep;
s = sqrt((a + b) ^ 2 - 4 * kpo * b); Kp = (a + b + s) / 2; Km = (a + b - s) / 2;
A = (Kp - (kep + kpe)) / (Kp - Km);
C = Fp * (A * ex(cp, Kp, del) + (1 - A) * ex(cp, Km, del));
end
