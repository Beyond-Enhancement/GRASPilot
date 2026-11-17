function [C, valid] = graspilot_signal2conc(S, S0, T10, TR, FA, r1)
%GRASPILOT_SIGNAL2CONC Contrast agent concentration from DCE signal: signal-ratio (RATIO) method.
%
%   [C, valid] = graspilot_signal2conc(S, S0, T10, TR, FA, r1)
%
%   S     [N x T] signal curves (one per row)
%   S0    [N x 1] pre-contrast signal (mean of the baseline frames)
%   T10   pre-contrast T1 (s): scalar (global T10) or [N x 1]
%   TR    repetition time (s)
%   FA    flip angle (degrees)
%   r1    longitudinal relaxivity of the contrast agent (1/(mM s))
%   C     [N x T] concentration (mM); NaN where the signal ratio cannot be
%         produced by the spoiled gradient echo equation (S/S0 above the
%         fully relaxed limit) or S0 <= 0
%   valid [N x T] logical, false where C is NaN
%
%   The spoiled gradient echo (SPGR / FLASH) signal is
%       S = M0 sin(a) (1 - E1) / (1 - cos(a) E1),   E1 = exp(-TR R1).
%   Dividing by the baseline signal removes M0 (and any receiver or
%   reconstruction scaling, including the global-maximum normalization of
%   the GRASP reconstruction):
%       S/S0 = (1 - E1)(1 - cos(a) E10) / ((1 - E10)(1 - cos(a) E1)),
%   so with A = (S/S0) (1 - E10) / (1 - cos(a) E10):
%       E1 = (1 - A) / (1 - A cos(a)),   R1 = -ln(E1) / TR,
%       C  = (R1 - 1/T10) / r1.
%   Only T10 (here one global value, no T1 map), TR, the flip angle and r1
%   are needed (e.g. Schabel & Parker, Phys Med Biol 2008;53:2345).

S = double(S);
[N, T] = size(S);
S0 = double(S0(:));
T10 = double(T10(:));
if isscalar(T10), T10 = repmat(T10, N, 1); end
a = FA * pi / 180;
E10 = exp(-TR ./ T10);
k = (1 - E10) ./ (1 - cos(a) * E10);
A = (S ./ repmat(S0, 1, T)) .* repmat(k, 1, T);
E1 = (1 - A) ./ (1 - A * cos(a));
valid = A < 1 & A >= 0 & isfinite(A) & repmat(S0 > 0, 1, T);
E1(~valid) = NaN;
R1 = -log(E1) / TR;
C = (R1 - repmat(1 ./ T10, 1, T)) / r1;
C(~valid) = NaN;
end
