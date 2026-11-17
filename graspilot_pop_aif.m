function Cb = graspilot_pop_aif(t)
%GRASPILOT_POP_AIF Population-averaged arterial input function (Parker et al. 2006).
%
%   Cb = graspilot_pop_aif(t)
%
%   t    times (s) after the contrast injection; Cb = 0 for t < 0
%   Cb   whole-blood gadolinium concentration (mM) for a standard dose of
%        0.1 mmol/kg; the plasma concentration is Cb / (1 - Hct)
%
%   Two Gaussians (first and recirculation peak) plus an exponential
%   modulated by a sigmoid (Parker GJM et al., Magn Reson Med
%   2006;56:993-1000), with the published population parameters. The time
%   origin of the model is the injection; its first peak is at ~10 s.

tm = double(t) / 60;                         % minutes
A1 = 0.809; A2 = 0.330;                      % mmol min
T1 = 0.17046; T2 = 0.365;                    % min
s1 = 0.0563; s2 = 0.132;                     % min
alpha = 1.050; beta = 0.1685;                % mM, 1/min
s = 38.078; tau = 0.483;                     % 1/min, min
Cb = A1 / (s1 * sqrt(2 * pi)) * exp(-(tm - T1) .^ 2 / (2 * s1 ^ 2)) + ...
    A2 / (s2 * sqrt(2 * pi)) * exp(-(tm - T2) .^ 2 / (2 * s2 ^ 2)) + ...
    alpha * exp(-beta * tm) ./ (1 + exp(-s * (tm - tau)));
Cb(tm < 0) = 0;
end
