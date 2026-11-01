function res= TempPCASub(Phi)
% 
% implements a temporal FFT along the dim dimenison
%

res.adjoint = 0;
res.Phi=Phi;
res = class(res,'TempPCASub');
