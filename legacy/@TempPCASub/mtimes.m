function res = mtimes(a,b)

if isa(a,'TempPCASub') == 0
    error('In  A*B only A can be TempPCASub operator');
end
if a.adjoint
    [nx,ny,nt]=size(b);
    Data_Seq=reshape(b,nx*ny,nt);
    Data_Seq=Data_Seq';
    Data_PCA=(a.Phi*Data_Seq)';
    res=reshape(Data_PCA,nx,ny,size(a.Phi,1));
else
    [nx,ny,nt]=size(b);
    Data_Seq=reshape(b,nx*ny,nt);
    Data_PCA = (a.Phi' * Data_Seq')';
    res=reshape(Data_PCA,nx,ny,size(a.Phi,2));
end