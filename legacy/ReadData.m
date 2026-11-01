
clear data1
Recon_Allslice=permute(vk_all,[2 3 4 1]);
for ii=1:size(basis_all,1)
    
    Phi=single(squeeze(squeeze(basis_all(ii,:,:))));
    PCA=TempPCASub(Phi);
    data1(:,:,:,ii)=abs(PCA'*Recon_Allslice(:,:,:,ii));
end
data1=CropImg(data1/max(data1(:)),size(vk_all,2)/2,size(vk_all,2)/2);
% figure,imagescn(abs(data1(:,:,:,60)),[0 .6],[],[],3),axis image

% size of data1: nx,ny, time/frame, slices,