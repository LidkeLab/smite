function [OTFVector]=rescaleOTF(PSF,Data)
%rescaleOTF OTF rescaling.
%
% INPUTS:
%    PSF         PSF image stack
%    Data
%
% OUTPUTS:
%    OFTVector   Optical Transfer Function

P.SZ=size(PSF,1);

[XGridS,YGridS]=meshgrid((-P.SZ/2:P.SZ/2-1),(-P.SZ/2:P.SZ/2-1));
    
Theta0=[1 2];
Optim=optimset();
OTFVector=fminsearch(@findCC,Theta0,Optim,PSF,Data,XGridS,YGridS)

end

function NegCC=findCC(Theta,PSF,Data,XGridS,YGridS)

SmoothKer=normpdf(gpuArray(XGridS),0,Theta(2)).*normpdf(gpuArray(YGridS),0,Theta(1));

A=gpuArray(PSF);
for zz=1:size(A,3)
    A(:,:,zz)=conv2(PSF(:,:,zz),SmoothKer,'same');
end

A=A-mean(A(:));
B=Data;
B=B-mean(B(:));
A=A/sqrt(sum(A(:).^2));
B=B/sqrt(sum(B(:).^2));
NegCC=-gather(sum(A(:).*B(:)));


end
