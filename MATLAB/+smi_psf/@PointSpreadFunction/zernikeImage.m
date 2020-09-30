function [Image]=zernikeImage(NollIndex,SZ,Radius,R,Theta,Mask)
%zernikeImage Generate a Zernike Polynomical Image from Noll Coefficient
%
%
%   Noll ordering is used as a linear index.
%
%   More about formula and conventions.
%
% INPUTS:
%   NollIndex:      Noll Linear Index
%   SZ:             Size of returned image. Image is square.
%   Radius:         Circle Edge Radius (Pixels)
%   R:              Radius Image (Recalculated by Default)
%   Theta:          Theta Image (Recalculated by Default)
%   Mask:           Mask to apply to image (Default: pixels with R <= 1)
%
% OUTPUTS:
%   Image:          Zernike Polynomial Image
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%

OSZ=SZ;

if nargin <5
    [XGrid,YGrid]=meshgrid((-OSZ/2:OSZ/2-1),(-OSZ/2:OSZ/2-1));
    R=sqrt(gpuArray(XGrid.^2+YGrid.^2));
    R=R/Radius;
    Mask=R<=1;
    R=R.*Mask;
    Theta=(gpuArray(atan2(YGrid,XGrid))); %CHECK!
end

%convert to NM
[N, M] = smi_psf.Zernike.zNoll2NM(NollIndex);

if M<0
    Image=Mask.*sin(M*Theta);
else
    Image=Mask.*cos(M*Theta);
end


N=abs(N);
M=abs(M);

%calculate the radial polynomial
RP=0;
for kk=0:(N-M)/2
    Poly=R.^(N-2*kk);
    %Coef = (-1).^kk*factorial(N-kk)/...
    %    (factorial(kk)*factorial((N+M)/2-kk)*factorial((N-M)/2-kk));
    Coef = (-1)^kk * prod((N - M)/2 - kk + 1 : N - kk) / ...
                     (factorial(kk) * factorial((N + M)/2 - kk));
    RP=RP+Coef*Poly;
end

Image=Mask.*Image.*RP;
