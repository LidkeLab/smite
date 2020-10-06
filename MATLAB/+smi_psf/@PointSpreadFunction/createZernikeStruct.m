function  [ZStruct]=createZernikeStruct(SZ,Radius,NMax)
%createZernikeStruct Creates a createZernikeStruct.
%
%   ZStruct contains pre-computed images for quick expansion and sum
%
%   It contains the following fields:
%
%   SZ:             Image size (Default = 256)
%   Radius:         Radius of Pupil (Pixels) (Default = 64)
%   NMax:           Number of Zernike Images (Default = 21)
%   ZImages
%   R
%   Theta
%   NMax 
%   Norms
%
% OUTPUTS:
%   ZStruct:        ZStruct with all fields set to default
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%

%Check for GPU, if not present, overwrite gpuArray
% try
%     gpuDevice
% catch
%     gpuArray=@(x)x;
% end


if nargin <1
    SZ=256;
end

if nargin <2
    Radius=64;
end

if nargin <1
    NMax=21;
end


[XGrid,YGrid]=meshgrid((-SZ/2:SZ/2-1),(-SZ/2:SZ/2-1));
R=sqrt(gpuArray(XGrid.^2+YGrid.^2));
R=R/Radius;
Mask=R<=1;
R=R.*Mask;
Theta=(gpuArray(atan2(YGrid,XGrid)));
ZStruct.ZImages=gpuArray(zeros(NMax,SZ*SZ));

for nn=1:NMax
    Im=smi_psf.PointSpreadFunction.zernikeImage(nn,SZ,Radius,R,Theta,Mask);
    ZStruct.ZImages(nn,:)=Im(:);
end

A=1/sum(Mask(:));

for nn=1:NMax
    [N,M]=smi_psf.Zernike.zNoll2NM(nn);
    if M==0
        Em=2;
    else
        Em=1;
    end
    ZStruct.Norms(nn)=(2+2*N)/(Em)*A;   
end

ZStruct.SZ=SZ;
ZStruct.Radius=Radius;
ZStruct.NMax=NMax;
ZStruct.R=R;
ZStruct.Theta=Theta;
ZStruct.NMax=NMax;
