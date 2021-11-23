function [PSFStruct]=phaseRetrievalEM(PSFStruct,Data)
%phaseRetrievalEM Phase retrieval plus EM to optimize NA, Lambda.
%
%   The PSF is generated using the following model:
%   PSF = |F[OTF]|^2
%   % more...
%
%   PSF is normalized such that the integral over all space = 1
%
%   Input PSFStruct does not require Pupil field.
%
% INPUTS:
%   PSFStruct:  PSF Structure.
%   Data:       PSF Stack
%
% OUTPUTS:
%   PSFStruct:  PSF Image Stack
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%
% CITATION:
%   Hanser, ....

Data=gpuArray(single(Data));
Data=deTilt(Data);

Theta0=[PSFStruct.NA];
Theta_Found=fminsearch(@SSE,Theta0,optimset,Data,PSFStruct);

PSFStruct.NA=Theta_Found(1);
%PSFStruct.Lambda=Theta_Found(2);

end

function [fval,Model]=SSE(Theta,Data,PSFStruct)

NZF=49;

PSFStruct.NA=Theta(1);
%PSFStruct.Lambda=Theta(2);
[P]=smi_psf.PointSpreadFunction.phaseRetrieval(PSFStruct,Data);

KPixelSize=1/(size(P.Pupil,1)*P.PixelSize);
PupilRadius=(P.NA/P.Lambda)/KPixelSize;
ZS=smi_psf.PointSpreadFunction.createZernikeStruct(P.OSZ,PupilRadius,NZF);
P.ZC_Phase=gather(smi_psf.PointSpreadFunction.zernikeExpansion(P.Pupil(:,:,2),ZS));
P.ZC_Mag=gather(smi_psf.PointSpreadFunction.zernikeExpansion(P.Pupil(:,:,1),ZS));
[Model]=smi_psf.PointSpreadFunction.scalarPSFZernike(P);

fval=gather(sum( (Data(:)-Model(:)).^2));
L=size(Data,3);
M=[];
D=[];
for nn=1:L
    M=cat(1,M,stretch(gather(Model(:,:,nn))));
    D=cat(1,D,stretch(gather(Data(:,:,nn))));
end

%dipshow(1234,cat(2,D,M))
pause(.1)

end

function Out=deTilt(In)
%Remove offset from image. Important!

Out=In;
SZ=size(In,1);
PlaneFun=@(Theta,X,Y,Data)double(gather(mean((Theta(1)*X+Theta(2)*Y+Theta(3)-Data).^2)));
[XGrid,YGrid]=meshgrid((1:SZ),(1:SZ));

for ii=1:size(In,3)
    Left = In(:,1,ii);
    Right = In(:,end,ii);
    Top = In(1,:,ii);
    Bottom = In(end,:,ii);
    Data=cat(1,Left(:),Right(:),Top(:),Bottom(:));
    
    V=(1:SZ)';
    V1=ones(SZ,1);
    
    X=cat(1,V1,SZ*V1,V,V);
    Y=cat(1,V,V,V1,V1*SZ);
    
    Theta=fminsearch(PlaneFun,[1 1 1],optimset(),X,Y,Data);
    Out(:,:,ii)=In(:,:,ii)-XGrid*Theta(1)-YGrid*Theta(2)-Theta(3);
    Out(:,:,ii)=Out(:,:,ii)/sum(sum(Out(:,:,ii))); %Normalize
end

end
