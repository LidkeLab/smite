function [PSFStruct,PSF]=phaseRetrieval_Spiral(PSFStruct,Data,MaxZCMag,MaxZCPhase,ZInfo,XYsubSample)
%phaseRetrieval_Spiral Phase retrieval using GS Algorithm.
%
%   The PSF is generated using the following model:
%   PSF = |F[OTF]|^2
%   % more....
%
%   PSF is normalized such that the integral over all space = 1
%
%   Input PSFStruct does not require Pupil field.
%
% INPUTS:
%   PSFStruct:  PSF Structure.
%   Data:       PSF Stack
%   MaxZCMag:   Limit for Zernike Expansion smoothing magnitude (Default=22)
%   MaxZCphase: Limit for Zernike Expansion smoothing phase (Default=81)
%               0 for no expansion. 
%   ZInfo:      Used to contruct PSFStruct.Z
%   XYsubSample: Used to normalize PSFStruct.PixelSize
%   
%
% OUTPUTS:
%   PSFStruct:  PSF Image Stack
%   PSF:        Point Spread Function
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%
% CITATION:
%   Hanser, ....

%Check for GPU, if not present, overwrite gpuArray
% try
%     gpuDevice
% catch
%     gpuArray=@(x)x;
% end

NMax_Mag = MaxZCMag;  %Zernike Expansion of Magnitude
NMax_Phase=MaxZCPhase;
NIterations=1000;
P=PSFStruct;
ImageSize=size(Data,1);
P.OSZ=ImageSize*3; %Factor of 3 from tiling

%Setup the arrays needed
KPixelSize=1/(P.OSZ*P.PixelSize);
[XGrid,YGrid]=meshgrid((-P.OSZ/2:P.OSZ/2-1),(-P.OSZ/2:P.OSZ/2-1));
R=sqrt(gpuArray(XGrid.^2+YGrid.^2));
PupilRadius=(P.NA/P.Lambda)/KPixelSize;
Mask=gpuArray(R<PupilRadius);

ZStruct_Mag=smi_psf.PointSpreadFunction.createZernikeStruct(P.OSZ,PupilRadius,NMax_Mag);
ZStruct_Phase=smi_psf.PointSpreadFunction.createZernikeStruct(P.OSZ,PupilRadius,NMax_Phase);

Kr_Image=KPixelSize.*Mask.*R;

Pupil_Mag=Mask;
Pupil_Phase=Mask*0;

OTFA_Stack=gpuArray(zeros(P.OSZ,P.OSZ,length(P.Z),'single'));
DefocusKernel=2*pi*sqrt(complex((P.N/P.Lambda)^2-Kr_Image.^2));

%Prepare the magnitude term from data
Data=gpuArray(single(Data));
Data=deTilt(Data);
PData=pad(Data); %Pad to OTFA Size

for ii=1:size(Data,3) %Inverse FFT Shift
    PSFA_Mag(:,:,ii) =ifftshift(sqrt(abs(PData(:,:,ii))));
end

L=2;
[~,PSFStruct_Spiral]=smi_psf.PointSpreadFunction.scalarPSFPrasadZone(P,L);
PupilSpiral = PSFStruct_Spiral.Pupil(:,:,2);

% Main GS Loop
for ii=1:NIterations
    ii
    %Apply Defocus
    for zz=1:length(P.Z)
        Defocus=P.Z(zz);
        PhaseIm=DefocusKernel.*Defocus+Pupil_Phase+PupilSpiral;
        OTFA_Stack(:,:,zz)=Mask.*Pupil_Mag.*exp(1i*PhaseIm);
        
        %Parseval Normalization
        OTFA_Stack(:,:,zz)=OTFA_Stack(:,:,zz)/sqrt(sum(sum(abs(OTFA_Stack(:,:,zz)).^2)))/ImageSize;
    end
    
    %Transform to real space and replace magnitude
    PSFA_Stack=fft2(OTFA_Stack);
    PSFA_Phase = angle(PSFA_Stack);
    
    for zz=1:length(P.Z)
        PSFA_Stack(:,:,zz)=PSFA_Mag(:,:,zz).*exp(1i*PSFA_Phase(:,:,zz));
    end
    
    %Transform back
    OTFA_Stack=ifft2(PSFA_Stack);
    
    %Remove defocus
    for zz=1:length(P.Z)
        Defocus=P.Z(zz);
        PhaseIm=DefocusKernel.*Defocus+PupilSpiral;
        OTFA_Stack(:,:,zz)=OTFA_Stack(:,:,zz).*exp(-1i*PhaseIm);
    end
    
    %Average stack
    Pupil=Mask.*mean(OTFA_Stack,3);
    
    %Constrain magnitude with zernike expansion
     [PSFStruct.ZC_Mag,Pupil_Mag]=smi_psf.PointSpreadFunction.zernikeExpansion(abs(Pupil),ZStruct_Mag);
    
    %Constrain phase with zernike expansion
    Pupil_Phase=angle(Pupil);
    if NMax_Phase>0
        [PSFStruct.ZC_Phase,Pupil_Phase]=smi_psf.PointSpreadFunction.zernikeExpansion(Pupil_Phase,ZStruct_Phase);
    end
    
    %Apply OTF Retriction
    Pupil_Mag=Mask.*Pupil_Mag;
    Pupil_Phase=Mask.*Pupil_Phase;
    
end

PSFStruct.Pupil=gpuArray(zeros(P.OSZ,P.OSZ,2,'single'));
PSFStruct.Pupil(:,:,1)=Pupil_Mag;
PSFStruct.Pupil(:,:,2)=Pupil_Phase+PupilSpiral;
PSFStruct.OSZ=P.OSZ;

if nargin > 4
    PSFStruct.Z = ZInfo(1):ZInfo(2):ZInfo(3);
end
if nargin > 5
    PSFStruct.PixelSize = PSFStruct.PixelSize / XYsubSample;
end

Model=smi_psf.PointSpreadFunction.scalarPSFPupil(PSFStruct);
PSFStruct.OTFSigma=smi_psf.PointSpreadFunction.rescaleOTF(Model,Data)*PSFStruct.PixelSize;
PSF=smi_psf.PointSpreadFunction.scalarPSFPupil(PSFStruct);

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


function Out=pad(In)
% Pad image by flipping and windowing
% 
% try
%     gpuDevice
% catch
%     gpuArray=@(x)x;
% end

L=8; %window length
INSZ=size(In,1); %Assumed square

[XGrid,YGrid]=meshgrid((0:INSZ-1),(0:INSZ-1));
WR=(gpuArray((cos(pi/L*XGrid))+1).*(XGrid<L))/2;
WB=(gpuArray((cos(pi/L*YGrid))+1).*(YGrid<L))/2;
WB(WB<0)=0;

WL=flip(WR,2);
WT=flip(WB,1);

TSZ=INSZ*3;
%Flip image to make 3x3
Out=gpuArray(zeros([TSZ,TSZ,size(In,3)],'single'));

Im1=flip(In,1);
Im2=flip(In,2);
Im12=flip(flip(In,1),2);

% loop over z
for zz = 1 : size(In,3)

%Left Column
Out(1:INSZ,1:INSZ,zz)=Im12(:,:,zz).*WL.*WT;
Out(INSZ+1:INSZ+INSZ,1:INSZ,zz)=Im2(:,:,zz).*WL;
Out(2*INSZ+1:2*INSZ+INSZ,1:INSZ,zz)=Im12(:,:,zz).*WL.*WB;

%Middle Column
Out(1:INSZ,INSZ+1:INSZ+INSZ,zz)=Im1(:,:,zz).*WT;
Out(INSZ+1:INSZ+INSZ,INSZ+1:INSZ+INSZ,zz)=In(:,:,zz);
Out(2*INSZ+1:2*INSZ+INSZ,INSZ+1:INSZ+INSZ,zz)=Im1(:,:,zz).*WB;

%Right Column
Out(1:INSZ,2*INSZ+1:2*INSZ+INSZ,zz)=Im12(:,:,zz).*WT.*WR;
Out(INSZ+1:INSZ+INSZ,2*INSZ+1:2*INSZ+INSZ,zz)=Im2(:,:,zz).*WR;
Out(2*INSZ+1:2*INSZ+INSZ,2*INSZ+1:2*INSZ+INSZ,zz)=Im12(:,:,zz).*WB.*WR;

end

end
