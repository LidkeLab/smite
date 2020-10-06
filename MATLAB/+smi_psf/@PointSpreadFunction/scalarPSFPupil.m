function [PSF,PSFStruct]=scalarPSFPupil(PSFStruct)
%scalerPSFPupil PSF stack from Pupil based on a scalar model with OTF Rescaling.
%   The PSF is generated using the following model:
%   PSF = |F[OTF]|^2
%   % more...
%
%   PSF is normalized such that the integral over all space = 1
%   
% INPUTS:
%   PSFStruct:  PSF Structure.  (Default = createPSFStruct())
%
% OUTPUTS:
%   PSF:        PSF Image Stack
%   PSFStruct:  Updated PSF Structure.
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

% OTF:   Optical Transfer Function
% Pupil: Pupil Magnitude and Phase (SZ x SZ x 2)

%Set Defaults
if nargin <1
    P=smi_psf.PointSpreadFunction.createPSFStruct();
else
    P=PSFStruct;
end

if max(PSFStruct.OTFSigma)>0 %build smoothing kernel
    [XGridS,YGridS]=meshgrid((-P.SZ/2:P.SZ/2-1),(-P.SZ/2:P.SZ/2-1));
    SmoothKer=normpdf(gpuArray(XGridS),0,P.OTFSigma(2)/P.PixelSize).*...
        normpdf(gpuArray(YGridS),0,P.OTFSigma(1)/P.PixelSize);
end


[XGrid,YGrid]=meshgrid((-P.OSZ/2:P.OSZ/2-1),(-P.OSZ/2:P.OSZ/2-1));

R=sqrt(gpuArray(XGrid.^2+YGrid.^2));

KPixelSize=1/(P.OSZ*P.PixelSize);
PupilRadius=(P.NA/P.Lambda)/KPixelSize;
Mask=R<PupilRadius;
Kr_Image=KPixelSize.*Mask.*R;

Pupil_Mag=gpuArray(single(P.Pupil(:,:,1)));
Pupil_Phase=gpuArray(single(P.Pupil(:,:,2)));

%Sometimes this is too big, default to CPU
try
    OTFA_Stack=gpuArray(complex(zeros(P.OSZ,P.OSZ,length(P.Z),'single')));
catch
    OTFA_Stack=complex(zeros(P.OSZ,P.OSZ,length(P.Z),'single')); 
end

DefocusKernel=2*pi*sqrt(complex((P.N/P.Lambda)^2-Kr_Image.^2));

for zz=1:length(P.Z)
    Defocus=P.Z(zz);  %in microns
    PhaseIm=DefocusKernel.*Defocus+Pupil_Phase;
    
    OTFA=Mask.*Pupil_Mag.*exp(1i*PhaseIm);
    
    %Parseval Normalizaton
    Norm=sqrt(sum(sum(abs(OTFA).^2)))*P.OSZ;
    
    if isa(OTFA_Stack,'gpuArray') %Deal with large CPU arrays
        OTFA_Stack(:,:,zz)=OTFA/Norm;
    else
        OTFA_Stack(:,:,zz)=gather(OTFA/Norm);
    end
end

%Transform to real space
try %sometimes this is too big for some reason
    PSFA_Stack=fft2(OTFA_Stack);
catch
    PSFA_Stack=(fft2(gather(OTFA_Stack)));
    if exist('SmoothKer')
        SmoothKer=gather(SmoothKer);
    end
end

try
    SimDataStack=gpuArray(zeros(P.OSZ,P.OSZ,length(P.Z),'single'));
catch
    SimDataStack=(zeros(P.OSZ,P.OSZ,length(P.Z),'single'));
end


for zz=1:length(P.Z)
    PSFA_Stack(:,:,zz)=fftshift(PSFA_Stack(:,:,zz));
    if max(P.OTFSigma)>0 %apply smoothing in real space
        SimDataStack(:,:,zz)=conv2(abs(PSFA_Stack(:,:,zz)).^2,SmoothKer,'same');
    else
        SimDataStack(:,:,zz)=abs(PSFA_Stack(:,:,zz)).^2;
    end
end


% now cut to desired size: (only for even sizes) TODO
CenterPixel=P.OSZ/2+1;
PSF=SimDataStack(CenterPixel-P.SZ/2:CenterPixel+P.SZ/2-1,CenterPixel-P.SZ/2:CenterPixel+P.SZ/2-1,:);


end
