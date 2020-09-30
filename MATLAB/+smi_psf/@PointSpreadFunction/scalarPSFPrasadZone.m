function [PSF,PSFStruct]=scalarPSFPrasadZone(PSFStruct,L)
%scalarPSFPrasadZone PSF stack based on a scalar model Prasad Zones 
%   The PSF is generated using the following model:
%   PSF = |F[OTF]|^2
%   % more...
%
%   PSF is normalized such that the integral over all space = 1
%
%   Input PSFStruct does not require Pupil field.
%
% INPUTS:
%   PSFStruct:  PSF Structure.  (Default = createPSFStruct())
%   L:          expansion order
%
% OUTPUTS:
%   PSF:        PSF Image Stack
%   PSFStruct:  updated PSF Structure
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%

% We are really just setting up a pupil to send to the scalarPSFPupil()
% method

%Check for GPU, if not present, overwrite gpuArray
% try
%     gpuDevice
% catch
%     gpuArray=@(x)x;
% end


if nargin <1
    P=smi_psf.PointSpreadFunction.createPSFStruct();
else %recalc Pupil
    P=PSFStruct;
    [XGrid,YGrid]=meshgrid((-P.OSZ/2:P.OSZ/2-1),(-P.OSZ/2:P.OSZ/2-1));
    R=sqrt(gpuArray(XGrid.^2+YGrid.^2));
    KPixelSize=1/(P.OSZ*P.PixelSize);
    PupilRadius=(P.NA/P.Lambda)/KPixelSize;
    Mask=R<PupilRadius;
    R=R/PupilRadius;
    Pupil_Mag=Mask;
    Pupil_Phase=gpuArray(zeros(P.OSZ,P.OSZ,'single'));
    Theta=(gpuArray(atan2(YGrid,XGrid))); %CHECK!
    
    Alpha=1/2;
    for ll=1:L
        M=(R>=((ll-1)/L).^Alpha)&(R<(ll/L).^Alpha);
        Pupil_Phase(M)=mod((2*(ll-1)+1)*Theta(M),2*pi); %make hole
    end
    %dipshow(gather(Pupil_Phase))
    P.Pupil=cat(3,Pupil_Mag,Pupil_Phase);
        

end



[PSF]=smi_psf.PointSpreadFunction.scalarPSFPupil(P);

PSFStruct=P;


end
