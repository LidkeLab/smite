function [PSF,P] =scalarPSFZernike(PSFStruct)
%scalarPSFZernike PSF stack from Zernike Coefficients based on a scalar model with OTF Rescaling.
%   The PSF is generated using the following model:
%   PSF = |F[OTF]|^2
%   % more...
%
%   Noll ordering is used as a linear index.    
%
%   PSF is normalized such that the integral over all space = 1
%
% INPUTS:
%   PSFStruct:  PSF Structure.  (Default = createPSFStruct())
%
% OUTPUTS:
%   PSF:        PSF Image Stack
%   P:          An updated PSFStruct 
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%
    
if nargin <1
    P=smi_psf.PointSpreadFunction.createPSFStruct();
else
    P=PSFStruct;
end


KPixelSize=1/(P.OSZ*P.PixelSize);
PupilRadius=(P.NA/P.Lambda)/KPixelSize;

NMax=max(length(P.ZC_Phase),length(P.ZC_Mag));
ZStruct=smi_psf.PointSpreadFunction.createZernikeStruct(P.OSZ,PupilRadius,NMax);

%magnitude
[Pupil_Mag]=smi_psf.PointSpreadFunction.zernikeSum(P.ZC_Mag,ZStruct);

%phase
[Pupil_Phase]=smi_psf.PointSpreadFunction.zernikeSum(P.ZC_Phase,ZStruct);

P.Pupil=cat(3,Pupil_Mag,Pupil_Phase);

[PSF,OTF]=smi_psf.PointSpreadFunction.scalarPSFPupil(P);


end
