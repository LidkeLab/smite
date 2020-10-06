function [PSF,PSFStruct] = scalarPSF(PSFStruct)
%scalarPSF PSF stack based on a scalar model with OTF rescaling.
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
%
% OUTPUTS:
%   PSF:        PSF Image Stack
%   OTF:        Optical Transfer Function
%   Pupil:      Pupil Magnitude and Phase (SZ x SZ x 2)
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%

% We are really just setting up a pupil to send to the scalarPSFPupil()
% method

if nargin <1
    P=smi_psf.PointSpreadFunction.createPSFStruct();
end

if nargin <1
    P=smi_psf.PointSpreadFunction.createPSFStruct();
else %recalc Pupil
    P=PSFStruct;
    [XGrid,YGrid]=meshgrid((-P.OSZ/2:P.OSZ/2-1),(-P.OSZ/2:P.OSZ/2-1));
    R=sqrt(gpuArray(XGrid.^2+YGrid.^2));
    KPixelSize=1/(P.OSZ*P.PixelSize);
    PupilRadius=(P.NA/P.Lambda)/KPixelSize;
    Mask=R<PupilRadius;
    Pupil_Phase=gpuArray(zeros(P.OSZ,P.OSZ,'single'));
    Pupil_Mag=Mask;
    P.Pupil=cat(3,Pupil_Mag,Pupil_Phase);
end

[PSF,PSFStruct]=smi_psf.PointSpreadFunction.scalarPSFPupil(P);

end
