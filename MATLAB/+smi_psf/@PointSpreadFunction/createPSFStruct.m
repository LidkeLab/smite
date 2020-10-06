function  [PSFStruct]=createPSFStruct()
%createPSFStruct Creates a default PSF Structure.
%
%   PSFStruct is a compact representation of a PSF
%
%   It contains the following fields: 
%   
%   Pupil:          Pupil image [Magnitude;Phase]. (SZ x SZ x 2)
%   Z:              Vector of Z positions (microns) (Default = (-2:.1:2))
%   Lambda:         Emission Wavelength (microns) (Default = .69)
%   NA:             Objective Numerical Aperture (Default = 1.35)
%   N:              Index of Refraction (Default = 1.4)
%   PixelSize:      Lateral pixel size (microns) (Default = 0.1)
%   SZ:             PSF lateral size (Pixels) (Default = 64)
%   OTFSigma:       OTF Rescaling factor Sigma [Y X]. (micron) (Default = [0 0])
%   OSZ:            OTF size (Pixels). (Default = 256)
%   ZC_Mag:         Zernike Coef Magnitude (Noll Order) (Default = [1])
%   ZC_Phase:       Zernike Coef Phase (Noll Order) (Default = [0])

% OUTPUTS:
%   PSFStruct:     	PSFStruct will all fields set to default
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

PSFStruct.Z=(-2:.1:2);
PSFStruct.Lambda=.69;
PSFStruct.NA=1.35;
PSFStruct.N=1.4;
PSFStruct.PixelSize=.1;
PSFStruct.SZ=64;
PSFStruct.OTFSigma=[0 0];
PSFStruct.OSZ=256;
PSFStruct.ZC_Mag=[1];
PSFStruct.ZC_Phase=[0];

[XGrid,YGrid]=meshgrid((-PSFStruct.OSZ/2:PSFStruct.OSZ/2-1),(-PSFStruct.OSZ/2:PSFStruct.OSZ/2-1));

R=sqrt(gpuArray(XGrid.^2+YGrid.^2));

KPixelSize=1/(PSFStruct.OSZ*PSFStruct.PixelSize);
PupilRadius=(PSFStruct.NA/PSFStruct.Lambda)/KPixelSize;
Mask=R<PupilRadius;
Pupil_Mag=Mask;
Pupil_Phase=gpuArray(zeros(PSFStruct.OSZ,PSFStruct.OSZ,'single'));
PSFStruct.Pupil=cat(3,Pupil_Mag,Pupil_Phase);



