function  [CRLB]=crlbPSFPupil_unitTest()
%crlbPSFPupil_unitTest Tests crlbPSFPupil functionality.
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%



%% Astigmatism PSF
PSFStruct=smi_psf.PointSpreadFunction.createPSFStruct()
PSFStruct.Z=(-1:.05:1);
PSFStruct.ZC_Phase=0;
PSFStruct.ZC_Phase(6)=1;
[PSF,PSFStruct]=smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct); 
imshow(gather(PSF))
[CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil(PSFStruct)
title('Astigmatism')

%% Prasad Spiral
L=2;
PSFStruct=smi_psf.PointSpreadFunction.createPSFStruct()
PSFStruct.Z=(-1:.05:1);
[PSF,PSFStruct]=smi_psf.PointSpreadFunction.scalarPSFPrasadZone(PSFStruct,L);
imshow(gather(PSF))
[CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil(PSFStruct)
title('Prasad Spiral L=1')

%% Tetrapod
PSFStruct=smi_psf.PointSpreadFunction.createPSFStruct()
PSFStruct.Z=(-1:.05:1);
PSFStruct.ZC_Phase(6)=1;
PSFStruct.ZC_Phase(12)=-2;
PSFStruct.Z=(-1:.05:1);
[PSF,PSFStruct]=smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct); 
imshow(gather(PSF))
[CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil(PSFStruct)
title('Tetrapod')

%% Spherical Aberration
PSFStruct=smi_psf.PointSpreadFunction.createPSFStruct()
PSFStruct.ZC_Phase(11)=1;
PSFStruct.Z=(-1:.05:1);
[PSF,PSFStruct]=smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct); 
imshow(gather(PSF))
[CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil(PSFStruct)
title('Spherical Abberation')
