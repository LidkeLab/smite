function  [CRLB]=crlbPSFPupil_unitTest()
%crlbPSFPupil_unitTest Tests crlbPSFPupil functionality.
%
% REQUIRES:
%   Parallel Procesing Toolbox
%   NVidia GPU
%

SaveDir = smi_helpers.mkSMITETmpDir('unitTest', 'crlbPSFPupil');

%% Astigmatism PSF
PSFStruct=smi_psf.PointSpreadFunction.createPSFStruct()
PSFStruct.Z=(-1:.05:1);
PSFStruct.ZC_Phase=0;
PSFStruct.ZC_Phase(6)=1;
[PSF,PSFStruct]=smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct); 
sliceViewer(gather(PSF))
[CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil(PSFStruct)
title('Astigmatism')
saveas(gcf, fullfile(SaveDir, 'cPP1.png'));

%% Prasad Spiral
L=2;
PSFStruct=smi_psf.PointSpreadFunction.createPSFStruct()
PSFStruct.Z=(-1:.05:1);
[PSF,PSFStruct]=smi_psf.PointSpreadFunction.scalarPSFPrasadZone(PSFStruct,L);
sliceViewer(gather(PSF))
[CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil(PSFStruct)
title('Prasad Spiral L=1')
saveas(gcf, fullfile(SaveDir, 'cPP2.png'));

%% Tetrapod
PSFStruct=smi_psf.PointSpreadFunction.createPSFStruct()
PSFStruct.Z=(-1:.05:1);
PSFStruct.ZC_Phase(6)=1;
PSFStruct.ZC_Phase(12)=-2;
PSFStruct.Z=(-1:.05:1);
[PSF,PSFStruct]=smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct); 
sliceViewer(gather(PSF))
[CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil(PSFStruct)
title('Tetrapod')
saveas(gcf, fullfile(SaveDir, 'cPP3.png'));

%% Spherical Aberration
PSFStruct=smi_psf.PointSpreadFunction.createPSFStruct()
PSFStruct.ZC_Phase(11)=1;
PSFStruct.Z=(-1:.05:1);
[PSF,PSFStruct]=smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct); 
sliceViewer(gather(PSF))
[CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil(PSFStruct)
title('Spherical Abberation')
saveas(gcf, fullfile(SaveDir, 'cPP4.png'));
