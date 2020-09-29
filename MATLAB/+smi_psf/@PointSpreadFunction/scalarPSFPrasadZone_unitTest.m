function [Report]=scalarPSFPrasadZone_unitTest()
%scalarPSFPrasadZone_unitTest tests scalarPSFPrasadZone functionality.

%%
L=5;
Photons=500;
Bg=5;

close all;clc;
P=smi_psf.PointSpreadFunction.createPSFStruct();
P.Z=(-10:.5:10);
P.SZ=256;
P.OSZ=512;
[PSF,P]=smi_psf.PointSpreadFunction.scalarPSFPrasadZone(P,L);

dipshow(gather(PSF));
dipshow(gather(P.Pupil(:,:,2)));
smi_psf.PointSpreadFunction.crlbPSFPupil(P,Photons,Bg);
smi_psf.PointSpreadFunction.crlbPSFPupil(smi_psf.PointSpreadFunction.createPSFStruct(),Photons,Bg)

end
