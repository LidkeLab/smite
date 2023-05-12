function [Report]=scalarPSFPrasadZone_unitTest()
%scalarPSFPrasadZone_unitTest Test scalarPSFPrasadZone functionality.

Report = 0;

SaveDir = smi_helpers.mkSMITETmpDir('unitTest', 'scalarPSFPrasadZone');

L=5;
Photons=500;
Bg=5;

%close all;clc;
P=smi_psf.PointSpreadFunction.createPSFStruct();
P.Z=(-10:.5:10);
P.SZ=256;
P.OSZ=512;
[PSF,P]=smi_psf.PointSpreadFunction.scalarPSFPrasadZone(P,L);

figure; sliceViewer(gather(PSF));
saveas(gcf, fullfile(SaveDir, 'sPPZ1.png'));
figure; imshow(gather(P.Pupil(:,:,2)));
saveas(gcf, fullfile(SaveDir, 'sPPZ2.png'));
smi_psf.PointSpreadFunction.crlbPSFPupil(P,Photons,Bg);
saveas(gcf, fullfile(SaveDir, 'sPPZ3.png'));
smi_psf.PointSpreadFunction.crlbPSFPupil(smi_psf.PointSpreadFunction.createPSFStruct(),Photons,Bg)
saveas(gcf, fullfile(SaveDir, 'sPPZ4.png'));

Report = 1;

end
