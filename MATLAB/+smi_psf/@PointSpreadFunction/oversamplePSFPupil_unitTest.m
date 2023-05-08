function [Report] = oversamplePFSPupil_unitTest(PSFStruct,Sampling)
%oversamplePFSPupil_unitTest Test and Demonstrate oversamplePFSPupil.
%
% REQUIRES:
%    DIPimage (https://diplib.org/DIPimage.html)

Report = 0;

SaveDir = fullfile(tempdir, 'smite', 'unitTest', 'oversamplePSFPupil');
if ~isfolder(SaveDir)
   mkdir(fullfile(tempdir, 'smite'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest', 'oversamplePSFPupil'));
end

%%
%clc; close all
%Create PSFStruct
P=smi_psf.PointSpreadFunction.createPSFStruct();
P.Z=(-2:.02:2)
P.OTFSigma=[.1,.1];
%Show unsampled PSF
[PSF]=smi_psf.PointSpreadFunction.scalarPSFPupil(P);
figure;
sliceViewer(gather(PSF));
colormap('hot')
saveas(gcf, fullfile(SaveDir, 'osPP1.png'));

[PFS_OS,P_OS]=smi_psf.PointSpreadFunction.oversamplePSFPupil(P,4);
%Show upsampled PSF
figure;
sliceViewer(gather(PFS_OS));
colormap('hot')
saveas(gcf, fullfile(SaveDir, 'osPP2.png'));

Report = 1;

end
