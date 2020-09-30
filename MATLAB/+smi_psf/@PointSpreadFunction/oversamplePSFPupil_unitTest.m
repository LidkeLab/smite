function [Report] = oversamplePFSPupil_unitTest(PSFStruct,Sampling)
%oversamplePFSPupil_unitTest Test and Demonstrate oversamplePFSPupil

Report = 0;

%%
%clc; close all
%Create PSFStruct
P=smi_psf.PointSpreadFunction.createPSFStruct();
P.Z=(-2:.02:2)
P.OTFSigma=[.1,.1];
%Show unsampled PSF
[PSF]=smi_psf.PointSpreadFunction.scalarPSFPupil(P);
dipshow(gather(PSF))
colormap('hot')

[PFS_OS,P_OS]=smi_psf.PointSpreadFunction.oversamplePSFPupil(P,4);
%Show upsampled PSF
dipshow(gather(PFS_OS))
colormap('hot')

Report = 1;

end
