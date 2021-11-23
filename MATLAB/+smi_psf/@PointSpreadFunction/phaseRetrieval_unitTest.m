function [PSFStruct]=phaseRetrieval_unitTest()
%phaseRetrieval_unitTest Tests phaseRetrieval using GS Algorithm.

%%
%close all

%Create a noise defocus stack with offset
P=smi_psf.PointSpreadFunction.createPSFStruct();
P.OSZ=64*3; %Factor 3 to match padding done in PR
P.ZC_Phase(6)=6;
P.ZC_Phase(12)=-2;
P.ZC_Mag(4)=-.2
P.Z=(-1:.2:1);
[PSF,P]=smi_psf.PointSpreadFunction.scalarPSFZernike(P);

Data=poissrnd(gather(PSF)*100000+10)+100;
%Data=Data+shift(Data,[10 0 0])


%Do phase retrieval and compare data stack
P_Out=smi_psf.PointSpreadFunction.phaseRetrieval(P,Data)
%dipshow(gather(P_Out.Pupil))
Model=gather(smi_psf.PointSpreadFunction.scalarPSFPupil(P_Out));

%dipshow(gather(Data))
%dipshow(gather(Model))
joinchannels('RGB',stretch(gather(Data)),stretch(Model))


%Expand into Zernike modes and see if we get correct phase back
KPixelSize=1/(size(P_Out.Pupil,1)*P_Out.PixelSize);
PupilRadius=(P_Out.NA/P_Out.Lambda)/KPixelSize;

ZS=smi_psf.PointSpreadFunction.createZernikeStruct(P_Out.OSZ,PupilRadius,21);
P_Out.ZC_Phase=gather(smi_psf.PointSpreadFunction.zernikeExpansion(P.Pupil(:,:,2),ZS));
P_Out.ZC_Phase


%% Trying with some real data

close all
load('C:\Data\Data_FastVarLaser2-2016-8-3-12-23-8.mat')
SZ=128;
%%dipshow(1,CleanData);
CNR=[141 107];
Data=cut(CleanData,[SZ SZ size(CleanData,3)],[CNR-SZ/2 0]);

P=smi_psf.PointSpreadFunction.createPSFStruct();
P.SZ=SZ;
P.NA=1.4;
P.N=1.52;
P.Lambda=.69;
P.Z=ZVector;
DataIn=Data;
P_Out=smi_psf.PointSpreadFunction.phaseRetrieval(P,DataIn);

%dipshow(gather(P_Out.Pupil))
Model=gather(smi_psf.PointSpreadFunction.scalarPSFPupil(P_Out));

KPixelSize=1/(size(P_Out.Pupil,1)*P_Out.PixelSize);
PupilRadius=(P_Out.NA/P_Out.Lambda)/KPixelSize;
ZS=smi_psf.PointSpreadFunction.createZernikeStruct(P_Out.OSZ,PupilRadius,49);
P_Out.ZC_Phase=gather(smi_psf.PointSpreadFunction.zernikeExpansion(P_Out.Pupil(:,:,2),ZS));
P_Out.ZC_Mag=gather(smi_psf.PointSpreadFunction.zernikeExpansion(P_Out.Pupil(:,:,1),ZS));
[PSFZ,P_OutZ]=smi_psf.PointSpreadFunction.scalarPSFZernike(P_Out);
%dipshow(gather(P_OutZ.Pupil))
%dipshow(gather(DataIn))
%dipshow(gather(PSFZ))
P_OutZ.ZC_Phase

%%


close all
load('C:\Data\Data_FastVarLaser2-2016-8-3-12-23-8.mat')
SZ=128;
%%dipshow(1,CleanData);
CNR=[141 107];
Data=cut(CleanData,[SZ SZ size(CleanData,3)],[CNR-SZ/2 0]);

P=smi_psf.PointSpreadFunction.createPSFStruct();
P.SZ=SZ;
P.NA=1.48;
P.N=1.52;
P.Lambda=.69;
P.Z=ZVector([1,7]);
DataIn=Data(:,:,[1,6]);
P_Out=smi_psf.PointSpreadFunction.phaseRetrievalEM(P,DataIn);
P_Out=smi_psf.PointSpreadFunction.phaseRetrieval(P_Out,DataIn);

%dipshow(gather(P_Out.Pupil))
Model=gather(smi_psf.PointSpreadFunction.scalarPSFPupil(P_Out));

KPixelSize=1/(size(P_Out.Pupil,1)*P_Out.PixelSize);
PupilRadius=(P_Out.NA/P_Out.Lambda)/KPixelSize;
ZS=smi_psf.PointSpreadFunction.createZernikeStruct(P_Out.OSZ,PupilRadius,49);
P_Out.ZC_Phase=gather(smi_psf.PointSpreadFunction.zernikeExpansion(P_Out.Pupil(:,:,2),ZS));
P_Out.ZC_Mag=gather(smi_psf.PointSpreadFunction.zernikeExpansion(P_Out.Pupil(:,:,1),ZS));
[PSFZ,P_OutZ]=smi_psf.PointSpreadFunction.scalarPSFZernike(P_Out);
%dipshow(gather(P_OutZ.Pupil))
%dipshow(gather(DataIn))
%dipshow(gather(PSFZ))
P_OutZ.ZC_Phase

PSFStruct = P;
