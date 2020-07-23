
%% Make Data
PSFSigma=1.3;
SZ=[128 256 2000];
Photons=1000;
Data=poissrnd(Photons*single(gaussf(randn(SZ)>2.5,PSFSigma*[1 1 0])));
dipshow(Data)

%% Test Find ROI
SMF=smi_core.SingleMoleculeFitting.createSMF()
SMF.Fitting.PSFSigma=PSFSigma;
SMF.BoxFinding.BoxOverlap=2;
SMF.BoxFinding.BoxSize=5;

FB=smi_core.FindROI(SMF,Data);
tic
[ROIStack,SMD]=FB.findROI();
toc
FB.showOverlay()



%% Test gaussMLE


