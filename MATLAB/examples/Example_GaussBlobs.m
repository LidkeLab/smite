%% Demonstrate uses of GaussBlobs Class

%% smi_sim.GaussBlobs.gaussBlobImage
% Generate an image stack of randomly located blobs

Rho = .001; %blobs per pixel
SZ=256; 
Photons=500;
PSFSigma=1.3;
Bg=0;

%setup our SMD and SMF
SMD=smi_core.SingleMoleculeData.createSMD()
SMF=smi_core.SingleMoleculeFitting()

SMF.Data.DataROI=[1 1 SZ SZ];

%create some random data
SMD.NFrames=1000;
SMD.NDatasets=1;

for nn=1:SMD.NFrames
    N=poissrnd(Rho*SZ*SZ);
    SMD.FrameNum=cat(1,SMD.FrameNum,nn*ones(N,1));
    SMD.X=cat(1,SMD.X,SZ*rand(N,1));
    SMD.Y=cat(1,SMD.Y,SZ*rand(N,1));
    SMD.Photons=cat(1,SMD.Photons,Photons*rand(N,1));
    SMD.Bg=cat(1,SMD.Bg,Bg*ones(N,1));
    SMD.DatasetNum=cat(1,SMD.DatasetNum,nn*ones(N,1));
    SMD.PSFSigma=cat(1,SMD.PSFSigma,PSFSigma*ones(N,1));
end

%call image stack generator
im=smi_sim.GaussBlobs.gaussBlobImage(SMD,SMF);

sliceViewer(im/max(im(:)))

%% smi_sim.GaussBlobs.gaussBlobImage.genRandomBlobImage
% quick generation of a random blob image

B=smi_sim.GaussBlobs.genRandomBlobImage();
sliceViewer(B/max(B(:)))








