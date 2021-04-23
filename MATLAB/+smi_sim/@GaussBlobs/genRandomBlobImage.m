function BlobStack=genRandomBlobImage(SZ,NFrames,Rho,Photons,PSFSigma,Bg)
%genRandomBlobImage Generate a stack of images with randomly placed blobs
%
% INPUTS: 
%   SZ:         Y,X Size of image. Scalar or [1 x 2]. (Default=256)
%   NFrames:    Number of frames. (Default=1000)
%   Rho:        Density of blobs (blobs/pixel) (Default = .001) 
%   Photons:    Photons per blob (Default = 1000)
%   PSFSigma:   Blob 2D Sigma (Default = 1)
%   Bg:         Background (photons/pixel) (Default = 0)
%
% OUTPUTS: 
%   BlobStack:  Output image stack. SZ x SZ x NFrames

if nargin < 1
    SZ=256; 
end
if nargin < 2
    NFrames=1000;
end
if nargin < 3
    Rho=.001;
end
if nargin < 4
    Photons=1000;
end
if nargin < 5
    PSFSigma=1;
end
if nargin < 6
    Bg=0;
end

%setup our SMD and SMF
SMD=smi_core.SingleMoleculeData.createSMD();
SMF=smi_core.SingleMoleculeFitting();

if length(SZ)==2
    SMF.Data.DataROI=[1 1 SZ];
else   
    SMF.Data.DataROI=[1 1 SZ SZ];
end

%create some random data
SMD.NFrames=NFrames;
SMD.NDatasets=1;

for nn=1:SMD.NFrames
    N=poissrnd(Rho*SZ*SZ);
    SMD.FrameNum=cat(1,SMD.FrameNum,nn*ones(N,1));
    SMD.X=cat(1,SMD.X,SZ*rand(N,1));
    SMD.Y=cat(1,SMD.Y,SZ*rand(N,1));
    SMD.Photons=cat(1,SMD.Photons,Photons*ones(N,1));
    SMD.Bg=cat(1,SMD.Bg,Bg*ones(N,1));
    SMD.DatasetNum=cat(1,SMD.DatasetNum,nn*ones(N,1));
    SMD.PSFSigma=cat(1,SMD.PSFSigma,PSFSigma*ones(N,1));
end

BlobStack=smi_sim.GaussBlobs.gaussBlobImage(SMD,SMF);


end