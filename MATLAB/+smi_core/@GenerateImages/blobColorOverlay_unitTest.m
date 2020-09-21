function success = blobColorOverlay_unitTest()
%histogramImageUnitTest tests all functionality of smi_core.GenerateImages.blobColorOverlay

success = 0;
fprintf('\nTesting smi_core.GenerateImages.blobColorOverlay...\n');

% create sequence and SMD
SZ=256;
NFrames = 100;
NBlobs = 10000;
Background = 5;

SMD.X = (SZ-1)*rand(NBlobs,1);
SMD.Y = (SZ-1)*rand(NBlobs,1);
SMD.Photons=1000+100*randn(NBlobs,1);
SMD.Photons(SMD.Photons<=0) = 0;
SMD.Background=zeros(NBlobs,1); % no background for individual emitters
SMD.PSFSigma=1.3+0.1*randn(NBlobs,1);
SMD.FrameNum = 1+round((NFrames-1)*rand(NBlobs,1));
[~,Sequence] = smi_sim.GaussBlobs.gaussBlobImage(SZ,NFrames,SMD);
Sequence = Sequence+Background; % add background to the whole sequence
% remove particles from SMD to simulated not all got fit
SMD.X = SMD.X(1:NBlobs/2);
SMD.Y = SMD.Y(1:NBlobs/2);
SMD.Photons=SMD.Photons(1:NBlobs/2);
SMD.PSFSigma=SMD.PSFSigma(1:NBlobs/2);
SMD.FrameNum = SMD.FrameNum(1:NBlobs/2);
SMD.Background=5*ones(NBlobs/2,1); % now add background to SMD for overlay

% test with no output
fprintf('Testing with no output...\n');
smi_core.GenerateImages.blobColorOverlay(Sequence,SMD);
fprintf('Check sequence if you like, it will close in 20 sec...\n');
pause(20)
close all

% test with output
fprintf('Testing with output...\n');
[~] = smi_core.GenerateImages.blobColorOverlay(Sequence,SMD);
close all

% finish
fprintf('Done, test successful!\n\n');
success = 1;
end
