function success = blobColorOverlay_unitTest()
%blobColorOverlay_unitTest tests all functionality of blobColorOverlay.

success = 0;
fprintf('\nTesting smi_vis.GenerateImages.blobColorOverlay...\n');

SaveDir = smi_helpers.mkSMITETmpDir('unitTest', 'blobColorOverlay');

% create sequence and SMD
SZ=256;
NFrames = 100;
NBlobs = 10000;
Background = 5;

SMD.X = (SZ-1)*rand(NBlobs,1);
SMD.Y = (SZ-1)*rand(NBlobs,1);
SMD.Z = [];
SMD.Photons=1000+100*randn(NBlobs,1);
SMD.Photons(SMD.Photons<=0) = 0;
SMD.Background=zeros(NBlobs,1); % no background for individual emitters
SMD.PSFSigma=1.3+0.1*randn(NBlobs,1);
SMD.FrameNum = 1+round((NFrames-1)*rand(NBlobs,1));
SMD.Bg=5*ones(NBlobs,1); % now add background to SMD for overlay
SMD.XSize = SZ;
SMD.YSize = SZ;
SMD.NFrames = NFrames;
[~,Sequence] = smi_sim.GaussBlobs.gaussBlobImage(SMD);

Sequence = Sequence+Background; % add background to the whole sequence
% remove particles from SMD to simulate not all getting fit
SMD.X = SMD.X(1:NBlobs/2);
SMD.Y = SMD.Y(1:NBlobs/2);
SMD.Photons=SMD.Photons(1:NBlobs/2);
SMD.PSFSigma=SMD.PSFSigma(1:NBlobs/2);
SMD.FrameNum = SMD.FrameNum(1:NBlobs/2);
SMD.Bg=5*ones(NBlobs/2,1); % now add background to SMD for overlay

% test with no output
%fprintf('Testing with no output...\n');
%smi_vis.GenerateImages.blobColorOverlay(Sequence,SMD);
%fprintf('Check sequence if you like, it will close in 20 sec...\n');
%pause(20)
%close all

% test with output
fprintf('Testing with output...\n');
[BlobImage] = smi_vis.GenerateImages.blobColorOverlay(Sequence,SMD);
fprintf('Check sequence if you like, it will close in 20 sec...\n');
sliceViewer(BlobImage);
saveas(gcf, fullfile(SaveDir, 'BlobImage.png'));
pause(20)
close all

% finish
fprintf('Done, test successful!\n\n');
success = 1;
end
