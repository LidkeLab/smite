function success = gaussianImage_unitTest()
%gaussianImage_unitTest tests all functionality of gaussianImage.

success = 0;
fprintf('\nTesting smi_vis.GenerateImages.gaussianImage...\n');

% setting display options
TrueSize = dipgetpref('TrueSize');
dipsetpref('TrueSize',true)

% create random input
fprintf('Creating data...\n');
rng('default')
SMR.X = 1 + 98*rand(100000,1);
SMR.Y = 1 + 48*rand(100000,1);
SMR.Z = [];
SMR.X_SE = 0.20*rand(100000,1);
SMR.Y_SE = 0.20*rand(100000,1);
SMR.XSize = 100;
SMR.YSize = 50;
SMR.FrameNum = ones(100000,1);
SMR.Photons = 1000*rand(100000,1);
SMR.Bg = 100*rand(100000,1);
SMR.PixelSize = 100;
SRImageZoom = 4;

% test with no output
fprintf('Testing with no output...\n');
smi_vis.GenerateImages.gaussianImage(SMR,SRImageZoom);
pause(3)
close all
% test with output
fprintf('Testing with output and all input...\n');
[gaussIm] = smi_vis.GenerateImages.gaussianImage(SMR,SRImageZoom);
dipshow(gaussIm)
pause(3)
close all

% setting display options back
dipsetpref('TrueSize',TrueSize)

% finish
fprintf('Done, test successful!\n\n');
success = 1;
end
