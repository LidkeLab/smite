function success = gaussianImage_unitTest()
%gaussianImage_unitTest tests all functionality of gaussianImage.

success = 0;
fprintf('\nTesting smi_vis.GenerateImages.gaussianImage...\n');

SaveDir = fullfile(tempdir, 'smite', 'unitTest', 'gaussianImage');
if ~isfolder(SaveDir)
   mkdir(fullfile(tempdir, 'smite'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest', 'gaussianImage'));
end

% setting display options
%TrueSize = dipgetpref('TrueSize');
%dipsetpref('TrueSize',true)

% create random input
fprintf('Creating data...\n');
rng('default')
SMD.X = 1 + 98*rand(100000,1);
SMD.Y = 1 + 48*rand(100000,1);
SMD.Z = [];
SMD.X_SE = 0.20*rand(100000,1);
SMD.Y_SE = 0.20*rand(100000,1);
SMD.XSize = 100;
SMD.YSize = 50;
SMD.FrameNum = ones(100000,1);
SMD.Photons = 1000*rand(100000,1);
SMD.Bg = 100*rand(100000,1);
SMD.PixelSize = 0.1;
SRImageZoom = 4;

% test with no output
fprintf('Testing with no output...\n');
smi_vis.GenerateImages.gaussianImage(SMD,SRImageZoom);
saveas(gcf, fullfile(SaveDir, 'gI1.png'));
pause(3)
close all
% test with output
fprintf('Testing with output and all input...\n');
[gaussIm] = smi_vis.GenerateImages.gaussianImage(SMD,SRImageZoom);
imshow(gaussIm)
saveas(gcf, fullfile(SaveDir, 'gI2.png'));
pause(3)
close all

% setting display options back
%dipsetpref('TrueSize',TrueSize)

% finish
fprintf('Done, test successful!\n\n');
success = 1;
end
