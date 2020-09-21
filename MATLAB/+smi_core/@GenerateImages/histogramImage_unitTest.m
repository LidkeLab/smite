function success = histogramImage_unitTest()
%histogramImageUnitTest tests all functionality of smi_core.GenerateImages.histogramImage

success = 0;
fprintf('\nTesting smi_core.GenerateImages.histogramImage...\n');

% setting display options
TrueSize = dipgetpref('TrueSize');
dipsetpref('TrueSize',true)

% create random input
fprintf('Creating data...\n');
rng('default')
SMR.X = 1 + 98*rand(100000,1);
SMR.Y = 1 + 48*rand(100000,1);
SMR.XSize = 100;
SMR.YSize = 50;
SRImageZoom = 4;
ColorMap = 'jet';

% test with no output
fprintf('Testing with no output...\n');
smi_core.GenerateImages.histogramImage(SMR,SRImageZoom,ColorMap);
pause(3)
close all
% test with output
fprintf('Testing with output and all input...\n');
[histIm,RGBim] = smi_core.GenerateImages.histogramImage(SMR,SRImageZoom,ColorMap);
dipshow(histIm)
h = dipshow(RGBim);
pos = h.Position;
pos(2) = pos(2)-300;
h.Position = pos;
pause(3)
close all
% test without colormap (default should be hot)
fprintf('Testing with output and no colormap input...\n');
[histIm,RGBim] = smi_core.GenerateImages.histogramImage(SMR,SRImageZoom);
dipshow(histIm)
h = dipshow(RGBim);
pos = h.Position;
pos(2) = pos(2)-300;
h.Position = pos;
pause(3)
close all

% setting display options back
dipsetpref('TrueSize',TrueSize)

% finish
fprintf('Done, test successful!\n\n');
success = 1;
end
