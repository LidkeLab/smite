function success = histogramImage_unitTest()
%histogramImage_unitTest tests all functionality of histogramImage.

success = 0;
fprintf('\nTesting smi_vis.GenerateImages.histogramImage...\n');

SaveDir = smi_helpers.mkSMITETmpDir('unitTest', 'histogramImage');

% setting display options
%TrueSize = dipgetpref('TrueSize');
%dipsetpref('TrueSize',true)

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
smi_vis.GenerateImages.histogramImage(SMR,SRImageZoom,ColorMap);
saveas(gcf, fullfile(SaveDir, 'hI1.png'));
pause(3)
close all
% test with output
fprintf('Testing with output and all input...\n');
[histIm,RGBim] = smi_vis.GenerateImages.histogramImage(SMR,SRImageZoom,ColorMap);
imshow(histIm)
saveas(gcf, fullfile(SaveDir, 'hI2.png'));
h = imshow(RGBim);
saveas(gcf, fullfile(SaveDir, 'hI3.png'));
% DIPimage:
%pos = h.Position;
%pos(2) = pos(2)-300;
%h.Position = pos;
pause(3)
close all
% test without colormap (default should be hot)
fprintf('Testing with output and no colormap input...\n');
[histIm,RGBim] = smi_vis.GenerateImages.histogramImage(SMR,SRImageZoom);
imshow(histIm)
saveas(gcf, fullfile(SaveDir, 'hI4.png'));
h = imshow(RGBim);
saveas(gcf, fullfile(SaveDir, 'hI5.png'));
% DIPimage:
%pos = h.Position;
%pos(2) = pos(2)-300;
%h.Position = pos;
pause(3)
close all

% setting display options back
%dipsetpref('TrueSize',TrueSize)

% finish
fprintf('Done, test successful!\n\n');
success = 1;
end
