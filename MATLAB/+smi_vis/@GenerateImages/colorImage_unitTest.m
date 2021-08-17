function success = colorImage_unitTest()
%colorImage_unitTest tests all functionality of colorImage.

success = 0;
fprintf('\nTesting smi_vis.GenerateImages.colorImage...\n');
fprintf('Creating image...\n');

% create image
SZ = 256;
NFrames = 1;
Data = smi_sim.GaussBlobs.genRandomBlobImage(SZ,NFrames);

% inputs
ColorMap = parula(256);
sortedVals = sort(Data(:));
NumEl = numel(sortedVals);
MinMax = [sortedVals(round(NumEl*0.1)),sortedVals(round(NumEl*0.9))];

% test with only image input
fprintf('Testing with image as only input...\n');
[RGBim]=smi_vis.GenerateImages.colorImage(Data);
imshow(RGBim);
pause(3);
close all
% test with image and colormap input
fprintf('Testing with image and colormap input...\n');
[RGBim]=smi_vis.GenerateImages.colorImage(Data,ColorMap);
imshow(RGBim)
pause(3)
close all
% test without colormap (default should be hot)
fprintf('Testing with image, colormap and minmax input...\n');
[RGBim]=smi_vis.GenerateImages.colorImage(Data,ColorMap,MinMax);
imshow(RGBim)
pause(3)
close all

% finish
fprintf('Done, test successful!\n\n');
success = 1;
end
