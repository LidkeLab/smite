function success = driftImage_unitTest()
%driftImage_unitTest tests all functionality of driftImage.

success = 0;
fprintf('\nTesting smi_core.GenerateImages.driftImage...\n');

% create random input
fprintf('Creating data...\n');
rng('default')
pointsPerPoint = 400;
X = 1 + 37*rand(pointsPerPoint,1);
Y = 1 + 27*rand(pointsPerPoint,1);
% simulate drift
SMR.X = zeros(10*pointsPerPoint,1);
SMR.Y = zeros(10*pointsPerPoint,1);
SMR.XSize = 40;
SMR.YSize = 30;
SMR.DatasetNum = zeros(10*pointsPerPoint,1);
SMR.FrameNum = zeros(size(SMR.X));
for ii = 1 : 10
    SMR.X(((ii-1)*pointsPerPoint)+1 : ii*pointsPerPoint) = X + (ii-1)/10;
    SMR.Y(((ii-1)*pointsPerPoint)+1 : ii*pointsPerPoint) = Y + (ii-1)/10;
    SMR.DatasetNum(((ii-1)*pointsPerPoint)+1 : ii*pointsPerPoint) = ii;
    SMR.FrameNum(((ii-1)*pointsPerPoint)+1 : ii*pointsPerPoint) = 1:1:pointsPerPoint;
end
% input parameters
SRImageZoom = 4;

% test with no output
fprintf('Testing with no output...\n');
smi_core.GenerateImages.driftImage(SMR,SRImageZoom);
diptruesize(gcf,200);
pause(3);
close gcf
% test with single output variable
fprintf('Testing with single output variable...\n');
[driftIm] = smi_core.GenerateImages.driftImage(SMR,SRImageZoom);
h = dipshow(driftIm);
diptruesize(h,200);
pause(3)
close(h)
% test with 2 output variables
fprintf('Testing with two output variables...\n');
[driftIm,driftImRGB] = smi_core.GenerateImages.driftImage(SMR,SRImageZoom);
h1 = dipshow(driftIm);
diptruesize(h1,200);
h2 = dipshow(driftImRGB);
diptruesize(h2,200);
pos = h2.Position;
pos(2) = pos(2)-300;
h2.Position = pos;
pause(3)
close(h1,h2)

% finish
fprintf('Done, test successful!\n\n');
success = 1;
end
