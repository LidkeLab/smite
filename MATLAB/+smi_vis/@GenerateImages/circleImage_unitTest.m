function [Success] =  circleImage_unitTest()
%circleImage_unitTest tests the class method circleImage().
%
% OUTPUTS:
%   Success: A boolean to indicate the success (1) or failure (0) of this
%            unit test.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Seed the random number generator so that the simulated SMD is predictable
% NOTE: If this is changed, there will almost certainly be entries of
%       Success that are 0.
rng(1234)

% Simulate several emitter positions.
% NOTE: Most of these parameters were chosen arbitrarily, HOWEVER, there
%       may be some checks below that might be affected by changing these
%       parameters (i.e., causing some elements of Success to be 0 even if
%       things worked correctly).
FrameSize = 64;
NEmitters = 100;
Coordinates = FrameSize * rand(NEmitters, 2);

% Construct two SMD structures, each with half of the emitters.
NEmittersPerChannel = round(NEmitters / 2);
SMD1 = smi_core.SingleMoleculeData.createSMD();
SMD1.XSize = FrameSize;
SMD1.YSize = FrameSize;
SMD1.X = Coordinates(1:NEmittersPerChannel, 1);
SMD1.Y = Coordinates(1:NEmittersPerChannel, 2);
SMD1.X_SE = 0.1*rand(NEmittersPerChannel, 1) + 0.1;
SMD1.Y_SE = SMD1.X_SE;
SMD2 = smi_core.SingleMoleculeData.createSMD();
SMD2.XSize = FrameSize;
SMD2.YSize = FrameSize;
SMD2.X = Coordinates(NEmittersPerChannel+1:end, 1);
SMD2.Y = Coordinates(NEmittersPerChannel+1:end, 2);
SMD2.X_SE = 0.1*rand(NEmittersPerChannel, 1) + 0.1;
SMD2.Y_SE = SMD2.X_SE;

% Concatenate the two SMDs and try to make a two-color image (I'm doing
% this for the sake of demonstration in case anybody uses this unit test to
% guide their usage of circleImage()).
SMDAll = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2);
ColorMapGreen = repmat([0, 1, 0], NEmittersPerChannel, 1);
ColorMapMagenta = repmat([1, 0, 1], NEmittersPerChannel, 1);
ColorMap = [ColorMapGreen; ColorMapMagenta];
[CircleImage, CircleImageRGB, SRImageZoom] = ...
    smi_vis.GenerateImages.circleImage(SMDAll, ColorMap, [], 16);
imshow(CircleImage);
imshow(CircleImageRGB);
Success = 1;


end
