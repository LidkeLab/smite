function [Params] = defineCropROI(TR, Params)
%defineCropOI defines the ROI surrounding the data for an auto-cropping.
% This method defines a ROI encompassing the data in TR based on some of
% the parameters in Params (e.g., enforcing a minimum ROI size, boundary
% padding, etc.).
%
% INPUTS:
%   TR: Tracking Results structure containing trajectories.
%   Params: Structure of movie parameters (see
%           smi_vis.GenerateMovies.prepDefaults())
%
% OUTPUTS:
%   Params: Input 'Params' with XPixels, YPixels, and ZFrames (possibly)
%           being updated to define an auto-cropping ROI around the data in
%           'TR'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Enforce the temporal cropping and padding.
AllFrames = cell2mat({TR.FrameNum}.');
MinFrame = smi_helpers.arrayMUX({min(AllFrames), 1}, isempty(AllFrames));
MaxFrame = smi_helpers.arrayMUX({max(AllFrames), TR(1).NFrames}, ...
    isempty(AllFrames));
Params.ZFrames = [max(1, MinFrame-Params.NPadFrames), ...
    min(TR(1).NFrames, MaxFrame+Params.NPadFrames)];

% Enforce the min. XY range and the XY padding.
IsDisplayed = ...
    ismember(AllFrames, Params.ZFrames(1):Params.ZFrames(2));
AllX = cell2mat({TR.X}.');
MinX = smi_helpers.arrayMUX({min(AllX(IsDisplayed)), 1}, isempty(AllX));
MaxX = smi_helpers.arrayMUX({max(AllX(IsDisplayed)), TR(1).XSize}, ...
    isempty(AllX));
XPixels = [max(1, MinX-Params.NPadPixels), ...
    min(TR(1).XSize, MaxX+Params.NPadPixels)];
XWidth = max(diff(XPixels), Params.MinXYRange);
XCenterIdeal = mean(XPixels);
XStart = max(1, floor(min(XCenterIdeal-XWidth/2, TR(1).XSize-XWidth)));
XEnd = min(TR(1).XSize, ceil(XStart+XWidth));
Params.XPixels = [XStart, XEnd];
AllY = cell2mat({TR.Y}.');
MinY = smi_helpers.arrayMUX({min(AllY(IsDisplayed)), 1}, isempty(AllY));
MaxY = smi_helpers.arrayMUX({max(AllY(IsDisplayed)), TR(1).YSize}, ...
    isempty(AllY));
YPixels = [max(1, MinY-Params.NPadPixels), ...
    min(TR(1).YSize, MaxY+Params.NPadPixels)];
YWidth = max(diff(YPixels), Params.MinXYRange);
YCenterIdeal = mean(YPixels);
YStart = max(1, floor(min(YCenterIdeal-YWidth/2, TR(1).YSize-YWidth)));
YEnd = min(TR(1).YSize, ceil(YStart+YWidth));
Params.YPixels = [YStart, YEnd];


end