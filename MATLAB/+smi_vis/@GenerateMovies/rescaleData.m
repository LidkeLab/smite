function [] = rescaleData(obj)
%rescaleData crops and rescales the movie data for display purposes.
% This method crops and rescales obj.RawData based on parameters in
% obj.Params.  Furthermore, this method will also revise some fields of
% obj.Params to ensure obj.Params.AutoClip is being enforced.
%
% REQUIRES:
%   Image Processing Toolbox (to use padarray())

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Grab some parameters from 'TR' which may be empty.  If they are empty,
% we'll need to come up with reasonable guesses for these values.
NFrames = obj.TR(1).NFrames;
XSize = obj.TR(1).XSize;
YSize = obj.TR(1).YSize;
DataSize = size(obj.RawData, 1:3);
if isempty(NFrames)
    NFrames = DataSize(3);
end
if isempty(XSize)
    XSize = DataSize(2);
end
if isempty(YSize)
    YSize = DataSize(1);
end

% Make sure other vital parameters are set in obj.Params.
obj.setVitalParams()

% Enforce obj.Params.AutoClip if needed.
if obj.Params.AutoClip
    % Enforce the temporal clipping and padding.
    AllFrames = cell2mat({obj.TR.FrameNum}.');
    MinFrame = isempty(AllFrames) + ~isempty(AllFrames)*min(AllFrames);
    MaxFrame = isempty(AllFrames)*NFrames ...
        + ~isempty(AllFrames)*max(AllFrames);
    obj.Params.ZFrames = [max(1, MinFrame-obj.Params.NPadFrames), ...
        min(NFrames, MaxFrame+obj.Params.NPadFrames)];
    
    % Enforce the min. XY range and the XY padding.
    IsDisplayed = ...
        ismember(AllFrames, obj.Params.ZFrames(1):obj.Params.ZFrames(2));
    AllX = cell2mat({obj.TR.X}.');
    MinX = isempty(AllX) + ~isempty(AllX)*min(AllX(IsDisplayed));
    MaxX = isempty(AllX)*XSize + ~isempty(AllX)*max(AllX(IsDisplayed));
    XPixels = [max(1, MinX-obj.Params.NPadPixels), ...
        min(XSize, MaxX+obj.Params.NPadPixels)];
    XWidth = max(diff(XPixels), obj.Params.MinXYRange);
    XCenterIdeal = mean(XWidth);
    XStart = max(1, floor(min(XCenterIdeal-XWidth/2, XSize-XWidth)));
    XEnd = min(XSize, ceil(XStart+XWidth));
    obj.Params.XPixels = [XStart, XEnd];
    AllY = cell2mat({obj.TR.Y}.');
    MinY = isempty(AllY) + ~isempty(AllY)*min(AllY(IsDisplayed));
    MaxY = isempty(AllY)*YSize + ~isempty(AllY)*max(AllY(IsDisplayed));
    YPixels = [max(1, MinY-obj.Params.NPadPixels), ...
        min(YSize, MaxY+obj.Params.NPadPixels)];
    YWidth = max(diff(YPixels), obj.Params.MinXYRange);
    YCenterIdeal = mean(YWidth);
    YStart = max(1, floor(min(YCenterIdeal-YWidth/2, YSize-YWidth)));
    YEnd = min(YSize, ceil(YStart+YWidth));
    obj.Params.YPixels = [YStart, YEnd];
end

% Rescale the raw data after isolating the portion that will be displayed.
RawData = ...
    padarray(obj.RawData, [max(0, YSize-DataSize(1)), ...
    max(0, XSize-DataSize(2)), ...
    max(0, NFrames-DataSize(3))], 'post');
if ~(isempty(obj.Params.XPixels) && isempty(obj.Params.YPixels) ...
        && isempty(obj.Params.ZFrames))
    RawData = RawData(obj.Params.YPixels(1):obj.Params.YPixels(2), ...
        obj.Params.XPixels(1):obj.Params.XPixels(2), :);
end
obj.ScaledData = smi_vis.contrastStretch(single(RawData), [0; 1], ...
    obj.Params.PercentileCeiling, obj.Params.MinScaleIntensity);


end