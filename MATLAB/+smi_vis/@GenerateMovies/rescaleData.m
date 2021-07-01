function [] = rescaleData(obj)
%rescaleData rescales the raw data in obj.RawData for display purposes.
% This method crops and rescales obj.RawData based on parameters in
% obj.Params.
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
DataSize = size(obj.RawData);
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

% Rescale the raw data after isolating the portion that will be displayed.
RawData = ...
    padarray(obj.RawData, [min(0, YSize-DataSize(1)), ...
    min(0, XSize-DataSize(2)), ...
    min(0, NFrames-DataSize(3))], 'post');
if ~(isempty(obj.Params.XPixels) && isempty(obj.Params.YPixels) ...
        && isempty(obj.Params.ZFrames))
    RawData = RawData(obj.Params.YPixels(1):obj.Params.YPixels(2), ...
        obj.Params.XPixels(1):obj.Params.XPixels(2), ...
        obj.Params.ZFrames(1):obj.Params.ZFrames(2));
end
obj.ScaledData = smi_vis.contrastStretch(single(RawData), [0; 1], ...
    obj.Params.PercentileCeiling, obj.Params.MinScaleIntensity);


end