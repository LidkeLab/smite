function [] = prepRawData(obj)
%prepRawData crops and rescales the movie data for display purposes.
% This method crops and rescales obj.RawData based on parameters in
% obj.Params.  Furthermore, this method will also revise some fields of
% obj.Params to ensure obj.Params.AutoClip is being enforced.
%
% REQUIRES:
%   Image Processing Toolbox (to use padarray())

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Check if the flag obj.DataIsScaled is set.  If it is, we can exit early.
if obj.DataIsPrepped
    return
end

% If needed, reshape the raw data.
DataSize = size(obj.RawData, 1:4);
if ((DataSize(4)==1) && (DataSize(3)~=3))
    obj.RawData = reshape(obj.RawData, [DataSize(1:2), 1, DataSize(3)]);
    obj.RawData = repmat(obj.RawData, [1, 1, 3, 1]);
    DataSize = size(obj.RawData);
end
if isempty(obj.RawData)
    obj.setVitalParams()
    obj.RawData = zeros(obj.Params.YPixels(2), ...
        obj.Params.XPixels(2), ...
        3, ...
        obj.Params.ZFrames(2));
    obj.ScaledData = zeros(obj.Params.YPixels(2), ...
        obj.Params.XPixels(2), ...
        3, ...
        obj.Params.ZFrames(2));
    return
end

% Grab some parameters from 'TRInternal' which may be empty.  If they are 
% empty, we'll need to come up with reasonable guesses for these values.
NFrames = obj.TRInternal(1).NFrames;
XSize = obj.TRInternal(1).XSize;
YSize = obj.TRInternal(1).YSize;
if isempty(NFrames)
    NFrames = DataSize(4);
end
if isempty(XSize)
    XSize = DataSize(2);
end
if isempty(YSize)
    YSize = DataSize(1);
end
[obj.TRInternal.NFrames] = deal(NFrames);
[obj.TRInternal.XSize] = deal(XSize);
[obj.TRInternal.YSize] = deal(YSize);

% Make sure other vital parameters are set in obj.Params.
obj.setVitalParams()

% Enforce obj.Params.AutoCrop if needed.
if obj.Params.AutoCrop
    obj.Params = obj.defineCropROI(obj.TRInternal, obj.Params);
end

% Pad the raw data to ensure the cropping region is valid.
RawData = single(padarray(obj.RawData, ...
    [max(0, YSize-DataSize(1)), ...
    max(0, XSize-DataSize(2)), ...
    0, ...
    max(0, NFrames-DataSize(4))], 'post'));

% Crop the raw data to the specified ROI.
RawData = obj.cropRawData(RawData, obj.Params);

% Rescale the cropped data, treating each color channel independently.
NColorChannels = 3;
obj.ScaledData = zeros(size(RawData));
for cc = 1:NColorChannels
    obj.ScaledData(:, :, cc, :) = ...
        smi_vis.contrastStretch(RawData(:, :, cc, :), [0; 1], ...
        obj.Params.PercentileCeiling, obj.Params.PercentileFloor, ...
        obj.Params.MinScaleIntensity);
end
obj.DataIsPrepped = true;


end