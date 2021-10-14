function [] = setVitalParams(obj)
%setVitalParams ensures vital parameters in obj.Params are set.
% This method will check obj.Params and ensure that the "vital" parameters
% are given useful values (i.e., values that'll let the user create a movie
% without crashing the code).

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set default values for the vital parameters.
DataProvided = ~isempty(obj.RawData);
TrajProvided = ~isempty(obj.TR.FrameNum);
assert(DataProvided || TrajProvided, ...
    ['You must set GenerateMovies.RawData or ', ...
    'GenerateMovies.TR before calling this method.'])
DataSize = size(obj.RawData, 1:4);
if isempty(obj.Params.ZFrames)
    NFramesData = (DataSize(4)>1)*DataSize(4) ...
        + (DataSize(4)==1)*DataSize(3);
    if (TrajProvided && DataProvided)
        obj.Params.ZFrames = ...
            [1, max(obj.TR(1).NFrames, NFramesData)];
    elseif TrajProvided
        obj.Params.ZFrames = [1, obj.TR(1).NFrames];
    elseif DataProvided
        obj.Params.ZFrames = [1, NFramesData];
    end
end
if isempty(obj.Params.XPixels)
    if (TrajProvided && DataProvided)
        obj.Params.XPixels = ...
            [1, max(obj.TR(1).XSize, size(obj.RawData, 2))];
    elseif TrajProvided
        obj.Params.XPixels = [1, obj.TR(1).XSize];
    elseif DataProvided
        obj.Params.XPixels = [1, size(obj.RawData, 2)];
    end
end
if isempty(obj.Params.YPixels)
    if (TrajProvided && DataProvided)
        obj.Params.YPixels = ...
            [1, max(obj.TR(1).YSize, size(obj.RawData, 1))];
    elseif TrajProvided
        obj.Params.YPixels = [1, obj.TR(1).YSize];
    elseif DataProvided
        obj.Params.YPixels = [1, size(obj.RawData, 1)];
    end
end
if isempty(obj.Params.TrajColor)
    obj.Params.TrajColor = lines(numel(obj.TR));
else
    % Make sure TrajColor has as many rows as the number of trajectories in
    % obj.TR.
    NTraj = numel(obj.TR);
    NColors = size(obj.Params.TrajColor, 1);
    obj.Params.TrajColor = cat(1, obj.Params.TrajColor, ...
        repmat(obj.Params.TrajColor(end, :), ...
        [max(NTraj-NColors, 1), 1, 1])); 
end
if isempty(obj.Params.SMDColor)
    obj.Params.SMDColor = lines(numel(obj.SMD.FrameNum));
else
    % Make sure SMDColor has as many rows as the number of datapoints in
    % obj.SMD
    NData = numel(obj.SMD);
    NColors = size(obj.Params.SMDColor, 1);
    obj.Params.SMDColor = cat(1, obj.Params.SMDColor, ...
        repmat(obj.Params.SMDColor(end, :), ...
        [max(NData-NColors, 1), 1, 1])); 
end


end