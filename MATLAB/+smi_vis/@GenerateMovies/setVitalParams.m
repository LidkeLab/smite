function [] = setVitalParams(obj)
%setVitalParams ensures vital parameters in obj.Params are set.
% This method will check obj.Params and ensure that the "vital" parameters
% are given useful values (i.e., values that'll let the user create a movie
% without crashing the code).

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set default values for the vital parameters.
DataProvided = ~isempty(obj.RawData);
TrajProvided = ~isempty(obj.TR);
assert(DataProvided || TrajProvided, ...
    ['You must set GenerateMovies.RawData or ', ...
    'GenerateMovies.TR before calling this method.'])
if isempty(obj.Params.ZFrames)
    if (TrajProvided && DataProvided)
        obj.Params.ZFrames = ...
            [1, max(obj.TR(1).NFrames, size(obj.RawData, 3))];
    elseif TrajProvided
        obj.Params.ZFrames = [1, obj.TR(1).NFrames];
    elseif DataProvided
        obj.Params.ZFrames = [1, size(obj.RawData, 3)];
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
end


end