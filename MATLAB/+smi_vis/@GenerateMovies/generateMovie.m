function generateMovie(obj)
%generateMovie generates a movie of 2D trajectories over time.
% This method is intended to be a wrapper around playMovie() which will add
% some "decorations" (e.g., axis labels), prepare saved movies, etc.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Ensure some necessary parameters are given meaningful values based on the
% provided data/trajectories.
if isempty(obj.Params.ZFrames)
    obj.Params.ZFrames = ceil([1, max(cell2mat({obj.TR.FrameNum}.'))]);
end
if isempty(obj.Params.XPixels)
    obj.Params.XPixels = ceil([1, max(cell2mat({obj.TR.X}.'))]);
end
if isempty(obj.Params.YPixels)
    obj.Params.YPixels = ceil([1, max(cell2mat({obj.TR.Y}.'))]);
end
if isempty(obj.Params.TrajColor)
    obj.Params.TrajColor = lines(numel(obj.TR));
end

% Add some decorations to the movie axes.
xlabel(obj.MovieAxes, sprintf('X (%s)', obj.LengthUnitString), ...
    'Interpreter', 'Latex')
ylabel(obj.MovieAxes, sprintf('Y (%s)', obj.LengthUnitString), ...
    'Interpreter', 'Latex')
zlabel(obj.MovieAxes, ...
    sprintf('%s (%s)', obj.TimeDimensionString, obj.TimeUnitString), ...
    'Interpreter', 'Latex')

% Play the movie.
obj.playMovie(obj.MovieAxes, obj.TR, obj.RawData, ...
    obj.Params, obj.SMF, obj.SMD, obj.VideoObject)


end