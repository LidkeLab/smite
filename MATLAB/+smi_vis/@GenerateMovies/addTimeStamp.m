function [] = addTimeStamp(PlotAxes, Frame, FrameRate, Params)
%addTimeStamp places a time stamp in a corner of the PlotAxes.
%
% INPUTS:
%   PlotAxes: Axes in which the timestamp is added.
%   Frame: Current frame of the movie.
%   FrameRate: FrameRate used to convert Frame to units of seconds.
%   Params: Structure of several parameters (see
%           smi_vis.GenerateMovies.prepDefaults())

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Add the timestamp, changing units/appearance based on the unit flag.
Time = smi_helpers.arrayMUX({Frame, (Frame-1) / FrameRate}, ...
    Params.UnitFlag);
TimeUnitString = smi_helpers.arrayMUX({'frame', 's'}, Params.UnitFlag);
Timestamp = smi_helpers.arrayMUX(...
    {sprintf('%s %i', TimeUnitString, Time), ...
    sprintf('%.2f %s', Time, TimeUnitString)}, Params.UnitFlag);
text(PlotAxes, double(Params.XPixels(1)), double(Params.YPixels(1)), ...
    double(Params.ZFrames(2)), ...
    Timestamp, 'Color', [1, 1, 1], ...
    'VerticalAlignment', 'top', ...
    'FontSize', 12*(Params.XPixels(2)-Params.XPixels(1))/128)


end