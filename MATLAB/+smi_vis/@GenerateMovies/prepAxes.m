function [] = prepAxes(obj, PlotAxes)
%prepAxes prepares some features of the movie axes.
% This method does a few things to the axes 'PlotAxes' based on fields 
% in obj.Params (e.g., setting the axis limits).
%
% INPUTS:
%   PlotAxes: Axes that will be modified. (Default = obj.MovieAxes)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Before proceeding, check if the axes are in need of an update.
if obj.AxesPrepped
    return
end

% Set defaults if needed.
if (~exist('PlotAxes', 'var') || isempty(PlotAxes) || ~isvalid(PlotAxes))
    PlotAxes = obj.MovieAxes;
end

% Revise several properties of the axes.
hold(PlotAxes, 'on')
PlotAxes.ActivePositionProperty = 'position';
PlotAxes.DataAspectRatio = [1, 1, ...
    max(obj.Params.ZFrames)/max([obj.Params.XPixels, obj.Params.YPixels])];
PlotAxes.YDir = 'reverse';
PlotAxes.XLimMode = 'manual';
PlotAxes.YLimMode = 'manual';
PlotAxes.ZLimMode = 'manual';
PlotAxes.XLim = obj.Params.XPixels + [0, 1];
PlotAxes.YLim = obj.Params.YPixels + [0, 1];
PlotAxes.ZLim = obj.Params.ZFrames;
view(PlotAxes, obj.Params.LineOfSite(1, :));
colormap(PlotAxes, 'gray')

% Revise axes ticks based on the unit flag.
obj.addAxesTicks(PlotAxes, 5, '%.1f')

% Add some decorations to the movie axes.
xlabel(PlotAxes, sprintf('X (%s)', obj.LengthUnitString), ...
    'Interpreter', 'Latex')
ylabel(PlotAxes, sprintf('Y (%s)', obj.LengthUnitString), ...
    'Interpreter', 'Latex')
zlabel(PlotAxes, ...
    sprintf('%s (%s)', obj.TimeDimensionString, obj.TimeUnitString), ...
    'Interpreter', 'Latex')

% Update the flag indicating the axes have been prepped.
obj.AxesPrepped = true;


end