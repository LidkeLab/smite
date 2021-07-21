function [] = prepAxes(obj, PlotAxes)
%prepAxes prepares some features of the movie axes.
% This method does a few things to the axes 'PlotAxes' based on fields 
% in obj.Params (e.g., setting the axis limits).
%
% INPUTS:
%   PlotAxes: Axes that will be modified. (Default = obj.MovieAxes)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


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
PlotAxes.XTickMode = 'manual';
PlotAxes.YTickMode = 'manual';
PlotAxes.ZTickMode = 'manual';
PlotAxes.XTick = linspace(PlotAxes.XLim(1), PlotAxes.XLim(2), 5);
PlotAxes.YTick = linspace(PlotAxes.YLim(1), PlotAxes.YLim(2), 5);
PlotAxes.ZTick = linspace(PlotAxes.ZLim(1), PlotAxes.ZLim(2), 5);
if obj.Params.UnitFlag
    PlotAxes.XTickLabelMode = 'manual';
    PlotAxes.YTickLabelMode = 'manual';
    PlotAxes.ZTickLabelMode = 'manual';
    xtickformat(PlotAxes, '%.1f')
    ytickformat(PlotAxes, '%.1f')
    ztickformat(PlotAxes, '%.1f')
    XTicks = (PlotAxes.XTick-1) * obj.SMF.Data.PixelSize;
    YTicks = (PlotAxes.YTick-1) * obj.SMF.Data.PixelSize;
    ZTicks = (PlotAxes.ZTick-1) / obj.SMF.Data.FrameRate;
    PlotAxes.XTickLabel = num2str(XTicks.', '%.1f');
    PlotAxes.YTickLabel = num2str(YTicks.', '%.1f');
    PlotAxes.ZTickLabel = num2str(ZTicks.', '%.1f');
else
    xtickformat(PlotAxes, '%i')
    ytickformat(PlotAxes, '%i')
    ztickformat(PlotAxes, '%i')
end

% Add some decorations to the movie axes.
xlabel(PlotAxes, sprintf('X (%s)', obj.LengthUnitString), ...
    'Interpreter', 'Latex')
ylabel(PlotAxes, sprintf('Y (%s)', obj.LengthUnitString), ...
    'Interpreter', 'Latex')
zlabel(PlotAxes, ...
    sprintf('%s (%s)', obj.TimeDimensionString, obj.TimeUnitString), ...
    'Interpreter', 'Latex')


end