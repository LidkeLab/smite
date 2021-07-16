function [] = prepAxes(obj)
%prepAxes prepares some features of the movie axes.
% This method does a few things to the axes obj.MovieAxes based on fields 
% in obj.Params (e.g., setting the axis limits).

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Revise several properties of the axes.
hold(obj.MovieAxes, 'on')
obj.MovieAxes.ActivePositionProperty = 'position';
obj.MovieAxes.DataAspectRatio = [1, 1, ...
    max(obj.Params.ZFrames)/max([obj.Params.XPixels, obj.Params.YPixels])];
obj.MovieAxes.YDir = 'reverse';
obj.MovieAxes.XLimMode = 'manual';
obj.MovieAxes.YLimMode = 'manual';
obj.MovieAxes.ZLimMode = 'manual';
obj.MovieAxes.XLim = obj.Params.XPixels + [0, 1];
obj.MovieAxes.YLim = obj.Params.YPixels + [0, 1];
obj.MovieAxes.ZLim = obj.Params.ZFrames;
view(obj.MovieAxes, obj.Params.LineOfSite(1, :));
colormap(obj.MovieAxes, 'gray')

% Revise axes ticks based on the unit flag.
obj.MovieAxes.XTickMode = 'manual';
obj.MovieAxes.YTickMode = 'manual';
obj.MovieAxes.ZTickMode = 'manual';
obj.MovieAxes.XTick = linspace(obj.MovieAxes.XLim(1), obj.MovieAxes.XLim(2), 5);
obj.MovieAxes.YTick = linspace(obj.MovieAxes.YLim(1), obj.MovieAxes.YLim(2), 5);
obj.MovieAxes.ZTick = linspace(obj.MovieAxes.ZLim(1), obj.MovieAxes.ZLim(2), 5);
if obj.Params.UnitFlag
    obj.MovieAxes.XTickLabelMode = 'manual';
    obj.MovieAxes.YTickLabelMode = 'manual';
    obj.MovieAxes.ZTickLabelMode = 'manual';
    xtickformat(obj.MovieAxes, '%.1f')
    ytickformat(obj.MovieAxes, '%.1f')
    ztickformat(obj.MovieAxes, '%.1f')
    XTicks = (obj.MovieAxes.XTick-1) * obj.SMF.Data.PixelSize;
    YTicks = (obj.MovieAxes.YTick-1) * obj.SMF.Data.PixelSize;
    ZTicks = (obj.MovieAxes.ZTick-1) / obj.SMF.Data.FrameRate;
    obj.MovieAxes.XTickLabel = num2str(XTicks.', '%.1f');
    obj.MovieAxes.YTickLabel = num2str(YTicks.', '%.1f');
    obj.MovieAxes.ZTickLabel = num2str(ZTicks.', '%.1f');
else
    xtickformat(obj.MovieAxes, '%i')
    ytickformat(obj.MovieAxes, '%i')
    ztickformat(obj.MovieAxes, '%i')
end

% Add some decorations to the movie axes.
xlabel(obj.MovieAxes, sprintf('X (%s)', obj.LengthUnitString), ...
    'Interpreter', 'Latex')
ylabel(obj.MovieAxes, sprintf('Y (%s)', obj.LengthUnitString), ...
    'Interpreter', 'Latex')
zlabel(obj.MovieAxes, ...
    sprintf('%s (%s)', obj.TimeDimensionString, obj.TimeUnitString), ...
    'Interpreter', 'Latex')


end