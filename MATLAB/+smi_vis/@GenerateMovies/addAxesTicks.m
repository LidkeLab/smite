function [] = addAxesTicks(obj, PlotAxes, NTicks, FormatSpec)
%addAxesTicks prepares/displays tick labels on the axes.
% This method attempts to make nice choices for the ticks/tick labels on
% the provided axes.
%
% INPUTS:
%   PlotAxes: Axes that will be modified. (Default = obj.MovieAxes)
%   NTicks: Number of ticks to be displayed in the visible axes. 
%           (Default = 5)
%   FormatSpec: Tick number formatting specifiers (see usage of sprintf()
%               below when obj.UnitFlag is true)(Default = '%.1f')

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Set defaults if needed.
if (~exist('PlotAxes', 'var') || isempty(PlotAxes) || ~isvalid(PlotAxes))
    PlotAxes = obj.MovieAxes;
end
if (~exist('NTicks', 'var') || isempty(NTicks))
    NTicks = 5;
end
if (~exist('FormatSpec', 'var') || isempty(FormatSpec))
    FormatSpec = '%.1f';
end

% Revise axes ticks based on the unit flag.
PlotAxes.XTickMode = 'manual';
PlotAxes.YTickMode = 'manual';
PlotAxes.ZTickMode = 'manual';
PlotAxes.XTick = linspace(PlotAxes.XLim(1), PlotAxes.XLim(2), NTicks);
PlotAxes.YTick = linspace(PlotAxes.YLim(1), PlotAxes.YLim(2), NTicks);
PlotAxes.ZTick = linspace(PlotAxes.ZLim(1), PlotAxes.ZLim(2), NTicks);
if obj.Params.UnitFlag
    PlotAxes.XTickLabelMode = 'manual';
    PlotAxes.YTickLabelMode = 'manual';
    PlotAxes.ZTickLabelMode = 'manual';
    xtickformat(PlotAxes, FormatSpec)
    ytickformat(PlotAxes, FormatSpec)
    ztickformat(PlotAxes, FormatSpec)
    XTicks = (PlotAxes.XTick-1) * obj.SMF.Data.PixelSize;
    YTicks = (PlotAxes.YTick-1) * obj.SMF.Data.PixelSize;
    ZTicks = (PlotAxes.ZTick-1) / obj.SMF.Data.FrameRate;
    PlotAxes.XTickLabel = num2str(XTicks.', FormatSpec);
    PlotAxes.YTickLabel = num2str(YTicks.', FormatSpec);
    PlotAxes.ZTickLabel = num2str(ZTicks.', FormatSpec);
else
    xtickformat(PlotAxes, '%i')
    ytickformat(PlotAxes, '%i')
    ztickformat(PlotAxes, '%i')
end


end