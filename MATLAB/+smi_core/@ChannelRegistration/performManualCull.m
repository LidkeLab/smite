function [CulledCoordinates] = performManualCull(RawData, Coordinates)
%performManualCull performs the interactive (graphical) culling process
% This method will plot RawData and points in Coordinates on the same
% axes.  The user can click points of interest to remove them (the 
% intention is that this method can be used to remove "bad" points so they
% don't contribute the the final transform computed from a set of SMDs).
%
% INPUTS:
%   RawData: Stack of images that will be plotted and have coordinates 
%            overlain. (numeric array, MxNx2)
%   Coordinates: Stack of paired coordinates. (numeric array, PxQx2)
%
% OUTPUTS:
%   CulledCoordinates: Subset of the input Coordinates containing only
%                      those points which weren't culled.
%                      (numeric array, RxSx2, R,S<=P,Q)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Prepare a figure and axes for the plotting.
PlotFigure = figure();
PlotAxes = axes(PlotFigure);
axis(PlotAxes, 'tight')
PlotAxes.YDir = 'reverse';
title(PlotAxes, {'Click localizations to cull points', ...
    ['Right click in figure window (outside of axes) to undo ', ...
    'most recent cull'], 'Close the figure when finished'});
hold(PlotAxes, 'on');

% Plot the RawData and the localizations in Coordinates.  If there are more
% than 2 fiducials included, we'll do the culling process in a loop,
% treating the first image in RawData (and first element of Coordinates) as 
% the reference.
NPoints = size(Coordinates, 1);
[~, LineHandles] = smi_core.ChannelRegistration.plotCoordsOnData(...
    PlotAxes, RawData, {Coordinates(:, :, 1), Coordinates(:, :, 2)});

% Add the 'ButtonDownFcn' to each of LineHandles and the PlotFigure.
PlotFigure.ButtonDownFcn = @rightClickUndoCull;
MostRecentCull = [];
KeepPointsBool = ones(NPoints, 1, 'logical');
for ii = 1:NPoints
    LineHandles{1}(ii).ButtonDownFcn = {@localizationClicked, ii};
    LineHandles{2}(ii).ButtonDownFcn = {@localizationClicked, ii};
end

% Wait for the figure to be closed and then remove the points
% marked for removal.
waitfor(PlotFigure);
CulledCoordinates = Coordinates(KeepPointsBool, :, :);
clear('LineHandles')

    function localizationClicked(~, ~, ii)
        % This is the callback function called when a user clicks a
        % localization which they wish to remove from the FixedPoints and
        % MovingPoints arrays.        
        
        % Hide the culled points in the plot.
        LineHandles{1}(ii).Visible = 'off';
        LineHandles{2}(ii).Visible = 'off';
        
        % Mark the culled points in the point arrays that are used to
        % find the transform.
        KeepPointsBool(ii) = 0;
        MostRecentCull = ii;
    end

    function rightClickUndoCull(Source, ~)
        % This is the callback function called when a user right clicks the
        % figure, which will bring back the most recently culled pair of
        % localizations.
        
        % Check if this was a right click event, exiting this function if
        % it wasn't.
        if (~strcmpi(Source.SelectionType, 'alt') ...
                || isempty(MostRecentCull))
            return
        end
        
        % Undo the most recent cull.
        LineHandles{1}(MostRecentCull).Visible = 'on';
        LineHandles{2}(MostRecentCull).Visible = 'on';
        KeepPointsBool(MostRecentCull) = 1;
    end


end