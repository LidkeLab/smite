function [PlotAxes, LineHandles] = plotCoordsOnData(PlotAxes, ...
    RawData, Coordinates)
%plotCoordsOnData plots coordinates in Coordinates on top of ScaledData.
% This method will take coordinates in Coordinates and plot them on top of
% ScaledData.  The intention is that RawData will contain a stack of
% fiducial images, which will be displayed in some projection so they are
% both visible at once.  Coordinates will be a cell array of numeric arrays
% corresponding to localizations in each image of ScaledData.
%
% NOTE: Input PlotAxes can be set to an empty array to use the default.
% NOTE: This method is best used to compare only two images/localization
%       sets.  If numel(Coordinates) > 2, the results probably won't look
%       very nice.
%
% INPUTS:
%   PlotAxes: Axes in which we will make the plots. (Default = gca())
%   RawData: Stack of images that will be plotted and have coordinates 
%            overlain. (numeric array, MxNxP with P>1)
%   Coordinates: Cell array of coordinate arrays, where numel(Coordinates)
%                must match the third dimension of RawData. 
%                (Px1 cell array of numeric arrays)
%
% OUTPUTS:
%   PlotAxes: Axes containing the plots.
%   LineHandles: Handles to the points in Coordinates added to the plot.
%                This is a cell array such that LineHandles{2} contains an
%                array of line handles for the points in Coordinates{2}.
%                (This is useful if you intend to modify the plot outside
%                of this method.) (cell array of line handles)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set defaults/modify inputs if needed.
if (~exist('PlotAxes', 'var') || isempty(PlotAxes) || ~isvalid(PlotAxes))
    PlotAxes = gca();
end
if ~iscell(Coordinates)
    Coordinates = {Coordinates};
end

% Plot the images in ScaledData.
SumImage = sum(RawData, 3);
SumImage = (SumImage-min(SumImage(:))) ...
    ./ max(max(SumImage-min(SumImage(:))));
surface(PlotAxes, ...
    [1, size(SumImage, 2)], [1, size(SumImage, 1)], ...
    [0, 0; 0, 0], repmat(SumImage, [1, 1, 3]), ...
    'facecolor', 'texturemap');

% Plot the localizations in FiducialSMD.
NPoints = cellfun(@(X) size(X, 1), Coordinates);
NDatasets = numel(Coordinates);
LineHandles = cell(NDatasets, 1);
MarkerOptions = repelem('h', NDatasets);
MarkerOptions(1:14) = ['o', '+', '*', '.', 'x', '_', '|', 's', 'd', ...
    '^', 'v', '>', '<', 'p'];
LineColors = colormap(PlotAxes, lines(max(NPoints)));
for ii = 1:NDatasets
    % Plot each of the localizations for each datset in a different color
    % and with different markers if possible. 
    LineHandles{ii} = gobjects(NPoints(ii), 1);
    for jj = 1:NPoints(ii)
        LineHandles{ii}(jj) = line(PlotAxes, ...
            Coordinates{ii}(jj, 1)-0.5, Coordinates{ii}(jj, 2)-0.5, ...
            'Marker', MarkerOptions(ii), 'Color', LineColors(jj, :));
    end
end


end