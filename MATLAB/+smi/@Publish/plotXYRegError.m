function [PlotAxes, RegError] = plotXYRegError(PlotAxes, SMD)
%plotXYRegError plots the x,y registration error in a scatter plot
% This method will produce a parametric plot for the intra-dataset
% registration errors in x and y.
%
% INPUT:
%    SMD: Single Molecule Data structure containing the results from
%         a previous single molecule analysis.
%
% OUTPUT:
%    PlotAxes: Figure handle for the figure containing the parametric
%               plots.
%    RegError: A two column array with the first column being the
%              registration error in x, second column in y
%
% REQUIRES: 
%   Statistics and Machine Learning Toolbox (to use pdist() and 
%       squareform()). 
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Check inputs/set defaults.
if ~exist('SMD', 'var')
    error(['You must enter an SMR structure', ...
        'to plot the registration errors.'])
end
if (~exist('PlotAxes', 'var') || isempty(PlotAxes))
    PlotFigure = figure();
    PlotAxes = axes(PlotFigure);
end

% Extract relevant parameters from the SMD structure.
DriftX = SMD.DriftX;
DriftY = SMD.DriftY;
NDatasets = SMD.NDatasets;

% Compute the intra-dataset registration errors (the displacement of the
% starting point of the drift model from [0; 0]), defining the initial
% registration error to be [0; 0].
RegError = zeros(NDatasets, 2);
for ii = 2:NDatasets
    RegError(ii, :) = [DriftX(1, ii), DriftY(1, ii)];
end

% Create the scatter plot.
plot(PlotAxes, RegError(:, 1), RegError(:, 2), '.')
hold(PlotAxes, 'on')

% Plot a circle with diameter equal to the maximum pairwise distance
% between any two points, centered on the line connecting those two points.
PairwiseDistance = triu(squareform(pdist(RegError)));
[MaxDist, Index] = max(PairwiseDistance(:));
[Point1Row, Point2Row] = ind2sub(size(PairwiseDistance), Index);
MaxCircleCenter = (RegError(Point1Row, :)+RegError(Point2Row, :)) / 2;
Theta = linspace(0, 2*pi, 1000);
PwMaxCircle = (MaxDist/2)*[cos(Theta); sin(Theta)] + MaxCircleCenter.';
plot(PlotAxes, PwMaxCircle(1, :), PwMaxCircle(2, :), 'r:')

% Plot a circle centered at the center of mass of the errors with radius
% equal to the standard deviation of the distance of the errors from [0; 0]
RegErrorDistance = sqrt(sum(RegError.^2, 2));
RegErrorStd = std(RegErrorDistance);
CenterOfMass = (1/NDatasets) * sum(RegError);
StdCircle = RegErrorStd*[cos(Theta); sin(Theta)] + CenterOfMass.';
plot(PlotAxes, StdCircle(1, :), StdCircle(2, :), 'k-.')

% Modify the appearance of the plot to improve readability.
PlotAxes.XLim = [-1, 1] * sqrt(2) * MaxDist;
PlotAxes.YLim = [-1, 1] * sqrt(2) * MaxDist;
PlotAxes.PlotBoxAspectRatioMode = 'manual';
title(PlotAxes, 'X and Y Registration Errors')
xlabel(PlotAxes, 'X Registration Error (pixels)')
ylabel(PlotAxes, 'Y Registration Error (pixels)')
legend(PlotAxes, {'Registration Error', ...
    sprintf('max. pairwise dist. = %.3f pixels', MaxDist), ...
    sprintf('st. dev. = %.3f pixels', RegErrorStd)}, 'Location', 'best');


end