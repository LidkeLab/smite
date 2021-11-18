function [varargout] = plotDriftAndReg(PlotAxes, SMD, SMF)
%plotDriftAndReg plots drift correction and brightfield corrections.
% This method plots drift correction results (both inter- and intra-drift
% correction) as well as brightfield registration results in the same plot.
% For now, this method only works for 2D results.
%
% INPUTS:
%   PlotAxes: Axes in which the plot will be made. (Default = gca())
%   SMD: SingleMoleculeData structure containing the fields DriftX and
%        DriftY.
%   SMF: SingleMoleculeFitting structure defining the raw data file
%        containing brightfield registration info., with SMF.Data.FileDir,
%        SMF.Data.FileName, and SMF.Data.PixelSize all being populated.
%
% OUTPUTS:
%   varargout: varargout contains PlotAxes if an output is requested.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('PlotAxes', 'var') || isempty(PlotAxes) || ~isvalid(PlotAxes))
    PlotAxes = gca();
end

% Extract the brightfield registration corrections made during an
% experiment.
PxToNm = 1e3 * SMF.Data.PixelSize;
RegCorrection = -PxToNm * smi.Publish.computeRegCorrection(SMF);

% Separate intra- and inter-drift correction (even if not necessary, this
% might help clarify the code below).
InterX = PxToNm * SMD.DriftX(1, :);
InterY = PxToNm * SMD.DriftY(1, :);
IntraX = PxToNm*SMD.DriftX - InterX;
IntraY = PxToNm*SMD.DriftY - InterY;

% Define the size of points in our drift plot (lower frames -> larger
% points, higher frames -> smaller points, to help indicate time ordering)
% as well as their color (blue -> red corresponds to increasing time).
NFrames = SMD.NFrames;
ColorMap = parula(NFrames * SMD.NDatasets);
Mag = 20;
MarkerSize = (NFrames:-1:1) / Mag;

% Prepare the axes.
hold(PlotAxes, 'on');
axis(PlotAxes, 'equal')
xlabel(PlotAxes, 'X (nm)')
ylabel(PlotAxes, 'Y (nm)')
title(PlotAxes, 'Corrected drift and brightfield registration')

% For each dataset, plot the drift correction results and the brightfield
% registration results (ignoring the first dataset brightfield registration
% results, which are just initial alignments and might not be very useful).
CumulativeInter = cumsum([InterX.', InterY.'], 1);
RegCorrection(1, :) = [0, 0, 0];
CumulativeReg = [0, 0; cumsum(RegCorrection(2:SMD.NDatasets, 1:2), 1)];
BaseCoords = CumulativeReg + CumulativeInter;
for ii = 1:SMD.NDatasets
    % Plot points for each intra-drift correction, starting from 'BaseCoords'.
    scatter(PlotAxes, BaseCoords(ii, 1)+IntraX(:, ii), ...
        BaseCoords(ii, 2)+IntraY(:, ii), ...
        'Marker', '.', 'SizeData', MarkerSize, ...
        'CData', ColorMap((ii-1)*NFrames + (1:NFrames), :))
    
    % Add an arrow to indicate the contribution from inter-dataset DC.
    quiver(PlotAxes, ...
        BaseCoords(ii, 1)-InterX(ii), BaseCoords(ii, 2)-InterY(ii), ...
        InterX(ii), InterY(ii), 'Color', [1, 0, 1], 'LineWidth', 2, ...
        'AutoScale', 'off')
    
    % Add an arrow indicating the brightfield registration correction.
    if (ii > 1)
        quiver(PlotAxes, ...
            BaseCoords(ii, 1)-InterX(ii)-RegCorrection(ii, 1), ...
            BaseCoords(ii, 2)-InterY(ii)-RegCorrection(ii, 2), ...
            RegCorrection(ii, 1), RegCorrection(ii, 2), ...
            'Color', [0, 0, 0], 'LineWidth', 2, ...
            'AutoScale', 'off', 'MaxHeadSize', 0.2)
    end
end

% Add some extra decorations (legend, colorbar, ...).
DummyPoints(1) = line(PlotAxes, NaN, NaN, 'Color', [0, 0, 0]);
DummyPoints(2) = line(PlotAxes, NaN, NaN, 'Color', [1, 0, 1]);
legend(PlotAxes, DummyPoints, {'Brightfield', 'Inter-DS DC'}, ...
    'Location', 'best')
ColorBar = colorbar(PlotAxes);
ColorBar.Label.String = sprintf('fraction of total frames (= %d x %d)', ...
    SMD.NDatasets, SMD.NFrames);

% Return the plot axes if needed.
if (nargout > 0)
    varargout{1} = PlotAxes;
end


end