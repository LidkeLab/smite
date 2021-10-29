function [FigureHandle, DisplayParams] = ...
    createSummaryPlot(FigureHandle, TRArray, SMF, DisplayParams, UnitFlag)
%createSummaryPlot creates a multi-panel summary plot of the HMM analysis.
% This method creates a multi-panel plot with several sub-plots
% corresponding to a single dimer event identified by the HMM.  The
% purpose of this multi-panel plot is to contain as much useful information
% related to the HMM analysis results as possible, so that the user can
% look at the resulting figure and quickly determine if the dimer event was
% correctly identified/worth saving.
%
% INPUTS:
%   FigureHandle: Figure in which we should plot stuff.
%                 (Default = figure())
%   TRArray: A structure array of TR structures, where the constituent TR
%            structures correspond to the relevant segments of dimer
%            trajectories.  TRArray(1) will contain information about a
%            trajectory from TR1 that was found to have dimerized with
%            TRArray(2).
%   SMF: Single Molecule Fitting structure (for now, only
%        SMF.Data.PixelSize and SMF.Data.FrameRate are  used).
%        (Default = smi_core.SingleMoleculeFitting)
%   DisplayParams: A structure of display parameters for the movie.
%                  ChannelNames: Name labels for each channel.
%                                (Default = {'Channel 1'; 'Channel 2'})
%                  FileNames: List of filenames corresponding to
%                             trajectories in TRArray (e.g.,
%                             FileNames{2} is the filename corresponding to
%                             TRArray(2)). (Default = {''; ''})
%                  RegFileName: Registration filename applied to the
%                               "moving" channel of TRArray. (Default = '')
%                  MaxYDisplaySep: Maximum of y axis for any plots of
%                                  separation vs time (pixels)(Default =
%                                  max(cell2mat({TRArray.Separation}))
%                  MaxYDisplayState: Maximum of y axis for plots of the
%                                    state (i.e., states 1, 2, 3).
%                                    (Default = 1.1 * MaxStateNumber)
%                  MinXYRange: Minimum XY display range of the movie
%                              (pixels)(Default = 10)
%                  PairNumber: The pair number identifying the dimer pair
%                              being plotted. (Default = 1)
%                  StateNames: see SMA_HMM.m
%                              (Default = {'Dimer'; 'Free'})
%                  StateColormap: A colormap corresponding to the states
%                                 given in ModelSpecifier/StateNames
%                                 (NStates x 3 array)
%                                 (Default = colormap(lines(NStates)))
%   UnitFlag: Flag indicating we should use physical units. 
%             (Default = false)
%
% OUTPUTS:
%   FigureHandle: Figure in which we plotted stuff.
%   DisplayParams: Set of parameters used to produce plots, including the
%                  defaults set internally.
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Prepare some defaults if needed.
if (~exist('FigureHandle', 'var') || isempty(FigureHandle) ...
        || ~isvalid(FigureHandle))
    % Generate the figure.
    FigureHandle = figure('Units', 'inches', ...
        'Position', [0, 0, 8.5, 11], ...
        'PaperSize', [8.5, 11]);
end
if (~exist('SMF', 'var') || isempty(SMF))
    SMF = smi_core.SingleMoleculeFitting;
end
if (~exist('DisplayParams', 'var') || isempty(DisplayParams))
    DisplayParams = struct();
end
if (~exist('UnitFlag', 'var') || isempty(UnitFlag))
    UnitFlag = false;
end

% Set default parameter values where needed.
StateSequence = cell2mat({TRArray.StateSequence}.');
NStates = max([numel(unique(StateSequence(~isnan(StateSequence)))), ...
    size(TRArray(1, 1).EmissionProbabilities, 2), ...
    max(StateSequence)]);
DefaultDisplayParams.ChannelNames = {'Channel 1'; 'Channel 2'};
DefaultDisplayParams.FileNames = {''; ''};
DefaultDisplayParams.RegFileName = '';
DefaultDisplayParams.MaxYDisplaySep = ...
    max(cell2mat({TRArray.Separation}.'));
DefaultDisplayParams.MaxYDisplayState = ...
    1.1 * numel(unique(StateSequence(~isnan(StateSequence))));
DefaultDisplayParams.MinXYRange = 10;
DefaultDisplayParams.PairNumber = 1;
DefaultDisplayParams.StateNames = {'Dimer'; 'Free'};
DefaultDisplayParams.StateColormap = lines(NStates);
DisplayParams = smi_helpers.padStruct(DisplayParams, DefaultDisplayParams);

% Define some more parameters.
LengthUnitString = ...
    smi_helpers.arrayMUX({'pixel', '$\mu m$'}, UnitFlag);
TimeDimensionString = ...
    smi_helpers.arrayMUX({'Frame number', 'Time'}, UnitFlag);
TimeUnitString = ...
    smi_helpers.arrayMUX({'frame', 'second'}, UnitFlag);

% Add a title to the summary plot.
TrajectoryID1 = TRArray(1).ConnectID;
TrajectoryID2 = TRArray(2).ConnectID;
TitleString = sprintf('Dimer pair %i, %s ID: %i, %s ID: %i\n', ...
    DisplayParams.PairNumber, ...
    DisplayParams.ChannelNames{1}, TrajectoryID1, ...
    DisplayParams.ChannelNames{2}, TrajectoryID2);
if ~(isempty(DisplayParams.FileNames{1}) ...
        ||isempty(DisplayParams.FileNames{2}))
    TitleString = [TitleString, ...
        sprintf('%s file: %s\n%s file: %s\n', ...
        DisplayParams.ChannelNames{1}, DisplayParams.FileNames{1}, ...
        DisplayParams.ChannelNames{2}, DisplayParams.FileNames{2})];
end
if (((~isempty(TRArray(2).IsTransformed)&&TRArray(2).IsTransformed) ...
        || (~isempty(TRArray(1).IsTransformed)&&TRArray(1).IsTransformed)) ...
        && ~isempty(DisplayParams.RegFileName))
    TitleString = [TitleString, ...
        sprintf('registration file: %s\n', DisplayParams.RegFileName)];
end
sgtitle(FigureHandle, TitleString, 'FontSize', 8, 'Interpreter', 'none');

% Generate the state sequence + separation plot.
PlotAxes = subplot(4, 2, 1:2, 'Parent', FigureHandle);
hold(PlotAxes, 'on')
smi_stat.HMM.plotDimerPairInfo(PlotAxes, ...
    'ViterbiSequence', TRArray, SMF, DisplayParams, UnitFlag);
xlabel(PlotAxes, ...
    sprintf('%s (%s)', TimeDimensionString, TimeUnitString), ...
    'Interpreter', 'Latex')
yyaxis(PlotAxes, 'left')
ylabel(PlotAxes, sprintf('Separation (%s)', LengthUnitString), ...
    'Interpreter', 'Latex')
yyaxis(PlotAxes, 'right')
ylabel(PlotAxes, 'State', 'Interpreter', 'Latex')

% Plot the emission probabilities and the state sequence.
PlotAxes = subplot(4, 2, 3:4, 'Parent', FigureHandle);
hold(PlotAxes, 'on');
smi_stat.HMM.plotDimerPairInfo(PlotAxes, ...
    'EmissionProbability', TRArray, SMF, DisplayParams, UnitFlag);
xlabel(PlotAxes, ...
    sprintf('%s (%s)', TimeDimensionString, TimeUnitString), ...
    'Interpreter', 'Latex')
yyaxis(PlotAxes, 'left')
ylabel(PlotAxes, 'Emission probability density', ...
    'Interpreter', 'Latex')
yyaxis(PlotAxes, 'right')
ylabel(PlotAxes, 'State', 'Interpreter', 'Latex')
legend(PlotAxes, ...
    {DisplayParams.StateNames{1}, DisplayParams.StateNames{2}, ...
    'State Sequence'}, 'Location', 'best')

% Plot the entire track, color coding the dimer state.
PlotAxes = subplot(4, 2, 5, 'Parent', FigureHandle);
hold(PlotAxes, 'on');
smi_stat.HMM.plotDimerPairInfo(PlotAxes, ...
    'TrajectoryPlot3D', TRArray, SMF, DisplayParams, UnitFlag);
xlabel(PlotAxes, sprintf('X (%s)', LengthUnitString), ...
    'Interpreter', 'Latex')
ylabel(PlotAxes, sprintf('Y (%s)', LengthUnitString), ...
    'Interpreter', 'Latex')
zlabel(PlotAxes, ...
    sprintf('%s (%s)', TimeDimensionString, TimeUnitString), ...
    'Interpreter', 'Latex')
legend(PlotAxes, ...
    {DisplayParams.ChannelNames{1}, DisplayParams.ChannelNames{2}, ...
    'Dimer State'}, 'Location', 'best');

% Plot the entire track in 2D to help identify edge effects.
PlotAxes = subplot(4, 2, 6, 'Parent', FigureHandle);
hold(PlotAxes, 'on');
smi_stat.HMM.plotDimerPairInfo(PlotAxes, ...
    'TrajectoryPlot2D', TRArray, SMF, DisplayParams, UnitFlag);
xlabel(PlotAxes, sprintf('X (%s)', LengthUnitString), ...
    'Interpreter', 'Latex')
ylabel(PlotAxes, sprintf('Y (%s)', LengthUnitString), ...
    'Interpreter', 'Latex')
legend(PlotAxes, ...
    {DisplayParams.ChannelNames{1}, DisplayParams.ChannelNames{2}, ...
    'Dimer State'}, 'Location', 'best');

% Plot the x and y separations in a scatterplot.
PlotAxes = subplot(4, 2, 7, 'Parent', FigureHandle);
hold(PlotAxes, 'on');
smi_stat.HMM.plotDimerPairInfo(PlotAxes, ...
    'XYSeparationPlot', TRArray, SMF, DisplayParams, UnitFlag);
xlabel(PlotAxes, sprintf('X Separation (%s)', LengthUnitString), ...
    'Interpreter', 'Latex')
ylabel(PlotAxes, sprintf('Y Separation (%s)', LengthUnitString), ...
    'Interpreter', 'Latex')

% Plot the registration corrections over time.
if (isfield(TRArray, 'XRegCorrection') ...
        && isfield(TRArray, 'YRegCorrection'))
    PlotAxes = subplot(4, 2, 8, 'Parent', FigureHandle);
    hold(PlotAxes, 'on');
    smi_stat.HMM.plotDimerPairInfo(PlotAxes, ...
        'RegistrationPlot', TRArray, SMF, DisplayParams, UnitFlag);
    xlabel(PlotAxes, ...
        sprintf('%s (%s)', TimeDimensionString, TimeUnitString), ...
        'Interpreter', 'Latex')
    ylabel(PlotAxes, sprintf('Reg. Correction (%s)', LengthUnitString), ...
        'Interpreter', 'Latex')
    zlabel(PlotAxes, ...
        sprintf('%s (%s)', TimeDimensionString, TimeUnitString), ...
        'Interpreter', 'Latex')
    legend(PlotAxes, {'X Correction', 'YCorrection'}, 'Location', 'best');
end


end