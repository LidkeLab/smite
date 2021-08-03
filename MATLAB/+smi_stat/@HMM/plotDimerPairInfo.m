function [PlotAxes, DisplayParams] = ...
    plotDimerPairInfo(PlotAxes, TRArray, SMF, DisplayParams, PlotType)
%plotDimerPairInfo creates various plots related to dimer pair trajectories
% This method organizes several plotting codes meant to plot information
% related to a pair of dimer trajectories.
%
% NOTE: We should split this up into a few simpler methods.  Specifically,
%       each "case" of the switch/case below should be a simpler call to a
%       different method.  For now, we'll leave this as is.
% 
% INPUTS:
%   TRArray: A structure array of TR structures, where the constituent TR
%            structures correspond to the relevant segments of dimer
%            trajectories.  TRArray(1) will contain information about a 
%            trajectory from TR1 that was found to have dimerized with
%            TRArray(2).  TRArray can additionally contain a "PairIDs" 
%            cell array, which contains the trajectory ID's of trajectories
%            which were part of an oligomer in each frame (this information
%            will be plotted when PlotType = 'ViterbiSequence') as might be 
%            known for a simulated TRArray.  For example, if
%            TRArray(1).PairIDs(23) = [5; 8], this means that the
%            trajectory TRArray(1) was part of an oligomer with
%            trajectories 5 and 8 in frame 23.
%   SMF: Single Molecule Fitting structure (for now, only
%        SMF.Data.PixelSize and SMF.Data.FrameRate are  used).
%        (Default = smi_core.SingleMoleculeFitting)
%   DisplayParams: A structure of display parameters for the plots.
%                  StateColormap: A colormap corresponding to the states
%                                 given in ModelSpecifier/StateNames
%                                 (NStates x 3 array)
%                                 (Default = colormap(lines(NStates)))
%                  UnitFlag: 0 for camera units (pixels, frames) 
%                            1 for physical units (micrometers, seconds)
%                            (Default = 1)
%                  MinXYRange: Minimum XY display range of the movie
%                              (pixels)(Default = 10) 
%                  MaxYDisplaySep: Maximum of y axis for any plots of
%                                  separation vs time (pixels)(Default =
%                                  max(cell2mat({TRArray.Separation}))
%                  MaxYDisplayState: Maximum of y axis for state plots.
%                                    (Default = 1.1*NStates)
%                  PairNumber: The pair number identifying the dimer pair 
%                              being plotted. (Default = 1)
%                  StateNames: see smi_stat.HMM 
%                              (Default = {'Dimer'; 'Free'})
%                  StateToMark: Numeric value corresponding to a specific 
%                               state that should be marked when
%                               PlotType = 'TrajectoryPlot'
%                               (Default = 1)
%   PlotType: Type of plot that will be created. 
%             ViterbiSequence: Plot of the separation vs. time with state
%                              as identified by Viterbi algorithm overlain.
%             EmissionProbability: Plot of the emission probability vs.
%                                  time with the state sequence overlain.
%             TrajectoryPlot3D: 3D plot of the two trajectories over time,
%                               with the state specified by StateToIndicate
%                               being marked in some way.
%             TrajectoryPlot2D: 2D plot of the two trajectories over time,
%                               with the state specified by StateToIndicate
%                               being marked in some way.  The frame will
%                               be forced to the full frame size in this
%                               plot (the goal of this plot is to help
%                               identify edge effects, so we want to see
%                               how close to the ROI edge the trajectories
%                               were).
%             RegistrationPlot: Plot of the registration corrections made
%                               over time.
%             XYSeparationPlot: Plot the separation between the two
%                               trajectories separately in x, y as a
%                               scatterplot.
%             (Default = 'ViterbiSequence')
%   PlotAxes: Axes in which we should plot stuff. (Default = gca())
% 
% OUTPUTS:
%   PlotAxes: Axes object in which we've plotted stuff.
%   DisplayParams: The same structure as the input structure, but with the
%                  default parameters added in.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)

% Set default values if needed.
if (~exist('PlotAxes', 'var') || isempty(PlotAxes))
    PlotAxes = gca();
end
if (~exist('SMF', 'var') || isempty(SMF))
    SMF = smi_core.SingleMoleculeFitting;
end
if (~exist('PlotType', 'var') || isempty(PlotType))
    PlotType = 'ViterbiSequence';
end

% Define some parameters that we'll need for setting some defaults.
DimerCandidateBool1 = TRArray(1).DimerCandidateBool;
DimerCandidateBool2 = TRArray(2).DimerCandidateBool;
StateSequence = TRArray(1).StateSequence(DimerCandidateBool1);
EmissionProbabilities = TRArray(1).EmissionProbabilities(...
    DimerCandidateBool1, :);
NStates = max([numel(unique(StateSequence)), ...
    size(EmissionProbabilities, 2), ...
    max(StateSequence)]);

% Set default display parameters where needed.
if (~exist('DisplayParams', 'var') || isempty(DisplayParams))
    DisplayParams = struct();
end
DefaultDisplayParams.StateNames = {'Dimer'; 'Free'};
DefaultDisplayParams.StateColormap = lines(NStates);
DefaultDisplayParams.UnitFlag = 1;
DefaultDisplayParams.MinXYRange = 10;
DefaultDisplayParams.MaxYDisplaySep = ...
    max(cell2mat({TRArray.Separation}.'));
DefaultDisplayParams.MaxYDisplayState = ...
    1.1 * numel(unique(StateSequence(~isnan(StateSequence))));
DefaultDisplayParams.PairNumber = 1;
DefaultDisplayParams.StateToMark = 1;
DefaultParameterNames = fieldnames(DefaultDisplayParams);
InputParameterNames = fieldnames(DisplayParams);
for ii = 1:numel(DefaultParameterNames)
    if ~any(ismember(DefaultParameterNames{ii}, InputParameterNames))
        % The field DefaultParameterNames{ii} is not present in the
        % DisplayParams structure and so the default must be added.
        DisplayParams.(DefaultParameterNames{ii}) = ...
            DefaultDisplayParams.(DefaultParameterNames{ii});
    end
end

% Convert misc. arrays to the desired units, isolate segments of
% interest of certain arrays, and define misc. parameters.
% NOTE: Some of the default parameters are dependent on these so I do this
%       first.
FrameRate = SMF.Data.FrameRate;
PixelSize = SMF.Data.PixelSize;
MaxYDisplaySepConverted = DisplayParams.MaxYDisplaySep ...
    * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
FrameNum = TRArray(1).FrameNum(DimerCandidateBool1) ...
    * (DisplayParams.UnitFlag/FrameRate + ~DisplayParams.UnitFlag);
Separation = TRArray(1).Separation(DimerCandidateBool1) ...
    * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
if isfield(TRArray(1), 'PairIDs')
    PairIDs = TRArray(1).PairIDs(DimerCandidateBool1);
    IsDimerGroundTruth = cellfun(...
        @(X) any(TRArray(2).TrajectoryID == X), PairIDs);
else
    IsDimerGroundTruth = zeros(sum(DimerCandidateBool1), 1);
end

% Prepare the figure/axes for the plots.
axes(PlotAxes);
hold(PlotAxes, 'on');

% Proceed to plot the user selected plot.
switch PlotType
    case 'ViterbiSequence'
        % Prepare the figure/axes for the Viterbi plot.
        axis(PlotAxes, 'tight');
        yyaxis(PlotAxes, 'left');
        
        % Define misc. parameters for later use and define parameter
        % dependent defaults.
        NObservations = numel(StateSequence);
        StatePlotValues = 1:NStates;
        MinYValue = 0;
        
        % Create an RGBA color from our selected RGB color (if we just add
        % a 4-th element to the color for the alpha value, .fig file isn't
        % saving that correctly as of R2020a).  I found this conversion
        %  here: http://marcodiiga.github.io/rgba-to-rgb-conversion
        AlphaValue = 0.5;
        StateColormap = (1-AlphaValue) * repmat(PlotAxes.Color, ...
            [size(DisplayParams.StateColormap, 1), 1]) ...
            + AlphaValue*DisplayParams.StateColormap(:, 1:3);

        % Plot the state color code rectangle for the first observation.
        CurrentStateColor = StateColormap(...
            StatePlotValues == StateSequence(1), :);
        RectangleStartFrame = zeros(NObservations, 1);
        RectangleStartFrame(1) = FrameNum(1);
        if (NObservations > 1)
            RectangleXWidth = ((FrameNum(1)+FrameNum(2))/2) ...
                - RectangleStartFrame(1);
        else
            RectangleXWidth = FrameNum(1) - RectangleStartFrame(1) + 0.5;
        end
        rectangle(PlotAxes, 'Position',  ...
            [RectangleStartFrame(1), MinYValue, RectangleXWidth, ...
            MaxYDisplaySepConverted], ...
            'FaceColor', CurrentStateColor, 'EdgeColor', 'none');
        if IsDimerGroundTruth(1)
            rectangle(PlotAxes, 'Position',  ...
                [RectangleStartFrame(1), ...
                MinYValue+MaxYDisplaySepConverted/2, ...
                RectangleXWidth, MaxYDisplaySepConverted/2], ...
                'FaceColor', CurrentStateColor, 'EdgeColor', 'none');
        end
        
        % Plot the state color code rectangles for the remaining
        % observations.
        for ii = 2:(NObservations-1)
            % Plot a rectangle as a background color behind each state.
            CurrentStateColor = StateColormap(...
                StatePlotValues == StateSequence(ii), :);
            RectangleStartFrame(ii) = (FrameNum(ii-1) + FrameNum(ii)) / 2;
            RectangleXWidth = (FrameNum(ii+1)+FrameNum(ii))/2 ...
                - RectangleStartFrame(ii);
            rectangle(PlotAxes, 'Position',  ...
                [RectangleStartFrame(ii), MinYValue, ...
                RectangleXWidth, MaxYDisplaySepConverted], ...
                'FaceColor', CurrentStateColor, 'EdgeColor', 'none');
            if IsDimerGroundTruth(ii)
                rectangle(PlotAxes, 'Position',  ...
                    [RectangleStartFrame(ii), ...
                    MinYValue+MaxYDisplaySepConverted/2, ...
                    RectangleXWidth, MaxYDisplaySepConverted/2], ...
                    'FaceColor', CurrentStateColor, 'EdgeColor', 'none');
            end
        end
        CurrentStateColor = StateColormap(...
            StatePlotValues == StateSequence(NObservations), :);
        if (NObservations > 1)
            RectangleStartFrame(NObservations) = ...
                (FrameNum(NObservations-1) + FrameNum(NObservations)) / 2;
        else
            RectangleStartFrame(1) = FrameNum(1) - (1/2);
        end
        RectangleXWidth = FrameNum(NObservations) ...
            - RectangleStartFrame(NObservations);
        rectangle(PlotAxes, 'Position',  ...
            [RectangleStartFrame(NObservations), MinYValue, ...
            RectangleXWidth, MaxYDisplaySepConverted], ...
            'FaceColor', CurrentStateColor, 'EdgeColor', 'none');
        if IsDimerGroundTruth(NObservations)
            rectangle(PlotAxes, 'Position',  ...
                [RectangleStartFrame(NObservations), ...
                MinYValue + MaxYDisplaySepConverted/2, ...
                RectangleXWidth, MaxYDisplaySepConverted/2], ...
                'FaceColor', CurrentStateColor, 'EdgeColor', 'none');
        end
        line(PlotAxes, FrameNum, Separation, 'Color', 'k', ...
            'LineWidth', 1, 'LineStyle', '-')
        PlotAxes.YColor = [0, 0, 0];
        PlotAxes.YLim = [MinYValue, MaxYDisplaySepConverted];
        yyaxis(PlotAxes, 'right');
        stairs(PlotAxes, ...
            [RectangleStartFrame; ...
            RectangleStartFrame(end) + RectangleXWidth], ...
            [StateSequence; StateSequence(end)], ...
            'Color', 'k', 'LineWidth', 2, 'LineStyle', ':')
        PlotAxes.YColor = [0, 0, 0];
        PlotAxes.YLim = [0, DisplayParams.MaxYDisplayState];
        PlotAxes.YTick = min(StateSequence):max(StateSequence);
        PlotAxes.YTickLabel = DisplayParams.StateNames;
        PlotAxes.YColor = [0, 0, 0];
    case 'EmissionProbability'
        % Generate the emission probability vs. time plot with the Viterbi
        % state sequence shown on top.
        yyaxis(PlotAxes, 'left')
        for jj = 1:NStates
            plot(PlotAxes, FrameNum, ...
                EmissionProbabilities(:, jj), ...
                'Color', DisplayParams.StateColormap(jj, :), ...
                'LineWidth', 1.5, 'LineStyle', '-');
        end
        yyaxis(PlotAxes, 'right')
        stairs(PlotAxes, [FrameNum(1); ...
            (FrameNum(1:(end-1))+FrameNum(2:end))/2; ...
            FrameNum(end)], ...
            [StateSequence; StateSequence(end)], ...
            'k:', 'LineWidth', 2)
        PlotAxes.XLim = [min(FrameNum), max(FrameNum)];
        PlotAxes.YLim = [0, DisplayParams.MaxYDisplayState];
        PlotAxes.YTick = min(StateSequence):max(StateSequence);
        PlotAxes.YTickLabel = DisplayParams.StateNames;
        PlotAxes.YColor = [0, 0, 0];
    case 'TrajectoryPlot3D'
        % Plot the two trajectories over time in 3D.
        view(PlotAxes, -45, 5) % force corner view
        axis(PlotAxes, 'tight'); % will be overridden in x, y below
        XCurrent1 = TRArray(1).X(DimerCandidateBool1) ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        YCurrent1 = TRArray(1).Y(DimerCandidateBool1) ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        XCurrent2 = TRArray(2).X(DimerCandidateBool2) ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        YCurrent2 = TRArray(2).Y(DimerCandidateBool2) ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        line(PlotAxes, XCurrent1, YCurrent1, FrameNum, ...
            'Color', [0, 0, 0], 'LineWidth', 1.5, 'LineStyle', '--');
        line(PlotAxes, XCurrent2, YCurrent2, FrameNum, ...
            'Color', [0, 0, 0], 'LineWidth', 1.5, 'LineStyle', ':');
        
        % Mark the appropriate state (if requested) in the plot.
        MarkedStateBool = (StateSequence == DisplayParams.StateToMark);
        EventChanges = [MarkedStateBool(1); diff(MarkedStateBool)];
        StartIndices = find(EventChanges == 1);
        EndIndices = find(EventChanges == -1) - 1;
        if (numel(EndIndices) ~= numel(StartIndices))
            EndIndices = [EndIndices; numel(MarkedStateBool)];
        end
        for kk = 1:numel(StartIndices)
            % Loop through all instances of the jj-th state event and
            % plot them (this is done to make the lines continuous, not
            % jumping all over when switching to another state).
            EventIndices = StartIndices(kk):EndIndices(kk);
            line(PlotAxes, ...
                XCurrent1(EventIndices), ...
                YCurrent1(EventIndices), ...
                FrameNum(EventIndices), ...
                'Color', [DisplayParams.StateColormap(1, :), 0.5], ...
                'LineWidth', 2, 'LineStyle', '-');
            line(PlotAxes, ...
                XCurrent2(EventIndices), ...
                YCurrent2(EventIndices), ...
                FrameNum(EventIndices), ...
                'Color', [DisplayParams.StateColormap(1, :), 0.5], ...
                'LineWidth', 2, 'LineStyle', '-');
        end
        
        % Modify the x, y plot limits to improve appearance.
        XRange = [min(min(XCurrent1), min(XCurrent2)), ...
            max(max(XCurrent1), max(XCurrent2))]; % x range of data
        XWidth = max(XRange(2) - XRange(1), DisplayParams.MinXYRange ...
            * (DisplayParams.UnitFlag*PixelSize ...
            + ~DisplayParams.UnitFlag));
        YRange = [min(min(YCurrent1), min(YCurrent2)), ...
            max(max(YCurrent1), max(YCurrent2))];
        YWidth = max(YRange(2)-YRange(1), DisplayParams.MinXYRange ...
            * (DisplayParams.UnitFlag*PixelSize ...
            + ~DisplayParams.UnitFlag));
        PlotAxes.XLim = mean(XRange) + [-XWidth, XWidth]/2;
        PlotAxes.YLim = mean(YRange) + [-YWidth, YWidth]/2;
    case 'TrajectoryPlot2D'
        % Plot the two trajectories in 2D
        view(PlotAxes, 2)
        axis(PlotAxes, 'equal')
        XCurrent1 = TRArray(1).X(DimerCandidateBool1) ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        YCurrent1 = TRArray(1).Y(DimerCandidateBool1) ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        XCurrent2 = TRArray(2).X(DimerCandidateBool2) ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        YCurrent2 = TRArray(2).Y(DimerCandidateBool2) ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        line(PlotAxes, XCurrent1, YCurrent1, ...
            'Color', [0, 0, 0], 'LineWidth', 1.5, 'LineStyle', '--');
        line(PlotAxes, XCurrent2, YCurrent2, ...
            'Color', [0, 0, 0], 'LineWidth', 1.5, 'LineStyle', ':');
        
        % Mark the appropriate state (if requested) in the plot.
        MarkedStateBool = (StateSequence == DisplayParams.StateToMark);
        EventChanges = [MarkedStateBool(1); diff(MarkedStateBool)];
        StartIndices = find(EventChanges == 1);
        EndIndices = find(EventChanges == -1) - 1;
        if (numel(EndIndices) ~= numel(StartIndices))
            EndIndices = [EndIndices; numel(MarkedStateBool)];
        end
        for kk = 1:numel(StartIndices)
            % Loop through all instances of the jj-th state event and
            % plot them (this is done to make the lines continuous, not
            % jumping all over when switching to another state).
            EventIndices = StartIndices(kk):EndIndices(kk);
            line(PlotAxes, ...
                XCurrent1(EventIndices), ...
                YCurrent1(EventIndices), ...
                'Color', [DisplayParams.StateColormap(1, :), 0.5], ...
                'LineWidth', 2, 'LineStyle', '-');
            line(PlotAxes, ...
                XCurrent2(EventIndices), ...
                YCurrent2(EventIndices), ...
                'Color', [DisplayParams.StateColormap(1, :), 0.5], ...
                'LineWidth', 2, 'LineStyle', '-');
        end
        
        % Modify the x, y plot limits to improve appearance.
        PlotAxes.XLim = [0, TRArray(1).XSize] ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        PlotAxes.YLim = [0, TRArray(1).YSize] ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
    case 'RegistrationPlot'
        XRegCorrection = sum(...
            [cell2mat({TRArray(1).XRegCorrection(DimerCandidateBool1)}),...
            cell2mat({TRArray(2).XRegCorrection(DimerCandidateBool2)})], 2);
        XRegCorrection = XRegCorrection ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        YRegCorrection = sum(...
            [cell2mat({TRArray(1).YRegCorrection(DimerCandidateBool1)}),...
            cell2mat({TRArray(2).YRegCorrection(DimerCandidateBool2)})], 2);
        YRegCorrection = YRegCorrection ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        plot(PlotAxes, FrameNum, XRegCorrection, 'kx')
        plot(PlotAxes, FrameNum, YRegCorrection, 'ko')
    case 'XYSeparationPlot'
        % Plot the separations in x and y over time.
        XSeparation = cell2mat({TRArray(2).X(DimerCandidateBool2)}.') ...
            - cell2mat({TRArray(1).X(DimerCandidateBool1)}.');
        XSeparation = XSeparation ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        YSeparation = cell2mat({TRArray(2).Y(DimerCandidateBool2)}.') ...
            - cell2mat({TRArray(1).Y(DimerCandidateBool1)}.');
        YSeparation = YSeparation ...
            * (DisplayParams.UnitFlag*PixelSize + ~DisplayParams.UnitFlag);
        plot(PlotAxes, XSeparation, YSeparation, 'k-')
        
        % Mark the appropriate state (if requested) in the plot.
        MarkedStateBool = (StateSequence == DisplayParams.StateToMark);
        EventChanges = [MarkedStateBool(1); diff(MarkedStateBool)];
        StartIndices = find(EventChanges == 1);
        EndIndices = find(EventChanges == -1) - 1;
        if (numel(EndIndices) ~= numel(StartIndices))
            EndIndices = [EndIndices; numel(MarkedStateBool)];
        end
        for kk = 1:numel(StartIndices)
            % Loop through all instances of the jj-th state event and
            % plot them (this is done to make the lines continuous, not
            % jumping all over when switching to another state).
            EventIndices = StartIndices(kk):EndIndices(kk);
            line(PlotAxes, ...
                XSeparation(EventIndices), ...
                YSeparation(EventIndices), ...
                FrameNum(EventIndices), ...
                'Color', [DisplayParams.StateColormap(1, :), 0.5], ...
                'LineWidth', 2, 'LineStyle', '-');
        end
        
        % Modify the x and y limits to improve plot appearance.
        XRange = [min(XSeparation), max(XSeparation)];
        XWidth = max(2*max(abs(XRange)), DisplayParams.MinXYRange ...
            * (DisplayParams.UnitFlag*PixelSize ...
            + ~DisplayParams.UnitFlag));
        YRange = [min(YSeparation), max(YSeparation)];
        YWidth = max(2*max(abs(YRange)), DisplayParams.MinXYRange ...
            * (DisplayParams.UnitFlag*PixelSize ...
            + ~DisplayParams.UnitFlag));
        PlotAxes.XLim = [-XWidth, XWidth] / 2;
        PlotAxes.YLim = [-YWidth, YWidth] / 2;
        plot(PlotAxes, [0, 0], PlotAxes.YLim, 'm:')
        plot(PlotAxes, PlotAxes.XLim, [0, 0], 'm:')
end


end