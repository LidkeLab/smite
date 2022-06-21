function [LineHandles] = plotTrajectories(PlotAxes, ...
    Params, TR, FrameRange, Color, varargin)
%plotTrajectories plots trajectories in the specified axes.
% This method plots trajectories in 'TR'.  This method is intended to
% remain "lightweight" for speed purposes, so minimal input checks/default
% settings should be implemented here.
%
% INPUTS:
%   PlotAxes: Axes in which the trajectories will be plotted.
%   Params: Structure of parameters (see usage below and settings in
%           obj.prepOligoDefaults())
%   TR: Tracking Results structure (see smi_core.TrackingResults).
%   FrameRange: Range of frames over which trajectories should be plotted.
%   Color: Color of the trajectory lines. (NTrajectoriesx3(4) float array)
%   varargin: Additional inputs that can be provided as keyword arguments
%             for the MATLAB line() method (see Line Properties from line()
%             documentation for options).  For example, you can change
%             the Marker property to 'x' as plotTrajectories(PlotAxes, TR,
%             FrameRange, MaxTrajLength, Color, 'Marker', 'x')
%
% OUTPUTS:
%   LineHandles: Array of line handles corresponding to the trajectories
%                in 'TR'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Loop through trajectories present in 'TR' and plot them, enforcing the
% requested 'FrameRange' and 'MaxTrajLength'.
NTraj = numel(TR);
LineHandles = gobjects(NTraj, 1);
IndicateDimer = Params.IndicateDimer;
IndicateDimerCandidate = Params.IndicateDimerCandidate;
AddDimerDecorations = ...
    ((IndicateDimerCandidate&&isfield(TR, 'DimerCandidateBool')) ...
    || (IndicateDimer&&isfield(TR, 'StateSequence')));
for nn = 1:NTraj
    % Determine which datapoints we should plot.
    KeepBool = ((TR(nn).FrameNum>=FrameRange(1)) ...
        & (TR(nn).FrameNum<=FrameRange(2)));
    if ~any(KeepBool)
        continue
    end
    X = TR(nn).X(KeepBool);
    Y = TR(nn).Y(KeepBool);
    FrameNum = TR(nn).FrameNum(KeepBool);

    % Plot the trajectories, adding extra decorations for dimerization if
    % needed.
    if AddDimerDecorations
        % Add a dotted line for the entire trajectory (this should get
        % covered up by the other lines where relevant, otherwise it'll
        % stay as a dotted line to indicate parts that weren't dimer
        % candidates).
        line(PlotAxes, X, Y, FrameNum, ...
            'Color', Color(nn, :), 'LineStyle', ':');

        % Indicate dimer candidate frames, making sure to not have overlaps
        % (e.g., if there's a frame gap between candidate frames, don't
        % connect a line across the frame gap!).
        if isfield(TR, 'DimerCandidateBool')
            DimerCandidateBool = TR(nn).DimerCandidateBool(KeepBool);
            if any(DimerCandidateBool)
                % Find and plot each segment of dimer candidates.
                [CandidateStartInd, CandidateEndInd] = ...
                    smi_helpers.findStartEndInds(DimerCandidateBool);
                for ii = 1:numel(CandidateStartInd)
                    line(PlotAxes, ...
                        X(CandidateStartInd(ii):CandidateEndInd(ii)), ...
                        Y(CandidateStartInd(ii):CandidateEndInd(ii)), ...
                        FrameNum(CandidateStartInd(ii):CandidateEndInd(ii)), ...
                        'Color', Color(nn, :), ...
                        'LineStyle', '-', 'LineWidth', 2)
                end
            end
        else
            DimerCandidateBool = zeros(size(FrameNum), 'logical');
        end

        % Indicate dimer frames.
        if isfield(TR, 'StateSequence')
            DimerBool = (TR(nn).StateSequence(KeepBool) == 1);
            if any(DimerBool)
                % Find and plot each segment of dimer candidates.
                [DimerStartInd, DimerEndInd] = ...
                    smi_helpers.findStartEndInds(DimerBool);
                for ii = 1:numel(DimerStartInd)
                    line(PlotAxes, ...
                        X(DimerStartInd(ii):DimerEndInd(ii)), ...
                        Y(DimerStartInd(ii):DimerEndInd(ii)), ...
                        FrameNum(DimerStartInd(ii):DimerEndInd(ii)), ...
                        'Color', [0, 0, 1], ...
                        'LineStyle', '-', 'LineWidth', 2)
                end
            end
        else
            DimerBool = zeros(size(FrameNum), 'logical');
        end

        % Plot an invisible line over the whole trajectory for use as the
        % linehandle (which will allow us to click the trajectory in a
        % GUI).
        LineHandles(nn) = line(PlotAxes, ...
            X, Y, FrameNum, ...
            'Color', [0, 0, 0, 0], 'LineWidth', 2, ...
            'PickableParts', 'all');
    else
        % Plot the full trajectory as usual.
        LineHandles(nn) = line(PlotAxes, ...
            X, Y, FrameNum, ...
            'Color', Color(nn, :), varargin{:});
    end
end


end