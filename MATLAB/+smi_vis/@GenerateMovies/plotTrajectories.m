function [LineHandles] = plotTrajectories(PlotAxes, ...
    TR, FrameRange, MaxTrajLength, Color, varargin)
%plotTrajectories plots trajectories in the specified axes.
% This method plots trajectories in 'TR'.  This method is intended to
% remain "lightweight" for speed purposes, so minimal input checks/default
% settings should be implemented here.
%
% INPUTS:
%   PlotAxes: Axes in which the trajectories will be plotted.
%   TR: Tracking Results structure (see smi_core.TrackingResults).
%   FrameRange: Range of frames over which trajectories should be plotted.
%   MaxTrajLength: Maximum number of points to be plotted for each
%                  trajectory.
%   Color: Color of the trajectory lines. (NTrajectoriesx3(4) float array)
%   varargin: Additional inputs that can be provided as keyword arguments
%             for the MATLAB line() method (see Line Properties from line()
%             documentation for options).  For example, you can change
%             the Marker property to 'x' as plotTrajectories(PlotAxes, TR,
%             FrameRange, MaxTrajLength, Color, 'Marker', 'x')
%   Marker: Marker for each localization in the trajectory. (see Line
%           Properties for MATLAB method line())
%   LineStyle: Style of line for each trajectory. (see Line Properties
%              for MATLAB method line())
%
% OUTPUTS:
%   LineHandles: Array of line handles corresponding to the trajectories
%                in 'TR'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Make sure the input 'Color' is long enough (this is a convenient input
% check that doesn't add much overhead).
NTraj = numel(TR);
Color = repmat(Color, [min(NTraj-size(Color, 1)+1, NTraj), 1]);

% Loop through trajectories present in 'TR' and plot them, enforcing the
% requested 'FrameRange' and 'MaxTrajLength'.
LineHandles = gobjects(NTraj, 1);
for nn = 1:NTraj
    KeepBool = ((TR(nn).FrameNum>=FrameRange(1)) ...
        & (TR(nn).FrameNum<=FrameRange(2)));
    if ~any(KeepBool)
        continue
    end
    X = TR(nn).X(KeepBool);
    Y = TR(nn).Y(KeepBool);
    FrameNum = TR(nn).FrameNum(KeepBool);
    NLoc = numel(FrameNum);
    IndexArray = (1:min(NLoc, MaxTrajLength));
    LineHandles(nn) = line(PlotAxes, ...
        X(IndexArray), Y(IndexArray), FrameNum(IndexArray), ...
        'Color', Color(nn, :), varargin{:});
end


end