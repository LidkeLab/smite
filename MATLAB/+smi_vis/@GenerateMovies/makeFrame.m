function [LineHandles] = makeFrame(PlotAxes, ...
    TR, ScaledData, Params, SMF, SMD, Frame)
%makeFrame plots a single frame for the trajectory movie.
% This method plots a single frame of the movie being prepared.
%
% INPUTS:
%   PlotAxes: Axes in which the trajectories will be plotted.
%   TR: Tracking Results structure (see smi_core.TrackingResults).
%   ScaledData: Individual image corresponding to the trajectories in
%               'TR' for frame 'Frame'.
%   Params: Structure of display parameters that will be applied to
%           the movie (see smi_vis.GenerateMovies.prepDefaults()).
%   SMF: Single Molecule Fitting structure (used for Data.PixelSize and
%        Data.FrameRate).
%   SMD: Single Molecule Data structure containing additional localizations
%        that we want to show in the movie (e.g., this might be the results
%        from a box finding algorithm, where the localizations aren't
%        associated as trajectories yet).
%   Frame: Frame of the movie that will be plotted.
%
% OUTPUTS:
%   LineHandles: Array of line handles corresponding to the trajectories
%                in 'TR'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Clear the axes to make sure deleted objects aren't accumulating (which 
% can slow things down when trying to capture the frame by, e.g.,
% getframe()).
cla(PlotAxes)

% Display the raw data in the axes.
surface(PlotAxes, Params.XPixels + [0, 1], Params.YPixels + [0, 1], ...
    repmat(min(Params.ZFrames), [2, 2]), ...
    ScaledData, 'facecolor', 'texturemap')

% Plot the trajectories.
[LineHandles] = smi_vis.GenerateMovies.plotTrajectories(PlotAxes, ...
    TR, [Frame-Params.MaxTrajLength, Frame], Params.MaxTrajLength, ...
    Params.TrajColor, 'Marker', Params.PlotMarker);

% Plot the points in the optional input 'SMD'.
smi_vis.GenerateMovies.plotTrajectories(PlotAxes, ...
    SMD, [Frame-Params.MaxTrajLength, Frame], Params.MaxTrajLength, ...
    Params.SMDColor, 'Marker', Params.SMDMarker, 'LineStyle', 'none');

% Add a timestamp.
if Params.AddTimeStamp
    smi_vis.GenerateMovies.addTimeStamp(PlotAxes, Frame, ...
        SMF.Data.FrameRate, Params);
end


end