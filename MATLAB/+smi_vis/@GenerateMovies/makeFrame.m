function [] = makeFrame(PlotAxes, TR, ScaledData, Params, SMD, Frame)
%makeFrame plots a single frame for the trajectory movie.
% This method plots a single frame of the movie being prepared.  This
% method is meant to be lightweight without unnecessary code/default
% settings.
%
% INPUTS:
%   PlotAxes: Axes in which the trajectories will be plotted.
%   TR: Tracking Results structure (see smi_core.TrackingResults).
%   ScaledData: Individual image corresponding to the trajectories in
%               'TR' for frame 'Frame'.
%   Params: Structure of display parameters that will be applied to
%           the movie (see smi_vis.GenerateMovies.prepDefaults()).
%   SMD: Single Molecule Data structure containing additional localizations
%        that we want to show in the movie (e.g., this might be the results
%        from a box finding algorithm, where the localizations aren't
%        associated as trajectories yet).
%   Frame: Frame of the movie that will be plotted.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Define a few other parameters that'll be used in this method.
% NOTE: XRange and YRange are defined to resolve some coordinate
%       differences between raw data plots and the real data.
XRange = Params.XPixels + [0, 1];
YRange = Params.YPixels + [0, 1];
ZPosition = repmat(min(Params.ZFrames), [2, 2]);

% Display the raw data in the axes.
ScaledData = repmat(ScaledData, [1, 1, 3]);
surface(PlotAxes, XRange, YRange, ZPosition, ...
    ScaledData, 'facecolor', 'texturemap')

% Plot the trajectories.
smi_vis.GenerateMovies.plotTrajectories(PlotAxes, ...
    TR, [Frame-Params.MaxTrajLength, Frame], Params.MaxTrajLength, ...
    Params.TrajColor, 'Marker', Params.PlotMarker);

% Plot the points in the optional input 'SMD'.
smi_vis.GenerateMovies.plotTrajectories(PlotAxes, ...
    SMD, [Frame-Params.MaxTrajLength, Frame], Params.MaxTrajLength, ...
    Params.TrajColor, 'Marker', Params.SMDMarker);


end