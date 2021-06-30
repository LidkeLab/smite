function playMovie(PlotAxes, TR, RawData, Params, SMD, VideoObject)
%playMovie prepares a movie of the trajectories in TR.
% This method prepares a movie of single-particle tracking trajectories for
% 2D tracking.
%
% INPUTS:
%   PlotAxes: Axes in which the trajectories will be plotted.
%             (Default = gca())
%   TR: Tracking Results structure (see smi_core.TrackingResults).
%   RawData: Raw data corresponding to the trajectories in 'TR'.  The third
%            dimension index is assumed to correspond directly with the
%            frame number.
%            (Default = zeros(TR(1).YSize, TR(1).XSize, TR(1).NFrames))
%   Params: Structure of display parameters that will be applied to
%           the movie (see smi_vis.GenerateMovies.prepDefaults()).
%   SMD: Single Molecule Data structure containing additional localizations
%        that we want to show in the movie (e.g., this might be the results
%        from a box finding algorithm, where the localizations aren't
%        associated as trajectories yet). (Default is an empty SMD)
%   VideoObject: Video writer object defining the movie that will be saved
%                while preparing this movie (see MATLAB VideoWriter object)
%                (Default = [] and the movie isn't saved)
%
% REQUIRES:
%   Image Processing Toolbox (to use padarray())

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Grab some parameters from 'TR' which may be empty.  If they are empty,
% we'll need to come up with reasonable guesses for these values.
NFrames = TR(1).NFrames;
XSize = TR(1).XSize;
YSize = TR(1).YSize;
if isempty(NFrames)
    NFrames = max(cell2mat({TR.FrameNum}.'));
end
if isempty(XSize)
    XSize = max(cell2mat({TR.X}.'));
end
if isempty(YSize)
    YSize = max(cell2mat({TR.Y}.'));
end

% Set default parameter values where needed.
if (~exist('Params', 'var') || isempty(Params))
    Params = struct();
end
if (~exist('PlotAxes', 'var') || isempty(PlotAxes))
    PlotAxes = gca();
end
if ~exist('RawData', 'var')
    RawData = zeros(YSize, XSize, NFrames);
end
if (~exist('SMD', 'var') || isempty(SMD))
    SMD = smi_core.SingleMoleculeData.createSMD();
end
if (~exist('VideoWriter', 'var') || isempty(VideoObject))
    VideoObject = [];
end
DefaultParams = smi_vis.GenerateMovies.prepDefaults();
NTraj = numel(TR);
DefaultParams.ZFrames = [1, NFrames];
DefaultParams.XPixels = [1, XSize];
DefaultParams.YPixels = [1, YSize];
DefaultParams.TrajColor = lines(NTraj);
Params = smi_helpers.padParams(Params, DefaultParams);

% Define a few other parameters that'll be used in this method.
% NOTE: XRange and YRange are defined to resolve some coordinate
%       differences between raw data plots and the real data.
XRange = Params.XPixels + [0, 1]; 
YRange = Params.YPixels + [0, 1];
ZPosition = repmat(min(Params.ZFrames), [2, 2]);
IsRotating = (size(Params.LineOfSite, 1) > 1);

% Rescale the raw data after isolating the portion that will be displayed.
RawData = ...
    padarray(RawData, [0, 0, min(0, NFrames-size(RawData, 3))], 'post');
if ~(isempty(Params.XPixels) && isempty(Params.YPixels) ...
        && isempty(Params.ZFrames))
    RawData = RawData(Params.XPixels(1):Params.XPixels(2), ...
        Params.YPixels(1):Params.YPixels(2), ...
        Params.ZFrames(1):Params.ZFrames(2));
end
RawData = smi_vis.contrastStretch(single(RawData), [0; 1], ...
    Params.PercentileCeiling, Params.MinScaleIntensity);

% Prepare the designated axes.
PlotAxes.ActivePositionProperty = 'position';
PlotAxes.DataAspectRatio = ...
    [1, 1, max(Params.ZFrames)/max([Params.XPixels, Params.YPixels])];
PlotAxes.YDir = 'reverse';
PlotAxes.XLim = XRange;
PlotAxes.YLim = YRange;
PlotAxes.ZLim = Params.ZFrames;
view(PlotAxes, Params.LineOfSite(1, :));
colormap(PlotAxes, 'gray')

% Loop through the frames of raw data and prepare the movie.
for ff = 1:NFrames
    % Clear the axes to make sure deleted objects aren't accumulating
    % (which slows things down a lot!).
    cla(PlotAxes)
    
    % Display the raw data in the axes.
    CurrentData = repmat(RawData(:, :, ff), [1, 1, 3]);
    surface(PlotAxes, XRange, YRange, ZPosition, ...
        CurrentData, 'facecolor', 'texturemap')
    
    % Plot the trajectories.
    smi_vis.GenerateMovies.plotTrajectories(PlotAxes, ...
        TR, [ff-Params.MaxTrajLength, ff], Params.MaxTrajLength, ...
        Params.TrajColor, 'Marker', Params.PlotMarker);
    
    % Plot the points in the optional input 'SMD'.
    smi_vis.GenerateMovies.plotTrajectories(PlotAxes, ...
        SMD, [ff-Params.MaxTrajLength, ff], Params.MaxTrajLength, ...
        Params.TrajColor, 'Marker', Params.SMDMarker);
    
    % Update the line of site if needed.
    if IsRotating
        view(PlotAxes, Params.LineOfSite(ff, :))
    end
    
    % Update the axes to ensure all new objects and changes are present.
    drawnow()
    
    % If needed, write this movie frame.
    if isempty(VideoObject)
        % This movie isn't being saved: add a pause so the movie doesn't go
        % too fast.
        pause(1 / Params.FrameRate);
    else
        FrameData = getframe(PlotAxes);
        VideoObject.writeVideo(FrameData.cdata);
    end
end


end