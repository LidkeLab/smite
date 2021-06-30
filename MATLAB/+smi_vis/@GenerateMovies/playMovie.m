function playMovie(PlotAxes, TR, RawData, Params, SMF, SMD, VideoObject)
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
%   SMF: Single Molecule Fitting structure with SMF.Data.PixelSize and
%        SMF.Data.FrameRate populated.
%        (Default = smi_core.SingleMoleculeFitting)
%   SMD: Single Molecule Data structure containing additional localizations
%        that we want to show in the movie (e.g., this might be the results
%        from a box finding algorithm, where the localizations aren't
%        associated as trajectories yet). (Default is an empty SMD)
%   VideoObject: Video writer object defining the movie that will be saved
%                while preparing this movie (see MATLAB VideoWriter
%                object).  This object should be opened and closed outside
%                of this method for proper usage.
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
if (~exist('RawData', 'var') || isempty(RawData))
    RawData = zeros(YSize, XSize, NFrames);
end
if (~exist('SMF', 'var') || isempty(SMF))
    SMF = smi_core.SingleMoleculeFitting;
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
HiRes = (Params.Resolution ~= 0);
ResolutionString = sprintf('-r%i', Params.Resolution);

% Rescale the raw data after isolating the portion that will be displayed.
RawData = ...
    padarray(RawData, [min(0, YSize-size(RawData, 1)), ...
    min(0, XSize-size(RawData, 2)), ...
    min(0, NFrames-size(RawData, 3))], 'post');
if ~(isempty(Params.XPixels) && isempty(Params.YPixels) ...
        && isempty(Params.ZFrames))
    RawData = RawData(Params.YPixels(1):Params.YPixels(2), ...
        Params.XPixels(1):Params.XPixels(2), ...
        Params.ZFrames(1):Params.ZFrames(2));
end
RawData = smi_vis.contrastStretch(single(RawData), [0; 1], ...
    Params.PercentileCeiling, Params.MinScaleIntensity);

% Prepare the designated axes.
PlotAxes.ActivePositionProperty = 'position';
PlotAxes.DataAspectRatio = [1, 1, ...
    max(Params.ZFrames)/max([Params.XPixels, Params.YPixels])];
PlotAxes.YDir = 'reverse';
PlotAxes.XLimMode = 'manual';
PlotAxes.YLimMode = 'manual';
PlotAxes.ZLimMode = 'manual';
PlotAxes.XLim = XRange;
PlotAxes.YLim = YRange;
PlotAxes.ZLim = Params.ZFrames;
view(PlotAxes, Params.LineOfSite(1, :));
colormap(PlotAxes, 'gray')

% Revise axes ticks based on the unit flag.
PlotAxes.XTickMode = 'manual';
PlotAxes.YTickMode = 'manual';
PlotAxes.ZTickMode = 'manual';
PlotAxes.XTick = linspace(PlotAxes.XLim(1), PlotAxes.XLim(2), 5);
PlotAxes.YTick = linspace(PlotAxes.YLim(1), PlotAxes.YLim(2), 5);
PlotAxes.ZTick = linspace(PlotAxes.ZLim(1), PlotAxes.ZLim(2), 5);
if Params.UnitFlag
    PlotAxes.XTickLabelMode = 'manual';
    PlotAxes.YTickLabelMode = 'manual';
    PlotAxes.ZTickLabelMode = 'manual';
    XTicks = (PlotAxes.XTick-1) * SMF.Data.PixelSize;
    YTicks = (PlotAxes.YTick-1) * SMF.Data.PixelSize;
    ZTicks = (PlotAxes.ZTick-1) / SMF.Data.FrameRate;
    PlotAxes.XTickLabel = num2str(XTicks, '%.1f');
    PlotAxes.YTickLabel = num2str(YTicks, '%.1f');
    PlotAxes.ZTickLabel = num2str(ZTicks, '%.1f');
    xtickformat(PlotAxes, '%.1f')
    ytickformat(PlotAxes, '%.1f')
    ztickformat(PlotAxes, '%.1f')
else
    xtickformat(PlotAxes, '%i')
    ytickformat(PlotAxes, '%i')
    ztickformat(PlotAxes, '%i')
end

% Loop through the frames of raw data and prepare the movie.
for ff = Params.ZFrames(1):Params.ZFrames(2)
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
        if HiRes
            FrameData = print(PlotFigure, ...
                '-RGBImage', '-opengl', ResolutionString);
        else
            FrameData = getframe(PlotAxes);
        end
        VideoObject.writeVideo(FrameData);
    end
end


end