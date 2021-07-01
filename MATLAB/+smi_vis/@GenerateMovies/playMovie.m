function [] = playMovie(PlotAxes, TR, ScaledData, Params, SMF, SMD, ...
    VideoObject)
%playMovie prepares a movie of the trajectories in TR.
% This method prepares a movie of single-particle tracking trajectories for
% 2D tracking.
%
% INPUTS:
%   PlotAxes: Axes in which the trajectories will be plotted.
%             (Default = gca())
%   TR: Tracking Results structure (see smi_core.TrackingResults).
%   ScaledData: Data corresponding to the trajectories in 'TR'.  The third
%               dimension index is assumed to correspond directly with the
%               frame number.
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

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set default parameter values where needed.
if (~exist('Params', 'var') || isempty(Params))
    Params = struct();
end
if (~exist('PlotAxes', 'var') || isempty(PlotAxes))
    PlotAxes = gca();
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
Params = smi_helpers.padParams(Params, DefaultParams);

% Define a few other parameters that'll be used in this method.
% NOTE: XRange and YRange are defined to resolve some coordinate
%       differences between raw data plots and the real data.
XRange = Params.XPixels + [0, 1]; 
YRange = Params.YPixels + [0, 1]; 
IsRotating = (size(Params.LineOfSite, 1) > 1);
CustomRes = (Params.Resolution ~= 0);
ResolutionString = sprintf('-r%i', Params.Resolution);

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
    xtickformat(PlotAxes, '%.1f')
    ytickformat(PlotAxes, '%.1f')
    ztickformat(PlotAxes, '%.1f')
    XTicks = (PlotAxes.XTick-1) * SMF.Data.PixelSize;
    YTicks = (PlotAxes.YTick-1) * SMF.Data.PixelSize;
    ZTicks = (PlotAxes.ZTick-1) / SMF.Data.FrameRate;
    PlotAxes.XTickLabel = num2str(XTicks.', '%.1f');
    PlotAxes.YTickLabel = num2str(YTicks.', '%.1f');
    PlotAxes.ZTickLabel = num2str(ZTicks.', '%.1f');
else
    xtickformat(PlotAxes, '%i')
    ytickformat(PlotAxes, '%i')
    ztickformat(PlotAxes, '%i')
end

% Loop through the frames of raw data and prepare the movie.
AxesParent = PlotAxes.Parent;
for ff = Params.ZFrames(1):Params.ZFrames(2)
    % Clear the axes to make sure deleted objects aren't accumulating
    % (which slows things down a lot!).
    cla(PlotAxes)
            
    % Add the data and trajectories to the axes.
    smi_vis.GenerateMovies.makeFrame(PlotAxes, TR, ScaledData(:, :, ff), ...
        Params, SMD, ff)
        
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
        if CustomRes
            FrameData = print(AxesParent, ...
                '-RGBImage', '-opengl', ResolutionString);
        else
            % If the resolution is set to 0 (which uses a default screen
            % resolution in print() above), we should just use getframe(),
            % which does something similar but seems to be faster.
            FrameData = getframe(AxesParent);
        end
        VideoObject.writeVideo(FrameData);
    end
end


end