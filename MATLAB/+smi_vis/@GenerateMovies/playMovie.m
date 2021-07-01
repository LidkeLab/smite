function [] = playMovie(PlotAxes, TR, ScaledData, Params, SMD, VideoObject)
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
IsRotating = (size(Params.LineOfSite, 1) > 1);
CustomRes = (Params.Resolution ~= 0);
ResolutionString = sprintf('-r%i', Params.Resolution);

% Loop through the frames of raw data and prepare the movie.
AxesParent = PlotAxes.Parent;
for ff = Params.ZFrames(1):Params.ZFrames(2)           
    % Make the current frame of the movie.
    smi_vis.GenerateMovies.makeFrame(PlotAxes, TR, ScaledData(:, :, ff), ...
        Params, SMD, ff);
        
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