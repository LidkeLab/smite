function [Params] = prepDefaults()
%prepDefaults prepares a default 'Params' structure.
% This method prepares a 'Params' structure with all fields present with
% default values.
%
% OUTPUTS:
%   Params: Structure of parameters that will be applied to movies.
%           UnitFlag: Boolean flag to specify the units of the movie                            false for camera units (pixels, frames)
%                     true for physical units (micrometers, seconds)
%                     (boolean)(Default = false)
%           MaxTrajLength: Maximum number of points of a given
%                          trajectory that will be plotted.
%                          (frames)(Default = inf)
%           XPixels: Range of X values shown in movie.
%                    ([MinPixel, MaxPixel])(pixels)
%                    (Default = [1, XSize])
%           YPixels: Range of Y values shown in movie.
%                    ([MinPixel, MaxPixel])(pixels)
%                    (Default = [1, YSize])
%           ZFrames: Range of frames shown in movie.
%                    ([MinFrame, MaxFrame])(frames)
%                    (Default = [1, NFrames])
%           MinScaleIntensity: Minimum value that we will use to scale the 
%                              raw data. If the raw data ranges from 
%                              [0, MaxIntensity], we scale as
%                              RawData / MaxIntensity, but if raw data is 
%                              too noisy this just brightens noise, so we
%                              set MaxIntensity = max(MaxIntensity,
%                              MinScaleIntensity).
%                              (numeric scalar)(Default = 1)
%           PercentileCeiling: Percentile ceiling of pixel values in the 
%                              raw data above which values are clipped.
%                               (percentage)(Default = 100)
%           PlotMarker: Marker to be placed on top of each localization. 
%                       Options: MATLAB markers, e.g., 'x' (see Line 
%                       Properties in MATLAB documentation) or 'none'.
%                        (char. array)(Default = 'none')
%           SMDMarker: Marker to be placed on top of each
%                      localization in SMD. The color is currently
%                      hard coded to be red. (Default = '.')
%           LineOfSite: Line of site vector for view() (rotates
%                       "camera" orientation for 3D movie.  When
%                       Is2D = 0, LineOfSite can be given as an
%                       NFramex2 array and the movie will rotate as
%                       prescribed by these angles.
%                       (1x2 array or NFramesx2 array)(degrees)
%                       (Default = [0, 90])
%           Is2D: false if movie should be 3D (LineOfSite~=[0, 90])
%                 true if movie should be 2D (might run faster for long
%                 movies)(Default = true)
%           Resolution: Resolution of the movie that gets saved to
%                       FilePath. This should be set to 0 when the movie is
%                       generated within a "decorated" figure (i.e., within
%                       a GUI), otherwise the decorations will be stored in
%                       the saved movie (see usage of getframe() vs print()
%                       in, e.g., smi_vis.GenerateMovies.playMovie()).
%                       (scalar)(dpi)(Default = 0, which sets to screen
%                       resolution)
%           FrameRate: Approximate playback framerate used when the 
%                      input 'VideoWriter' is not provided. 
%                      (Default = 10 fps)
%           TrajColor: Color of each trajectory in TR.
%                      (NTrajx3 numeric array)
%                      (Default = lines(NTraj))

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set some default parameters.
Params.UnitFlag = false;
Params.MaxTrajLength = inf;
Params.MinXYRange = 20;
Params.NPadPixels = 5;
Params.NPadFrames = 10;
Params.MinScaleIntensity = 1;
Params.PercentileCeiling = 100;
Params.PlotMarker = 'none';
Params.SMDMarker = '.';
Params.LineOfSite = [0, 90];
Params.Is2D = true;
Params.Resolution = 0;
Params.FrameRate = 10;
Params.ZFrames = [1, 1];
Params.XPixels = [1, 1];
Params.YPixels = [1, 1];
Params.TrajColor = [0, 0, 1];


end