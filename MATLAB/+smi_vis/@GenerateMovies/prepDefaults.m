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
%                    ([MinPixel, MaxPixel])(pixels)(Default = [])
%                    NOTE: XPixels, YPixels, and ZFrames are kept as
%                          separate parameters (i.e., instead of lumping
%                          into one ROI parameter) to allow one to be set
%                          without having to set the others.
%           YPixels: Range of Y values shown in movie.
%                    ([MinPixel, MaxPixel])(pixels)
%                    (Default = [])
%           ZFrames: Range of frames shown in movie.
%                    ([MinFrame, MaxFrame])(frames)
%                    (Default = [])
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
%           PercentileFloor: Percentile floor of pixel values in the raw
%                            data below which values are clipped.
%                            (percentage)(Default = 0)
%           PlotMarker: Marker to be placed on top of each localization. 
%                       Options: MATLAB markers, e.g., 'x' (see Line 
%                       Properties in MATLAB documentation) or 'none'.
%                        (char. array)(Default = 'none')
%           SMDMarker: Marker to be placed on top of each
%                      localization in SMD. (Default = 'none')
%           SMDColor: Color of the SMD markers. (Default = [1, 0, 0])
%           LineOfSite: Line of site vector for view() (rotates
%                       "camera" orientation for 3D movie.  When
%                       Is2D = 0, LineOfSite can be given as an
%                       NFramex2 array and the movie will rotate as
%                       prescribed by these angles.
%                       (1x2 array or NFramesx2 array)(degrees)
%                       (Default = [0, 90])
%           Resolution: Resolution of the movie that gets saved to
%                       FilePath.
%                       (scalar)(dpi)(Default = 0, which sets to screen
%                       resolution)
%           FrameRate: Approximate playback framerate used when not saving
%                      a movie (i.e., if the movie is being displayed but
%                      not saved, an additional pause(1/FrameRate) is added
%                      between each frame). (Default = 10 fps)
%           AutoCrop: Flag to request cropping to data present in TR.  This
%                     is enforced in GenerateMovies.prepRawData() and 
%                     nowhere else as of this writing. (Default = false)
%           MinXYRange: Minimum XY range allowed when AutoClip = true,
%                       otherwise this is not used. (pixels)(Default = 20)
%           NPadPixels: Number of padding pixels added to edges of XY data
%                       when AutoClip = true, otherwise this is not used.
%                       (pixels)(Default = 5)
%           NPadFrames: Number of padding frames added around frames of
%                       trajectories in TR when AutoClip = true, otherwise
%                       this is not used. (pixels)(Default = 10)
%           AddTimeStamp: Add a frame/time stamp in a corner of the movie.
%                         (Default = false)
%           ChannelNames: Label names for each channel. 
%                         (Default = {'Channel 1'; 'Channel 2'})
%           RawDataColors: Color of raw data in each channel.
%                          (NChannelsx3 numeric array, [RGB])
%                          (Default = [0, 1, 0; 1, 0, 1])
%           TrajColor: Color of each trajectory in TR.
%                      (NTrajx3 numeric array)
%                      (Default = [], values set automatically elsewhere)
%           CropToDimerCandidates: 0 for normal behavior.
%                                  1 if 'AutoCrop' should apply temporally
%                                  to only those frames in which the
%                                  trajectory was a dimer candidate.
%                                  (Default = 0)
%           IndicateDimer: 0 if you don't want special dimer marker
%                          1 if you want to indicate dimer events
%                          If the field TR.StateSequence doesn't exist/is
%                          empty, this property doesn't do anything.
%                          (Default = 0)
%           IndicateDimerCandidate: 0 if you don't want to indicate the 
%                                       section of data considered as a
%                                       dimer candidate.
%                                   1 if you want to indicate the section
%                                     of data considered as a dimer
%                                     candidate.
%                                   If the field TR.DimerCandidateBool 
%                                   doesn't exist/is empty, this property 
%                                   doesn't do anything.
%                                   (Default = 0)

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Set some default parameters.
Params.UnitFlag = false;
Params.MaxTrajLength = inf;
Params.ZFrames = [];
Params.XPixels = [];
Params.YPixels = [];
Params.MinScaleIntensity = 1;
Params.PercentileCeiling = 100;
Params.PercentileFloor = 0;
Params.PlotMarker = 'none';
Params.SMDMarker = 'none';
Params.SMDColor = [1, 0, 0];
Params.LineOfSite = [0, 90];
Params.Resolution = 0;
Params.FrameRate = 10;
Params.AutoCrop = false;
Params.MinXYRange = 20;
Params.NPadPixels = 5;
Params.NPadFrames = 10;
Params.AddTimeStamp = false;
Params.RawDataColors = [0, 1, 0; 1, 0, 1];
Params.TrajColor = [];
Params.ChannelNames = {'Channel 1'; 'Channel 2'};
Params.CropToDimerCandidates = false;
Params.IndicateDimer = false;
Params.IndicateDimerCandidate = false;


end