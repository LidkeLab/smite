% Produce a movie from SPT results.
%
TrackingResults = 'NEEDS_TO_BE_SET!';
RawData = 'NEEDS_TO_BE_SET!'

% Load SMITE processed tracking results: SMF, SMD, TR.
load(TrackingResults);
% Load raw blinking data.
load(RawData);

%sliceViewer(Channel_1);

GM = smi_vis.GenerateMovies();
GM.SMF = SMF;
GM.SMD = SMD;
GM.TR  = TR;
GM.RawData = Channel_1;

% Set some parameters; see smi_vis.GenerateMovies.prepDefaults for a full list,
%
% Place the indicated marker on each localization given in GM.SMD
%GM.Params.SMDMarker = 'x';
% Maximum number of points of a given trajectory that will be plotted. 
GM.Params.MaxTrajLength = 5;
% Choose a ROI to examine (default is entire image):
%raw_size = size(GM.RawData);
%GM.Params.XPixels = [0, raw_size(1) - 1];
%GM.Params.YPixels = [0, raw_size(2) - 1];
GM.Params.XPixels = [64, 127];
GM.Params.YPixels = [64, 127];

%sliceViewer(Channel_1);

GM = smi_vis.GenerateMovies();
GM.SMF = SMF;
GM.SMD = SMD;
GM.TR  = TR;
GM.RawData = Channel_1;
% Set some parameters; see smi_vis,GenerateMovies.prepDefaults for a full list,
%
% Place the indicated marker on each localization given in GM.SMD
%GM.Params.SMDMarker = 'x';
% Maximum number of points of a given trajectory that will be plotted. 
GM.Params.MaxTrajLength = 5;
% Choose a ROI to examine (default is entire image):
%raw_size = size(GM.RawData);
%GM.Params.XPixels = [0, raw_size(1) - 1];
%GM.Params.YPixels = [0, raw_size(2) - 1];
GM.Params.XPixels = [64, 127];
GM.Params.YPixels = [64, 127];
% Choose from which range of frames to display trajectories (default is all).
%GM.Params.ZFrames = [1, raw_size(3)];
GM.Params.ZFrames = [1, 100];
GM.gui();
