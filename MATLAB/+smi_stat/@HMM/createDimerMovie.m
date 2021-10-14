function [MovieParams] = createDimerMovie(MovieAxes, ...
    TRArray, RawDataChannel1, RawDataChannel2, MovieParams, VideoObject)
%createDimerMovie creates a movie of trajectories superimposed on raw data.
% This method will create a movie of a dimer event contained in TRArray
% superimposed on the images present in RawData.
%
% INPUTS:
%   MovieAxes: Axes object in which we've plotted stuff.
%              (Default = axes(figure()))
%   TRArray: A structure array of two TR structures, with each TR structure
%            consisting of one trajectory.
%   RawDataChannel1: 3D matrix containing the raw data of the first
%                    channel corresponding to the trajectory TRArray(1).
%   RawDataChannel2: 3D matrix containing the raw data of the second
%                    channel corresponding to the trajectory TRArray(2).
%   MovieParams: A structure of display parameters for the movie. See
%                smi_vis.GenerateMovies.prepDefaults() for options. Note
%                that some of the defaults are changed to be better suited
%                for dimer movies. Some additional fields relevant to dimer
%                data are also added to this structure:
%                CreateGrayscaleMovie: 0 to create color overlay movie.
%                                      1 to create side-by-side movie
%                                      (Default = 0)
%                ColorMap: ColorMap used to display trajectories.
%                          (Default = [0, 1, 0; 1, 0, 1], i.e., green
%                          channel 1 and magenta channel 2)
%                ChannelNames: name labels for each channel.
%                              (Default = {'Channel 1', 'Channel 2'})
%                IndicateDimer: When true, indicate dimer events with a
%                               special marking. (Default = true)
%                IndicateDimerCandidate: When true, indicate which the data
%                                        which was considered to be a dimer
%                                        candidate in pre-processsing.
%                                        (Default = true)
%   VideoObject: Video writer object defining the movie that will be saved
%                while preparing this movie (see MATLAB VideoWriter
%                object).  This object should be opened and closed outside
%                of this method for proper usage.
%                (Default = [] and the movie isn't saved)
%
% OUTPUTS:
%   DisplayParams: The same structure as the input structure, but with the
%                  default parameters added in.
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set default parameter values where needed.
if (~exist('MovieParams', 'var') || isempty(MovieParams))
    MovieParams = struct();
end
if (~exist('MovieAxes', 'var') || isempty(MovieAxes))
    MovieAxes = axes(figure());
end
if (~exist('SMD', 'var') || isempty(SMD))
    SMD = smi_core.SingleMoleculeData.createSMD();
end
if (~exist('VideoObject', 'var') || isempty(VideoObject))
    VideoObject = [];
end
DefaultParams = smi_vis.GenerateMovies.prepDefaults();
DefaultParams.CreateGrayscaleMovie = false;
DefaultParams.ChannelNames = {'Channel 1', 'Channel 2'};
DefaultParams.ColorMap = [0, 1, 0; 1, 0, 1];
DefaultParams.TrajColor = [0, 1, 0; 1, 0, 1];
DefaultParams.IndicateDimerCandidate = true;
DefaultParams.IndicateDimer = true;
DefaultParams.AutoClip = true;
MovieParams = smi_helpers.padStruct(MovieParams, DefaultParams);
MovieParams.ColorMap = smi_vis.contrastStretch(MovieParams.ColorMap);

% Prepare the raw data.
DataSize = size(RawDataChannel1);
ScaledDataCh1 = reshape(smi_vis.contrastStretch(RawDataChannel1), ...
    [DataSize(1:2), 1 ,DataSize(3)]);
ScaledDataCh2 = reshape(smi_vis.contrastStretch(RawDataChannel2), ...
    [DataSize(1:2), 1 ,DataSize(3)]);
if DefaultParams.CreateGrayscaleMovie
    % For grayscale movies, we'll concatenate the data side-by-side into
    % one large set of grayscale images.
    ScaledData = [ScaledDataCh1, ScaledDataCh2];
    ScaledData = repmat(ScaledData, 1, 1, 3, 1);
else
    ScaledDataCh1Color = cat(3, ...
        ScaledDataCh1 * MovieParams.ColorMap(1, 1), ...
        ScaledDataCh1 * MovieParams.ColorMap(1, 2), ...
        ScaledDataCh1 * MovieParams.ColorMap(1, 3));
    ScaledDataCh2Color = cat(3, ...
        ScaledDataCh2 * MovieParams.ColorMap(2, 1), ...
        ScaledDataCh2 * MovieParams.ColorMap(2, 2), ...
        ScaledDataCh2 * MovieParams.ColorMap(2, 3));
    ScaledData = ScaledDataCh1Color + ScaledDataCh2Color;
end

% Prepare the movie generator.
MovieGenerator = smi_vis.GenerateMovies;
MovieGenerator.TR = TRArray;
MovieGenerator.Params = MovieParams;
MovieGenerator.RawData = ScaledData;
MovieGenerator.prepRawData()
MovieGenerator.prepAxes(MovieAxes)

% Define a few other parameters that'll be used in this method.
% NOTE: XRange and YRange are defined to resolve some coordinate
%       differences between raw data plots and the real data.
IsRotating = (size(MovieParams.LineOfSite, 1) > 1);
ResolutionString = sprintf('-r%i', MovieParams.Resolution);

% If the VideoObject is non-empty, open it.
if ~isempty(VideoObject)
    open(VideoObject)
end

% Define some parameters used in the display of dimer data.
DimerParams = MovieParams;
DimerParams.TrajColor = [0, 0, 1; 0, 0, 1];
TRDimerCh1 = smi_core.SingleMoleculeData.isolateSubSMD(TRArray(1), ...
    TRArray(1).StateSequence == 1);
TRDimerCh2 = smi_core.SingleMoleculeData.isolateSubSMD(TRArray(2), ...
    TRArray(2).StateSequence == 1);
TRDimer = smi_core.TrackingResults.catTR(TRDimerCh1, TRDimerCh2);
DimerCandParams = MovieParams;
DimerCandParams.LineStyle = '-';
TRDimerCandCh1 = smi_core.SingleMoleculeData.isolateSubSMD(TRArray(1), ...
    TRArray(1).DimerCandidateBool);
TRDimerCandCh2 = smi_core.SingleMoleculeData.isolateSubSMD(TRArray(2), ...
    TRArray(2).DimerCandidateBool);
TRDimerCand = smi_core.TrackingResults.catTR(TRDimerCandCh1, TRDimerCandCh2);
RemainderParams = MovieParams;
RemainderParams.LineStyle = ':';
TRRemainderCh1 = smi_core.SingleMoleculeData.isolateSubSMD(TRArray(1), ...
    ~(TRArray(1).DimerCandidateBool | (TRArray(1).StateSequence==1)));
TRRemainderCh2 = smi_core.SingleMoleculeData.isolateSubSMD(TRArray(2), ...
    ~(TRArray(2).DimerCandidateBool | (TRArray(2).StateSequence==1)));
TRRemainder = smi_core.TrackingResults.catTR(TRRemainderCh1, TRRemainderCh2);

% Loop through the frames of raw data and prepare the movie.
ScaledData = MovieGenerator.ScaledData;
MovieFigure = MovieAxes.Parent;
for ff = MovieGenerator.Params.ZFrames(1):MovieGenerator.Params.ZFrames(2)
    % Make the current frame of the movie.
    smi_vis.GenerateMovies.makeFrame(MovieAxes, ...
        [], ScaledData(:, :, :, ff), MovieGenerator.Params, SMD, ff);
    smi_vis.GenerateMovies.plotTrajectories(MovieAxes, ...
        TRRemainder, [ff-RemainderParams.MaxTrajLength, ff], ...
        RemainderParams.MaxTrajLength, ...
        RemainderParams.TrajColor, 'LineStyle', RemainderParams.LineStyle);
    smi_vis.GenerateMovies.plotTrajectories(MovieAxes, ...
        TRDimerCand, [ff-DimerCandParams.MaxTrajLength, ff], ...
        DimerCandParams.MaxTrajLength, ...
        DimerCandParams.TrajColor, 'LineStyle', DimerCandParams.LineStyle);
    smi_vis.GenerateMovies.plotTrajectories(MovieAxes, ...
        TRDimer, [ff-DimerParams.MaxTrajLength, ff], ...
        DimerParams.MaxTrajLength, ...
        DimerParams.TrajColor, 'Marker', DimerParams.PlotMarker);
    
    % Update the line of site if needed.
    if IsRotating
        view(MovieAxes, MovieGenerator.Params.LineOfSite(ff, :))
    end
    
    % Update the axes to ensure all new objects and changes are present.
    drawnow()
    
    % If needed, write this movie frame.
    if isempty(VideoObject)
        % This movie isn't being saved: add a pause so the movie doesn't go
        % too fast.
        pause(1 / MovieGenerator.Params.FrameRate);
    else
        FrameData = print(MovieFigure, ...
            '-RGBImage', '-opengl', ResolutionString);
        VideoObject.writeVideo(FrameData);
    end
end

% Close the VideoObject.
if ~isempty(VideoObject)
    close(VideoObject)
end


end