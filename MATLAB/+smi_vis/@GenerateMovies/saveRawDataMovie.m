function [] = saveRawDataMovie(RawData, FilePath, Params, FrameRate)
%saveRawDataMovie generates and saves a movie of raw data.
% This method generates and saves a movie of the raw data provided in
% 'RawData', with some parameters from 'Params' applied to scale the data.
%
% INPUTS:
%   RawData: Stack of raw data used to generate the movie.
%            (YxXx1xFrames or YxXx3xFrames, see MATLAB's writeVideo() input
%            'img' for more details)
%   FilePath: Path to the file in which the movie will be saved. The file
%             can have either the .mp4 extension or the .avi extension.
%             (Default saved into current directory as .mp4).
%   Params: Structure of parameters that may be used to scale the data.
%           (see GenerateMovies.prepDefaults() for the full set of
%           parameters, though aren't used in this method).
%   FrameRate: Physical frame-rate that should be passed along if
%              Params.UnitFlag is true. (Default = 1)(frames per second)

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Set defaults.
if (~exist('FilePath', 'var') || isempty(FilePath))
    FilePath = fullfile(pwd(), [smi_helpers.genTimeString(), '_movie.mp4']);
end
if (~exist('Params', 'var') || isempty(Params))
    Params = smi_vis.GenerateMovies.prepDefaults();
else
    Params = smi_helpers.padStruct(Params, smi_vis.GenerateMovies.prepDefaults());
end
if (~exist('FrameRate', 'var') || isempty(FrameRate))
    FrameRate = 1;
end

% Prepare the raw data.
assert(ismember(size(RawData, 3), [1; 3]), ...
    'Raw data must be provided as a YxXx1xFrames or YxXx3xFrames array.')
RawData = smi_vis.contrastStretch(single(RawData), [0; 1], ...
    Params.PercentileCeiling, Params.MinScaleIntensity);

% Create the video writer.
[Path, File, FileExtension] = fileparts(FilePath);
switch FileExtension
    case '.mp4'
        VideoObject = ...
            VideoWriter(fullfile(Path, File), 'MPEG-4');
        VideoObject.Quality = 100;
    case '.avi'
        VideoObject = ...
            VideoWriter(fullfile(Path, File), 'Uncompressed AVI');
end
VideoObject.FrameRate = Params.FrameRate;

% Generate the movie.
open(VideoObject)
PlotFigure = figure('Resize', 'off');
PlotAxes = axes(PlotFigure);
PlotAxes.YDir = 'reverse';
ResolutionString = sprintf('-r%i', Params.Resolution);
Params.XPixels = [1, size(RawData, 2)];
Params.YPixels = [1, size(RawData, 1)];
Params.ZFrames = [1, size(RawData, 4)];
for ff = 1:size(RawData, 4)
    imshow(RawData(:, :, :, ff), [0, 1], ...
        'Parent', PlotAxes, 'Border', 'tight')
    if Params.AddTimeStamp
        smi_vis.GenerateMovies.addTimeStamp(PlotAxes, ff, FrameRate, Params)
    end
    PlotAxes.Visible = 'off'; % must be w/in loop as of R2021b
    FrameData = print(PlotFigure, '-RGBImage', '-opengl', ResolutionString);
    VideoObject.writeVideo(FrameData);
end
close(VideoObject)
close(PlotFigure)


end