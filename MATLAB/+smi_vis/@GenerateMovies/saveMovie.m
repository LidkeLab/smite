function [] = saveMovie(obj, SavePath)
%saveMovie generates and saves a movie of 2D trajectories.
%
% INPUTS:
%   SavePath: Complete path to the movie file in which the movie will be
%             saved (Default is a file selection prompt).
%   Profile: Movie profile option sent to the video writer. See input
%            'profile' to MATLAB built in VideoWriter.
%            (Default = 'Motion JPEG AVI')

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('SavePath', 'var') || isempty(SavePath))
    [File, Path] = uiputfile({'*.mp4', 'MPEG-4 (*.mp4)'; ...
        '*.avi', 'Uncompressed AVI (*.avi)'});
    SavePath = fullfile(Path, File);
end

% Create the video writer.
[Path, File, FileExtension] = fileparts(SavePath);
switch FileExtension
    case '.mp4'
        obj.VideoObject = ...
            VideoWriter(fullfile(Path, File), 'MPEG-4');
        obj.VideoObject.Quality = 100;
    case '.avi'
        obj.VideoObject = ...
            VideoWriter(fullfile(Path, File), 'Uncompressed AVI');
end
obj.VideoObject.FrameRate = obj.Params.FrameRate;

% Create a fresh figure in which we'll prepare the movie for saving (this
% is done due to the usage of print() in obj.playMovie(), which will print
% the whole figure, and thus GUI elements if we prepare the movie within
% the GUI).
PlotFigure = figure();
OriginalAxes = obj.MovieAxes;
obj.MovieAxes = axes(PlotFigure);
obj.MovieAxes.Toolbar.Visible = 'off';

% Generate and save the movie.
obj.generateMovie()
close(PlotFigure)
obj.MovieAxes = OriginalAxes;


end