function [] = generateMovie(obj, SavePath)
%generateMovie generates a movie of 2D trajectories over time.
% This method is intended to be a wrapper around playMovie() which will add
% some "decorations" (e.g., axis labels), prepare saved movies, etc.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Ensure some necessary parameters are given meaningful values based on the
% provided data/trajectories.
obj.setVitalParams()

% Rescale the raw data to improve the display.
obj.rescaleData()

% Prepare the axes based on settings in obj.Params.
obj.prepAxes()

% If a save path has been defined, prepare a video writer object.
if exist('SavePath', 'var')
    obj.VideoObject = VideoWriter(SavePath, 'MPEG-4');
    obj.VideoObject.Quality = 100;
end

% Play the movie.
obj.playMovie(obj.MovieAxes, obj.TR, obj.ScaledData, ...
    obj.Params, obj.SMF, obj.SMD, obj.VideoObject)


end