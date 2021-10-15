function [] = generateMovie(obj)
%generateMovie generates a movie of 2D trajectories over time.
% This method is intended to be a wrapper around playMovie() which will
% do a couple of useful things before playing the movie.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Rescale/crop the raw data to improve the display.
obj.prepRawData()

% Prepare the axes based on settings in obj.Params.
obj.prepAxes()

% Play the movie.
obj.playMovie(obj.MovieAxes, obj.TR, obj.ScaledData, ...
    obj.Params, obj.SMF, obj.SMD, obj.VideoObject)


end