% c_GenColorChannels C++ function for applying colormap to greyscale image
% 
% SYNTAX
%   [R,G,B] = c_GenColorChannels(Image,ColorMap,MinValue,MaxValue);
% INPUTS
%   Image - 2D array (single floats)
%   ColorMap - N by 3 array (double floats) describing color values, 
%              columns represent Red, Green and Blue values for color 
%              levels (see help colormap)
%   MinValue (optional) - minimum gray scale value for displaying colormap
%                         (default min(Image)) 
%   MaxValue (optional) - maximum grey scale value for displaying colormap
%                         (default max(Image)) 
% OUTPUTS
%   R - 2D array representing Red values
%   G - 2D array representing Green values
%   B - 2D array representing Blue values
%
% Dipimage's joinchannel function can be used to combine R, G and B into an
% RGB image
%

% Marjolein Meddens, Lidke Lab 2017