function [RGBimage]=colorImage(Image,ColorMap,MinMax)
%colorImage generates RGB image with certain colormap from grey scale image
%
%   INPUT
%      Image - 2D array containing gray scale values
%      ColorMap (optional) - 3 by N matrix, interpreted as Red Green and
%                            Blue values, default is hot(256)
%      MinMax (optional) - [min max] limits for scaling colormap, default
%                          is min and max values of Image
%
%   OUTPUT
%      If no output args are given histogram image will be displayed with
%           hot colormap
%      RGBimage - RBG version of image with colormap
%
%   REQUIRES
%      Matlab 2014b or higher
%      Dipimage toolbox (http://www.diplib.org/)

% Created by
%    Marjolein Meddens, Lidke Lab 2017

% check input
if nargin < 1
     error('smi_vis:GenerateImages:colorImage:notEnoughInputArgs','Need at least Image input');
end
if nargin < 2
    ColorMap = hot(256);
end
if nargin < 3
    MinMax = [min(Image(:)) max(Image(:))];
end

% combine channels
Image=single(Image);
[r, g, b]=c_GenColorChannels(Image,ColorMap,MinMax(1),MinMax(2));
RGBimage=smi_vis.GenerateImages.rgbImage(r,g,b);
end
