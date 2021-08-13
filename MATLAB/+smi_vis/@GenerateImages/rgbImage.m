function [RGBimage]=rgbImage(R,G,B)
%rgbImage generates RGB image with certain colormap from grey scale image
%
% If no output is requested, the image is displayed in a figure
%
%   INPUT
%      R:   2D Array for Red Image
%      G:   2D Array for Red Image
%      B:   2D Array for Red Image
%   OUTPUT
%      RGBimage - RBG array scaled from 0 to 1
%
%   REQUIRES
%      Matlab 2014b or higher
%      Dipimage toolbox (http://www.diplib.org/)

% check input
if nargin < 3
    error('smi_vis.GenerateImages.rgbImage:: Must have 3 inputs: R,G,B')
end

RGBimage=rescale(cat(3,R,G,B));

if nargout<1
    figure
    imagesc(RGBimage)
    set(gca,'visible','off')
    axis equal
    axis tight
end


end

