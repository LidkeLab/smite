function [HistIm,RgbHistIm] = histogramImage(SMR, SRImageZoom, ColorMap)
%histogramImage generates 2D histogram image from SR data.
%
%   INPUT
%      SMR - structure containing at least fields 'X' and 'Y' having the x
%            and y coordinates of localizations to be plotted
%      SRImageZoom - zoom factor for histogram image (default = 10)
%      ColorMap (optional) - colormap for RBGimage 
%
%   OUTPUT
%      If no output args are given histogram image will be displayed with
%           hot colormap
%      HistIm - gray value histogram image
%      RgbHistIm - RBG histogram image with colormap
%
%   REQUIRES
%      Matlab 2014b or higher
%      Dipimage toolbox (http://www.diplib.org/)

% Created by:
%    Marjolein Meddens, Lidke Lab 2017

% check input
if nargin <1
    error('smi_vis:GenerateImages::histogramImage:notEnoughInputArgs','Need input of atleast SMR');
end

if nargin<2
    disp('No input for SRImageZoom, default value of 10 assigned');
    SRImageZoom=10;
end

if ~isfield(SMR,'X') || ~isfield(SMR,'Y')
    error('smi_vis:GenerateImages::histogramImage:noXYfields','First input should be structure with fields X and Y');
end
if nargin < 3
    ColorMap = 'hot';
end
% calculate size of histogram image
Xsize=(SMR.XSize*SRImageZoom);
Ysize=(SMR.YSize*SRImageZoom);
% transfer coordinates to match histogram image size
X=single(SMR.X*SRImageZoom);
Y=single(SMR.Y*SRImageZoom);
% generate histogram image (using c-function)
HistIm = c_HistRecon(Ysize,Xsize,Y,X,0);
% display image if no output args are given
[~,b]=hist(single(HistIm(HistIm>0)));
imagemax=b(1);
if nargout==0
    h = dipshow(HistIm);
    colormap(ColorMap);
    dipmapping(h,[0 imagemax]);
end
% create RGB image
if nargout > 1
    cm = eval([ColorMap '(256)']);
    MinMax = [0,imagemax];
    RgbHistIm = smi_vis.GenerateImages.colorImage(HistIm,cm,MinMax);
end

end

