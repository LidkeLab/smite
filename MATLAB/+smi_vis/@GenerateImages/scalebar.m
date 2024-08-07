function [ImageOut,Image] = scalebar(Image,PixelSize,Length,Location)
%scalebar puts a scalebar of a desired length on the input image.
%   
%   Creates a scale bar of desired length on the input image. 
%   
% INPUTS:
%   Image:      Input image
%   PixelSize:  Pixel size of input image (in um)
%   Length:     Desired length of the scalebar (in um)
%   Location:   Desired location of scalebar on image. Options are
%   {topleft, topright, bottomleft, bottomright}. Default is set as
%   'bottomright'.
%
% OUTPUTS:
%   Imageout:   Output image with scalebar superimposed
%   Image:      The same as input image
%
% REQUIRES:
%   Matlab
%   Statistics Toolbox

% Created by
%   Sandeep Pallikkuth, Lidke Lab, 2017       

if nargin<3
    msg1='Not enough input. Please input image, pixel size, scalebar length and scale bar location !';
    error(msg1);
elseif nargin<4
    Location='bottomright';
end

% Processing input image
sizein=size(Image);     % determine size of image
channels=length(sizein);    
in=Image;

% Scalebar Color
sbcolor=[255 255 255];    %white scalebar

% Calculate scale bar length
if PixelSize>Length
    disp('Input scalebar length smaller than pixel size. Scalebar length re-sized to 10 pixels');
    barlength=10;
else
    barlength=round(Length/PixelSize);
end

%fprintf( ...
% 'scalebar: size = %d, PixelSize = %g um, Length = %d um, barlength = %d\n',...
%        sizein(1), PixelSize, Length, barlength);


% Calculate scale bar width
barwidth=round(sizein(2)*0.015);

%
if channels==3
    switch Location
        case 'bottomright'
            in(round(0.9*sizein(1))-barwidth:round(0.9*sizein(1)),round(0.9*sizein(2))-barlength:round(0.9*sizein(2)),1)=sbcolor(1);
            in(round(0.9*sizein(1))-barwidth:round(0.9*sizein(1)),round(0.9*sizein(2))-barlength:round(0.9*sizein(2)),2)=sbcolor(2);
            in(round(0.9*sizein(1))-barwidth:round(0.9*sizein(1)),round(0.9*sizein(2))-barlength:round(0.9*sizein(2)),3)=sbcolor(3);
        case 'bottomleft'
            in(round(0.9*sizein(1))-barwidth:round(0.9*sizein(1)),round(0.1*sizein(2)):round(0.1*sizein(2))+barlength,1)=sbcolor(1);
            in(round(0.9*sizein(1))-barwidth:round(0.9*sizein(1)),round(0.1*sizein(2)):round(0.1*sizein(2))+barlength,2)=sbcolor(2);
            in(round(0.9*sizein(1))-barwidth:round(0.9*sizein(1)),round(0.1*sizein(2)):round(0.1*sizein(2))+barlength,3)=sbcolor(3);
        case 'topright'
            in(round(0.1*sizein(1))-barwidth:round(0.1*sizein(1)),round(0.9*sizein(2))-barlength:round(0.9*sizein(2)),1)=sbcolor(1);
            in(round(0.1*sizein(1))-barwidth:round(0.1*sizein(1)),round(0.9*sizein(2))-barlength:round(0.9*sizein(2)),2)=sbcolor(2);
            in(round(0.1*sizein(1))-barwidth:round(0.1*sizein(1)),round(0.9*sizein(2))-barlength:round(0.9*sizein(2)),3)=sbcolor(3);
        case 'topleft'
            in(round(0.1*sizein(1))-barwidth:round(0.1*sizein(1)),round(0.1*sizein(2)):round(0.1*sizein(2))+barlength,1)=sbcolor(1);
            in(round(0.1*sizein(1))-barwidth:round(0.1*sizein(1)),round(0.1*sizein(2)):round(0.1*sizein(2))+barlength,2)=sbcolor(2);
            in(round(0.1*sizein(1))-barwidth:round(0.1*sizein(1)),round(0.1*sizein(2)):round(0.1*sizein(2))+barlength,3)=sbcolor(3);
    end
else
    switch Location
        case 'bottomright'
            in(round(0.9*sizein(1))-barwidth:round(0.9*sizein(1)),round(0.9*sizein(2))-barlength:round(0.9*sizein(2)),1)=sbcolor(1);
        case 'bottomleft'
            in(round(0.9*sizein(1))-barwidth:round(0.9*sizein(1)),round(0.1*sizein(2)):round(0.1*sizein(2))+barlength,1)=sbcolor(1);
        case 'topright'
            in(round(0.1*sizein(1))-barwidth:round(0.1*sizein(1)),round(0.9*sizein(2))-barlength:round(0.9*sizein(2)),1)=sbcolor(1);
        case 'topleft'
            in(round(0.1*sizein(1))-barwidth:round(0.1*sizein(1)),round(0.1*sizein(2)):round(0.1*sizein(2))+barlength,1)=sbcolor(1);
    end
end    

ImageOut=in;

if nargout == 0
    imshow(ImageOut)
end

end

