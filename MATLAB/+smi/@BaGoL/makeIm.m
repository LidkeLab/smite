function [SRIm]=makeIm(SMD,SZ,PixSize,XStart,YStart,BoxSize)
%makeIm: Produces a Gaussian blob image from the list of input coordinates.
% [SRIm]=BaGoL.makeIm(SMD,SZ,PixSize,XStart,YStart)
%
% Each input coordinate and standard error is used to add a normalized, 
% 2D Gaussian blob to an intially empty image. The size and image region
% is given by the input parameters. 
%
% For computational speed, this operation is done by calculating all blobs 
% within small boxes of the same size in one operation and then copying the 
% results into the output image. 
% 
% INPUTS:
%   SMD:     SMD structure with fields:
%       X:        Vector of X SR-localization positions (nm)(Nx1)
%       Y:        Vector of Y SR-localization positions (nm)(Nx1)
%       X_SE:     Vector of X SR-localization standard errors (nm)(Nx1)
%       Y_SE:     Vector of Y SR-localization standard errors (nm)(Nx1)
%   SZ:      Size of the generated image (square) (nm) 
%   PixSize: Pixel size of the output image (nm)
%   XStart:  Start of generated image in the X-axis (nm), (Default = 0)
%   XStart:  Start of generated image in the Y-axis (nm), (Default = 0)
%   BoxSize: Box size for single Gaussian blobs (nm) (OPTIONAL) (Default=10*median(X_SE)) 
%
% OUTPUTS: 
%   SRIm:    Output image.
%

% Created by:
%    Mohamadreza Fazel (LidkeLab, 2019)

%% SR-Image

if nargin < 5
    MinX = 0; 
    MinY = 0;
else
    MinX = -XStart;
    MinY = -YStart;
end

if nargin < 6
    BoxSize = floor(15*median([SMD.X_SE;SMD.Y_SE])/PixSize); %Size of box in pixels
else
    BoxSize = round(BoxSize/PixSize); 
end
if floor(BoxSize/2)~=BoxSize/2
    BoxSize = BoxSize + 1; 
end

XBox = single(floor((SMD.X+MinX)/PixSize)-floor(BoxSize/2));
YBox = single(floor((SMD.Y+MinY)/PixSize)-floor(BoxSize/2));
X = single((SMD.X+MinX) - PixSize*floor((SMD.X+MinX)/PixSize) + PixSize*floor(BoxSize/2)-1);
Y = single((SMD.Y+MinY) - PixSize*floor((SMD.Y+MinY)/PixSize) + PixSize*floor(BoxSize/2)-1);

[Xg,Yg,~]=meshgrid((PixSize/2:PixSize:BoxSize*PixSize-PixSize/2),...
    (PixSize/2:PixSize:BoxSize*PixSize-PixSize/2),(1:length(SMD.X)));
Xg = single(Xg);
Yg = single(Yg);
MuX = ones(size(Xg),'single');
SigX = ones(size(Xg),'single');
MuY = MuX;
SigY = SigX;
for ii = 1:length(SMD.X)
     SigX(:,:,ii) = SMD.X_SE(ii);
     SigY(:,:,ii) = SMD.Y_SE(ii);
     MuX(:,:,ii) = X(ii);
     MuY(:,:,ii) = Y(ii);
end
Im = normpdf(Xg,MuX,SigX).*normpdf(Yg,MuY,SigY);
ExtendSZ=100;  %imae size extension
SRImT = zeros(ceil(SZ/PixSize)+ExtendSZ,'single');
for ii = 1:length(SMD.X) 
    try
        SRImT(YBox(ii)+(ExtendSZ/2+1):YBox(ii)+BoxSize+(ExtendSZ/2),...
            XBox(ii)+(ExtendSZ/2+1):XBox(ii)+BoxSize+(ExtendSZ/2)) = ...
            SRImT(YBox(ii)+(ExtendSZ/2+1):YBox(ii)+BoxSize+(ExtendSZ/2),...
            XBox(ii)+(ExtendSZ/2+1):XBox(ii)+BoxSize+(ExtendSZ/2)) + Im(:,:,ii);
    catch
    end
end
SRIm = SRImT((ExtendSZ/2+1):end-(ExtendSZ/2),(ExtendSZ/2+1):end-(ExtendSZ/2));

end
