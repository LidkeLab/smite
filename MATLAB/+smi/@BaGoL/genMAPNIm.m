function [SRIm]=genMAPNIm(obj,ImFlag)
%makeIm: Produces a Gaussian blob image from either SMD or MAPN.
% [SRIm]=obj.genMAPNIm(ImFlag)
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
%   ImFlag: Type of image to make. (Default = MAPN Image)
%      1: Retrieves the MAPN coordinates from obj and makes the MAPN image
%      2: Retrieves the SMD coordinates from obj and makes the SR image
%
% OUTPUTS: 
%   SRIm:    Output image.
%

% Created by:
%    Mohamadreza Fazel (LidkeLab, 2020)

%% SR-Image

MinX = -obj.XStart;
MinY = -obj.YStart;
if nargin==1||ImFlag==1
    SMR = obj.MAPN;
elseif ImFlag == 2
    SMR = obj.SMD;
else
    error('ImFlag must be either 1 or 2.')
end
BoxSize = floor(15*median([SMR.X_SE;SMR.Y_SE])/obj.PixelSize); %Size of box in pixels
if floor(BoxSize/2)~=BoxSize/2
    BoxSize = BoxSize + 1; 
end

XBox = single(floor((SMR.X+MinX)/obj.PixelSize)-floor(BoxSize/2));
YBox = single(floor((SMR.Y+MinY)/obj.PixelSize)-floor(BoxSize/2));
X = single((SMR.X+MinX) - obj.PixelSize*floor((SMR.X+MinX)/obj.PixelSize) + obj.PixelSize*floor(BoxSize/2)-1);
Y = single((SMR.Y+MinY) - obj.PixelSize*floor((SMR.Y+MinY)/obj.PixelSize) + obj.PixelSize*floor(BoxSize/2)-1);

[Xg,Yg,~]=meshgrid((obj.PixelSize/2:obj.PixelSize:BoxSize*obj.PixelSize-obj.PixelSize/2),...
    (obj.PixelSize/2:obj.PixelSize:BoxSize*obj.PixelSize-obj.PixelSize/2),(1:length(SMR.X)));
Xg = single(Xg);
Yg = single(Yg);
MuX = ones(size(Xg),'single');
SigX = ones(size(Xg),'single');
MuY = MuX;
SigY = SigX;
for ii = 1:length(SMR.X)
     SigX(:,:,ii) = SMR.X_SE(ii);
     SigY(:,:,ii) = SMR.Y_SE(ii);
     MuX(:,:,ii) = X(ii);
     MuY(:,:,ii) = Y(ii);
end
Im = normpdf(Xg,MuX,SigX).*normpdf(Yg,MuY,SigY);
ExtendSZ=100;  %imae size extension
SRImT = zeros(ceil(obj.PImageSize/obj.PixelSize)+ExtendSZ,'single');
for ii = 1:length(SMR.X) 
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
