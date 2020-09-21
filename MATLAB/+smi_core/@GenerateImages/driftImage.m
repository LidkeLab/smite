function [DriftIm, DriftImRGB] = driftImage(SMR, SRImageZoom)
%driftImage generates 2D histogram image from SR data
%   INPUT
%      SMR - structure containing at least fields 'X', 'Y', 'XSize', 'YSize','FileNum', and 
%            'FrameNum' having the x and y coordinates, rawdata x and y sizes, file number and 
%            frame number of localizations to be plotted
%      SRImageZoom - zoom factor for drift image 
%   OUTPUT
%      If no output args are given drift image will be displayed
%      DriftIm - gray value drift image
%      DriftImRGB - RBG drift image
%   REQUIRES
%      Matlab 2014b or higher
%      Dipimage toolbox (http://www.diplib.org/)

% Marjolein Meddens, Lidke Lab 2017

% check input
if nargin <1
    error('smi_core:GenerateImages:driftImage: Not enough input arguments. Please input atleast SMR');
end

if nargin<2
    disp('No input on SRImageZoom, default value of 10 assigned');
    SRImageZoom=10;
end


RawDataSize=[SMR.YSize SMR.XSize];

% calculate size of drift image
xsize=(RawDataSize(2)*SRImageZoom);
ysize=(RawDataSize(1)*SRImageZoom);
% transform coordinates to match image size
x=single(SMR.X*SRImageZoom);
y=single(SMR.Y*SRImageZoom);
% create time stamp of each localization
NperFrame=single(max(SMR.FrameNum)+1);
t=single(SMR.FrameNum)+NperFrame*single(SMR.DatasetNum-1);
t=single(t);
% create drift image
DriftIm = c_HistImTime(ysize,xsize,y,x,t);
% create output
cm=jet(256);
cm(1,:)=[0 0 0];
if nargout==0
    dipshow(DriftIm);
    colormap(cm);
elseif nargout == 2
    [DriftImRGB]=smi_core.GenerateImages.colorImage(DriftIm,cm);
end

end
