function plotBox(SMD,Data,Frame,BoxSize)
%plotBox plots the found boxes in the given Frame of the Data/SMD structure.
%
%plotBox takes the data and the structure containing the box info and the
%number of the frame that you are interested in and plots the found boxes.
%This routine is only really needed for debugging.
%
% INPUTS
%    SMD             SMD data structure with fields:
%       FrameNum     frame numbers of localizations
%       XBoxCorner   X coordinates of top right box corners (pixels)
%       YBoxCorner   Y coordinates of top right box corners (pixels)
%    Data            raw image data (XSize x YSize x NFrames)
%    Frame           user specified frame number of interest
%    BoxSize         linear box size for fitting (pixels)

% Created by
%    Mohamad Fazel

ParticleFrame = find(SMD.FrameNum==Frame);
XCorners = double(SMD.XBoxCorner(ParticleFrame));
YCorners = double(SMD.YBoxCorner(ParticleFrame));
PlotAxes = axes(figure());
imshow(Data(:,:,Frame), [], 'Parent', PlotAxes)
%PlotAxes.Visible = 'off';
hold(PlotAxes, 'on')
for nn = 1:length(ParticleFrame)
   plot(PlotAxes, XCorners(nn)*ones(1,BoxSize),YCorners(nn):YCorners(nn)+BoxSize-1,'y') 
   plot(PlotAxes, (XCorners(nn)+BoxSize-1)*ones(1,BoxSize),YCorners(nn):YCorners(nn)+BoxSize-1,'y')
   plot(PlotAxes, XCorners(nn):XCorners(nn)+BoxSize-1,YCorners(nn)*ones(1,BoxSize),'y') 
   plot(PlotAxes, XCorners(nn):XCorners(nn)+BoxSize-1,(YCorners(nn)+BoxSize-1)*ones(1,BoxSize),'y')
end
hold off

end
