function plotBox(SMD,Data,Frame,BoxSize)
%plotBox plots the found boxes in the given Frame of the Data/SMD structure.
%
%plotBox takes the data and the structure containing the box info and the
%number of the frame that you are interested in and plot the found boxes.

ParticleFrame = find(SMD.FrameNum==Frame);
XCorners = double(SMD.XBoxCorner(ParticleFrame));
YCorners = double(SMD.YBoxCorner(ParticleFrame));
dipshow(Data(:,:,Frame))
hold
for nn = 1:length(ParticleFrame)
   plot(XCorners(nn)*ones(1,BoxSize),YCorners(nn):YCorners(nn)+BoxSize-1,'r') 
   plot((XCorners(nn)+BoxSize-1)*ones(1,BoxSize),YCorners(nn):YCorners(nn)+BoxSize-1,'r')
   plot(XCorners(nn):XCorners(nn)+BoxSize-1,YCorners(nn)*ones(1,BoxSize),'r') 
   plot(XCorners(nn):XCorners(nn)+BoxSize-1,(YCorners(nn)+BoxSize-1)*ones(1,BoxSize),'r')
end

end
