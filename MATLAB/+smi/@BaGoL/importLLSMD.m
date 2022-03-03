function importLLSMD(obj,SMD,PixelSize)
%importLLSMD Imports a Lidke Lab SMD and converts to BaGoL format
%
% The Lidke Lab standard has SMD in units of pixels.  This converts 
% SMD fields from Pixels to nm. The BaGoL SMD property is updated with the
% result. 
%
% INPUTS:
%   SMD:    SMD Structure with fields: 
%       X:      (Pixels)
%       Y:      (Pixels)
%       X_SE:   (Pixels)
%       Y_SE:   (Pixels)
%       FrameNum: (Frames)
%   PixelSize:  Camera back-projected pixel size (nm/Pixel)
% 

obj.SMD=[];
obj.SMD.X=SMD.X*PixelSize;
obj.SMD.Y=SMD.Y*PixelSize;
obj.SMD.X_SE=SMD.X_SE*PixelSize;
obj.SMD.Y_SE=SMD.Y_SE*PixelSize;
obj.SMD.FrameNum=SMD.FrameNum; 


end