classdef GenerateImages

%GenerateImages contains static methods for all visualization functions
%   of single molecule super-resolution data
    
% Lidke Lab 2017
    
properties
end
    
methods(Static)
   [OverlayImage] = blobColorOverlay(Sequence, SMD)
   [RGBimage] = colorImage(Image, ColorMap, MinMax)
   dispIm()
   [DriftIm, DriftImRGB] = driftImage(SMR, RawDataSize, SRImageZoom)
   [GaussIm] = gaussianImage(SMR, SRImageZoom)
   [HistIm, RgbHistIm] = histogramImage(SMR, RawDataSize, SRImageZoom, ColorMap)
   plotBox(SMD, Data, Frame, BoxSize)
   [ImageOut, Image] = scalebar(Image, PixelSize, Length, Location)

   % Unit tests.
   [success] = blobColorOverlay_unitTest()
   [success] = colorImage_unitTest()
   [success] = driftImage_unitTest()
   [success] = histogramImage_unitTest()

   function [success] = unitTest()
      [success(1)] = blobColorOverlay_unitTest();
      [success(2)] = colorImage_unitTest();
      [success(3)] = driftImage_unitTest();
      [success(4)] = histogramImage_unitTest();
   end
end % methods(Static)

end % classdef GenerateImages
