classdef GenerateImages

%GenerateImages contains static methods for all visualization functions
%   of single molecule super-resolution data.
%
% REQUIRES:
%   Statistics Toolbox
%   Parallel Procesing Toolbox
%   Dipimage toolbox
%   NVidia GPU

% Lidke Lab 2017, 2020
    
properties
end
    
methods(Static)
   [OverlayImage] = blobColorOverlay(Sequence, SMD)
   [RGBimage] = colorImage(Image, ColorMap, MinMax)
   dispIm()
   [DriftIm, DriftImRGB] = driftImage(SMR, SRImageZoom)
   [GaussIm] = gaussianImage(SMR, SRImageZoom)
   [HistIm, RgbHistIm] = histogramImage(SMR, SRImageZoom, ColorMap)
   plotBox(SMD, Data, Frame, BoxSize)
   [ImageOut, Image] = scalebar(Image, PixelSize, Length, Location)

   % Unit tests.
   [success] = blobColorOverlay_unitTest()
   [success] = colorImage_unitTest()
   [success] = driftImage_unitTest()
   [success] = gaussianImage_unitTest()
   [success] = histogramImage_unitTest()

   function [success] = unitTest()
      [success(1)] = smi_core.GenerateImages.blobColorOverlay_unitTest();
      [success(2)] = smi_core.GenerateImages.colorImage_unitTest();
      [success(3)] = smi_core.GenerateImages.driftImage_unitTest();
      [success(4)] = smi_core.GenerateImages.gaussianImage_unitTest();
      [success(5)] = smi_core.GenerateImages.histogramImage_unitTest();
   end
end % methods(Static)

end % classdef GenerateImages
