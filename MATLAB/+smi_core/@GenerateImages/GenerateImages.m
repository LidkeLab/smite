classdef GenerateImages

%GenerateImages contains static methods for all visualization functions
%   of single molecule super-resolution data.
%
% blobColorOverlay creates color overlay of fitted emitters (green) onto data
%   (red). 
% colorImage generates RGB image with given colormap from a grey scale image
% dispIm produces a flexible GUI for multiple images.
% driftImage generates 2D histogram image from SR data.
% gaussianImage creates gaussian-blob image from SR data.
% histogramImage generates 2D histogram image from SR data.
% scaleBar creates a scale bar of desired length on the input image.
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
   [ImageOut, Image] = scalebar(Image, PixelSize, Length, Location)

   % Unit tests.
   [success] = blobColorOverlay_unitTest()
   [success] = colorImage_unitTest()
   [success] = driftImage_unitTest()
   [success] = gaussianImage_unitTest()
   [success] = histogramImage_unitTest()

   function [success] = unitTest()
      %unitTest tests various class functionality.
      [success(1)] = smi_core.GenerateImages.blobColorOverlay_unitTest();
      [success(2)] = smi_core.GenerateImages.colorImage_unitTest();
      [success(3)] = smi_core.GenerateImages.driftImage_unitTest();
      [success(4)] = smi_core.GenerateImages.gaussianImage_unitTest();
      [success(5)] = smi_core.GenerateImages.histogramImage_unitTest();
   end
end % methods(Static)

end % classdef GenerateImages
