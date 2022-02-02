classdef GenerateImages
    
    %GenerateImages contains static methods for general visualization functions
    %   of single molecule super-resolution data.
    %
    % blobColorOverlay creates color overlay of fitted emitters (green) onto data
    %    (red).
    % colorImage generates RGB image with given colormap from a grey scale image
    % dispIm produces a flexible GUI for multiple images.
    % driftImage generates 2D histogram image from SR data.
    % gaussianImage creates gaussian-blob image from SR data.
    % histogramImage generates 2D histogram image from SR data.
    % plotHistogram creates a histogram of a field from an SMD structure.
    % scaleBar creates a scale bar of desired length on the input image.
    %
    % REQUIRES:
    %    Statistics Toolbox
    %    Parallel Procesing Toolbox
    %    NVidia GPU
    %
    % SEE ALSO (other visualization routines):
    %    smi_core.DriftCorrection.plotDriftCorrection
    %    smi_core.FindROI.plotBox
    %    smi_core.Threshold.rejectedLocalizations
    
    % Lidke Lab 2017, 2020
    
    properties
    end
    
    methods(Static)
        [OverlayImage] = blobColorOverlay(Sequence, SMD)
        [RGBimage] = colorImage(Image, ColorMap, MinMax)
        [RGBimage] = rgbImage(R,G,B)
        showim(Image)
        dispIm()
        [DriftIm, DriftImRGB] = driftImage(SMR, SRImageZoom)
        [GaussIm] = gaussianImage(SMR, SRImageZoom, ScalebarLength)
        [CircleImage, CircleImageRGB, SRImageZoom] = ...
            circleImage(SMR, ColorMap, ...
            SRImageZoom, MinPixelsPerCircle, SEScaleFactor);
        [CircleDriftImage, SRImageZoom] = ...
            circleDriftImage(SMR, SRImageZoom, ...
            MinPixelsPerCircle, SEScaleFactor);
        [HistIm, RgbHistIm] = histogramImage(SMR, SRImageZoom, ColorMap)
        [FigHandle] = plotHistogram(Vector_in, Hist_Name)
        [ImageOut, Image] = scalebar(Image, PixelSize, Length, Location)
        [OverlayImage, ColorOrderTag] = overlayNImages(ImageStack);
        
        % Unit tests.
        [success] = blobColorOverlay_unitTest()
        [success] = colorImage_unitTest()
        [success] = driftImage_unitTest()
        [success] = gaussianImage_unitTest()
        [success] = histogramImage_unitTest()
        [Success] =  circleImage_unitTest()
        
        function [success] = unitTest()
            %unitTest tests various class functionality.
            [success(1)] = smi_vis.GenerateImages.blobColorOverlay_unitTest();
            [success(2)] = smi_vis.GenerateImages.colorImage_unitTest();
            [success(3)] = smi_vis.GenerateImages.driftImage_unitTest();
            [success(4)] = smi_vis.GenerateImages.gaussianImage_unitTest();
            [success(5)] = smi_vis.GenerateImages.histogramImage_unitTest();
            [success(6)] = smi_vis.GenerateImages.circleImage_unitTest();
        end
    end % methods(Static)
    
end % classdef GenerateImages
