### +smi_vis/@GenerateImages
    
GenerateImages contains static methods for general visualization functions
of single molecule super-resolution data.

methods:
- **blobColorOverlay** creates color overlay of fitted emitters (green) onto
  data (red).
- **colorImage** generates RGB image with given colormap from a grey scale
  image
- **dispIm** produces a flexible GUI for multiple images.
- **driftImage** generates 2D histogram image from SR data.
- **gaussianImage** creates gaussian-blob image from SR data.
- **histogramImage** generates 2D histogram image from SR data.
- **plotHistogram** creates a histogram of a field from an SMD structure.
- **scaleBar** creates a scale bar of desired length on the input image.

REQUIRES:
- Statistics Toolbox
- Parallel Procesing Toolbox
- NVidia GPU

SEE ALSO (other visualization routines):
- smi_core.DriftCorrection.plotDriftCorrection
- smi_core.FindROI.plotBox
- smi_core.Threshold.rejectedLocalizations
