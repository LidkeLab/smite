### +smi_vis/@GenerateImages
    
GenerateImages contains static methods for general visualization functions
of single molecule super-resolution data.

REQUIRES:
- Statistics Toolbox
- Parallel Procesing Toolbox
- NVidia GPU

SEE ALSO (other visualization routines):
- smi_core.DriftCorrection.plotDriftCorrection
- smi_core.FindROI.plotBox
- smi_core.Threshold.rejectedLocalizations

---

methods:
- **[blobColorOverlay](blobColorOverlay.m)**:
  creates color overlay of fitted emitters (green) onto data (red).
- **[colorImage](colorImage.m)**:
  generates RGB image with given colormap from a grey scale image
- **[dispIm](dispIm.m)**:
  produces a flexible GUI for multiple images.
- **[driftImage](driftImage.m)**:
  generates 2D histogram image from SR data.
- **[gaussianImage](gaussianImage.m)**:
  creates gaussian-blob image from SR data.
- **[grayscaleImage](grayscaleImage.m)**:
  creates grayscale gaussian-blob image from SR data.
- **[histogramImage](histogramImage.m)**:
  generates 2D histogram image from SR data.
- **[plotHistogram](plotHistogram.m)**:
  creates a histogram of a field from an SMD structure.
- **[scaleBar](scaleBar.m)**:
  creates a scale bar of desired length on the input image.
