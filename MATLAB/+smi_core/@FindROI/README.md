### +smi_core/@FindROI

FindROI Finds and collates subregions from a stack of 2D images

Subregions are found by filtering with a difference of Gaussians filter
followed by finding the local maximum.  Local maxima are used as the
center of ROIs if the estimated single molecule intesity is greater than
MinPhotons.

Gaussian filtering is done using the recursive method of Young.
Local maximum finding is done using a novel separable filter (Lidke).

Filtering and maximum finding are implemented on GPU using CUDA
and the Parallel Processsing Toolbox.

REQUIRES:
- MATLAB 2014a or later versions
- Parallel Procesing Toolbox
- NVidia GPU
- smi_cuda_FindROI.ptx
- smi_cuda_FindROI.cu

CITATION:
  Ian T. Young, Lucas J. van Vliet,
  Recursive implementation of the Gaussian filter,
  Signal Processing, Volume 44, Issue 2, 1995

---

```
properties:
   Data            %Stack of 2D images
   BoxSize=7       %Linear box size for fitting (Pixels)(Default=7)
   BoxOverlap=2    %Overlap of boxes (Integer Pixels)(Default=2)
   MinPhotons=200  %Minimum number of photons from emitter (Default=200)
   PSFSigma=1.3    %Sigma of 2D Gaussian PSF (Pixels) (Default=1.3)
   ROIs            %Stack of subregions
   LocalMaxIm      %Binary Image showing local maxima above the threshold
   PlotBoxFrame=1  %If Verbose >= 3, plot boxes for this frame (Default=1)
   Verbose=1       %Verbosity level
   IsSCMOS         %Is a SCMOS Camera and will use Variance image 
   Varim           %Variance Image (photons^2)
```

---

methods:
- **[extractROIs](extractROIs.m)**:
  reloads raw data and extracts ROIs corresponding to SMD locs
- **[plotBox](plotBox.m)**:
  plots the found boxes in the given Frame of the Data/SMD structure
- **[plotBoxStack](plotBoxStack.m)**:
  plots the found boxes in Data/SMD structure
