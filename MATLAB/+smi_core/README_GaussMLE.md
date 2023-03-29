GaussMLE: Maximum Likelihood Estimate of 2D Gaussian blob parameters

GaussMLE fits a 2D Gaussian blob model to the each image in a stack of 2D
images.  Several fit types are available. The MLE is maximized with an
iterative Newton-Raphson method and implemented on GPU using NVIDIA
CUDA as described by Smith et al.

The input is assumed to be gain and offset corrected so that
pixel values indicate effective photon counts. Noise model is 'Poisson',
which is suitable for EMCCD cameras used at high gain, or 'SCMOS' which
includes a pixel-wise readnoise as described by Huang et al.

REQUIRES:
- MATLAB 2014a or later versions
- Parallel Procesing Toolbox
- NVidia GPU
- smi_cuda_FindROI.ptx
- smi_cuda_FindROI.cu

CITATION:
Smith, C., Joseph, N., Rieger, B. et al.
Fast, single-molecule localization that achieves theoretically minimum
uncertainty. Nat Methods 7, 373–375 (2010).
https://doi.org/10.1038/nmeth.1449

Huang, F., Hartwich, T., Rivera-Molina, F. et al. Video-rate nanoscopy
using sCMOS camera–specific single-molecule localization algorithms.
Nat Methods 10, 653–658 (2013).
https://doi.org/10.1038/nmeth.2488
