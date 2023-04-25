### +smi_sim/@GaussBlobs

GaussBlobs A collection of methods for generating 2D Gaussian Blob Images

This class consists of two static methods that can be used to generate
image stacks containing Gaussian blobs: 
- **[gaussBlobROIStack](gaussBlobROIStack.m)**
  generates a stack of images, each containing a single blob.
- **[gaussBlobImage](gaussBlobImage.m)**
  generates a stack of images with multiple blobs and
  internally uses gaussBlobROIStack.  This function is used in simulation
  and display of 2D single molecule data. 
- **[unitTest](GaussBlobs.m)**

REQUIRES:
- MATLAB 2014a or later versions
- Parallel Procesing Toolbox
- NVidia GPU
- smi_cuda_gaussBlobROIStack.ptx
- smi_cuda_gaussBlobROIStack.cu
