# ***smite***: Single Molecule Imaging Toolbox Extraordinaire

This MATLAB-based toolbox provides analysis tools for fluorescence single molecule imaging with an emphasis on single molecule localization microscopy (SMLM) and single particle tracking (SPT).

## Overview

### Workflow concept
***smite*** is designed around the concept that a parameter structure, the Single Molecule Fitting (SMF) structure, uniquely and completely defines the data analysis.  The results are completely contained in a Single Molecule Data (SMD) structure.  ***smite*** is desingend to make lowest-level tools just as easy to use the higher-level application-specific classes.  All tools make use of the SMF and SMD structures.   

### Code organization
***smite*** is oragnized into a set of namespaces that group similar tools and concepts. The namespace  `+smi`   containes the highest level tools that will be the most common entry point for processing SMLM and SPT data sets.  

### Image and Detector Model
Image arrays follow MATLAB's column-major format.  An image coordinate of (1,1) means the center of the top-left pixel, whereaas (2,1) would indicate the center of the pixel that is one down from the top, but in the left-most column.   

### Dependencies
For full functionality ***smite*** requires:
- Nvidia GPU with CUDA compute capability [supported by your version of MATLAB](https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html)
- MATLAB Parallel Processing Toolbox
- MATLAB Statistics and Machine Learning Toolbox

## Examples
### High level SMLM analysis

### Finding coordinates from a stack of images containing blobs








