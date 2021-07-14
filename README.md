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
### Working with SMF
SMF is implemented as a class to enable a gui and to provide useful helper methods.  However, the most common use will be as a structure with fixed fields.  

Create an SMF object:
```
  SMF=smi_core.SingleMoleculeFitting()
```
Get an SMF property
```
  B=SMF.BoxFinding.BoxOverlap
```
Set an SMF property
```
SMF.BoxFinding.BoxOverlap=0
```
Use the SMF GUI to interactive set values:
```
SMF.gui()
```

### Finding coordinates from a stack of images containing blobs

Create a test dataset and make it noisy
```
B=smi_sim.GaussBlobs.genRandomBlobImage();
B=poissrnd(B);
```
Create an `SMF` object with default values:
```
SMF=smi_core.SingleMoleculeFitting()
```
Create a `LocalizeData` object with our `SMF`
```
LD = smi_core.LocalizeData(B, SMF)
```
Localize
```
[SMD] = LD.genLocalizations();
```

Localize again with `Verbose=2` to show color overlay output
```
LD.Verbose=2
[SMD] = LD.genLocalizations();
```

### High level SMLM analysis

Create a SMLM object.  When there is no input aruments, it will open the GUI. 
```
SMLMobj=smi.SMLM()  
```
Use GUI to navigate to a test dataset such as this TIRF DNA-PAINT: 

Y:\Sandeep\20-11-2020-DNA_PAINT_Tubulin\Dock2-Cell1-2020-11-12-10-29-58.h5

Set SMF values from within GUI and run either a test dataset or analyze all datasets. 







