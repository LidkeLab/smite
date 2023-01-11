# ***smite***: Single Molecule Imaging Toolbox Extraordinaire

This MATLAB-based toolbox provides analysis tools for fluorescence single molecule imaging with an emphasis on single molecule localization microscopy (SMLM) and single particle tracking (SPT).

## Overview

### Workflow concept
***smite*** is designed around the concept that a parameter structure, the Single Molecule Fitting (SMF) structure, uniquely and completely defines the data analysis.  The results are completely contained in a Single Molecule Data (SMD) structure.  ***smite*** is designed to make lowest-level tools just as easy to use as the higher-level application-specific classes.  All tools make use of the SMF and SMD structures.   

### Code organization
***smite*** is organized into a set of namespaces that group similar tools and concepts.  The namespace  `+smi`  contains the highest level tools that will be the most common entry point for processing SMLM and SPT data sets.  The file [SMITEclasses.txt](SMITEclasses.txt) provides a short 1-line desscription of each class in the distribution.

### Image and Detector Model
Image arrays follow MATLAB's column-major format.  An image coordinate of (1,1) means the center of the top-left pixel, whereaas (2,1) would indicate the center of the pixel that is one down from the top, but in the left-most column.   

## Installation
Clone (MacOS/Linux example; similar for Windows) into ~/Documents/MATLAB the ***smite*** GitHub distribution (https://github.com/LidkeLab/smite.git).  Add to Documents/MATLAB/startup.m the following:
```
   addpath '~/Documents/MATLAB/smite/MATLAB'
   setupSMITE;
```

### Dependencies
For full functionality ***smite*** requires:
- MATLAB version R2019b or later
- Nvidia GPU with CUDA compute capability [supported by your version of MATLAB](https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html)
- MATLAB Curve Fitting Toolbox [smi_cluster, smi_core.FRC, smi_stat.DiffusionEstimator]
- MATLAB Image Processing Toolbox
- MATLAB Optimization Toolbox [smi_cluster.PairCorrelation, smi_stat.DiffusionEstimator]
- MATLAB Parallel Computing Toolbox
- MATLAB Signal Processing Toolbox [smi_core.FRC]
- MATLAB Statistics and Machine Learning Toolbox

## Examples
### Working with SMF
SMF is implemented as a class to enable a gui and to provide useful helper methods.  However, the most common use will be as a structure with fixed fields.  

Create an SMF object:
```
  SMF = smi_core.SingleMoleculeFitting()
```
Get an SMF property:
```
  B = SMF.BoxFinding.BoxOverlap
```
Set an SMF property:
```
  SMF.BoxFinding.BoxOverlap = 0
```
Use the SMF GUI to interactively set values:
```
  SMF.gui()
```

### Finding coordinates from a stack of images containing blobs

Create a test dataset and make it noisy:
```
  B = smi_sim.GaussBlobs.genRandomBlobImage();
  B = poissrnd(B);
```
Create an `SMF` object with default values:
```
  SMF = smi_core.SingleMoleculeFitting()
```
Create a `LocalizeData` object with our `SMF`:
```
  LD = smi_core.LocalizeData(B, SMF)
```
Localize:
```
  [SMD] = LD.genLocalizations();
```

Localize again with `Verbose = 3` to show color overlay output:
```
  LD.Verbose = 3;
  [SMD] = LD.genLocalizations();
```

### High level SMLM analysis

Create an SMLM object.  When there are no input aruments, it will open the GUI:
```
  SMLMobj = smi.SMLM()  
```
Use the GUI to navigate to a test dataset such as this TIRF DNA-PAINT: 

Y:\Sandeep\20-11-2020-DNA_PAINT_Tubulin\Dock2-Cell1-2020-11-12-10-29-58.h5

Set SMF values from within GUI and run either a test dataset or analyze all datasets. 

### Additional Examples
Additional ***smite*** examples can be found in the examples subdirectory of MATLAB as well as the unitTests for some of the classes.  Some of the examples generate and analyze their own data, while others provide a template for how to run the example given supplied data.

## Contributions/Support
Issues or problems with the software should be reported via the Issues tab of the ***smite*** GitHub repository.  People seeking support or wishing to contribute to ***smite*** should contact Keith Lidke (klidke@unm.edu).
