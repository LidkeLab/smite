# ***smite***: Single Molecule Imaging Toolbox Extraordinaire

This MATLAB-based toolbox provides analysis tools for fluorescence
single molecule imaging with an emphasis on single molecule
localization microscopy (SMLM) and single particle tracking (SPT).

## Overview

### Workflow concept
***smite*** is designed around the concept that a parameter structure,
the Single Molecule Fitting (SMF) structure, uniquely and completely
defines the data analysis.  The results are completely contained
in a Single Molecule Data (SMD) structure.  ***smite*** is designed
to make lowest-level tools just as easy to use as the higher-level
application-specific classes.  All tools make use of the SMF and
SMD structures.

### Code organization
***smite*** is organized into a set of namespaces that group similar
tools and concepts.  The namespace  `+smi`  contains the highest
level tools that will be the most common entry point for processing
SMLM and SPT data sets.  The file [SMITEclasses.md](doc/SMITEclasses.md)
provides a short 1-line description of each class in the distribution.

### Image and Detector Model
Image arrays follow MATLAB's column-major format.  An image coordinate
of (1,1) means the center of the top-left pixel, whereas (2,1) would
indicate the center of the pixel that is one down from the top, but
in the left-most column.

## Installation
Clone (MacOS/Linux example; similar for Windows) into ~/Documents/MATLAB
the ***smite*** GitHub distribution
(https://github.com/LidkeLab/smite.git).  Add to Documents/MATLAB/startup.m
the following:
```
   addpath '~/Documents/MATLAB/smite/MATLAB'
   setupSMITE;
```
Directories
MATLAB/source/C
cuda_Make
MATLAB/source/C
mex_Make

### Dependencies
For full functionality ***smite*** requires:
- Linux, MacOS or Windows
- MATLAB version R2021a or later
- Nvidia GPU with CUDA compute capability [supported by your version of MATLAB](https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html)
- MATLAB Curve Fitting Toolbox [smi_cluster, smi_core.FRC,
  smi_stat.DiffusionEstimator]
- MATLAB Image Processing Toolbox
- MATLAB Optimization Toolbox [smi_cluster.PairCorrelation,
  smi_stat.DiffusionEstimator]
- MATLAB Parallel Computing Toolbox
- MATLAB Signal Processing Toolbox [smi_core.FRC]
- MATLAB Statistics and Machine Learning Toolbox
- ffmpeg installed for Linux (https://ffmpeg.org)
  [smi_core.LocalizeData.genLocalizations for obj.Verbose >= 3]

## Examples
### Working with SMF
SMF is implemented as a class to enable a gui and to provide useful
helper methods.  However, the most common use will be as a structure
with fixed fields.

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
Use the GUI to navigate to a test dataset such as available from

- [Pallikkuth, S., Martin, C., Farzam, F., Edwards, J. S., Lakin,
  M. R., Lidke, D. S., & Lidke, K. A. (2018). Supporting data for
  Sequential Super-Resolution Imaging using DNA Strand Displacement
  [Data set]. University of New Mexico.
  https://doi.org/10.25827/CS2A-DH13](https://digitalrepository.unm.edu/physics_data/3/#attach_additional_files)
- [Wester, Michael J., Mazloom-Farsibaf, Hanieh, Farzam, Farzin,
  Fazel, Mohamadreza, Meddens, Marjolein B. M., & Lidke, Keith A.
  (2020), Comparing Lifeact and Phalloidin for super-resolution imaging
  of actin in fixed cells, Dryad, Dataset,
  https://doi.org/10.5061/dryad.xsj3tx9cn](https://datadryad.org/stash/dataset/doi:10.5061/dryad.xsj3tx9cn).

Set SMF values from within the GUI and run either a test dataset
or analyze all datasets.

### [Additional Examples](MATLAB/examples/README.md)
Additional ***smite*** examples can be found in the examples
subdirectory of MATLAB as well as the unitTests for some of the
classes (see [here](MATLAB/examples/README.md) for a summary).  Some
of the examples generate and analyze their own data, while others
provide a template for how to run the example given supplied data.

### [Additional Details](doc/DetailedOverview.md)
Additional details on ***smite*** can be found
[here](doc/DetailedOverview.md).

## Contributions/Support
Issues or problems and people seeking support with the software
should be reported via the Issues tab of the ***smite*** GitHub
repository.  Contributions to ***smite*** should be performed on
new branches which are then requested to merge with the main branch
via Pull Requests.

## Related Software
Please note the related software: MATLAB Instrument Control
(***MIC***), a collection of MATLAB classes for automated data
collection on complex, multi-component custom built microscopes.
This software can be obtained from the ***MIC*** GitHub distribution
(https://github.com/LidkeLab/matlab-instrument-control.git).
