# Typical UNM workflows using ***smite***

## 1-color: Publish -> BaGoL -> cluster -> compare conditions OR overlay images

- NOTE: moat of the scripts mentioned below are in ***smite*** examples.

- run **smi.SMLM**

  - Try *Test fit* from the GUI after setting various parameters, and iterate
    until the results are satisfactory.
  - *Export SMF* then and use it below.

- use **Example_Publish** or **Example_Publish_SMF**

  This script will call smi.Publish to generate localizations and miscellaneous
  results for a SMLM experiment on the sequential microscope.

- use **simpleROIcluster**
  - *Set important parameters*:
    ALWAYS run this first MATLAB section before running the later sections
    so that necessary parameters are set.
  - *Define the ROIs*:

    \_Results.mat -> \_ROIs.mat

    Each image [Results file] produces a separate ROIs file containing all the
    ROIs for that image; the idea is to do all the ROI selection early on so it
    need not be repeated---each image can have multiple ROIs; note that this
    deals with files in a single directory.  Many times lab members select the
    ROIs from the SR _Results.mat files.  To transfer these to BaGoL ROIs
    (necessary, as the ROI files are self-contained, meaning they contain the
    actual coordinates of the localizations in each ROI), it is necessary to
    transfer the ROI regions to BaGoL localizations.  See below.

- use **hierBaGoL_wrapper**

  Script to produce BaGoL results from SMITE SR *_Results.mat files.  The BaGoL
  results are placed in the subdirectory Results_BaGoL, assigned a name below,
  under the directory containing the *_Results.mat files.  hierBaGoL_analysis
  is called separately on each file in a parfor loop, assigning a worker for
  each dataset, so conducive to be run on a multi-core machine.

  - If the variable *ROIs* is true and only one filename is provided (see
    below), and assuming SMD files are of the form Cell_nn_Label_0n_Results.mat
    and ROI files are of the form Cell_nn_Label_01_Results_ROIs.mat (and are
    located in the subdirectory 'Analysis' of DataDir) as done when choosing
    ROIs using the scripts simpleROIcluster or simplePairCorr, this current
    script will automatically use the ROI information in the _ROIs.mat file to
    produce a series of individual analyses for each ROI which
    smi.BaGoL.hierBaGoL_run will parallelize via a parfor loop.  The set of
    analyses will take on names of the form Cell_nn_Label_0n_Results_ROI_mm.
    PairAnalysis.overlayBaGoLROIs is a useful function for plotting the ROIs
    produced by this process, both 1-color and 2-color.
  - If *ROIs* is true above, the example script **generateBaGoLScripts** may
    come in handy:
    Generate a series of self-contained scripts to run BaGoL over a number of
    cells, one script per cell.  This allows the BaGoL analysis of each cell to
    be run independently, which is useful when running on different servers
    sharing a common data space.

- use **simpleROIcluster**:
  - *Define the BaGoL ROIs from the previous ROIs and BaGoL results*

    \_ROIs.mat + MAPN\_*.mat -> \_BaGoL_ROIs.mat

    The two sets of files should be in corresponding order; the BaGoL
    coordinates replace the coordinates in the original ROI files.
  - *Possibly, combine individually processed BaGoL ROIs into a single
    _ROIs.mat file using the _ROIs.mat file that was used to define the ROIs
    originally from the SR data*

    The above combiner is used when BaGoL is run with *ROIs* set to true, in
    preparation for the clustering step below.

- use either **singleConditionDriver** (standalone script) or
  **simpleROIcluster**: *Statistics for a single condition*.
  **singleConditionDriver** tends to run a great deal faster as it runs in
  batch mode; it is adapted from the latter code.

  **singleConditionDriver** performs cluster analysis for comparison of
  experimental conditions.  This is a batch version of simpleROIcluster for
  performing parameter studies.  See pathname/files definitions in
  *Statistics for a single condition*.  pathname designates the path where the
  ROIs or BaGoL_ROIs are stashed for a single condition, while files reference
  the *_ROIs.mat or *_BaGoL_ROIs.mat files containing ROI definitions.
  *_BaGol_ROIs.mat files are converted from from *_ROIs.mat and MAPN files by
  the section (*Define the BaGoL ROIs from the previous ROIs and BaGoL
  results*) in **simpleROIcluster**.

- use
  - **simROIcluster**: *Combined statistics for one or more conditions* (simple
    combined plots) or *Combined statistics for multiple conditions and
    experiments* (more complex combined plots) to compare multiple single
    conditions whose results (.mat files) are collected together into a single
    directory.
  - **Example_plotROIDriver** to make overlays of Circle, Gaussian and
    Guassian-like (constant localization standard error) [option
    'GaussSEConst'] images with ROI and cluster boundaries.  Inputs can be
    single condition results, plus one of:
    (1) a SR Results file [option 'SR']; (2) a single BaGoL MAPN file
    containing multiple ROIs [option 'MAPN']; a series of individual BaGoL MAPN
    Results_ROI files containing one ROI per file [option 'MAPNResultsROI'].
    This last can handle missing ROIs as long as the naming convention is
    observed:

    filesB and filesC.files are assumed to have the string 'Cell_nn' embedded
    in their names, where nn is a 2-digit cell number like 01, 12, etc.  filesB
    is also assumed to have a 2-digit ROI number: 'ROI_nn' embedded in their
    names.

## 2-color: Publish -> BaGoL -> cluster -> compare conditions OR overlay images

- The analyses for 2 colors has much in common with 1-color analysis, except
  Label_01 and Label_02 files are treated separately once ROIs are chosen.
