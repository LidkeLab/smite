### hierBaGoL_wrapper Summary

[examples/hierBaGoL_wrapper](../../MATLAB/examples/hierBaGoL_wrapper.m)
is a script for processing multiple hierarchical BaGoL datasets or
splitting a single dataset into multiple previously defined ROIs.
It calls
[smi.BaGoL.hierBaGoL_run](../../MATLAB/+smi/@BaGoL/hierBaGoL_run.m)
which runs one or more BaGoL analyses, calling
[smi.BaGoL.hierBaGoL_analysis](../../MATLAB/+smi/@BaGoL/hierBaGoL_analysis.m)
on each individual dataset to be processed, so acting as a dispatch
intermediary.  A single dataset is run directly, while a set of datasets
(which could be a single dataset split into multiple ROIs) are run in
parallel using a parfor loop.  hierBaGoL_analysis.m is adapted from a
version of BaGoL_EGFR_dSTORM.  hierBaGol_wrapper was designed to be
the only routine that the user needs to interact with, but new features
should regard this flow.

A basic description is as follows:

Script to produce BaGoL results from SMITE \*\_Results.mat files.  The BaGoL
results are placed in the subdirectory Results_BaGoL, assigned a name below,
under the directory containing the \*\_Results.mat files.  hierBaGoL_analysis
is called separately on each file in a parfor loop, assigning a worker for
each dataset, so conducive to be run on a multi-core machine.
The wrapper sets parameters and lists files (full paths) to be analyzed and
optional ROIs to apply in the next to last section.

For \_Results.mat files with large numbers of localizations (> 300,000 or so),
hierBaGoL may crash (or partially crash), so should not be part of a parfor
loop as a crash will cause ALL of the non-finished parallel jobs to restart.
Such \_Results.mat files should be analyzed in separate MATLABs.

NOTE: MAPN\_\*.mat files are always produced containing simply the MAPN
coordinates.  See
[smi.BaGoL.hierBaGoL_analysis](../../MATLAB/+smi/@BaGoL/hierBaGoL_analysis.m)
for more details on files produced.

If the variable ROIs is true and only one filename is provided (see below),
and assuming SMD files are of the form Cell_nn_Label_0n_Results.mat and ROI
files are of the form Cell_nn_Label_01_Results_ROIs.mat (and are located
in the subdirectory 'Analysis' of DataDir) as done when choosing ROIs using
the scripts
[examples/simpleROIcluster](../../MATLAB/examples/simpleROIcluster.m) or
[examples/simplePairCorr](../../MATLAB/examples/simplePairCorr.m),
this current script will
automatically use the ROI information in the \_ROIs.mat file to produce a
series of individual analyses for each ROI which smi.BaGoL.hierBaGoL_run will
parallelize via a parfor loop.  The set of analyses will take on names of the
form Cell_nn_Label_0n_Results_ROI_mm.
[smi_cluster.PairAnalysis.overlayBaGoLROIs](../../MATLAB/+smi_cluster/@PairAnalysis/overlayBaGoLROIs.m)
is a useful function for plotting the ROIs produced by this process,
both 1-color and 2-color.

---

An example input file structure with `ROIs = true;` and 2 ROIs selected
previously:
```
DATA/
  ResultsStructs/
    Cell_02_Label_01_Results.mat'
    Analysis/
      Cell_02_Label_01_Results_ROIs.mat'
```
results in these files:
```
      BaGoL_Results_Cell_02_Label_01_Results_ROI_01_ResultsStruct.mat
      BaGoL_Results_Cell_02_Label_01_Results_ROI_02_ResultsStruct.mat
      MAPN_Cell_02_Label_01_Results_ROI_01.mat
      MAPN_Cell_02_Label_01_Results_ROI_02.mat
      Cell_02_Label_01_Results_ROI_01/
        BaGoL_X-SE.png
        BaGoL_Y-SE.png
        FULL.png
        LocsScatter-MAPN.fig
        MAPN.mat
        MAPN-Im.png
        MAPN_NmeanHist.png
        NND.png
        NNDScaledData.png
        NNDScaledRandom.png
        NND.txt
        Overlay_cPost_rMap.png
        Overlay_gSR_bPost_rMap.png
        Overlay_gSR_mMap.png
        Overlay_gSR_mPost.png
        Overlay_SR_Map_circle.png
        Post-Im.png
        prior.txt
        ROI.png
        SR-Im.png
        Xi.png
        XiChain.png
      Cell_02_Label_01_Results_ROI_02/
        ...
```
- The scale bar throughout is 500 nm.
- FULL.png and ROI.png are the SR localizations in the full cell and the ROI,
  respectively.
- Overlay_cPost_rMap: posterior localizations are colored cyan and MAPN red.
- Overlay_gSR_bPost_rMap: SR, posterior and MAPN localizations are colored
  green, blue and red, respectively.
- Overlay_gSR_mMap: SR localizations are colored green and MAPN magenta.
- Overlay_SR_Map_circle use the same color scheme as above (SR: green, MAPN:
  magenta) and plots circles instead of Gaussian blobs where the radii of the
  circles are proportional to the standard error (SE) of the localizations.
- prior.txt is the estimated prior for a second run of the data, which is not
  needed in hierarchical BaGoL.

---

Suggested pre-filtering actions (frame connection and NN not used for dSTORM
data) [NOTE that the middle four actions should be performed during the SMLM
analysis, e.g., Publish]:
```
SR data -> remove localizations with negative coordinates
           [smi_helpers.Filters.filterNonNeg called by hierBaGoL_analysis]
        -> intensity filter [SMF.Thresholding.InMeanMultiplier]
        -> inflate standard errors [SMF.Data.SEAdjust]
        -> frame connection, removing connections which involve only a
           specified number of frames [SMF.FrameConnection.MinNFrameConns]
        -> Nearest Neighbor filter (N_NN) --- Do not use on dSTORM data!
           [SMF.Thresholding.NNMedianMultiplier,
            SMF.Thresholding.MinNumNeighbors]
        -> BaGoL (via parfor calling hierBaGoL_analysis on each dataset)
```
The pre-filtering (except for removing negative coordinates) is all now in
SMLM, although SE_Adjust can be set here as well.

---

hierBaGoL_wrapper collects together important BaGoL parameters in the
structure `BaGoLParams`.  These (along with some important local
parameters) are:
```
% Output directory name.
Results_BaGoL = 'Results_BaGoLHier';

% Generic parameters.
BaGoLParams.ImageSize = 256;        % (pixel)
%BaGoLParams.PixelSize = 108.018;    % camera back projected size (nm) [TIRF]
BaGoLParams.PixelSize = 97.8;       % (nm) [sequential]
BaGoLParams.OutputPixelSize = 4;    %2; % pixel size for posterior images (nm)
BaGoLParams.N_Burnin = 32000;       % Length of Burn-in chain
BaGoLParams.N_Trials = 8000;        % Length of post-burn-in chain
%BaGoLParams.N_Burnin = 8000;        % Length of Burn-in chain
%BaGoLParams.N_Trials = 2000;        % Length of post-burn-in chain
BaGoLParams.NSamples = 10;          % Number of samples before sampling Xi
BaGoLParams.ClusterDrift = 0;       % Expected magnitude of drift (nm/frame)

% Y_Adjust is sometimes needed to deal with lower left versus upper left
% y-origin issues.  Lower left with y increasing upwards is the default,
% requiring no changes, while upper left with y increasing downwards can
% sometimes occur, so Y must be replaced by Y - Y_Adjust, where Y_Adjust is the
% image size (see below) [pixels].
%BaGoLParams.Y_Adjust = BaGoLParams.ImageSize;
BaGoLParams.Y_Adjust = [];

% SE_Adjust adds to X_SE and Y_SE, so inflates the precision.  For DNA_PAINT
% data, SE_Adjust = 1--2 nm, while for dSTORM, slightly bigger values should
% be used.  Note that this quantity can be specified as an array of length
% n_files if applied differently to each file.
BaGoLParams.SE_Adjust = 0;          % Precision inflation applied to SE (nm)
%BaGoLParams.SE_Adjust = [0, 0];     % Precision inflation applied to SE (nm)

% The values for ROIsz and OverLap directly below are good overall for much
% data, but note that the larger the ROIsz, the more the computational effort.
% Artifacts in dense data can come about if the ROIsz is too large.  The
% pre-clustering cutoff should be around the localization precision.
BaGoLParams.ROIsz = 500;            % ROI size for RJMCMC (nm)
BaGoLParams.OverLap = 50;           % Size of overlapping region (nm)
BaGoLParams.Cutoff = 25;            % Pre-clustering cutoff (nm)
%BaGoLParams.ROIsz = 100;            % ROI size for RJMCMC (nm)
%BaGoLParams.OverLap = 25;           % Size of overlapping region (nm)
%BaGoLParams.ROIsz = 50;             % ROI size for RJMCMC (nm)
%BaGoLParams.OverLap = 10;           % Size of overlapping region (nm)

% k and theta below are the shape and scale parameters for the Gamma
% probability distribution function.  If just one parameter is provided,
% a Poisson distribution is used.
BaGoLParams.Xi = [20, 1];           % [k, theta] parameters for gamma prior

% Note for batch runs, in which Files and DataROI are input by hand, please see
% ### comments below.
BaGoLParams.DataROI = [];           % [Xmin, Xmax, Ymin, Ymax] (pixel)
DataROI = [];

% If ROIs is true, the input file has ROIs already defined (\*\_ROIs.mat),
% so use them below if only one filename is provided.
ROIs = false;
```
Results files (produced by SMLM analyis) to be batch-processed.
```
D1 = 'DATA';
Files = {
fullfile(D1, 'Cell_02_Label_01_Results.mat');
fullfile(D1, 'Cell_03_Label_01_Results.mat');
};
```
DataROI is used automatically when ROIs is true, so should be left commented
out in this situation, otherwise it can be helpful in downsizing the
computational effort needed for very dense datasets.
```
% DataROI is defined when running BaGoL over only part of the image.
% If DataROI is empty, use the whole image.
% 
% Define a single region of interest for each dataset (units are pixels).
% [YStart, XStart, YEnd, XEnd] = [163, 385, 233, 455]
% [Xmin, Xmax, Ymin, Ymax] (pixel)
% [385, 455, 163, 233]

%DataROI = [
%[120, 136, 190, 206]
%[110, 126,  90, 106]
%];
```
