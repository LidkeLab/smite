---
title: "API Reference: +smi Namespace"
category: "api-reference"
level: "intermediate"
tags: ["api", "smi", "SMLM", "SPT", "BaGoL", "Publish", "workflows"]
prerequisites: ["../core-concepts/architecture.md", "../core-concepts/smf-structure.md", "../core-concepts/smd-structure.md"]
related: ["../workflows/smlm-analysis.md", "../workflows/spt-tracking.md", "../workflows/bagol-clustering.md", "../workflows/batch-processing.md"]
summary: "Complete API reference for the top-level +smi namespace classes: SMLM, SPT, BaGoL, and Publish"
estimated_time: "30 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# API Reference: +smi Namespace

## Purpose

This document provides a complete API reference for the four main classes in smite's top-level `+smi` namespace: `SMLM`, `SPT`, `BaGoL`, and `Publish`. These classes represent the primary entry points for complete analysis workflows in smite. Understanding their properties, methods, and usage patterns is essential for both GUI-based and script-based analyses.

## Prerequisites

- Understanding of [smite architecture](../core-concepts/architecture.md)
- Familiarity with [SMF structure](../core-concepts/smf-structure.md)
- Familiarity with [SMD structure](../core-concepts/smd-structure.md)
- Basic MATLAB programming knowledge

## Overview

The `+smi` namespace contains four high-level workflow classes that orchestrate complete analyses from raw data to results:

- **`smi.SMLM`**: Single Molecule Localization Microscopy - localize and analyze blinking emitters
- **`smi.SPT`**: Single Particle Tracking - track moving particles over time
- **`smi.BaGoL`**: Bayesian Grouping of Localizations - determine true emitter positions from multiple localizations
- **`smi.Publish`**: Batch processing of structured datasets with standardized directory organization

Each class provides:
- A complete analysis pipeline
- Optional GUI interface
- Configurable parameters via SMF structure
- Automatic results saving and visualization
- Both interactive and programmatic usage modes

---

## smi.SMLM

### Description

`smi.SMLM` performs complete Single Molecule Localization Microscopy analysis. It handles loading raw camera data, finding and fitting single molecules, filtering poor quality localizations, connecting localizations across frames (for blinking emitters), correcting for stage drift, and generating super-resolution images and diagnostic plots.

### Class Definition

```matlab
classdef SMLM < handle
```

### Constructor

```matlab
obj = smi.SMLM()              % Opens GUI
obj = smi.SMLM(SMF)           % Uses provided SMF, no GUI
obj = smi.SMLM(SMF, StartGUI) % Control GUI opening
```

**Parameters:**
- `SMF`: Single Molecule Fitting structure (optional) - defines all analysis parameters
- `StartGUI`: Boolean flag to open GUI (optional, default: true if no args, false if SMF provided)

**Returns:**
- `obj`: SMLM object instance

**Example:**
```matlab
% GUI mode
SMLMobj = smi.SMLM();  % GUI opens automatically

% Script mode
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'experiment.h5'};
SMLMobj = smi.SMLM(SMF);
```

### Properties

#### Public Properties

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `SMD` | struct | [] | Final Single Molecule Data results structure |
| `SMF` | struct | - | Single Molecule Fitting parameter structure |
| `SMDPreThresh` | struct | [] | SMD before thresholding (tracks rejected localizations) |
| `FitFramePreFC` | cell | {} | Fits per frame before frame connection |
| `PlotDo` | string array | [] | Specify which plots to generate (all by default) |
| `SRImageZoom` | numeric | 20 | Magnification factor for super-resolution images |
| `SRCircImZoom` | numeric | 25 | Magnification factor for circle images |
| `Verbose` | integer | 1 | Verbosity level (0=silent, 1=normal, 2+=debug) |
| `VerboseTest` | integer | 3 | Verbosity level for test fits |
| `CalledByGUI` | logical | false | Internal flag tracking GUI invocation |
| `FCCount` | integer | 0 | Number of frame-connected localizations eliminated |

#### Hidden Properties

| Property | Type | Description |
|----------|------|-------------|
| `DC` | object | DriftCorrection class object |
| `DCMethod` | string | Drift correction method ('DC-KNN' or 'DC-BF') |
| `ResultsDir` | string | Top-level results directory |
| `ResultsSubDir` | string | Subdirectory for specific analysis results |
| `ShowPlots` | logical | Whether to display plots interactively |

### Methods

#### Primary Workflow Methods

##### fullAnalysis

```matlab
obj.fullAnalysis()
```

Performs complete analysis on all datasets specified in `SMF.Data.DatasetList`. This is the main method for production analyses.

**Workflow:**
1. Creates output directories
2. For each dataset:
   - Loads raw data
   - Converts to photons
   - Localizes molecules
   - Applies thresholding
   - Performs frame connection (if enabled)
   - Corrects intra-dataset drift (if enabled)
3. Concatenates all datasets
4. Corrects inter-dataset drift (if enabled)
5. Saves results and generates plots

**Example:**
```matlab
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();
% Results saved to SMF.Data.ResultsDir/DatasetName/
```

##### testFit

```matlab
obj.testFit(DatasetIndex)
```

Performs detailed analysis with extensive diagnostics on a single dataset. Use this to verify parameters before running full analysis.

**Parameters:**
- `DatasetIndex`: Index of dataset to analyze (from `SMF.Data.DatasetList`)

**Features:**
- More verbose output (`VerboseTest` level)
- Shows plots interactively
- Saves diagnostic information
- Faster than full analysis

**Example:**
```matlab
% Test on first dataset
SMLMobj.testFit(1);

% Examine results
disp(SMLMobj.SMD);
fprintf('Found %d localizations\n', length(SMLMobj.SMD.X));
```

##### analyzeAll

```matlab
obj.analyzeAll()
obj.analyzeAll(DatasetList)
```

Loops over datasets and creates SMD results. This is an internal method called by `fullAnalysis` and `testFit`.

**Parameters:**
- `DatasetList`: Array of dataset indices to analyze (optional, uses `SMF.Data.DatasetList` if not provided)

**Internal flow:**
```
analyzeAll
  └─► for each dataset
       └─► analyzeDataset
            ├─► LoadData
            ├─► DataToPhotons
            ├─► LocalizeData
            ├─► FrameConnection
            └─► DriftCorrection (intra)
  └─► DriftCorrection (inter-dataset)
  └─► filterNN (nearest neighbor filter)
```

##### analyzeDataset

```matlab
SMD = obj.analyzeDataset(DatasetIndex, DatasetCount)
```

Loads and analyzes a single dataset. Internal method.

**Parameters:**
- `DatasetIndex`: Index of dataset in file list
- `DatasetCount`: Sequential count for numbering (optional)

**Returns:**
- `SMD`: Single Molecule Data structure for this dataset

#### Visualization and Output

##### generatePlots

```matlab
obj.generatePlots(PlotSaveDir1, PlotSaveDir2, AnalysisID, ShowPlots, PlotDo)
```

Creates all histograms, statistics, and super-resolution images.

**Parameters:**
- `PlotSaveDir1`: Directory for priority plots (e.g., GaussIm)
- `PlotSaveDir2`: Directory for diagnostic plots
- `AnalysisID`: String to append to filenames
- `ShowPlots`: Boolean - display plots interactively
- `PlotDo`: String array of plots to generate

**Available plots:**
- `"Photons"` - Photon count histogram
- `"Bg"` - Background histogram
- `"PSFSigma"` - PSF width histogram
- `"PValue"` - Fit quality histogram
- `"X_SE"`, `"Y_SE"`, `"Z_SE"` - Precision histograms
- `"NCombined"` - Frame connections histogram
- `"CumDrift"` - Cumulative drift plot
- `"Drift"` - Absolute drift plot
- `"FitFrame"` - Localizations per frame
- `"DriftIm"` - Color-coded drift image
- `"GaussIm"` - Gaussian rendering
- `"HistIm"` - Histogram rendering
- `"CircleIm"` - Localization circles
- `"CircleImDrift"` - Time-coded circles

**Example:**
```matlab
% Generate only key plots
SMLMobj.PlotDo = ["GaussIm", "Photons", "X_SE"];
SMLMobj.generatePlots(resultsDir, resultsDir, 'test', true, SMLMobj.PlotDo);
```

##### saveResults

```matlab
obj.saveResults()
```

Saves SMD, SMF structures and generates all plots. Called automatically by `fullAnalysis()` and `testFit()`.

**Saves:**
- `DatasetName_Results.mat` containing `SMD` and `SMF`
- All plots specified in `PlotDo`
- Rejection statistics in `reject.txt`

#### GUI Methods

##### gui

```matlab
obj.gui()
```

Opens the interactive SMLM GUI for parameter configuration and analysis control.

**Features:**
- Browse for data files
- Configure all SMF parameters
- Run test fits
- Execute full analysis
- View results

**Example:**
```matlab
SMLMobj = smi.SMLM();
% GUI opens automatically
% Or explicitly:
SMLMobj.gui();
```

#### Utility Methods

##### createDirectories

```matlab
obj.createDirectories()
```

Creates directory structure for saving results. Called internally by `analyzeAll()`.

**Directory structure:**
```
ResultsDir/
  └─ DatasetName_AnalysisID/
      ├─ Results.mat
      ├─ GaussImage.png
      ├─ TestFit/  (for test runs)
      └─ ... diagnostic plots ...
```

### Static Methods

##### fitsPerFrame

```matlab
FitFrame = smi.SMLM.fitsPerFrame(SMD, DatasetIndex)
```

Computes number of localizations in each frame.

**Parameters:**
- `SMD`: Single Molecule Data structure
- `DatasetIndex`: Dataset index (optional)

**Returns:**
- `FitFrame`: Vector of counts per frame

**Example:**
```matlab
counts = smi.SMLM.fitsPerFrame(SMD);
plot(counts);
xlabel('Frame'); ylabel('Localizations');
```

##### unitTest

```matlab
Success = smi.SMLM.unitTest()
```

Runs comprehensive unit tests for SMLM workflow.

**Returns:**
- `Success`: Boolean indicating test success

### Usage Patterns

#### Pattern 1: GUI-Based Analysis

```matlab
% Start with GUI
SMLMobj = smi.SMLM();

% Use GUI to:
% - Set file paths
% - Configure parameters
% - Run test fit
% - Execute full analysis

% Access results
SMD = SMLMobj.SMD;
```

#### Pattern 2: Script-Based Analysis

```matlab
% Configure SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/2024-01-10';
SMF.Data.FileName = {'Cell1.h5'};
SMF.Data.PixelSize = 0.108;  % 108 nm
SMF.Data.CameraGain = 2.5;
SMF.Data.CameraOffset = 100;

% Box finding
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 250;

% Fitting
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Thresholding
SMF.Thresholding.MaxXY_SE = 0.15;
SMF.Thresholding.MinPhotons = 150;

% Frame connection
SMF.FrameConnection.On = true;
SMF.FrameConnection.MaxSeparation = 1.0;

% Drift correction
SMF.DriftCorrection.On = true;

% Run analysis
SMLMobj = smi.SMLM(SMF);
SMLMobj.Verbose = 1;
SMLMobj.fullAnalysis();

% Results
SMD = SMLMobj.SMD;
fprintf('Detected %d localizations\n', length(SMD.X));
```

#### Pattern 3: Parameter Testing

```matlab
% Load existing SMF
load('previous_analysis_Results.mat', 'SMF');

% Modify specific parameter
SMF.Thresholding.MaxXY_SE = 0.10;  % Stricter threshold

% Test on subset
SMF.Data.DatasetList = int32(1);
SMLMobj = smi.SMLM(SMF);
SMLMobj.testFit(1);

% Examine effect
fprintf('Original: %d, New: %d localizations\n', ...
    length(previous_SMD.X), length(SMLMobj.SMD.X));
```

### See Also

- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete workflow guide
- [First Analysis](../getting-started/first-analysis.md) - Hands-on tutorial
- `smi_core.LocalizeData` - Core localization engine
- `smi_core.FrameConnection` - Frame connection algorithm
- `smi_core.DriftCorrection` - Drift correction methods

---

## smi.SPT

### Description

`smi.SPT` performs Single Particle Tracking analysis. It localizes particles in each frame, then connects localizations across frames to form trajectories representing particle motion. The class supports frame-to-frame linking, gap closing (bridging missing detections), diffusion analysis, and optional channel registration for multi-color tracking.

### Class Definition

```matlab
classdef SPT < handle
```

### Constructor

```matlab
obj = smi.SPT()               % Default SMF, opens GUI
obj = smi.SPT(SMF)            % Uses provided SMF, opens GUI
obj = smi.SPT(SMF, StartGUI)  % Control GUI opening
```

**Parameters:**
- `SMF`: Single Molecule Fitting structure (optional)
- `StartGUI`: Boolean flag to open GUI (optional, default: true)

**Returns:**
- `obj`: SPT object instance

**Example:**
```matlab
% GUI mode
SPTobj = smi.SPT();

% Script mode
SMF = smi_core.SingleMoleculeFitting();
SPTobj = smi.SPT(SMF, false);
```

### Properties

#### Public Properties

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `SMF` | struct | - | Single Molecule Fitting parameter structure |
| `TransformDir` | string | '' | Directory containing channel registration transforms |
| `TransformPattern` | string | 'RegistrationTransform*.mat' | Pattern to match transform files |
| `FindFiles` | logical | true | Search for files in batchTrack() |
| `FilePattern` | string | '*.mat' | Pattern to match data files |
| `DiffusionEstimator` | object | - | Instance of `smi_stat.DiffusionEstimator` |
| `MovieParams` | struct | - | Parameters for generating trajectory movies |
| `NonLinkMarker` | numeric | -1 | Marker for non-linkable entries in cost matrices |
| `GenerateMovies` | logical | true | Flag to generate trajectory movies |
| `GeneratePlots` | logical | true | Flag to generate diagnostic plots |
| `IsTestRun` | logical | false | Flag indicating test mode (no saving) |
| `UnitFlag` | logical | false | Output in physical units |
| `UseSparseMatrices` | logical | true | Use sparse matrices for gap closing |
| `Verbose` | integer | 1 | Verbosity level |

#### Protected Properties

| Property | Type | Description |
|----------|------|-------------|
| `SMD` | struct | Thresholded localizations |
| `SMDPreThresh` | struct | All localizations before thresholding |
| `SMDBatch` | struct | Concatenated SMD across batch |
| `TR` | struct array | Tracking Results - array of SMD, one per trajectory |
| `TRPreCR` | struct | TR before channel registration |
| `SMDPreCR` | struct | SMD before channel registration |
| `ScaledData` | array | Photon-converted image stack |
| `RhoOff` | numeric | Density of dark emitters |
| `SMLM` | object | Internal SMLM instance |

### Methods

#### Primary Workflow Methods

##### performFullAnalysis

```matlab
[TR, SMD, SMDPreThresh] = obj.performFullAnalysis()
```

Main analysis method that performs complete tracking workflow from raw data to trajectories.

**Workflow:**
1. Loads raw data
2. Performs gain/offset correction
3. Localizes particles in each frame
4. Creates trajectories via frame-to-frame linking
5. Performs gap closing (bridges missing detections)
6. Applies channel registration (if specified)
7. Filters short trajectories
8. Saves results

**Returns:**
- `TR`: Tracking Results structure array
- `SMD`: Single Molecule Data with tracking IDs
- `SMDPreThresh`: Pre-threshold SMD

**Example:**
```matlab
% Configure for tracking
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileName = {'particles.h5'};
SMF.Tracking.MaxDistFF = 2.0;  % pixels, frame-to-frame
SMF.Tracking.MaxDistGC = 4.0;  % pixels, gap closing

SPTobj = smi.SPT(SMF, false);
[TR, SMD] = SPTobj.performFullAnalysis();

% Analyze trajectories
fprintf('Found %d trajectories\n', length(TR));
fprintf('Mean trajectory length: %.1f frames\n', ...
    mean(arrayfun(@(x) length(x.FrameNum), TR)));
```

##### batchTrack

```matlab
[TR, SMD, SMDPreThresh, FileList, TransformList] = obj.batchTrack()
```

Performs tracking on multiple files in batch mode with optional channel registration.

**Workflow:**
1. Finds all files matching `FilePattern` in `SMF.Data.FileDir`
2. For each file:
   - Localizes particles
   - Generates trajectories
   - Applies channel registration (if transforms available)
3. Concatenates results across files

**Returns:**
- `TR`: Tracking Results (all files concatenated)
- `SMD`: Single Molecule Data (all files)
- `SMDPreThresh`: Pre-threshold SMD
- `FileList`: Cell array of processed files
- `TransformList`: Cell array of registration transforms used

**Example:**
```matlab
% Setup for batch tracking
SPTobj = smi.SPT(SMF, false);
SPTobj.FindFiles = true;
SPTobj.FilePattern = '*Cell*.mat';
SPTobj.TransformDir = '/path/to/transforms';

[TR, SMD] = SPTobj.batchTrack();
```

##### autoTrack

```matlab
obj.autoTrack()
```

Automatically generates trajectories from localizations in `obj.SMD`. This method orchestrates frame-to-frame linking and gap closing.

**Internal workflow:**
1. Estimates emitter densities
2. Performs frame-to-frame linking
3. Performs gap closing
4. Converts SMD to TR format
5. Estimates diffusion coefficients (if enabled)

**Called by:** `performFullAnalysis()` and `batchTrack()`

##### generateTrajectories

```matlab
obj.generateTrajectories()
```

Generates trajectories from localizations using current tracking parameters. Updates `obj.TR` and `obj.SMD`.

**Used when:** Reprocessing existing localizations with different tracking parameters

**Example:**
```matlab
% Load previous localizations
load('localizations.mat', 'SMD');
SPTobj.SMD = SMD;

% Try different tracking parameters
SPTobj.SMF.Tracking.MaxDistFF = 1.5;  % Stricter
SPTobj.generateTrajectories();
```

##### updateTrackingParams

```matlab
obj.updateTrackingParams(IsBatch)
```

Updates tracking parameters based on analysis results. Internal method called during iterative parameter refinement.

**Parameters:**
- `IsBatch`: Boolean indicating batch processing mode

#### Visualization and Output

##### saveResults

```matlab
obj.saveResults()
```

Saves tracking results, generates plots, and creates movies.

**Saves:**
- `DatasetName_TR.mat` - Tracking Results
- `DatasetName_SMD.mat` - SMD with trajectory IDs
- Trajectory movies (if `GenerateMovies = true`)
- Diagnostic plots (if `GeneratePlots = true`)

**Example:**
```matlab
SPTobj.GenerateMovies = true;
SPTobj.MovieParams.MIPScale = 30;  % Super-resolution scale
SPTobj.saveResults();
```

#### GUI Methods

##### gui

```matlab
obj.gui()
```

Opens interactive SPT GUI for parameter configuration and tracking control.

### Static Methods

##### genTrajFF

```matlab
SMD = smi.SPT.genTrajFF(SMD, SMF, RhoOff, NonLinkMarker)
```

Performs frame-to-frame trajectory linking using Linear Assignment Problem (LAP).

**Parameters:**
- `SMD`: Single Molecule Data structure
- `SMF`: Single Molecule Fitting structure
- `RhoOff`: Density of dark emitters
- `NonLinkMarker`: Value marking non-linkable entries

**Returns:**
- `SMD`: Updated with `ConnectID` field linking localizations

##### genTrajGC

```matlab
SMD = smi.SPT.genTrajGC(SMD, SMF, RhoOff, NonLinkMarker, UseSparseMatrices)
```

Performs gap-closing to bridge missing detections in trajectories.

**Parameters:**
- `SMD`: Single Molecule Data with frame-to-frame links
- `SMF`: Single Molecule Fitting structure
- `RhoOff`: Density of dark emitters
- `NonLinkMarker`: Non-link marker value
- `UseSparseMatrices`: Use sparse cost matrix

**Returns:**
- `SMD`: Updated with gap-closing connections

##### createCostMatrixFF

```matlab
CostMatrix = smi.SPT.createCostMatrixFF(SMD, SMF, RhoOff, FrameNumber, NonLinkMarker)
```

Creates cost matrix for frame-to-frame linking between two consecutive frames.

**Parameters:**
- `SMD`: Single Molecule Data
- `SMF`: Single Molecule Fitting structure
- `RhoOff`: Density of off emitters
- `FrameNumber`: Current frame number
- `NonLinkMarker`: Non-link marker

**Returns:**
- `CostMatrix`: Cost matrix for LAP solver

##### createCostMatrixGC

```matlab
[CostMatrix, StartEndIndices] = smi.SPT.createCostMatrixGC(SMD, SMF, RhoOff, NonLinkMarker, CreateSparseMatrix)
```

Creates cost matrix for gap closing.

**Parameters:**
- `SMD`: Single Molecule Data
- `SMF`: Single Molecule Fitting structure
- `RhoOff`: Density of off emitters
- `NonLinkMarker`: Non-link marker
- `CreateSparseMatrix`: Boolean for sparse format

**Returns:**
- `CostMatrix`: Gap-closing cost matrix
- `StartEndIndices`: Indices of trajectory starts/ends

##### solveLAP

```matlab
[Assign12, Cost12] = smi.SPT.solveLAP(CostMatrix, NonlinkMarker)
```

Solves Linear Assignment Problem to find optimal linking.

**Parameters:**
- `CostMatrix`: Cost matrix
- `NonlinkMarker`: Non-link marker value

**Returns:**
- `Assign12`: Assignment vector
- `Cost12`: Total cost

##### estimateDensities

```matlab
[RhoOff, RhoOn] = smi.SPT.estimateDensities(SMD, SMF)
```

Estimates density of dark and bright emitters for cost matrix calculation.

**Parameters:**
- `SMD`: Single Molecule Data
- `SMF`: Single Molecule Fitting structure

**Returns:**
- `RhoOff`: Density of dark emitters (per area)
- `RhoOn`: Density of bright emitters

##### unitTestFFGC

```matlab
Success = smi.SPT.unitTestFFGC()
```

Unit test for frame-to-frame and gap-closing algorithms.

**Returns:**
- `Success`: Boolean indicating test success

### Usage Patterns

#### Pattern 1: Basic Tracking

```matlab
% Configure SMF for tracking
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/tracking/';
SMF.Data.FileName = {'particles.h5'};
SMF.Data.PixelSize = 0.108;  % µm
SMF.Data.FrameRate = 100;    % Hz

% Tracking parameters
SMF.Tracking.MaxDistFF = 2.0;        % pixels
SMF.Tracking.MaxDistGC = 4.0;        % pixels
SMF.Tracking.MaxFrameGap = 2;        % frames
SMF.Tracking.MinTrackLength = 5;     % frames

% Run tracking
SPTobj = smi.SPT(SMF, false);
[TR, SMD] = SPTobj.performFullAnalysis();

% Analyze results
nTraj = length(TR);
lengths = arrayfun(@(x) length(x.FrameNum), TR);
fprintf('%d trajectories, mean length %.1f frames\n', ...
    nTraj, mean(lengths));
```

#### Pattern 2: Multi-File Batch Tracking

```matlab
% Setup batch parameters
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/experiment/';

SPTobj = smi.SPT(SMF, false);
SPTobj.FindFiles = true;
SPTobj.FilePattern = 'Cell*.h5';
SPTobj.Verbose = 2;

% Run batch
[TR, SMD, ~, FileList] = SPTobj.batchTrack();

fprintf('Processed %d files\n', length(FileList));
fprintf('Total trajectories: %d\n', length(TR));
```

#### Pattern 3: Tracking with Diffusion Analysis

```matlab
% Enable trajectory-wise diffusion estimation
SMF.Tracking.TrajwiseD = true;

SPTobj = smi.SPT(SMF, false);

% Configure diffusion estimator
SPTobj.DiffusionEstimator.FitTarget = 'LikelihoodOfJumps';
SPTobj.DiffusionEstimator.FitIndividualTrajectories = true;
SPTobj.DiffusionEstimator.FrameLagRange = [2, 2];

[TR, SMD] = SPTobj.performFullAnalysis();

% Access diffusion coefficients
D_values = arrayfun(@(x) x.D, TR);
histogram(D_values * 1e12);  % Convert to µm²/s
xlabel('Diffusion Coefficient (\mum^2/s)');
```

### See Also

- [SPT Workflow](../workflows/spt-tracking.md) - Complete tracking guide
- [Tracking Example](../examples/tracking-diffusion.md) - Hands-on example
- `smi_core.TrackingResults` - TR structure definition
- `smi_stat.DiffusionEstimator` - Diffusion analysis

---

## smi.BaGoL

### Description

`smi.BaGoL` implements Bayesian Grouping of Localizations, a method for determining true emitter positions when multiple localizations arise from the same blinking or binding emitter. BaGoL uses Reversible Jump Markov Chain Monte Carlo (RJMCMC) to explore possible numbers of emitters and their positions, producing both a posterior probability image and Maximum A Posteriori Number (MAPN) coordinates with improved precision.

### Class Definition

```matlab
classdef BaGoL < handle
```

### Constructor

```matlab
obj = smi.BaGoL()
```

Creates BaGoL object. Configure properties, set `SMD`, then call `analyze_all()`.

**Example:**
```matlab
B = smi.BaGoL();
B.SMD = mySMD;  % From SMLM analysis
B.ROIsize = 200;  % nm
B.Xi = 20;  % Expected localizations per emitter
B.analyze_all();
```

### Properties

#### Algorithm Parameters

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `Xi` | numeric | 20 | Localizations per emitter: scalar (Poisson) or [k, theta] (Gamma) |
| `Alpha_Xi` | numeric | 1 | Shape parameter for Xi gamma hyper-prior |
| `Beta_Xi` | numeric | 50 | Scale parameter for Xi gamma hyper-prior |
| `N_Burnin` | integer | 2000 | RJMCMC burn-in samples |
| `N_Trials` | integer | 3000 | RJMCMC post-burn-in samples |
| `NSamples` | integer | 10 | Samples before Xi update (hierarchical mode) |
| `P_Jumps` | array | [0.25, 0.25, 0.25, 0.25] | Proposal probabilities [Move, Allocate, Add, Remove] |
| `HierarchFlag` | logical | 0 | Learn Xi via hierarchical Bayes |
| `ChainFlag` | logical | 0 | Save RJMCMC chains |

#### Data Organization

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `ROIsize` | numeric | 200 | ROI size for RJMCMC (nm) |
| `Overlap` | numeric | 20 | Overlap between ROIs (nm) |
| `Cutoff` | numeric | ROIsize | Pre-clustering cutoff (nm) |
| `Drift` | numeric | 0 | Expected drift magnitude (nm/frame) |
| `SE_Adjust` | numeric | 0 | Localization precision adjustment (nm) |

#### Output Control

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `PImageFlag` | logical | 1 | Generate posterior image |
| `PixelSize` | numeric | 2 | Posterior image pixel size (nm) |
| `PImageSize` | numeric | [] | Posterior image size (nm) |
| `XStart` | numeric | [] | Posterior image X start (nm) |
| `YStart` | numeric | [] | Posterior image Y start (nm) |
| `SaveName` | string | 'BaGoL' | Base name for saved results |

#### Results

| Property | Type | Description |
|----------|------|-------------|
| `SMD` | struct | Input localizations with X, Y, X_SE, Y_SE, FrameNum, PixelSize |
| `ClusterSMD` | struct array | SMD for each pre-cluster subregion |
| `MAPN` | struct | Output emitter positions and uncertainties |
| `PImage` | array | Posterior probability image |
| `Chain` | cell array | RJMCMC chains (if ChainFlag=1) |
| `XiChain` | array | Chain of Xi samples (hierarchical mode) |

### Methods

#### Primary Analysis Method

##### analyze_all

```matlab
obj = obj.analyze_all()
```

Performs complete BaGoL analysis on the SMD structure.

**Workflow:**
1. Initializes posterior image parameters
2. Filters localizations with zero precision
3. Generates ROI grid
4. Pre-clusters localizations within ROIs
5. For each cluster:
   - Runs RJMCMC sampler
   - Generates posterior image contribution
   - Extracts MAPN coordinates
6. Collates results

**Example:**
```matlab
% Load SMLM results
load('SMLM_Results.mat', 'SMD', 'SMF');

% Configure BaGoL
B = smi.BaGoL();
B.SMD = SMD;
B.SMD.PixelSize = SMF.Data.PixelSize;
B.ROIsize = 200;
B.Overlap = 20;
B.Xi = 15;  % Expect ~15 localizations per emitter
B.N_Burnin = 2000;
B.N_Trials = 3000;
B.PImageFlag = 1;
B.PixelSize = 1;  % 1 nm pixel in output

% Run analysis
B.analyze_all();

% Access results
MAPN = B.MAPN;
fprintf('Found %d emitters from %d localizations\n', ...
    length(MAPN.X), length(SMD.X));

% Display posterior image
figure; imagesc(B.PImage); axis image;
colormap(hot); title('BaGoL Posterior Image');
```

#### ROI Management

##### genROIs

```matlab
ROIs = obj.genROIs()
```

Splits input coordinates into overlapping rectangular subregions.

**Returns:**
- `ROIs`: 2D array of SMD structures, one per subregion

**Example:**
```matlab
B = smi.BaGoL();
B.SMD = mySMD;
B.ROIsize = 200;
B.Overlap = 20;
ROIs = B.genROIs();
fprintf('Created %d x %d ROI grid\n', size(ROIs));
```

##### precluster

```matlab
obj.precluster(ROIs)
```

Performs hierarchical clustering within ROIs to create analysis subregions. Updates `obj.ClusterSMD`.

**Parameters:**
- `ROIs`: ROI array from `genROIs()`

##### assignROIs

```matlab
obj.assignROIs(ROIs)
```

Alternative ROI assignment method. Assigns localizations to ROIs without pre-clustering.

**Parameters:**
- `ROIs`: ROI array

#### Posterior Image Generation

##### genPosterior

```matlab
PostIm = obj.genPosterior(PostIm, SZ, TChain, ROIs, nn)
```

Generates posterior image from RJMCMC chain for one cluster.

**Parameters:**
- `PostIm`: Accumulator posterior image
- `SZ`: Image size
- `TChain`: RJMCMC chain for this cluster
- `ROIs`: ROI array
- `nn`: Cluster index

**Returns:**
- `PostIm`: Updated posterior image

##### genMAPN

```matlab
obj.genMAPN(TChain, ROIs, nn)
```

Extracts MAPN coordinates from RJMCMC chain.

**Parameters:**
- `TChain`: RJMCMC chain
- `ROIs`: ROI array
- `nn`: Cluster index

**Updates:** `obj.MAPN` by appending emitter coordinates

#### Initialization

##### initPostIm

```matlab
obj.initPostIm()
```

Automatically calculates posterior image size and starting coordinates from input data.

### Static Methods

##### BaGoL_RJMCMC

```matlab
Chain = smi.BaGoL.BaGoL_RJMCMC(SMD, Xi, MaxAlpha, PMove, NChain, NBurnin, DEBUG, MuX, MuY, AlphaX, AlphaY)
```

Core RJMCMC sampler for fixed Xi (non-hierarchical mode).

**Parameters:**
- `SMD`: Localization data for one cluster
- `Xi`: Localizations per emitter parameter
- `MaxAlpha`: Maximum drift velocity (nm/frame)
- `PMove`: Jump proposal probabilities [Move, Allocate, Add, Remove]
- `NChain`: Number of MCMC samples
- `NBurnin`: Burn-in samples
- `DEBUG`: Animation flag (0 or 1)
- `MuX`, `MuY`, `AlphaX`, `AlphaY`: Initial values (optional)

**Returns:**
- `Chain`: Structure array with fields N (num emitters), X, Y, AlphaX, AlphaY, ID

##### BaGoL_RJMCMC_Hierarchical

```matlab
[K, X, Y, AlphaX, AlphaY, ID] = smi.BaGoL.BaGoL_RJMCMC_Hierarchical(SMD, PDFgrid, AlphaSTD, P_Jumps, NSamples, Xi, MuX, MuY, AlphaX, AlphaY)
```

RJMCMC sampler for hierarchical Bayes mode (learning Xi).

**Returns:**
- `K`: Number of emitters
- `X`, `Y`: Emitter positions
- `AlphaX`, `AlphaY`: Drift velocities
- `ID`: Localization-to-emitter assignments

##### samplePoiss

```matlab
Xi = smi.BaGoL.samplePoiss(NPoints, K, Xi, Alpha, Beta)
```

Samples Poisson Xi parameter from posterior in hierarchical mode.

##### sampleGam

```matlab
Xi = smi.BaGoL.sampleGam(NPoints, K, Xi, Alpha, Beta)
```

Samples Gamma Xi parameters from posterior in hierarchical mode.

##### makeIm

```matlab
[SRIm, MapIm] = smi.BaGoL.makeIm(SMD, MAPN, SZ, PixSize, XStart, YStart)
```

Generates super-resolution images from SMD and MAPN coordinates.

**Parameters:**
- `SMD`: Input localizations
- `MAPN`: MAPN coordinates
- `SZ`: Image size (pixels)
- `PixSize`: Pixel size (nm)
- `XStart`, `YStart`: Image origin (nm)

**Returns:**
- `SRIm`: SR image from input SMD
- `MapIm`: SR image from MAPN coordinates

##### genSRMAPNOverlay

```matlab
Im = smi.BaGoL.genSRMAPNOverlay(SMD, MAPN, XSize, YSize, PixelSize, SaveDir, XStart, YStart, RadiusScale, ScaleBarLength)
```

Creates overlay image comparing input localizations and MAPN emitters.

**Parameters:**
- `SMD`, `MAPN`: Input and output data
- `XSize`, `YSize`: Image dimensions (pixels)
- `PixelSize`: Pixel size (nm)
- `SaveDir`: Directory to save image
- `XStart`, `YStart`: Origin (nm)
- `RadiusScale`: Circle size scaling
- `ScaleBarLength`: Scale bar length (nm)

**Returns:**
- `Im`: RGB overlay image

##### saveMAPN

```matlab
smi.BaGoL.saveMAPN(Directory, FileType, MAPN)
```

Saves MAPN coordinates to file.

**Parameters:**
- `Directory`: Save directory
- `FileType`: 'mat', 'csv', or 'both'
- `MAPN`: MAPN structure

##### frameConnect

```matlab
[SMD, SMD_combined] = smi.BaGoL.frameConnect(SMDin, LOS, MaxDistance, MaxFrameGap, FitType)
```

Frame connects localizations before BaGoL analysis.

**Parameters:**
- `SMDin`: Input SMD
- `LOS`: Line-of-sight frames
- `MaxDistance`: Max linking distance (pixels)
- `MaxFrameGap`: Max frame gap
- `FitType`: Fit model type

**Returns:**
- `SMD`: Frame-connected SMD
- `SMD_combined`: Combined SMD with improved precision

##### genCluster

```matlab
[SMD, SMC] = smi.BaGoL.genCluster(StructType, Scale, Ndist, PSF, MeanPhoton, Prec_Cutoff, DriftVec, PlotFlag)
```

Generates synthetic clustering data for testing.

##### loadPICASSOh5

```matlab
SMD = smi.BaGoL.loadPICASSOh5(DataDir, FileName)
```

Loads localizations from Picasso HDF5 format.

##### hierBaGoL_analysis

```matlab
BGL = smi.BaGoL.hierBaGoL_analysis(SMD, FileNameIn, SaveDir, BaGoLParams)
```

High-level hierarchical BaGoL analysis wrapper.

### Usage Patterns

#### Pattern 1: Basic BaGoL Analysis

```matlab
% Load SMLM results
load('SMLM_Results.mat', 'SMD', 'SMF');

% Setup BaGoL
B = smi.BaGoL();
B.SMD = SMD;
B.SMD.PixelSize = SMF.Data.PixelSize;

% Key parameters
B.ROIsize = 200;      % 200 nm ROIs
B.Xi = 20;            % ~20 localizations per emitter
B.N_Burnin = 2000;
B.N_Trials = 3000;
B.PImageFlag = 1;
B.PixelSize = 2;      % 2 nm/pixel in posterior

% Run
B.analyze_all();

% Results
MAPN = B.MAPN;
PImage = B.PImage;
fprintf('Input: %d locs → Output: %d emitters\n', ...
    length(SMD.X), length(MAPN.X));

% Precision improvement
input_prec = median(SMD.X_SE) * SMD.PixelSize;
output_prec = median(MAPN.X_SE);
fprintf('Precision: %.1f → %.1f nm\n', input_prec, output_prec);
```

#### Pattern 2: Hierarchical Bayes (Learning Xi)

```matlab
% When you don't know Xi ahead of time
B = smi.BaGoL();
B.SMD = SMD;
B.SMD.PixelSize = pixelSize;

% Hierarchical settings
B.HierarchFlag = 1;      % Learn Xi
B.Xi = 10;               % Initial guess
B.Alpha_Xi = 1;          % Hyper-prior shape
B.Beta_Xi = 50;          % Hyper-prior scale
B.NSamples = 10;         % Samples between Xi updates
B.ChainFlag = 1;         % Save chains
B.N_Trials = 5000;       % More samples for learning

B.analyze_all();

% Examine learned Xi
if length(B.XiChain(1,:)) == 1
    % Poisson
    final_Xi = mean(B.XiChain(end-500:end));
    fprintf('Learned Xi (Poisson): %.1f\n', final_Xi);
else
    % Gamma
    k = mean(B.XiChain(end-500:end, 1));
    theta = mean(B.XiChain(end-500:end, 2));
    fprintf('Learned Xi (Gamma): k=%.1f, theta=%.1f\n', k, theta);
end
```

#### Pattern 3: High-Density Data

```matlab
% For dense data with many emitters
B = smi.BaGoL();
B.SMD = SMD;
B.SMD.PixelSize = pixelSize;

% Smaller ROIs, more overlap
B.ROIsize = 150;        % Smaller ROIs
B.Overlap = 30;         % More overlap
B.Cutoff = 100;         % Tighter pre-clustering

% Adjust for density
B.Xi = 25;              % Expect more localizations per emitter
B.N_Trials = 4000;      % More samples

B.analyze_all();
```

### See Also

- [BaGoL Workflow](../workflows/bagol-clustering.md) - Complete guide
- [Clustering Example](../examples/clustering-analysis.md) - Hands-on example
- Publication: Fazel et al., Nature Communications 13:7152 (2022)

---

## smi.Publish

### Description

`smi.Publish` performs batch processing of super-resolution data organized in a standardized directory structure. It is designed for experiments with multiple cells and labels (e.g., sequential super-resolution imaging), automatically processing all matching datasets, applying consistent analysis parameters, and generating comparative statistics across cells and labels.

### Class Definition

```matlab
classdef Publish < handle
```

### Constructor

```matlab
obj = smi.Publish()
obj = smi.Publish(SMF, CoverslipDir)
```

**Parameters:**
- `SMF`: Single Molecule Fitting structure (optional)
- `CoverslipDir`: Top-level directory containing Cell*/Label*/Data*.h5 (optional)

**Example:**
```matlab
Pub = smi.Publish();
Pub.CoverslipDir = '/data/Experiment_2024-01-10';
Pub.SMF.Data.PixelSize = 0.108;
```

### Properties

#### Configuration

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `SMF` | struct | - | Single Molecule Fitting structure |
| `SMFperLabel` | cell array | {} | Per-label SMF overrides |
| `CoverslipDir` | string | '' | Top directory (Cell*/Label*/Data*.h5) |
| `SaveBaseDir` | string | CoverslipDir/Results | Base directory for results |
| `LogFilePath` | string | SaveBaseDir/Log_*.mat | Error log file path |

#### Analysis Control

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `LabelID` | array | [] | Specific labels to analyze ([] = all) |
| `CellList` | array | [] | Specific cells to analyze ([] = all) |
| `GenerateSR` | logical | true | Generate SR analysis results |
| `GenerateImagingStats` | logical | true | Generate imaging statistics |
| `GenerateOverlayStats` | logical | false | Generate multi-label overlay stats |
| `AnalyzeBleaching` | logical | false | Perform bleaching analysis |

#### Imaging Parameters

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `SRImageZoom` | numeric | 20 | SR image magnification |
| `SRCircleImageZoom` | numeric | 50 | Circle image magnification |

#### Drift Correction

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `UseBrightfieldDC` | logical | false | Apply brightfield drift correction |
| `MaxBrightfieldShift` | numeric | inf | Max shift for overlay masks (pixels) |
| `ShiftToReg` | logical | false | Shift localizations based on brightfield |
| `ShiftToRegViaDC` | logical | false | Shift via drift correction registration |

#### Verbosity

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `Verbose` | integer | 1 | Verbosity level |

#### Results

| Property | Type | Description |
|----------|------|-------------|
| `ResultsStruct` | struct | Collected analysis results across all cells/labels |

### Methods

#### Primary Workflow Methods

##### performFullAnalysis

```matlab
obj.performFullAnalysis()
```

Main analysis method that processes all cells and labels in `CoverslipDir`.

**Expected directory structure:**
```
CoverslipDir/
├── Cell_01/
│   ├── Label_01/
│   │   ├── Data_01.h5
│   │   ├── Data_02.h5
│   │   └── ...
│   ├── Label_02/
│   │   └── Data_01.h5
│   └── ...
├── Cell_02/
│   └── ...
└── Results/          (created automatically)
    ├── Cell_01/
    ├── Cell_02/
    ├── ResultsStruct.mat
    └── Log_*.mat
```

**Workflow:**
1. Creates results directory structure
2. Finds all Cell_* directories
3. For each cell:
   - Calls `processCell()`
4. Generates cross-cell statistics (if enabled)
5. Saves results and error log

**Example:**
```matlab
% Configure
Pub = smi.Publish();
Pub.CoverslipDir = '/data/Sequential_SR_Experiment/';
Pub.SMF.Data.PixelSize = 0.108;
Pub.SMF.Data.CameraGain = 2.5;
Pub.SMF.Data.CameraOffset = 100;

% Analysis parameters (apply to all cells/labels)
Pub.SMF.BoxFinding.BoxSize = 7;
Pub.SMF.Fitting.PSFSigma = 1.3;
Pub.SMF.Thresholding.MaxXY_SE = 0.15;

% Run
Pub.performFullAnalysis();
```

##### processCell

```matlab
obj.processCell(CellName)
```

Processes all labels within one cell directory.

**Parameters:**
- `CellName`: Name of cell directory (e.g., 'Cell_01')

**Workflow:**
1. Finds all Label_* directories in this cell
2. Loads brightfield images (if available)
3. For each label:
   - Calls `processLabel()`
4. Generates per-cell statistics
5. Generates overlay images (if multi-label)

**Called by:** `performFullAnalysis()`

##### processLabel

```matlab
obj.processLabel(CellName, LabelName)
```

Processes one label (one SR dataset) within a cell.

**Parameters:**
- `CellName`: Cell directory name
- `LabelName`: Label directory name (e.g., 'Label_01')

**Workflow:**
1. Finds all Data_*.h5 files
2. Configures SMF for this label
3. Runs SMLM analysis
4. Applies brightfield drift correction (if enabled)
5. Saves results
6. Generates images and plots

**Called by:** `processCell()`

#### Results and Visualization

##### saveResults

```matlab
obj.saveResults()
```

Saves the `ResultsStruct` to file.

**Saves:**
- `ResultsStruct.mat` in `SaveBaseDir`

##### genOverlayResults

```matlab
obj.genOverlayResults()
```

Generates multi-label overlay statistics and images across all cells.

**Called by:** `performFullAnalysis()` if `GenerateOverlayStats = true`

##### genAlignResults

```matlab
AlignResultsStruct = obj.genAlignResults(FilePath, SaveDir)
```

Generates alignment statistics from registration data.

**Parameters:**
- `FilePath`: Path to data file
- `SaveDir`: Directory for saving results

**Returns:**
- `AlignResultsStruct`: Alignment statistics

### Static Methods

##### genSROverlays

```matlab
smi.Publish.genSROverlays(ResultsCellDir, SaveDir, AnalysisID, Mask, MaskName)
```

Generates multi-channel SR overlay images.

**Parameters:**
- `ResultsCellDir`: Cell results directory
- `SaveDir`: Save directory
- `AnalysisID`: Analysis identifier
- `Mask`: Binary mask to apply
- `MaskName`: Mask name for filenames

##### overlayNImages

```matlab
[OverlayImage, ColorOrderTag] = smi.Publish.overlayNImages(ImageStack)
```

Creates RGB overlay from N grayscale images.

**Parameters:**
- `ImageStack`: Y × X × N array of images

**Returns:**
- `OverlayImage`: RGB overlay
- `ColorOrderTag`: Color assignment description

##### genAlignStats

```matlab
StatsStruct = smi.Publish.genAlignStats(AlignRegStruct, SMD, SaveDir)
```

Generates alignment statistics plots.

**Parameters:**
- `AlignRegStruct`: Alignment registration structure
- `SMD`: Single Molecule Data
- `SaveDir`: Save directory

**Returns:**
- `StatsStruct`: Statistics structure

##### genAlignXCorr

```matlab
XCorrStruct = smi.Publish.genAlignXCorr(AlignRegStruct, SaveDir)
```

Generates cross-correlation analysis of alignment.

**Parameters:**
- `AlignRegStruct`: Alignment structure
- `SaveDir`: Save directory

**Returns:**
- `XCorrStruct`: Cross-correlation results

##### computeRegCorrection

```matlab
RegCorrection = smi.Publish.computeRegCorrection(SMF)
```

Computes registration correction from SMF parameters.

**Parameters:**
- `SMF`: Single Molecule Fitting structure

**Returns:**
- `RegCorrection`: Registration correction parameters

##### genBFMask

```matlab
[Mask, Shifts, ImageROIs] = smi.Publish.genBFMask(FocusImageStruct, RefImage, MaxBrightfieldShift)
```

Generates mask based on brightfield image shifts.

**Parameters:**
- `FocusImageStruct`: Brightfield image structure
- `RefImage`: Reference image
- `MaxBrightfieldShift`: Maximum allowed shift (pixels)

**Returns:**
- `Mask`: Binary mask (1 = good, 0 = excessive shift)
- `Shifts`: Local shift magnitudes
- `ImageROIs`: ROI definitions

### Usage Patterns

#### Pattern 1: Standard Batch Processing

```matlab
% Setup
Pub = smi.Publish();
Pub.CoverslipDir = '/data/2024-01-10_Experiment/';

% Configure SMF (applied to all datasets)
Pub.SMF.Data.PixelSize = 0.108;
Pub.SMF.Data.CameraGain = 2.5;
Pub.SMF.Data.CameraOffset = 100;
Pub.SMF.BoxFinding.MinPhotons = 250;
Pub.SMF.Fitting.PSFSigma = 1.3;
Pub.SMF.FrameConnection.On = true;
Pub.SMF.DriftCorrection.On = true;

% Run analysis
Pub.performFullAnalysis();

% Results saved to: CoverslipDir/Results/
% - Cell_01/Label_01/Data_01_Results.mat
% - Cell_01/Label_02/Data_01_Results.mat
% - ...
% - ResultsStruct.mat
```

#### Pattern 2: Label-Specific Parameters

```matlab
% Base SMF
Pub = smi.Publish();
Pub.CoverslipDir = '/data/Sequential_Imaging/';
Pub.SMF.Data.PixelSize = 0.108;

% Label 1: DNA-PAINT (dim, dense)
SMF1 = smi_core.SingleMoleculeFitting();
SMF1.BoxFinding.MinPhotons = 150;
SMF1.Thresholding.MaxXY_SE = 0.20;

% Label 2: Antibody (bright, sparse)
SMF2 = smi_core.SingleMoleculeFitting();
SMF2.BoxFinding.MinPhotons = 500;
SMF2.Thresholding.MaxXY_SE = 0.10;

% Assign
Pub.SMFperLabel = {SMF1, SMF2};

% Common fields (must match across SMFs)
for ii = 1:2
    Pub.SMFperLabel{ii}.Data.FileDir = Pub.SMF.Data.FileDir;
    Pub.SMFperLabel{ii}.Data.PixelSize = 0.108;
    Pub.SMFperLabel{ii}.Data.CameraGain = 2.5;
end

Pub.performFullAnalysis();
```

#### Pattern 3: Selective Cell Processing

```matlab
% Process only specific cells
Pub = smi.Publish();
Pub.CoverslipDir = '/data/Experiment/';
Pub.CellList = [1, 3, 5];  % Only Cell_01, Cell_03, Cell_05
Pub.LabelID = [1];          % Only Label_01

Pub.performFullAnalysis();
% Processes only Cell_01/Label_01/, Cell_03/Label_01/, Cell_05/Label_01/
```

#### Pattern 4: With Brightfield Drift Correction

```matlab
% Use brightfield images for drift correction
Pub = smi.Publish();
Pub.CoverslipDir = '/data/Sequential_with_BF/';
Pub.UseBrightfieldDC = true;
Pub.MaxBrightfieldShift = 5;  % pixels

% Also perform standard drift correction on SR data
Pub.SMF.DriftCorrection.On = true;

Pub.performFullAnalysis();
```

### Directory Structure Requirements

Publish expects this structure:

```
CoverslipDir/
├── Cell_01/
│   ├── Label_01/
│   │   ├── Data_01.h5         (required: raw SR data)
│   │   ├── Data_02.h5
│   │   ├── PreSeqBF.mat       (optional: brightfield before)
│   │   └── PostSeqBF.mat      (optional: brightfield after)
│   └── Label_02/
│       └── Data_01.h5
├── Cell_02/
│   └── ...
└── Results/                    (created by Publish)
    ├── Cell_01/
    │   ├── Label_01/
    │   │   ├── Data_01_Results/
    │   │   │   ├── Data_01_Results.mat
    │   │   │   ├── Data_01_GaussImage.png
    │   │   │   └── ... diagnostic plots ...
    │   │   └── ...
    │   └── Label_02/
    │       └── ...
    ├── ResultsStruct.mat
    └── Log_*.mat
```

### See Also

- [Batch Processing Workflow](../workflows/batch-processing.md) - Complete guide
- `smi.SMLM` - SMLM analysis
- `smi_core.DriftCorrection` - Drift correction methods

---

## Comparison Summary

| Feature | SMLM | SPT | BaGoL | Publish |
|---------|------|-----|-------|---------|
| **Primary Use** | Localize blinking emitters | Track moving particles | Cluster localizations | Batch multi-cell processing |
| **Input** | Raw images | Raw images | SMD (localizations) | Directory of raw images |
| **Output** | SMD | TR (trajectories) | MAPN (refined positions) | Multiple SMD/Results |
| **Frame Connection** | Optional (same emitter) | Required (trajectory) | Required (pre-analysis) | Configurable |
| **Drift Correction** | Optional | No (tracking inherent) | Handled in algorithm | Optional + brightfield |
| **GUI Available** | Yes | Yes | No | No |
| **Typical Use Case** | dSTORM, PAINT | Particle tracking | PAINT, dSTORM | Multi-cell experiments |
| **Complexity** | Moderate | High | High | Moderate |

## Common Workflows

### SMLM → BaGoL

```matlab
% 1. Localize with SMLM
SMF = smi_core.SingleMoleculeFitting();
% ... configure SMF ...
SMF.FrameConnection.On = true;  % Important for BaGoL
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();
SMD = SMLMobj.SMD;

% 2. Refine with BaGoL
B = smi.BaGoL();
B.SMD = SMD;
B.Xi = 20;
B.analyze_all();
MAPN = B.MAPN;  % Refined positions

% Compare precision
fprintf('SMLM precision: %.1f nm\n', median(SMD.X_SE)*SMD.PixelSize*1000);
fprintf('BaGoL precision: %.1f nm\n', median(MAPN.X_SE));
```

### Batch SMLM with Publish

```matlab
% Process many cells with consistent parameters
Pub = smi.Publish();
Pub.CoverslipDir = '/data/Experiment/';
Pub.SMF = configuredSMF;
Pub.performFullAnalysis();

% Access aggregated results
load(fullfile(Pub.SaveBaseDir, 'ResultsStruct.mat'));
```

### Tracking with Diffusion Analysis

```matlab
% Track and analyze diffusion
SMF = smi_core.SingleMoleculeFitting();
SMF.Tracking.TrajwiseD = true;
SPTobj = smi.SPT(SMF, false);
[TR, SMD] = SPTobj.performFullAnalysis();

% Analyze diffusion
D_vals = arrayfun(@(x) x.D, TR);
figure; histogram(D_vals);
```

## See Also

- [Architecture Overview](../core-concepts/architecture.md) - Namespace organization
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete SMLM guide
- [SPT Workflow](../workflows/spt-tracking.md) - Tracking guide
- [BaGoL Workflow](../workflows/bagol-clustering.md) - Clustering guide
- [Batch Processing](../workflows/batch-processing.md) - Publish workflow
