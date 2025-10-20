---
title: "API Reference: +smi_core Namespace"
category: "api-reference"
level: "intermediate"
tags: ["api", "smi_core", "core-classes", "localization", "fitting", "data-structures"]
prerequisites: ["../core-concepts/architecture.md", "../core-concepts/smf-structure.md", "../core-concepts/smd-structure.md"]
related: ["../workflows/smlm-analysis.md", "../workflows/spt-tracking.md", "../how-to/localize-molecules.md"]
summary: "Complete API reference for the +smi_core namespace covering core processing classes, data structures, and fundamental algorithms"
estimated_time: "30 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# API Reference: +smi_core Namespace

## Purpose

The +smi_core namespace contains the foundational classes and algorithms that power all smite analyses. This reference documents the core processing pipeline components, data structure definitions, and essential algorithms for single molecule localization and tracking. Understanding these classes is crucial for customizing analyses, debugging pipelines, and extending smite functionality.

## Prerequisites

- Understanding of [smite architecture](../core-concepts/architecture.md)
- Familiarity with [SMF structure](../core-concepts/smf-structure.md)
- Knowledge of [SMD structure](../core-concepts/smd-structure.md)
- Basic MATLAB class usage

## Overview

The +smi_core namespace provides:

**Data Structures:**
- `SingleMoleculeFitting` - Complete analysis parameter specification (SMF)
- `SingleMoleculeData` - Localization results storage (SMD)
- `TrackingResults` - Trajectory-organized results (TR)

**Core Processing:**
- `LoadData` - File I/O for raw microscopy data
- `DataToPhotons` - Camera calibration and ADU conversion
- `LocalizeData` - End-to-end localization pipeline
- `FindROI` - Candidate molecule detection
- `Threshold` - Quality filtering of localizations

**Post-Processing:**
- `FrameConnection` - Temporal clustering of repeated localizations
- `DriftCorrection` - Stage drift correction (intra/inter-dataset)
- `ChannelRegistration` - Multi-channel alignment

**Algorithms:**
- `GaussMLE` - GPU-accelerated maximum likelihood fitting
- `FRC` - Fourier Ring Correlation resolution estimation

## Class Reference

### Data Structure Classes

#### SingleMoleculeFitting

**Purpose:** Defines and manages the SMF structure containing all analysis parameters.

**Key Concept:** The SMF structure is a "structure of structures" that completely specifies an analysis pipeline. Same SMF + same data = reproducible results.

**Class Definition:**
```matlab
classdef SingleMoleculeFitting < matlab.mixin.Copyable
```

**Properties:**
```matlab
Data            % File paths, camera settings, acquisition parameters
BoxFinding      % Molecule detection parameters
Fitting         % PSF model and fitting configuration
Thresholding    % Quality filtering thresholds
FrameConnection % Temporal clustering parameters
DriftCorrection % Drift correction settings
Tracking        % Particle tracking parameters
```

**Constructor:**
```matlab
SMF = smi_core.SingleMoleculeFitting()
```

Creates SMF with all default values. Returns a structure (not an object) organized as sub-structures.

**Key Methods:**

`gui()` - Interactive parameter editor
```matlab
SMF = smi_core.SingleMoleculeFitting();
SMF.gui();  % Opens graphical interface for parameter editing
```

`padSMF(SMF)` - Static method to ensure SMF has all required fields
```matlab
% Useful when loading old SMF structures
SMF = smi_core.SingleMoleculeFitting.padSMF(OldSMF);
```

`revertToDefaults()` - Reset to default values
```matlab
SMF = SMF.revertToDefaults();
```

**Usage Examples:**

Basic creation and modification:
```matlab
% Create with defaults
SMF = smi_core.SingleMoleculeFitting();

% Set file information
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'experiment1.h5'};
SMF.Data.PixelSize = 0.108;  % micrometers
SMF.Data.FrameRate = 100;    % Hz

% Configure fitting
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Set thresholds
SMF.Thresholding.MaxXY_SE = 0.2;  % pixels
SMF.Thresholding.MinPhotons = 100;

% Save configuration
save('analysis_params.mat', 'SMF');
```

**See Also:**
- [SMF Structure Guide](../core-concepts/smf-structure.md) - Complete field reference
- [SMLM Workflow](../workflows/smlm-analysis.md) - SMF in context

---

#### SingleMoleculeData

**Purpose:** Defines and manages the SMD structure containing localization results.

**Key Concept:** SMD is the primary results structure in smite. Each row represents one localization with position, photons, fit quality, frame number, and metadata.

**Class Definition:**
```matlab
classdef SingleMoleculeData
```

**Structure Fields:**

Scalar metadata:
```matlab
SMD.NFrames      % Total frames in raw data
SMD.NDatasets    % Number of datasets analyzed
SMD.FrameRate    % Acquisition rate (Hz)
SMD.PixelSize    % Pixel size (micrometers)
SMD.XSize        % Image width (pixels)
SMD.YSize        % Image height (pixels)
```

Per-localization arrays (N localizations):
```matlab
SMD.X            % X position (pixels)
SMD.Y            % Y position (pixels)
SMD.Z            % Z position (pixels, if 3D)
SMD.X_SE         % Standard error in X (pixels)
SMD.Y_SE         % Standard error in Y (pixels)
SMD.Z_SE         % Standard error in Z (pixels)
SMD.Photons      % Integrated photons
SMD.Photons_SE   % Standard error in photons
SMD.Bg           % Background (photons/pixel)
SMD.Bg_SE        % Standard error in background
SMD.PSFSigma     % PSF width (pixels, symmetric)
SMD.PSFSigmaX    % PSF width X (pixels, asymmetric)
SMD.PSFSigmaY    % PSF width Y (pixels, asymmetric)
SMD.FrameNum     % Frame number (1-indexed)
SMD.DatasetNum   % Dataset number (1-indexed)
SMD.PValue       % Fit p-value
SMD.LogLikelihood % Fit log-likelihood
SMD.ThreshFlag   % Thresholding flags (0 = passed)
SMD.ConnectID    % Frame-connection identifier
```

Drift correction (NFrames × NDatasets):
```matlab
SMD.DriftX       % X drift per frame (pixels)
SMD.DriftY       % Y drift per frame (pixels)
SMD.DriftZ       % Z drift per frame (pixels)
```

**Constructor:**
```matlab
SMD = smi_core.SingleMoleculeData.createSMD()
```

Creates empty SMD structure with all fields initialized to `[]`.

**Key Methods:**

`catSMD(SMD1, SMD2)` - Concatenate two SMD structures
```matlab
SMD = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2);
% Combines localizations from both structures
```

`isolateSubSMD(SMD, Indices)` - Extract subset of localizations
```matlab
% Get only first 1000 localizations
SMDSub = smi_core.SingleMoleculeData.isolateSubSMD(SMD, 1:1000);
```

`isolateSubROI(SMD, ROI)` - Extract localizations in spatial ROI
```matlab
% ROI = [XMin, XMax, YMin, YMax] in pixels
ROI = [10, 50, 20, 60];
SMDRegion = smi_core.SingleMoleculeData.isolateSubROI(SMD, ROI);
```

`maskSMD(SMD, Mask)` - Apply boolean mask to filter localizations
```matlab
% Keep only high photon localizations
Mask.Photons = [500, Inf];  % Min 500 photons
[SMDFiltered, SMDRemoved, KeepBool] = ...
    smi_core.SingleMoleculeData.maskSMD(SMD, Mask);
```

`extractDatasets(SMD, Datasets)` - Isolate specific datasets
```matlab
% Extract only datasets 1, 3, and 5
SMDSub = smi_core.SingleMoleculeData.extractDatasets(SMD, [1, 3, 5]);
```

`computeDensity(SMD)` - Compute localization density
```matlab
Density = smi_core.SingleMoleculeData.computeDensity(SMD);
% Returns localizations per square micrometer
```

**Usage Examples:**

Creating and populating SMD:
```matlab
% Create empty structure
SMD = smi_core.SingleMoleculeData.createSMD();

% Populate with data (normally done by LocalizeData)
SMD.X = [10.5; 20.3; 15.7];  % 3 localizations
SMD.Y = [12.1; 18.9; 25.4];
SMD.Photons = [500; 800; 600];
SMD.FrameNum = [1; 1; 2];
SMD.NFrames = 100;
SMD.PixelSize = 0.108;
```

Working with SMD results:
```matlab
% Load results
load('Results.mat', 'SMD');

% Basic statistics
fprintf('Found %d localizations\n', length(SMD.X));
fprintf('Median photons: %.1f\n', median(SMD.Photons));
fprintf('Mean precision: %.2f nm\n', mean(SMD.X_SE) * SMD.PixelSize * 1000);

% Filter by quality
GoodLocs = SMD.ThreshFlag == 0;
SMDGood = smi_core.SingleMoleculeData.isolateSubSMD(SMD, find(GoodLocs));

% Spatial filtering
ROI = [50, 150, 50, 150];  % Center region
SMDCenter = smi_core.SingleMoleculeData.isolateSubROI(SMD, ROI);

% Combine datasets
load('Results1.mat', 'SMD');
SMD1 = SMD;
load('Results2.mat', 'SMD');
SMD2 = SMD;
SMDCombined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2);
```

**See Also:**
- [SMD Structure Guide](../core-concepts/smd-structure.md) - Complete field reference
- [How to Save and Load Results](../how-to/save-and-load-smd.md)

---

#### TrackingResults

**Purpose:** Trajectory-organized variant of SMD for single particle tracking data.

**Key Concept:** TR is an array of SMD structures where each element represents one trajectory. TR(n) contains all localizations from trajectory n organized with the same fields as SMD.

**Class Definition:**
```matlab
classdef TrackingResults
```

**Structure Organization:**
```matlab
% TR is an array: TR(1), TR(2), ..., TR(NTrajectories)
% Each element has SMD fields:
TR(n).X          % X positions for trajectory n
TR(n).Y          % Y positions for trajectory n
TR(n).FrameNum   % Frame numbers
TR(n).ConnectID  % All have same ConnectID value
% ... all other SMD fields
```

**Constructor:**
```matlab
TR = smi_core.TrackingResults.createTR()
```

Creates empty TR structure (equivalent to empty SMD).

**Key Methods:**

`convertSMDToTR(SMD)` - Convert SMD to trajectory format
```matlab
% After tracking, convert to TR for analysis
TR = smi_core.TrackingResults.convertSMDToTR(SMD);
% Each unique ConnectID becomes one trajectory
```

`convertTRToSMD(TR)` - Convert TR back to SMD
```matlab
SMD = smi_core.TrackingResults.convertTRToSMD(TR);
% Concatenates all trajectories into single SMD
```

`computeTrajLengths(TR)` - Get number of localizations per trajectory
```matlab
Lengths = smi_core.TrackingResults.computeTrajLengths(TR);
% Returns array of trajectory lengths
```

`computeTrajDurations(TR)` - Get temporal duration of trajectories
```matlab
Durations = smi_core.TrackingResults.computeTrajDurations(TR);
% Returns duration in seconds
```

`computeTrajFidelity(TR)` - Compute trajectory fidelity metric
```matlab
Fidelity = smi_core.TrackingResults.computeTrajFidelity(TR);
% Ratio of observed to expected localizations
```

`threshTrajLength(TR, MinLength)` - Filter by minimum trajectory length
```matlab
% Keep only trajectories with ≥5 localizations
TRFiltered = smi_core.TrackingResults.threshTrajLength(TR, 5);
```

`joinTraj(TR, IDs)` - Manually join trajectories
```matlab
% Join trajectories 10 and 15
TR = smi_core.TrackingResults.joinTraj(TR, [10, 15]);
```

`windowTR(TR, MinFrame, MaxFrame)` - Temporal windowing
```matlab
% Extract trajectories in frames 100-200
TRWindow = smi_core.TrackingResults.windowTR(TR, 100, 200);
```

**Usage Examples:**

Converting and analyzing trajectories:
```matlab
% Load tracking results
load('TrackingResults.mat', 'SMD');

% Convert to trajectory format
TR = smi_core.TrackingResults.convertSMDToTR(SMD);
fprintf('Found %d trajectories\n', length(TR));

% Analyze trajectory lengths
Lengths = smi_core.TrackingResults.computeTrajLengths(TR);
fprintf('Median trajectory length: %d localizations\n', median(Lengths));

% Filter short trajectories
TR = smi_core.TrackingResults.threshTrajLength(TR, 10);

% Analyze individual trajectory
traj_id = 5;
figure;
plot(TR(traj_id).X, TR(traj_id).Y, 'o-');
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('Trajectory %d (%d localizations)', ...
    traj_id, length(TR(traj_id).X)));
```

Working with trajectory properties:
```matlab
% Calculate mean squared displacement for each trajectory
NTraj = length(TR);
MSD = zeros(NTraj, 1);
for ii = 1:NTraj
    dx = diff(TR(ii).X);
    dy = diff(TR(ii).Y);
    MSD(ii) = mean(dx.^2 + dy.^2);
end

% Convert to physical units
PixelSize = TR(1).PixelSize;  % micrometers
MSD_um2 = MSD * PixelSize^2;

% Estimate diffusion coefficients
FrameRate = TR(1).FrameRate;
dt = 1 / FrameRate;  % seconds
D = MSD_um2 / (4 * dt);  % micrometers^2/second

fprintf('Median diffusion: %.3f um^2/s\n', median(D));
```

**See Also:**
- [TR Structure Guide](../core-concepts/tr-structure.md)
- [SPT Workflow](../workflows/spt-tracking.md)

---

### Core Processing Classes

#### LoadData

**Purpose:** Loads raw microscopy data from .h5 and .mat files.

**Key Concept:** Handles file I/O, dataset selection, and ROI extraction. Supports both single files and multi-file analyses.

**Class Definition:**
```matlab
classdef LoadData < handle
```

**Properties:**
```matlab
FileType        % 'mat' or 'h5' (auto-detected)
FullFileName    % Complete file path
```

**Constructor:**
```matlab
LD = smi_core.LoadData()
```

**Key Methods:**

`loadRawData(SMF, varargin)` - Main loading method
```matlab
% For .mat files
LD = smi_core.LoadData();
[~, Data, SMF] = LD.loadRawData(SMF, DatasetNum);

% For .h5 files
[~, Data, SMF] = LD.loadRawData(SMF, DatasetIdx, ChannelIdx);
```

Static methods:

`loadDataMat(FilePath, DatasetNum, VarName)` - Load from .mat
```matlab
Data = smi_core.LoadData.loadDataMat(FilePath, 1, 'sequence');
```

`loadDataH5(FilePath, DatasetIdx, ChannelIdx)` - Load from .h5
```matlab
Data = smi_core.LoadData.loadDataH5(FilePath, 1, 1);
```

`seqH5Data(FilePath)` - Query .h5 file structure
```matlab
[Channel, Ver, NDatasets, TopList, InstrList] = ...
    smi_core.LoadData.seqH5Data(FilePath);
% Inspect .h5 file contents before loading
```

**Usage Examples:**

Basic data loading:
```matlab
% Setup SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'experiment.h5'};

% Load first dataset
LD = smi_core.LoadData();
[~, Data, SMF] = LD.loadRawData(SMF, 1);

fprintf('Loaded: %d × %d × %d (Y × X × Frames)\n', size(Data));
```

Inspecting .h5 files:
```matlab
FilePath = '/path/to/data/experiment.h5';
[Channel, Version, NDatasets] = smi_core.LoadData.seqH5Data(FilePath);
fprintf('Found %d datasets in file\n', NDatasets);
```

Loading with ROI:
```matlab
% Set spatial ROI before loading
SMF.Data.DataROI = [50, 100, 50, 100, 1, 1];
% [YStart, XStart, YEnd, XEnd, FrameStart, FrameStep]

LD = smi_core.LoadData();
[~, Data, SMF] = LD.loadRawData(SMF, 1);
% Data is cropped to specified ROI
```

**See Also:**
- [How to Load Data](../how-to/load-data.md)

---

#### DataToPhotons

**Purpose:** Converts raw camera data from ADU (Analog-to-Digital Units) to photons using gain and offset calibration.

**Key Concept:** Proper calibration is essential for accurate photon counting and optimal fitting. Supports both EMCCD (scalar calibration) and sCMOS (per-pixel calibration).

**Class Definition:**
```matlab
classdef DataToPhotons < handle
```

**Properties:**
```matlab
RawData              % Input data in ADU
CameraGain           % Gain (ADU/photon), scalar or image
CameraOffset         % Offset (ADU), scalar or image
CameraReadNoise      % Read noise (ADU), scalar or image
CorrectedData        % Output data in photons
CorrectedReadNoise   % Output read noise in photons^2
```

**Constructor:**
```matlab
DTP = smi_core.DataToPhotons(SMF, RawData, RawDataROI, CalibrationROI)
```

**Key Methods:**

`convertData()` - Perform conversion
```matlab
DTP = smi_core.DataToPhotons(SMF, RawData);
[Data, ReadNoise] = DTP.convertData();
% Data now in photons, ReadNoise in photons^2
```

**Usage Examples:**

EMCCD conversion:
```matlab
% Load raw data
LD = smi_core.LoadData();
[~, RawData, SMF] = LD.loadRawData(SMF, 1);

% Set EMCCD calibration
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;    % ADU/photon
SMF.Data.CameraOffset = 100;  % ADU
SMF.Data.CameraReadNoise = 8; % ADU

% Convert to photons
DTP = smi_core.DataToPhotons(SMF, RawData);
[Data, ~] = DTP.convertData();

% Check conversion
fprintf('Raw range: [%d, %d] ADU\n', min(RawData(:)), max(RawData(:)));
fprintf('Photon range: [%.1f, %.1f] photons\n', min(Data(:)), max(Data(:)));
```

sCMOS conversion with calibration file:
```matlab
% Calibration file contains gain, offset, noise as images
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = '/path/to/calibration.mat';

% Conversion automatically uses per-pixel calibration
DTP = smi_core.DataToPhotons(SMF, RawData);
[Data, ReadNoiseVar] = DTP.convertData();
```

**See Also:**
- [How to Load Data](../how-to/load-data.md)
- [Camera Calibration Guide](../how-to/camera-calibration.md)

---

#### LocalizeData

**Purpose:** Complete localization pipeline from photon-converted images to thresholded localizations.

**Key Concept:** This is the workhorse class that orchestrates FindROI, GaussMLE fitting, and Threshold filtering. One-stop localization.

**Class Definition:**
```matlab
classdef LocalizeData < handle
```

**Properties:**
```matlab
SMF              % Analysis parameters
ScaledData       % Input data in photons
Verbose          % Verbosity level (0-3)
ResultsDir       % Output directory for diagnostics
SMD              % Thresholded results
SMDPreThresh     % Unthresholded results
```

**Constructor:**
```matlab
LD = smi_core.LocalizeData(ScaledData, SMF, Verbose, ResultsDir)
```

**Key Methods:**

`genLocalizations()` - Run complete localization pipeline
```matlab
LD = smi_core.LocalizeData(Data, SMF);
[SMD, SMDPreThresh] = LD.genLocalizations();
```

Or auto-run in constructor:
```matlab
[~, SMD, SMDPreThresh] = smi_core.LocalizeData(Data, SMF, 1);
```

**Pipeline Steps:**
1. FindROI: Detect candidate molecules
2. GaussMLE: Fit PSF model to each candidate
3. Threshold: Filter by quality metrics

**Usage Examples:**

Basic localization:
```matlab
% Load and convert data to photons
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data';
SMF.Data.FileName = {'cell1.h5'};
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;
SMF.Data.PixelSize = 0.108;
SMF.Data.FrameRate = 100;

% Load
LD = smi_core.LoadData();
[~, RawData, SMF] = LD.loadRawData(SMF, 1);

% Convert to photons
DTP = smi_core.DataToPhotons(SMF, RawData);
[Data, ~] = DTP.convertData();

% Localize
LD = smi_core.LocalizeData(Data, SMF);
[SMD, SMDPreThresh] = LD.genLocalizations();

fprintf('Found %d localizations\n', length(SMD.X));
fprintf('Rejected %d poor fits\n', ...
    length(SMDPreThresh.X) - length(SMD.X));
```

With verbosity for diagnostics:
```matlab
% Verbose level 2 shows progress
LD = smi_core.LocalizeData(Data, SMF, 2);
[SMD, SMDPreThresh] = LD.genLocalizations();

% Verbose level 3 saves diagnostic images and videos
LD = smi_core.LocalizeData(Data, SMF, 3, '/results/diagnostics');
[SMD, SMDPreThresh] = LD.genLocalizations();
```

Parameter tuning workflow:
```matlab
% Try different fitting parameters
SMF.Fitting.PSFSigma = 1.3;
LD1 = smi_core.LocalizeData(Data, SMF);
[SMD1, ~] = LD1.genLocalizations();

SMF.Fitting.FitType = 'XYNBS';  % Let sigma vary
LD2 = smi_core.LocalizeData(Data, SMF);
[SMD2, ~] = LD2.genLocalizations();

fprintf('Fixed sigma: %d localizations\n', length(SMD1.X));
fprintf('Variable sigma: %d localizations\n', length(SMD2.X));
```

**See Also:**
- [How to Localize Molecules](../how-to/localize-molecules.md)
- [SMLM Workflow](../workflows/smlm-analysis.md)

---

#### FindROI

**Purpose:** Detects candidate molecule locations using difference-of-Gaussians filtering and finds local maxima.

**Key Concept:** GPU-accelerated detection using CUDA. Fast and sensitive detection is crucial for not missing dim emitters while avoiding false positives.

**Class Definition:**
```matlab
classdef FindROI < handle
```

**Properties:**
```matlab
Data            % Input image stack
BoxSize         % Fitting box size (pixels)
BoxOverlap      % Allowed overlap (pixels)
MinPhotons      % Detection threshold (photons)
PSFSigma        % PSF width for filtering (pixels)
ROIs            % Output ROI stack
LocalMaxIm      % Binary local maxima image
```

**Constructor:**
```matlab
FR = smi_core.FindROI(SMF, Data)
```

**Key Methods:**

`findROI(SMD)` - Detect molecules and extract ROIs
```matlab
FR = smi_core.FindROI(SMF, Data);
[ROIStack, SMD] = FR.findROI(SMD);
% ROIStack: 4D array (BoxSize × BoxSize × NFrames × NBoxes)
% SMD: Contains box positions
```

**Algorithm:**
1. Apply difference-of-Gaussians filter (band-pass)
2. Find local maxima
3. Threshold by estimated photon count
4. Extract boxes with allowed overlap

**Usage Examples:**

Standalone ROI finding:
```matlab
% Prepare data
SMF = smi_core.SingleMoleculeFitting();
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 200;
SMF.Fitting.PSFSigma = 1.3;

% Find ROIs
FR = smi_core.FindROI(SMF, Data);
SMD = smi_core.SingleMoleculeData.createSMD();
[ROIStack, SMD] = FR.findROI(SMD);

fprintf('Found %d candidates\n', size(ROIStack, 4));
```

Visualizing detections:
```matlab
% Find ROIs
FR = smi_core.FindROI(SMF, Data);
FR.PlotBoxFrame = 10;  % Visualize frame 10
FR.Verbose = 3;
[ROIStack, SMD] = FR.findROI(SMD);
% Saves visualization of detected boxes
```

Tuning detection sensitivity:
```matlab
% Conservative (fewer false positives)
SMF.BoxFinding.MinPhotons = 500;
FR = smi_core.FindROI(SMF, Data);
[ROIStack1, SMD1] = FR.findROI(SMD);

% Aggressive (more sensitivity)
SMF.BoxFinding.MinPhotons = 100;
FR = smi_core.FindROI(SMF, Data);
[ROIStack2, SMD2] = FR.findROI(SMD);

fprintf('Conservative: %d ROIs\n', size(ROIStack1, 4));
fprintf('Aggressive: %d ROIs\n', size(ROIStack2, 4));
```

**See Also:**
- [How to Localize Molecules](../how-to/localize-molecules.md)

---

#### Threshold

**Purpose:** Filters localizations based on quality metrics and spatial criteria.

**Key Concept:** Thresholding removes poor fits, outliers, and artifacts. The ThreshFlag field in SMD encodes which thresholds each localization failed.

**Class Definition:**
```matlab
classdef Threshold < handle
```

**Properties:**
```matlab
Verbose          % Verbosity level
```

**Constant Properties:**
```matlab
Fields           % Available thresholding fields
```

**Key Methods:**

`setThreshFlag(SMD, MinMax)` - Set threshold flags
```matlab
T = smi_core.Threshold();
MinMax.X_SE = [0, 0.2];      % Max precision 0.2 pixels
MinMax.Photons = [100, Inf]; % Min 100 photons
[SMD, TFlag] = T.setThreshFlag(SMD, MinMax);
% ThreshFlag bitfield indicates which thresholds failed
```

`applyThresh(SMD)` - Static method to apply thresholds
```matlab
SMD = smi_core.Threshold.applyThresh(SMD);
% Removes localizations where ThreshFlag ~= 0
```

`setMinMax(SMF)` - Static method to convert SMF to MinMax structure
```matlab
MinMax = smi_core.Threshold.setMinMax(SMF);
% Extracts thresholds from SMF.Thresholding
```

`translateThreshFlag(ThreshFlag)` - Decode threshold flags
```matlab
[Readable, HotBits] = smi_core.Threshold.translateThreshFlag(ThreshFlag);
% Readable: cell array of reason strings
% HotBits: which thresholds were violated
```

`rejectedLocalizations(SMD, Options, SaveDir)` - Visualize rejections
```matlab
T = smi_core.Threshold();
T.rejectedLocalizations(SMD, struct(), SaveDir);
% Creates plots showing which localizations failed which thresholds
```

**Usage Examples:**

Manual thresholding:
```matlab
% Load unthresholded data
load('Results.mat', 'SMDPreThresh');

% Define thresholds
MinMax.X_SE = [0, 0.2];
MinMax.Y_SE = [0, 0.2];
MinMax.Photons = [200, Inf];
MinMax.PValue = [0.01, 1];

% Apply thresholds
T = smi_core.Threshold();
[SMD, TFlag] = T.setThreshFlag(SMDPreThresh, MinMax);
SMD = smi_core.Threshold.applyThresh(SMD);

fprintf('Kept %d / %d localizations\n', ...
    length(SMD.X), length(SMDPreThresh.X));
```

Understanding rejections:
```matlab
% Which localizations failed thresholds?
Failed = SMDPreThresh.ThreshFlag ~= 0;
fprintf('%d failed thresholds\n', sum(Failed));

% Why did specific localization fail?
loc_idx = 100;
[Reasons, ~] = smi_core.Threshold.translateThreshFlag(...
    SMDPreThresh.ThreshFlag(loc_idx));
fprintf('Localization %d failed: %s\n', loc_idx, Reasons{1});
```

Visualizing thresholding results:
```matlab
% Create diagnostic plots
T = smi_core.Threshold();
T.Verbose = 2;
Options = struct();
T.rejectedLocalizations(SMDPreThresh, Options, '/results/threshold_plots');
% Saves histograms showing rejected vs accepted localizations
```

**See Also:**
- [How to Threshold Results](../how-to/threshold-results.md)
- [SMF Thresholding Fields](../core-concepts/smf-structure.md#thresholding-quality-filtering)

---

### Post-Processing Classes

#### FrameConnection

**Purpose:** Combines localizations from the same molecule across multiple frames to improve precision.

**Key Concept:** When molecules blink or photoactivate repeatedly, connecting temporal clusters yields better precision than any single localization. Uses Linear Assignment Problem (LAP) algorithm.

**Class Definition:**
```matlab
classdef FrameConnection < handle
```

**Properties:**
```matlab
SMF              % Analysis parameters
SMD              % Input SMD
SMDCombined      % Frame-connected output
Verbose          % Verbosity level
```

**Constructor:**
```matlab
FC = smi_core.FrameConnection(SMD, SMF, Verbose)
```

**Key Methods:**

`performFrameConnection()` - Main frame connection method
```matlab
FC = smi_core.FrameConnection(SMD, SMF);
[SMDCombined, SMD] = FC.performFrameConnection();
% SMDCombined: Combined localizations
% SMD: Updated with ConnectID field
```

Or auto-run:
```matlab
[~, SMDCombined, SMD] = smi_core.FrameConnection(SMD, SMF, 1, 1);
```

Static methods:

`lapFC(SMD, SMF)` - LAP-based frame connection
```matlab
[SMD, Params] = smi_core.FrameConnection.lapFC(SMD, SMF);
```

`combineLocalizations(SMD, SMF)` - Combine connected localizations
```matlab
SMDCombined = smi_core.FrameConnection.combineLocalizations(SMD, SMF);
```

`findConnected(SMR, SMD, ID)` - Find localizations with specific ID
```matlab
Indices = smi_core.FrameConnection.findConnected(SMD, SMD, 42);
% Returns indices of all localizations with ConnectID = 42
```

**Algorithm Overview:**
1. Cluster nearby localizations spatially
2. Estimate blinking parameters (KOn, KOff, bleaching)
3. Create cost matrix for temporal connections
4. Solve LAP to assign connections
5. Combine connected localizations with precision weighting

**Usage Examples:**

Basic frame connection:
```matlab
% After localization
load('LocalizedData.mat', 'SMD', 'SMF');

% Enable frame connection
SMF.FrameConnection.On = true;
SMF.FrameConnection.Method = 'LAP-FC';
SMF.FrameConnection.MaxSeparation = 1.0;  % pixels
SMF.FrameConnection.MaxFrameGap = 5;

% Perform connection
FC = smi_core.FrameConnection(SMD, SMF);
[SMDCombined, SMD] = FC.performFrameConnection();

fprintf('Original: %d localizations\n', length(SMD.X));
fprintf('Combined: %d unique emitters\n', length(SMDCombined.X));

% Check precision improvement
MeanPrecBefore = mean(sqrt(SMD.X_SE.^2 + SMD.Y_SE.^2));
MeanPrecAfter = mean(sqrt(SMDCombined.X_SE.^2 + SMDCombined.Y_SE.^2));
fprintf('Precision improved by %.1f%%\n', ...
    100 * (1 - MeanPrecAfter/MeanPrecBefore));
```

Analyzing frame connections:
```matlab
% How many localizations per emitter?
UniqueIDs = unique(SMD.ConnectID);
LocsPerEmitter = histcounts(SMD.ConnectID, max(UniqueIDs));

figure;
histogram(LocsPerEmitter);
xlabel('Localizations per emitter');
ylabel('Count');
title(sprintf('Median: %d localizations/emitter', median(LocsPerEmitter)));
```

Tuning frame connection parameters:
```matlab
% Test different max frame gaps
Gaps = [1, 3, 5, 10];
NCombined = zeros(size(Gaps));

for ii = 1:length(Gaps)
    SMF.FrameConnection.MaxFrameGap = Gaps(ii);
    FC = smi_core.FrameConnection(SMD, SMF, 0);  % Quiet
    [SMDCombined, ~] = FC.performFrameConnection();
    NCombined(ii) = length(SMDCombined.X);
end

figure;
plot(Gaps, NCombined, 'o-');
xlabel('Max Frame Gap');
ylabel('Number of Combined Localizations');
```

**See Also:**
- [SMLM Workflow](../workflows/smlm-analysis.md#frame-connection)
- [SMF FrameConnection Fields](../core-concepts/smf-structure.md#frameconnection-temporal-clustering)

**Citation:**
David J. Schodt and Keith A. Lidke, "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", Frontiers in Bioinformatics, 2021. https://doi.org/10.3389/fbinf.2021.724325

---

#### DriftCorrection

**Purpose:** Corrects for stage drift within and between datasets using k-nearest neighbor (KNN) optimization.

**Key Concept:** Stage drift during long acquisitions causes spatial distortions. KNN drift correction estimates and removes drift without fiducial markers.

**Class Definition:**
```matlab
classdef DriftCorrection < handle
```

**Properties:**
```matlab
L_intra          % Intra-dataset threshold (pixels)
L_inter          % Inter-dataset threshold (pixels)
PixelSizeZUnit   % Pixel size for 3D DC (micrometers)
PDegree          % Polynomial degree for drift rate
BFRegistration   % Was brightfield registration done?
Verbose          % Verbosity level
SMF              % Analysis parameters
```

**Constructor:**
```matlab
DC = smi_core.DriftCorrection(SMF, SMD)
```

**Key Methods:**

`driftCorrectKNN(SMD)` - Complete drift correction
```matlab
DC = smi_core.DriftCorrection(SMF, SMD);
[SMD, Statistics] = DC.driftCorrectKNN(SMD);
% Performs both intra- and inter-dataset correction
```

`driftCorrectKNNIntra(SMD, cDataset, iDataset)` - Intra-dataset only
```matlab
[SMD, Stats] = DC.driftCorrectKNNIntra(SMD, 1, 1);
```

`driftCorrectKNNInter(SMD)` - Inter-dataset only
```matlab
[SMD, Stats] = DC.driftCorrectKNNInter(SMD);
```

`plotDriftCorrection(SMD, Option)` - Visualize drift
```matlab
FigHandle = DC.plotDriftCorrection(SMD, 'cumulative');
```

Static methods:

`driftCorrectBF(SMD, SMF, RefImage)` - Brightfield-based drift correction
```matlab
[SMD, BFStruct] = smi_core.DriftCorrection.driftCorrectBF(...
    SMD, SMF, RefImage);
```

`plotCumDrift(SMD, FieldName)` - Plot cumulative drift
```matlab
Fig = smi_core.DriftCorrection.plotCumDrift(SMD, 'DriftX');
```

**Algorithm Overview:**

Intra-dataset:
1. For each frame, find K nearest neighbors
2. Optimize drift trajectory to minimize neighbor distances
3. Fit polynomial to drift trajectory
4. Apply correction to coordinates

Inter-dataset:
1. Find correspondences between datasets
2. Estimate relative drift between datasets
3. Align all datasets to reference (first dataset)

**Usage Examples:**

Basic drift correction:
```matlab
% After localization and frame connection
load('Results.mat', 'SMD', 'SMF');

% Configure drift correction
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.Method = 'DC-KNN';
SMF.DriftCorrection.L_intra = 1.0;  % pixels
SMF.DriftCorrection.L_inter = 2.0;  % pixels

% Perform correction
DC = smi_core.DriftCorrection(SMF, SMD);
[SMD, Statistics] = DC.driftCorrectKNN(SMD);

% Visualize drift
DC.plotDriftCorrection(SMD, 'cumulative');
```

Multi-dataset drift correction:
```matlab
% Process each dataset separately first
SMDIntra = [];
for ii = 1:NDatasets
    % Extract dataset
    SMDi = smi_core.SingleMoleculeData.extractDatasets(SMD, ii);

    % Intra-dataset correction
    DC = smi_core.DriftCorrection(SMF);
    [SMDi, ~] = DC.driftCorrectKNNIntra(SMDi, ii, ii);

    % Accumulate
    SMDIntra = smi_core.SingleMoleculeData.catSMD(SMDIntra, SMDi);
end

% Inter-dataset correction
DC = smi_core.DriftCorrection(SMF);
[SMDFinal, Statistics] = DC.driftCorrectKNNInter(SMDIntra);
```

Analyzing drift:
```matlab
% Plot drift trajectory
figure;
subplot(2,1,1);
plot(SMD.DriftX(:,1));
xlabel('Frame'); ylabel('X Drift (pixels)');
title('X Drift Trajectory');

subplot(2,1,2);
plot(SMD.DriftY(:,1));
xlabel('Frame'); ylabel('Y Drift (pixels)');
title('Y Drift Trajectory');

% Compute total drift
TotalDriftX = range(SMD.DriftX(:,1));
TotalDriftY = range(SMD.DriftY(:,1));
fprintf('Total drift: %.1f × %.1f pixels\n', TotalDriftX, TotalDriftY);
```

**See Also:**
- [SMLM Workflow](../workflows/smlm-analysis.md#drift-correction)
- [SMF DriftCorrection Fields](../core-concepts/smf-structure.md#driftcorrection-stage-drift)

---

#### ChannelRegistration

**Purpose:** Computes and applies geometric transformations to align multi-channel data.

**Key Concept:** Multi-color imaging requires precise registration between channels. This class computes transformations from fiducial markers or images and applies them to localizations.

**Class Definition:**
```matlab
classdef ChannelRegistration < handle
```

**Properties:**
```matlab
SMF                      % Analysis parameters
Coordinates              % Fiducial coordinates per channel
FiducialImages           % Fiducial marker images
FiducialROI              % ROI definitions for channels
TransformationBasis      % 'coordinates' or 'images'
TransformationType       % 'lwm', 'affine', 'projective', etc.
RegistrationTransform    % Computed transform objects
```

**Constructor:**
```matlab
CR = smi_core.ChannelRegistration()
```

**Key Methods:**

`findTransform()` - Compute registration transform
```matlab
CR = smi_core.ChannelRegistration();
CR.SMF = SMF;
CR.TransformationType = 'lwm';  % Local weighted mean
CR.findTransform();
% CR.RegistrationTransform now contains transform
```

`transformSMD(SMD, ChannelNum)` - Apply transform to SMD
```matlab
% Transform channel 2 to channel 1 coordinates
SMDTransformed = CR.transformSMD(SMD, 2);
```

`saveTransform(FilePath)` - Save transform for reuse
```matlab
CR.saveTransform('/path/to/registration.mat');
```

**Usage Examples:**

Computing registration from fiducials:
```matlab
% Setup for fiducial-based registration
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data';
SMF.Data.FileName = {'fiducials_ch1.h5', 'fiducials_ch2.h5'};

% Create registration object
CR = smi_core.ChannelRegistration();
CR.SMF = SMF;
CR.TransformationBasis = 'coordinates';
CR.TransformationType = 'lwm';

% Compute transform
CR.findTransform();

% Save for later use
CR.saveTransform('/data/registration.mat');
```

Applying registration to data:
```matlab
% Load registration
load('/data/registration.mat', 'RegistrationTransform');

% Load channel 2 data
load('/data/channel2_results.mat', 'SMD');

% Transform to channel 1 coordinates
CR = smi_core.ChannelRegistration();
CR.RegistrationTransform = RegistrationTransform;
SMDTransformed = CR.transformSMD(SMD, 2);

% Check registration error
load('/data/channel1_results.mat', 'SMD');
SMD1 = SMD;
SMD2 = SMDTransformed;

% Find nearest neighbors between channels
Tree = KDTreeSearcher([SMD1.X, SMD1.Y]);
[Idx, Dist] = knnsearch(Tree, [SMD2.X, SMD2.Y]);

fprintf('Registration error: %.2f nm RMS\n', ...
    rms(Dist) * SMD1.PixelSize * 1000);
```

**See Also:**
- [Multi-Channel Analysis Guide](../how-to/multi-channel-analysis.md)

---

### Algorithm Classes

#### GaussMLE

**Purpose:** GPU-accelerated maximum likelihood estimation of 2D Gaussian PSF parameters.

**Key Concept:** Fast, precise fitting is the heart of SMLM. This CUDA implementation achieves theoretically minimal uncertainty using Newton-Raphson optimization.

**Class Definition:**
Function, not a class. Called via `smi_core.GaussMLE()`.

**Function Signature:**
```matlab
[P, CRLB, LogL, PSFX, PSFY] = smi_core.GaussMLE(Data, PSFSigma, ...
    SMFitType, Iterations, Varim)
```

**Inputs:**
- `Data`: ROI stack (BoxSize × BoxSize × NROIs) in photons
- `PSFSigma`: Initial PSF width(s) in pixels
- `SMFitType`: Fit model type (see below)
- `Iterations`: Newton-Raphson iterations (typically 20)
- `Varim`: Read noise variance image (sCMOS only)

**Outputs:**
- `P`: Fit parameters (NParams × NROIs)
- `CRLB`: Cramér-Rao lower bounds (precision estimates)
- `LogL`: Log-likelihood values
- `PSFX`, `PSFY`: Fit PSF widths

**Fit Types:**

| SMFitType | Parameters | Use Case |
|-----------|------------|----------|
| `'XYNB'` | X, Y, Photons, Bg | Standard 2D |
| `'XYNBS'` | X, Y, N, B, Sigma | Variable PSF |
| `'XYNBSXSY'` | X, Y, N, B, SigmaX, SigmaY | Asymmetric PSF |
| `'XYZNB'` | X, Y, Z, N, B | 3D astigmatism |

**Usage Examples:**

Direct fitting:
```matlab
% Create test data
ROI = zeros(7, 7, 100);
for ii = 1:100
    % Simulate Gaussian blob
    [X, Y] = meshgrid(1:7, 1:7);
    ROI(:,:,ii) = 1000 * exp(-((X-3.7).^2 + (Y-4.2).^2) / (2*1.3^2)) + 10;
end

% Fit
PSFSigma = 1.3;
SMFitType = 'XYNB';
Iterations = 20;
[P, CRLB, LogL] = smi_core.GaussMLE(ROI, PSFSigma, SMFitType, Iterations);

% Extract parameters
X = P(1, :)';
Y = P(2, :)';
Photons = P(3, :)';
Bg = P(4, :)';

% Precision estimates
X_SE = sqrt(CRLB(1, :)');
Y_SE = sqrt(CRLB(2, :)');

fprintf('Mean position: (%.2f, %.2f) pixels\n', mean(X), mean(Y));
fprintf('Mean precision: %.3f pixels\n', mean(X_SE));
```

Comparing fit types:
```matlab
% Fixed PSF
[P1, CRLB1, LogL1] = smi_core.GaussMLE(ROI, 1.3, 'XYNB', 20);

% Variable PSF
[P2, CRLB2, LogL2] = smi_core.GaussMLE(ROI, 1.3, 'XYNBS', 20);

% Compare log-likelihoods
fprintf('Fixed PSF LogL: %.1f\n', mean(LogL1));
fprintf('Variable PSF LogL: %.1f\n', mean(LogL2));
```

**Citation:**
Smith, C., Joseph, N., Rieger, B. et al. Fast, single-molecule localization that achieves theoretically minimum uncertainty. Nat Methods 7, 373–375 (2010). https://doi.org/10.1038/nmeth.1449

**See Also:**
- [How to Localize Molecules](../how-to/localize-molecules.md)
- [GPU Setup Guide](../how-to/use-gpu.md)

---

#### FRC

**Purpose:** Computes Fourier Ring Correlation to estimate average image resolution.

**Key Concept:** FRC provides an objective measure of achievable resolution by correlating two statistically independent reconstructions in Fourier space.

**Class Definition:**
```matlab
classdef FRC < handle
```

**Properties:**
```matlab
PixelSize        % Pixel size (nm)
SRImageZoom      % Super-resolution magnification
Repeats          % Number of averaging repeats
```

**Constructor:**
```matlab
FR = smi_core.FRC(SMF)
```

**Key Methods:**

`posToResolution(SMR)` - Compute resolution from localizations
```matlab
FR = smi_core.FRC(SMF);
[Resolution, FRCCurve] = FR.posToResolution(SMD);
% Resolution in nm
```

`qCorrectionLocs(SMR)` - Resolution with repeat correction
```matlab
[ResCorrected, ResUncorrected, Q, FRCCurve] = FR.qCorrectionLocs(SMD);
% Corrects for repeated localizations
```

**Usage Examples:**

Basic resolution calculation:
```matlab
% Load results
load('Results.mat', 'SMD', 'SMF');

% Compute FRC resolution
FR = smi_core.FRC(SMF);
FR.PixelSize = SMF.Data.PixelSize * 1000;  % Convert to nm
FR.SRImageZoom = 10;

[Resolution, FRCCurve] = FR.posToResolution(SMD);
fprintf('Estimated resolution: %.1f nm\n', Resolution);

% Plot FRC curve
figure;
plot(FRCCurve.SpatialFreq, FRCCurve.Correlation);
hold on;
yline(1/7, 'r--', 'Threshold');
xlabel('Spatial Frequency (1/nm)');
ylabel('FRC');
title(sprintf('Resolution: %.1f nm', Resolution));
```

With overcounting correction:
```matlab
FR = smi_core.FRC(SMF);
[ResCorrected, ResRaw, Q] = FR.qCorrectionLocs(SMD);

fprintf('Raw resolution: %.1f nm\n', ResRaw);
fprintf('Corrected resolution: %.1f nm\n', ResCorrected);
fprintf('Mean localizations per emitter: %.1f\n', Q);
```

**Requirements:**
- DIPlib Image Resolution add-on
- Curve Fitting Toolbox
- Parallel Processing Toolbox

**Citation:**
R.P.J. Nieuwenhuizen, K.A. Lidke, M. Bates, D. Leyton Puig, D. Grünwald, S. Stallinga, B. Rieger, "Measuring Image Resolution in Optical Nanoscopy", Nature Methods, 10(6):557-562, 2013.

**See Also:**
- [Resolution Analysis Guide](../how-to/resolution-analysis.md)

---

## Common Usage Patterns

### Complete SMLM Pipeline

End-to-end localization with all corrections:

```matlab
% 1. Setup parameters
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data';
SMF.Data.FileName = {'cell1.h5'};
SMF.Data.PixelSize = 0.108;
SMF.Data.FrameRate = 100;
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;
SMF.Fitting.PSFSigma = 1.3;
SMF.FrameConnection.On = true;
SMF.DriftCorrection.On = true;

% 2. Load data
LD = smi_core.LoadData();
[~, RawData, SMF] = LD.loadRawData(SMF, 1);

% 3. Convert to photons
DTP = smi_core.DataToPhotons(SMF, RawData);
[Data, ~] = DTP.convertData();

% 4. Localize
LD = smi_core.LocalizeData(Data, SMF);
[SMD, ~] = LD.genLocalizations();

% 5. Frame connection
FC = smi_core.FrameConnection(SMD, SMF);
[SMDCombined, SMD] = FC.performFrameConnection();

% 6. Drift correction
DC = smi_core.DriftCorrection(SMF, SMDCombined);
[SMDFinal, ~] = DC.driftCorrectKNN(SMDCombined);

% 7. Save results
save('Results.mat', 'SMDFinal', 'SMF', 'SMD');
```

### Batch Processing Multiple Datasets

Process multiple datasets with same parameters:

```matlab
% Setup
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data';
FileNames = {'cell1.h5', 'cell2.h5', 'cell3.h5'};

% Process each
AllSMD = [];
for ii = 1:length(FileNames)
    fprintf('Processing %s...\n', FileNames{ii});

    % Load
    SMF.Data.FileName = FileNames(ii);
    LD = smi_core.LoadData();
    [~, RawData, SMF] = LD.loadRawData(SMF, 1);

    % Convert
    DTP = smi_core.DataToPhotons(SMF, RawData);
    [Data, ~] = DTP.convertData();

    % Localize
    LD = smi_core.LocalizeData(Data, SMF, 0);  % Quiet
    [SMD, ~] = LD.genLocalizations();

    % Add dataset marker
    SMD.DatasetNum(:) = ii;

    % Accumulate
    AllSMD = smi_core.SingleMoleculeData.catSMD(AllSMD, SMD);
end

% Combined drift correction
DC = smi_core.DriftCorrection(SMF, AllSMD);
[AllSMD, ~] = DC.driftCorrectKNN(AllSMD);
```

### Custom Thresholding Workflow

Apply custom quality filters:

```matlab
% Load raw localizations
load('Results.mat', 'SMDPreThresh');

% Define strict thresholds
MinMax.X_SE = [0, 0.15];      % Tighter precision
MinMax.Y_SE = [0, 0.15];
MinMax.Photons = [300, Inf];  % Higher photon requirement
MinMax.PValue = [0.05, 1];    % Stricter p-value
MinMax.PSFSigma = [1.0, 1.6]; % Narrower PSF range

% Apply
T = smi_core.Threshold();
[SMD, ~] = T.setThreshFlag(SMDPreThresh, MinMax);
SMD = smi_core.Threshold.applyThresh(SMD);

fprintf('Strict filtering: kept %d / %d (%.1f%%)\n', ...
    length(SMD.X), length(SMDPreThresh.X), ...
    100 * length(SMD.X) / length(SMDPreThresh.X));

% Visualize rejections
T.rejectedLocalizations(SMDPreThresh, struct(), '/results/threshold_plots');
```

### Tracking Pipeline

Convert localizations to trajectories:

```matlab
% After localization with tracking parameters
SMF.Tracking.D = 0.1;  % Diffusion constant
SMF.Tracking.MaxDistFF = 5;
SMF.Tracking.MaxDistGC = 10;
SMF.Tracking.MinTrackLength = 5;

% Perform tracking (using smi.SPT, not detailed here)
% Results in SMD with ConnectID set

% Convert to trajectory format
TR = smi_core.TrackingResults.convertSMDToTR(SMD);

% Filter short trajectories
TR = smi_core.TrackingResults.threshTrajLength(TR, 10);

% Analyze
Lengths = smi_core.TrackingResults.computeTrajLengths(TR);
Durations = smi_core.TrackingResults.computeTrajDurations(TR);

fprintf('Found %d trajectories\n', length(TR));
fprintf('Median length: %d localizations\n', median(Lengths));
fprintf('Median duration: %.2f seconds\n', median(Durations));
```

---

## Performance Considerations

### GPU Acceleration

FindROI and GaussMLE require CUDA-capable GPU:

```matlab
% Check GPU
gpuDevice

% If multiple GPUs, select one
gpuDevice(1)

% Monitor GPU memory
g = gpuDevice;
fprintf('GPU memory: %.1f / %.1f GB used\n', ...
    (g.TotalMemory - g.AvailableMemory) / 1e9, ...
    g.TotalMemory / 1e9);
```

### Memory Management

For large datasets:

```matlab
% Process in chunks
NFramesTotal = size(Data, 3);
ChunkSize = 1000;
NChunks = ceil(NFramesTotal / ChunkSize);

SMDTotal = smi_core.SingleMoleculeData.createSMD();

for ii = 1:NChunks
    % Extract chunk
    StartFrame = (ii-1) * ChunkSize + 1;
    EndFrame = min(ii * ChunkSize, NFramesTotal);
    DataChunk = Data(:, :, StartFrame:EndFrame);

    % Process
    LD = smi_core.LocalizeData(DataChunk, SMF, 0);
    [SMD, ~] = LD.genLocalizations();

    % Adjust frame numbers
    SMD.FrameNum = SMD.FrameNum + StartFrame - 1;

    % Accumulate
    SMDTotal = smi_core.SingleMoleculeData.catSMD(SMDTotal, SMD);

    % Clear GPU memory
    reset(gpuDevice);
end
```

### Verbosity Levels

Control output detail:

```matlab
Verbose = 0;  % Silent (batch processing)
Verbose = 1;  % Basic progress (default)
Verbose = 2;  % Detailed progress
Verbose = 3;  % Full diagnostics + saved visualizations

LD = smi_core.LocalizeData(Data, SMF, Verbose);
```

---

## Troubleshooting

### GPU Issues

**Problem:** "Invalid PTX file" or GPU errors

**Solution:**
```matlab
% Recompile CUDA files
cd(fullfile(smite_path, 'MATLAB', 'source', 'cuda'));
cuda_Make
```

### Memory Errors

**Problem:** "Out of memory" during localization

**Solution:** Process in smaller chunks (see Memory Management above)

### Poor Localization Quality

**Problem:** Many localizations with low p-values

**Solution:** Check data quality and adjust parameters
```matlab
% Visualize raw data
imagesc(mean(Data, 3));
colorbar;
title('Mean image');

% Check photon levels
fprintf('Mean intensity: %.1f photons/pixel\n', mean(Data(:)));

% Adjust detection threshold
SMF.BoxFinding.MinPhotons = median(Data(:)) * 10;  % Rule of thumb
```

### Frame Connection Issues

**Problem:** Too many or too few connections

**Solutions:**
```matlab
% Too many connections (over-connecting)
SMF.FrameConnection.MaxSeparation = 0.5;  % Tighter spatial threshold
SMF.FrameConnection.LoS = 0.001;  % Stricter statistical threshold

% Too few connections (under-connecting)
SMF.FrameConnection.MaxSeparation = 2.0;  % Looser spatial threshold
SMF.FrameConnection.MaxFrameGap = 10;     % Allow longer gaps
```

### Drift Correction Failures

**Problem:** Drift correction diverges or produces artifacts

**Solutions:**
```matlab
% Adjust convergence parameters
SMF.DriftCorrection.L_intra = 2.0;  % Looser threshold
SMF.DriftCorrection.PDegree = 2;    % Higher order polynomial

% Try brightfield registration instead
SMF.DriftCorrection.Method = 'DC-BF';
```

---

## See Also

### Core Concepts
- [Architecture Overview](../core-concepts/architecture.md)
- [SMF Structure Guide](../core-concepts/smf-structure.md)
- [SMD Structure Guide](../core-concepts/smd-structure.md)
- [TR Structure Guide](../core-concepts/tr-structure.md)

### Workflows
- [SMLM Analysis Workflow](../workflows/smlm-analysis.md)
- [SPT Tracking Workflow](../workflows/spt-tracking.md)

### How-To Guides
- [How to Load Data](../how-to/load-data.md)
- [How to Localize Molecules](../how-to/localize-molecules.md)
- [How to Threshold Results](../how-to/threshold-results.md)
- [How to Use GPU](../how-to/use-gpu.md)

### Higher-Level APIs
- [+smi Namespace](smi.md) - Main workflow classes
- [+smi_vis Namespace](smi-vis.md) - Visualization tools

---

## Summary

The +smi_core namespace provides the foundational building blocks for smite analyses:

**Data Structures:** SMF, SMD, and TR define the complete analysis specification and results format.

**Core Pipeline:** LoadData → DataToPhotons → LocalizeData orchestrates the path from raw files to localizations.

**Post-Processing:** FrameConnection and DriftCorrection improve precision and accuracy.

**Algorithms:** GPU-accelerated GaussMLE provides fast, optimal fitting; FRC estimates resolution.

These classes are designed to be used together in pipelines (via smi.SMLM) or individually for custom analyses. Understanding +smi_core is essential for advanced smite usage, troubleshooting, and extending functionality.
