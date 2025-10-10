---
title: "SMF Structure: Analysis Parameters"
category: "core-concepts"
level: "intermediate"
tags: ["smf", "parameters", "configuration", "data-structures"]
prerequisites: ["architecture.md"]
related: ["smd-structure.md", "../workflows/smlm-analysis.md"]
summary: "Complete reference to the Single Molecule Fitting (SMF) parameter structure that defines all smite analyses"
estimated_time: "20 minutes"
last_updated: "2025-10-10"
status: "complete"
---

# SMF Structure: Analysis Parameters

## Purpose

The SMF (Single Molecule Fitting) structure is the cornerstone of smite's design philosophy: **a single structure completely defines the entire analysis pipeline**. Understanding SMF is essential for configuring analyses, ensuring reproducibility, and troubleshooting results. This document provides a complete reference to all SMF fields with practical guidance on setting them.

## Prerequisites

- Understanding of [smite architecture](architecture.md)
- Basic knowledge of SMLM/SPT concepts
- Familiarity with MATLAB structures

## Overview

SMF is a structure of structures containing all parameters required to go from raw camera data to final results. It answers every question about an analysis:

- Where is the data? → `SMF.Data`
- How do we find molecules? → `SMF.BoxFinding`
- How do we fit them? → `SMF.Fitting`
- How do we filter results? → `SMF.Thresholding`
- How do we connect across frames? → `SMF.FrameConnection`
- How do we correct drift? → `SMF.DriftCorrection`
- How do we track particles? → `SMF.Tracking`

**Key principle**: Same SMF + same data = same results (reproducibility).

## Creating and Using SMF

### Basic Creation

```matlab
% Create SMF with all default values
SMF = smi_core.SingleMoleculeFitting();

% SMF is now a structure with sub-structures:
% SMF.Data, SMF.BoxFinding, SMF.Fitting, etc.
```

### Accessing Fields

```matlab
% Read a value
boxSize = SMF.BoxFinding.BoxSize;  % Returns 7 (default)

% Modify a value
SMF.Fitting.PSFSigma = 1.3;

% Access nested fields
fileName = SMF.Data.FileName{1};  % Note: FileName is a cell array
```

### Using the GUI

For interactive parameter setting:

```matlab
SMF = smi_core.SingleMoleculeFitting();
SMF.gui();  % Opens graphical interface
```

The GUI provides:
- Parameter descriptions
- Valid ranges and types
- Real-time validation
- Save/load functionality

### Saving and Loading

```matlab
% SMF is just a structure, so standard save/load works
save('my_parameters.mat', 'SMF');
load('my_parameters.mat', 'SMF');

% Results files always include the SMF used
load('Results.mat', 'SMF', 'SMD');  % Reload exact parameters
```

## SMF Sub-Structures Reference

### Data: File and Camera Information

Specifies where data comes from and camera characteristics.

```matlab
SMF.Data
```

**Fields:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `FileName` | cell array of char | `{}` | File name(s) to analyze |
| `FileDir` | char array | `''` | Directory containing files |
| `ResultsDir` | char array | `'FileDir/Results'` | Where to save results |
| `AnalysisID` | char array | `''` | ID tagged onto saved results |
| `FileType` | char array | auto-detected | 'h5', 'mat', or custom |
| `DataVariable` | char array | `'sequence'` | Variable name in .mat files |
| `DatasetList` | int32 array | `[]` | Which datasets to analyze (empty = all) |
| `DatasetMods` | cell {incl; excl} | `{[]; []}` | Include/exclude dataset modifiers |
| `CameraType` | char | `'EMCCD'` | 'EMCCD' or 'SCMOS' |
| `CameraGain` | scalar or image | `1` | Conversion gain (ADU/photon) |
| `CameraOffset` | scalar or image | `0` | Camera baseline (ADU) |
| `CameraNoise` | scalar or image | `0` | Read noise (electrons RMS) |
| `CalibrationFilePath` | char | `''` | Path to camera calibration file |
| `RegistrationFilePath` | char | `''` | Path to channel registration file |
| `DataROI` | 4-element vector | `[]` | [YStart, XStart, YEnd, XEnd] ROI |
| `FrameRate` | scalar | required | Acquisition rate (Hz) |
| `PixelSize` | scalar | required | Camera pixel size (micrometers) |
| `SEAdjust` | scalar | `0` | Standard error inflation (pixels) |

**Example:**

```matlab
% Typical EMCCD setup
SMF.Data.FileDir = '/data/experiment1';
SMF.Data.FileName = {'Cell1_PAINT.h5'};
SMF.Data.PixelSize = 0.108;  % 108 nm
SMF.Data.FrameRate = 100;    % 100 Hz
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;   % Typical EMCCD gain
SMF.Data.CameraOffset = 100; % Typical offset

% Analyze only datasets 1-10, excluding 5
SMF.Data.DatasetMods = {1:10; 5};
```

**Camera calibration**: For SCMOS cameras, use a calibration file:

```matlab
SMF.Data.CalibrationFilePath = '/path/to/calibration.mat';
% File should contain: CameraGain, CameraOffset, CameraNoise as images
```

### BoxFinding: Detecting Molecules

Parameters for finding regions of interest (ROIs) around molecules.

```matlab
SMF.BoxFinding
```

**Fields:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `BoxSize` | integer | `7` | Size of fitting box (pixels) |
| `BoxOverlap` | integer | `2` | Allowed overlap between boxes (pixels) |
| `MinPhotons` | scalar | `200` | Minimum photons to trigger detection |

**Choosing BoxSize:**

```matlab
% Rule of thumb: BoxSize ≈ 4-5 × PSFSigma
PSFSigma = 1.3;  % pixels
BoxSize = ceil(4 * PSFSigma);  % typically 5-10 pixels
SMF.BoxFinding.BoxSize = BoxSize;
```

**Typical values:** 5-10 pixels depending on PSF sigma.

**Choosing MinPhotons:**

- Too low: Noise creates false detections
- Too high: Dim emitters missed
- Typical range: 100-500 for SMLM, 50-200 for SPT

```matlab
% For high SNR PAINT data
SMF.BoxFinding.MinPhotons = 300;

% For dim tracking data
SMF.BoxFinding.MinPhotons = 100;
```

**BoxOverlap** allows nearby emitters to both be detected:

```matlab
% Prevent overlapping fits (conservative)
SMF.BoxFinding.BoxOverlap = 0;

% Allow some overlap (more aggressive fitting)
SMF.BoxFinding.BoxOverlap = 2;
```

### Fitting: PSF Model and Fitting Method

Controls how molecule positions are estimated from detected boxes.

```matlab
SMF.Fitting
```

**Fields:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `PSFSigma` | scalar | `1` | Initial/fixed PSF sigma (pixels) |
| `FitType` | char | `'XYNB'` | Which parameters to fit |
| `NParams` | integer | auto-set | Number of fit parameters |
| `Iterations` | integer | `20` | Newton-Raphson iterations |
| `ZFitStruct` | structure | — | Astigmatism calibration (3D) |

**FitType options:**

| FitType | Parameters Fit | Use Case |
|---------|----------------|----------|
| `'XYNB'` | X, Y, Photons, Background | Standard 2D SMLM |
| `'XYNBS'` | X, Y, N, B, Sigma | Fit PSF width (variable PSF) |
| `'XYNBSXSY'` | X, Y, N, B, SigmaX, SigmaY | Asymmetric PSF |
| `'XYZNB'` | X, Y, Z, N, B | 3D astigmatism |

**Example:**

```matlab
% Standard 2D fitting with fixed PSF
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Let PSF width vary (useful for varying focus)
SMF.Fitting.PSFSigma = 1.3;  % Initial guess
SMF.Fitting.FitType = 'XYNBS';

% 3D astigmatism fitting (requires calibration)
SMF.Fitting.FitType = 'XYZNB';
SMF.Fitting.ZFitStruct.Ax = 266.5;  % Calibration parameters
SMF.Fitting.ZFitStruct.Bx = -1533;
SMF.Fitting.ZFitStruct.Ay = 266.5;
SMF.Fitting.ZFitStruct.By = 1533;
SMF.Fitting.ZFitStruct.Gamma = 0.5;
SMF.Fitting.ZFitStruct.D = 0.5;
```

**Determining PSFSigma:**

The relationship between PSF sigma and optical parameters:

```
σ_pixels = 0.21 * λ / (NA * pixel_size_um)
```

For λ=670nm, NA=1.49, pixel=108nm: σ ≈ 1.3 pixels

Empirically, fit some data with `FitType='XYNBS'` and look at `median(SMD.PSFSigma)`.

### Thresholding: Quality Filtering

Filters out poor quality localizations.

```matlab
SMF.Thresholding
```

**Fields:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `On` | logical | `true` | Enable thresholding? |
| `MaxXY_SE` | scalar | `0.2` | Max XY precision (pixels) |
| `MaxZ_SE` | scalar | `0.05` | Max Z precision (micrometers) |
| `MinPValue` | scalar | `0.01` | Min fit p-value |
| `AutoThreshLogL` | logical | `false` | Auto-threshold on log-likelihood |
| `AutoThreshPrctile` | scalar | `1e-4` | Percentile for auto-threshold |
| `MinPSFSigma` | scalar | `0.5` | Min PSF sigma (pixels) |
| `MaxPSFSigma` | scalar | `2.0` | Max PSF sigma (pixels) |
| `MinPhotons` | scalar | `100` | Min photons after fit |
| `MaxBg` | scalar | `Inf` | Max background (photons/pixel) |
| `InMeanMultiplier` | scalar | `Inf` | Max intensity acceptance |
| `NNMedianMultiplier` | scalar | `3` | Nearest-neighbor filter |
| `MinNumNeighbors` | scalar | `0` | Min neighbors required |

**Most important thresholds:**

```matlab
% Precision threshold - adjust based on required localization accuracy
SMF.Thresholding.MaxXY_SE = 0.2;  % 20 nm for 100 nm pixels

% Photon threshold - remove dim fits
SMF.Thresholding.MinPhotons = 100;

% Fit quality - remove poor fits
SMF.Thresholding.MinPValue = 0.01;

% PSF sigma range - remove misfits
SMF.Thresholding.MinPSFSigma = 0.8;
SMF.Thresholding.MaxPSFSigma = 1.8;
```

**Automatic thresholding:**

```matlab
% Let smite determine log-likelihood threshold automatically
SMF.Thresholding.AutoThreshLogL = true;
% This overrides MinPValue
```

**Checking thresholding effects:**

After analysis, check `SMD.ThreshFlag` to see why localizations were rejected (0 = kept, >0 = rejected).

### FrameConnection: Linking Across Frames

Connects the same emitter appearing in multiple frames.

```matlab
SMF.FrameConnection
```

**Fields:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `On` | logical | `true` | Enable frame connection? |
| `Method` | char | `'LAP-FC'` | Connection algorithm |
| `MaxSeparation` | scalar | `1.0` | Max distance (pixels) |
| `LoS` | scalar | `0.01` | Level of significance |
| `MaxFrameGap` | integer | `5` | Max frames to bridge |
| `NSigmaDev` | scalar | `5` | SE multiplier for clustering |
| `NNearestClusters` | integer | `2` | Clusters for density |
| `NIterations` | integer | `1` | Iterative FC attempts |
| `MinNFrameConns` | integer | `1` | Min connections to keep |

**Purpose:** SMLM emitters blink on/off. Frame connection links appearances of the same emitter, improving precision.

**Example:**

```matlab
% Standard SMLM frame connection
SMF.FrameConnection.On = true;
SMF.FrameConnection.Method = 'LAP-FC';  % Linear assignment
SMF.FrameConnection.MaxSeparation = 1;  % 1 pixel max movement
SMF.FrameConnection.MaxFrameGap = 5;    % Bridge up to 5 frames

% Require multiple frames for stability
SMF.FrameConnection.MinNFrameConns = 2;  % Must appear in ≥2 frames
```

**MaxSeparation** should be larger than typical localization uncertainty but smaller than inter-emitter spacing.

**After frame connection:**
- `SMD.ConnectID`: Integer identifying same emitter across frames
- Connected localizations can be combined for improved precision

### DriftCorrection: Stage Drift Compensation

Corrects stage drift during acquisition.

```matlab
SMF.DriftCorrection
```

**Fields:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `On` | logical | `true` | Enable drift correction? |
| `Method` | char | `'DC-KNN'` | Correction algorithm |
| `BFRegistration` | logical | `true` | Use brightfield registration? |
| `L_intra` | scalar | `1.0` | Intra-dataset threshold (pixels) |
| `L_inter` | scalar | `2.0` | Inter-dataset threshold (pixels) |
| `PixelSizeZUnit` | scalar | `0.1` | XY pixel size for 3D (micrometers) |
| `PDegree` | integer | `1` | Polynomial degree for drift rate |

**Methods:**
- `'DC-KNN'`: K-nearest neighbors (default, robust)
- Other methods available for specific use cases

**Example:**

```matlab
% Standard drift correction
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.Method = 'DC-KNN';
SMF.DriftCorrection.L_intra = 1.0;   % Smooth within dataset
SMF.DriftCorrection.L_inter = 2.0;   % Align across datasets

% For high drift
SMF.DriftCorrection.PDegree = 2;  % Quadratic drift model
```

**Results:**
- `SMD.DriftX`, `SMD.DriftY`, `SMD.DriftZ`: Estimated drift per frame
- Positions in SMD are already drift-corrected

### Tracking: Single Particle Tracking

Parameters for SPT (tracking moving particles).

```matlab
SMF.Tracking
```

**Fields:**

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `Method` | char | `'CostMatrix'` | Tracking algorithm |
| `D` | scalar | `0.01` | Diffusion constant (pixels²/frame) |
| `TrajwiseD` | logical | `true` | Use per-trajectory D? |
| `K_on` | scalar | `0.9` | Off→On rate (frame⁻¹) |
| `K_off` | scalar | `0.1` | On→Off rate (frame⁻¹) |
| `MaxDistFF` | scalar | `5` | Max distance frame-to-frame (pixels) |
| `MaxDistGC` | scalar | `10` | Max distance gap closing (pixels) |
| `MaxFrameGap` | integer | `10` | Max gap to close (frames) |
| `MinTrackLength` | integer | `3` | Min trajectory length (frames) |
| `NIterMax` | integer | `5` | Max tracking iterations |
| `NIterMaxBatch` | integer | `5` | Max batch iterations |
| `MaxRelativeChange` | scalar | `1e-5` | Convergence threshold |
| `MaxZScoreDist` | scalar | `Inf` | Max z-score for jumps |
| `MaxZScorePhotons` | scalar | `Inf` | Max z-score for photons |
| `TryLowPValueLocs` | logical | `false` | Include low p-value locs? |

**Example:**

```matlab
% Fast diffusion tracking
SMF.Tracking.D = 0.1;  % pixels²/frame
SMF.Tracking.MaxDistFF = 5;  % Allow 5 pixel jumps
SMF.Tracking.MaxDistGC = 10;  % Bridge gaps up to 10 pixels
SMF.Tracking.MaxFrameGap = 5;  % Maximum 5-frame blinks
SMF.Tracking.MinTrackLength = 10;  % Keep trajectories ≥10 frames

% Trajectory-wise diffusion estimation
SMF.Tracking.TrajwiseD = true;  % Estimate D per trajectory
```

**Choosing parameters:**

Diffusion constant in pixels²/frame:
```
D_pixels = D_physical [μm²/s] / (pixel_size² [μm²] × frame_rate [Hz])
```

For D=0.1 μm²/s, pixel=0.1 μm, rate=100 Hz:
```
D_pixels = 0.1 / (0.01 × 100) = 0.1 pixels²/frame
```

## Complete Example

Comprehensive SMF configuration for typical SMLM analysis:

```matlab
% Create SMF
SMF = smi_core.SingleMoleculeFitting();

% Data location and camera
SMF.Data.FileDir = '/data/2024-01-10/Cell1';
SMF.Data.FileName = {'DNA_PAINT_647.h5'};
SMF.Data.PixelSize = 0.108;  % 108 nm
SMF.Data.FrameRate = 100;    % 100 Hz
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;
SMF.Data.CameraOffset = 100;
SMF.Data.ResultsDir = '/data/2024-01-10/Cell1/Results';
SMF.Data.AnalysisID = 'v2_tighter_threshold';

% Box finding
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 250;

% Fitting
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Thresholding
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.15;
SMF.Thresholding.MinPhotons = 150;
SMF.Thresholding.MinPSFSigma = 0.9;
SMF.Thresholding.MaxPSFSigma = 1.7;
SMF.Thresholding.MinPValue = 0.01;

% Frame connection
SMF.FrameConnection.On = true;
SMF.FrameConnection.MaxSeparation = 1.0;
SMF.FrameConnection.MaxFrameGap = 5;
SMF.FrameConnection.MinNFrameConns = 2;

% Drift correction
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.Method = 'DC-KNN';

% Save for reproducibility
save('SMF_parameters.mat', 'SMF');
```

## See Also

- [SMD Structure](smd-structure.md) - Results data format
- [Architecture Overview](architecture.md) - How SMF fits into smite
- [SMLM Workflow](../workflows/smlm-analysis.md) - Using SMF in practice
- [How to Load Data](../how-to/load-data.md) - Data configuration details
- doc/DataStructures/SMF.md (in repository) - Complete field reference
