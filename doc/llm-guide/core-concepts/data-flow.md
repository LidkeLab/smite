---
title: "Data Flow in smite: From Raw Images to Results"
category: "core-concepts"
level: "intermediate"
tags: ["data-flow", "pipeline", "smlm", "spt", "workflow", "architecture"]
prerequisites: ["architecture.md", "smf-structure.md", "smd-structure.md"]
related: ["../workflows/smlm-analysis.md", "../workflows/spt-tracking.md", "../how-to/localize-molecules.md"]
summary: "Comprehensive guide to how data flows through smite from raw camera images to final analysis results, covering all transformation steps and decision points"
estimated_time: "30 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# Data Flow in smite: From Raw Images to Results

## Purpose

Understanding data flow is fundamental to using smite effectively. This document traces the complete journey of data through smite's analysis pipeline, explaining what happens at each transformation step, what controls the flow, and when different processing options are applied. Whether you're troubleshooting an analysis, optimizing parameters, or building custom workflows, this knowledge is essential.

## Prerequisites

- Understanding of [smite architecture](architecture.md)
- Familiarity with [SMF structure](smf-structure.md)
- Familiarity with [SMD structure](smd-structure.md)
- Basic knowledge of SMLM/SPT concepts

## Overview

smite transforms raw camera images into scientific results through a series of well-defined stages. The pipeline is:

1. **Input**: Raw camera data (ADU values from .h5/.mat files)
2. **Loading**: Read data from disk into memory
3. **Calibration**: Convert ADU to photons using camera parameters
4. **Localization**: Find and fit molecules to extract positions
5. **Quality Control**: Filter low-quality localizations
6. **Frame Connection**: Link same molecule across frames (SMLM)
7. **Tracking**: Link moving molecules across frames (SPT)
8. **Drift Correction**: Compensate for stage drift
9. **Output**: SMD/TR results, super-resolution images, analysis plots

Each stage is controlled by SMF parameters and produces data in standardized formats. The beauty of this design: the same SMF parameters applied to the same raw data always produce identical results.

## Complete Data Flow Diagram

```
┌──────────────────────────────────────────────────────────────────┐
│                         RAW INPUT DATA                           │
│  .h5 or .mat files containing camera frames (Y × X × NFrames)   │
│  Values: ADU (Analog-to-Digital Units)                           │
└────────────────────────┬─────────────────────────────────────────┘
                         │
                         │ SMF.Data.FileDir
                         │ SMF.Data.FileName
                         │ SMF.Data.DatasetList
                         │
                    ┌────▼────┐
                    │LoadData │ smi_core.LoadData
                    └────┬────┘
                         │
                         │ Output: Image stack (Y × X × Frames)
                         │         Data type: uint16, int16, or double
                         │         Values: Still in ADU
                         │
┌────────────────────────▼─────────────────────────────────────────┐
│                    IMAGE STACK (ADU)                             │
│  3D array: rows (Y) × columns (X) × frames (time)                │
│  DataROI applied if specified                                    │
└────────────────────────┬─────────────────────────────────────────┘
                         │
                         │ SMF.Data.CameraGain
                         │ SMF.Data.CameraOffset
                         │ SMF.Data.CameraType
                         │ SMF.Data.CalibrationFilePath
                         │
                  ┌──────▼───────┐
                  │DataToPhotons │ smi_core.DataToPhotons
                  └──────┬───────┘
                         │
                         │ Photons = (ADU - Offset) / Gain
                         │ Pixel-wise for sCMOS, scalar for EMCCD
                         │
┌────────────────────────▼─────────────────────────────────────────┐
│                  IMAGE STACK (PHOTONS)                           │
│  Same dimensions, now in calibrated photon units                 │
│  Noise model properly defined                                    │
└────────────────────────┬─────────────────────────────────────────┘
                         │
                         │ SMF.BoxFinding.MinPhotons
                         │ SMF.BoxFinding.BoxSize
                         │ SMF.BoxFinding.BoxOverlap
                         │
                    ┌────▼────┐
                    │FindROI  │ smi_core.FindROI (called by LocalizeData)
                    └────┬────┘
                         │
                         │ For each frame:
                         │   - Detect bright spots
                         │   - Place boxes around candidates
                         │ Output: Box coordinates (N_boxes × 2)
                         │
                         │ SMF.Fitting.PSFSigma
                         │ SMF.Fitting.FitType
                         │ SMF.Fitting.Iterations
                         │
                   ┌─────▼──────┐
                   │GaussMLE    │ smi_core.GaussMLE (called by LocalizeData)
                   └─────┬──────┘
                         │
                         │ For each box:
                         │   - Fit 2D Gaussian PSF
                         │   - Extract: X, Y, N, Bg (and σ, Z if requested)
                         │   - Compute uncertainties (CRLB)
                         │   - Calculate fit quality (LogL, PValue)
                         │ Output: Raw fit parameters
                         │
                         │ SMF.Thresholding.MaxXY_SE
                         │ SMF.Thresholding.MinPhotons
                         │ SMF.Thresholding.MinPValue
                         │ SMF.Thresholding.MinPSFSigma
                         │ SMF.Thresholding.MaxPSFSigma
                         │
                ┌────────▼─────────┐
                │ThresholdFits     │ smi_core.Threshold (called by LocalizeData)
                └────────┬─────────┘
                         │
                         │ Apply quality filters
                         │ Set ThreshFlag for rejected fits
                         │
┌────────────────────────▼─────────────────────────────────────────┐
│              SMDPreThresh (All Localizations)                    │
│  Contains ALL fits (good and bad)                                │
│  ThreshFlag: 0 = passed, >0 = failed (bitwise flags)            │
│  Fields: X, Y, X_SE, Y_SE, Photons, Bg, FrameNum, etc.         │
└────────────────────────┬─────────────────────────────────────────┘
                         │
                         │ Keep only ThreshFlag == 0
                         │
┌────────────────────────▼─────────────────────────────────────────┐
│                 SMD (Filtered Localizations)                     │
│  Only localizations that passed quality filters                  │
│  Ready for downstream processing                                 │
└────────────────────────┬─────────────────────────────────────────┘
                         │
                    ┌────┴────┐
                    │ BRANCH  │ Workflow decision point
                    └────┬────┘
                         │
          ┌──────────────┴──────────────┐
          │                             │
    [SMLM PATH]                   [SPT PATH]
          │                             │
          │ SMF.FrameConnection         │ SMF.Tracking
          │     .On = true              │     .Method = 'LAP'
          │     .Method = 'LAP-FC'      │     .D, .K_on, .K_off
          │     .MaxSeparation          │     .MaxDistFF
          │     .MaxFrameGap            │     .MaxDistGC
          │                             │
    ┌─────▼──────┐              ┌──────▼──────┐
    │FrameConn   │              │  genTrajFF  │ smi.SPT
    └─────┬──────┘              └──────┬──────┘
          │                             │
          │ Links same emitter          │ Frame-to-frame
          │ across frames               │ tracking
          │ Assigns ConnectID           │
          │                             │
          │                      ┌──────▼──────┐
          │                      │  genTrajGC  │ smi.SPT
          │                      └──────┬──────┘
          │                             │
          │                             │ Gap closing
          │                             │ Assigns ConnectID
          │                             │
┌─────────▼─────────────┐     ┌─────────▼─────────────┐
│ SMD (Connected)       │     │ SMD (Tracked)         │
│ Same emitter linked   │     │ Trajectories formed   │
│ via ConnectID         │     │ via ConnectID         │
└─────────┬─────────────┘     └─────────┬─────────────┘
          │                             │
          └──────────────┬──────────────┘
                         │
                         │ Both paths converge
                         │
                         │ SMF.DriftCorrection.On
                         │ SMF.DriftCorrection.Method
                         │ SMF.DriftCorrection.L_intra
                         │
              ┌──────────▼──────────┐
              │DriftCorrection      │ smi_core.DriftCorrection
              │  (Intra-dataset)    │
              └──────────┬──────────┘
                         │
                         │ Estimate drift within each dataset
                         │ Apply corrections to X, Y, Z
                         │ Store drift: DriftX, DriftY, DriftZ
                         │
┌────────────────────────▼─────────────────────────────────────────┐
│             SMD (Drift-Corrected, Per Dataset)                   │
│  Positions corrected for intra-dataset drift                     │
│  One SMD per dataset                                             │
└────────────────────────┬─────────────────────────────────────────┘
                         │
                         │ Loop over all datasets
                         │ smi_core.SingleMoleculeData.catSMD()
                         │
┌────────────────────────▼─────────────────────────────────────────┐
│              SMD (Concatenated Datasets)                         │
│  All datasets combined                                           │
│  DatasetNum field tracks origin                                  │
└────────────────────────┬─────────────────────────────────────────┘
                         │
                         │ If multiple datasets
                         │ SMF.DriftCorrection.L_inter
                         │
              ┌──────────▼──────────┐
              │DriftCorrection      │ smi_core.DriftCorrection
              │  (Inter-dataset)    │
              └──────────┬──────────┘
                         │
                         │ Align datasets to common reference
                         │ Update X, Y, Z and drift arrays
                         │
┌────────────────────────▼─────────────────────────────────────────┐
│              SMD (Final, Fully Corrected)                        │
│  Complete results ready for analysis                             │
│  Fields: X, Y, Z, X_SE, Y_SE, Z_SE, Photons, Bg, PSFSigma,     │
│          FrameNum, DatasetNum, ConnectID, DriftX, DriftY, etc.  │
└────────────────────────┬─────────────────────────────────────────┘
                         │
                    ┌────┴────┐
                    │ OUTPUT  │ Multiple output paths
                    └────┬────┘
                         │
          ┌──────────────┼──────────────┬──────────────┐
          │              │              │              │
    ┌─────▼─────┐  ┌─────▼─────┐ ┌─────▼─────┐ ┌─────▼─────┐
    │Save SMD/  │  │Generate   │ │Generate   │ │Convert to │
    │SMF to     │  │SR Images  │ │Plots &    │ │TR for SPT │
    │.mat file  │  │           │ │Histograms │ │analysis   │
    └───────────┘  └───────────┘ └───────────┘ └───────────┘
         │              │              │              │
         │              │              │              │
    Results.mat    GaussImage     Precision      TR array
                   HistImage      Photons        (1 SMD per
                   CircleImage    Bg, PSFσ       trajectory)
                                  Drift plots
```

## Data Transformations Explained

### Stage 1: Loading Raw Data

**Input**: File path and dataset specifications
**Output**: 3D image array (Y × X × Frames) in ADU

```matlab
% What happens
LD = smi_core.LoadData();
[~, sequence, SMF] = LD.loadRawData(SMF, DatasetIndex);

% sequence dimensions
[nY, nX, nFrames] = size(sequence);  % e.g., 256 × 256 × 1000
```

**Controlled by**:
- `SMF.Data.FileDir` - Where to find files
- `SMF.Data.FileName` - Which file(s) to load
- `SMF.Data.DatasetList` - Which datasets within file
- `SMF.Data.DataROI` - Spatial crop region [YStart, XStart, YEnd, XEnd]

**Data format**:
- HDF5 (.h5): Each dataset is a 3D array
- MAT (.mat): Variable specified by `SMF.Data.DataVariable` (default: 'sequence')
- Values: Raw camera ADU (typically uint16)

**Example**:
```matlab
SMF.Data.FileDir = '/data/2024-01-10';
SMF.Data.FileName = {'Cell1_PAINT.h5'};
SMF.Data.DatasetList = int32([1, 2, 3]);  % Load first 3 datasets
SMF.Data.DataROI = [50, 50, 200, 200];    % 150×150 pixel ROI

LD = smi_core.LoadData();
[~, seq, SMF] = LD.loadRawData(SMF, 1);  % Load dataset 1
fprintf('Loaded: %d × %d × %d\n', size(seq));
```

### Stage 2: Converting ADU to Photons

**Input**: Image stack in ADU
**Output**: Image stack in calibrated photons

```matlab
% What happens internally
DTP = smi_core.DataToPhotons(SMF, sequence);
photons = DTP.convertData();

% The conversion
if scalar_calibration
    photons = (sequence - CameraOffset) / CameraGain;
else  % sCMOS pixel-wise
    photons = (sequence - CameraOffset) ./ CameraGain;  % Element-wise
end
```

**Controlled by**:
- `SMF.Data.CameraType` - 'EMCCD' or 'SCMOS'
- `SMF.Data.CameraGain` - Scalar (EMCCD) or image (sCMOS)
- `SMF.Data.CameraOffset` - Scalar or image
- `SMF.Data.CalibrationFilePath` - For sCMOS calibration

**Camera types**:

**EMCCD** (scalar calibration):
```matlab
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;    % ADU per photon
SMF.Data.CameraOffset = 100;  % ADU
% All pixels use same gain/offset
```

**sCMOS** (pixel-wise calibration):
```matlab
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = '/path/to/calibration.mat';
% File contains: CameraGain, CameraOffset, CameraNoise (all Y×X images)
% Each pixel has its own gain/offset/noise
```

**Why this matters**:
- Proper noise model for fitting (Poisson for photons)
- Correct uncertainties in CRLB
- Fair comparison across pixels (sCMOS)

### Stage 3: Detecting Molecules (Box Finding)

**Input**: Image stack in photons
**Output**: List of box coordinates for each frame

```matlab
% What happens (inside LocalizeData)
FR = smi_core.FindROI();
for frame = 1:nFrames
    [boxes, ~] = FR.findROI(photons(:,:,frame), SMF.BoxFinding);
    % boxes is N_detections × 2 array of [Y, X] coordinates
end
```

**Controlled by**:
- `SMF.BoxFinding.BoxSize` - Size of fitting region (pixels)
- `SMF.BoxFinding.MinPhotons` - Detection threshold
- `SMF.BoxFinding.BoxOverlap` - Allowed box overlap

**Algorithm**:
1. Apply Gaussian smoothing to reduce noise
2. Find local maxima in smoothed image
3. For each maximum:
   - Estimate photon count in neighborhood
   - If photons > MinPhotons, create box
4. Place BoxSize × BoxSize box centered on detection
5. Handle overlapping boxes based on BoxOverlap setting

**Key parameters**:

```matlab
% BoxSize: Should contain ~95% of PSF
% Rule: BoxSize ≈ 4-5 × PSFSigma
SMF.BoxFinding.BoxSize = 7;  % For PSFSigma ≈ 1.3-1.5

% MinPhotons: Threshold for detection
% Too low → many false detections from noise
% Too high → miss dim emitters
SMF.BoxFinding.MinPhotons = 250;  % Typical SMLM

% BoxOverlap: Prevents multiple fits to same emitter
SMF.BoxFinding.BoxOverlap = 2;  % Allow 2-pixel overlap
```

**Typical output**: 10-100 boxes per frame for SMLM, 1-10 for SPT

### Stage 4: Fitting PSF Models

**Input**: Image boxes (BoxSize × BoxSize × N_boxes)
**Output**: Fit parameters and uncertainties for each box

```matlab
% What happens (inside LocalizeData)
GM = smi_core.GaussMLE();
for box = 1:N_boxes
    box_data = photons(y:y+BoxSize-1, x:x+BoxSize-1, frame);
    [params, CRLB, LogL] = GM.gaussMLE(box_data, SMF.Fitting);
    % params = [X, Y, Photons, Bg] or more depending on FitType
    % CRLB = Cramér-Rao Lower Bound (uncertainties)
end
```

**Controlled by**:
- `SMF.Fitting.PSFSigma` - Initial/fixed PSF width
- `SMF.Fitting.FitType` - Which parameters to fit
- `SMF.Fitting.Iterations` - Max Newton-Raphson iterations

**Fit types and parameters**:

| FitType | Parameters | Use Case |
|---------|------------|----------|
| XYNB | X, Y, Photons, Bg | Standard 2D SMLM (fixed σ) |
| XYNBS | X, Y, N, Bg, Sigma | Variable PSF width |
| XYNBSXSY | X, Y, N, Bg, σ_x, σ_y | Asymmetric PSF |
| XYZNB | X, Y, Z, N, Bg | 3D astigmatism |

**Algorithm** (simplified):
1. Initialize: X₀, Y₀, N₀, B₀ from image data
2. For each iteration:
   - Compute expected image from current parameters
   - Compute log-likelihood: L = Σ [data·log(model) - model]
   - Compute gradient: ∂L/∂θ and Hessian: ∂²L/∂θ²
   - Update: θ_new = θ_old - H⁻¹·∇L (Newton-Raphson)
3. Check convergence (gradient magnitude)
4. Compute CRLB: uncertainties = √(diag(H⁻¹))

**Example**:
```matlab
SMF.Fitting.PSFSigma = 1.3;  % pixels (from optics or empirical)
SMF.Fitting.FitType = 'XYNB';  % Standard 2D fitting
SMF.Fitting.Iterations = 20;   % Usually converges in 5-10

% For 3D astigmatism
SMF.Fitting.FitType = 'XYZNB';
SMF.Fitting.ZFitStruct.Ax = 266.5;   % From calibration
SMF.Fitting.ZFitStruct.Bx = -1533;
% ... other calibration parameters
```

**Requirements**:
- CUDA GPU (NVIDIA, compute capability ≥ 5.0)
- Fitting runs on GPU for speed (~1000× faster than CPU)

**Output fields populated**:
- `X`, `Y` - Position (pixels, sub-pixel precision)
- `X_SE`, `Y_SE` - Position uncertainties (pixels)
- `Photons` - Integrated photon count
- `Bg` - Background (photons/pixel)
- `PSFSigma` - PSF width if fitted
- `LogLikelihood` - Fit quality metric
- `PValue` - Goodness of fit (0-1)

### Stage 5: Quality Filtering (Thresholding)

**Input**: All fits (raw SMD)
**Output**: SMDPreThresh (all fits with flags) and SMD (filtered)

```matlab
% What happens
THR = smi_core.Threshold();
[SMD, SMDPreThresh] = THR.thresholdSMD(SMD_raw, SMF.Thresholding);

% SMDPreThresh contains ALL fits
% SMD contains only fits where ThreshFlag == 0
```

**Controlled by**:
```matlab
SMF.Thresholding.On              % Enable filtering
SMF.Thresholding.MaxXY_SE        % Max position error
SMF.Thresholding.MinPhotons      % Min photon count
SMF.Thresholding.MinPValue       % Min fit quality
SMF.Thresholding.MinPSFSigma     % PSF width range
SMF.Thresholding.MaxPSFSigma
```

**Filters applied** (ThreshFlag bit assignments):

| Bit | Filter | Condition | Typical Value |
|-----|--------|-----------|---------------|
| 1 | Precision | X_SE or Y_SE > MaxXY_SE | 0.2 pixels |
| 2 | Photons | Photons < MinPhotons | 100-200 |
| 4 | Fit quality | PValue < MinPValue | 0.01 |
| 8 | PSF sigma | σ outside range | 0.8-1.8 pixels |
| 16 | Background | Bg > MaxBg | Usually Inf |

**Example**:
```matlab
% Standard quality filtering
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.2;      % 20 nm for 100 nm pixels
SMF.Thresholding.MinPhotons = 150;    % Post-fit photons
SMF.Thresholding.MinPValue = 0.01;    % 99% confidence
SMF.Thresholding.MinPSFSigma = 0.8;   % Reasonable PSF range
SMF.Thresholding.MaxPSFSigma = 1.8;

% After filtering
load('Results.mat', 'SMDPreThresh');
passed = sum(SMDPreThresh.ThreshFlag == 0);
total = length(SMDPreThresh.ThreshFlag);
fprintf('Passed: %d / %d (%.1f%%)\n', passed, total, 100*passed/total);

% See what failed
failed_precision = bitand(SMDPreThresh.ThreshFlag, 1) > 0;
fprintf('Failed precision: %d\n', sum(failed_precision));
```

**Automatic thresholding**:
```matlab
% Let smite determine optimal threshold
SMF.Thresholding.AutoThreshLogL = true;
SMF.Thresholding.AutoThreshPrctile = 1e-4;  % Keep top 99.99%
% Overrides MinPValue, uses log-likelihood distribution
```

### Stage 6A: Frame Connection (SMLM Path)

**Input**: SMD with independent localizations
**Output**: SMD with ConnectID linking same emitter

```matlab
% What happens
FC = smi_core.FrameConnection(SMD, SMF);
SMD_connected = FC.performFrameConnection();

% SMD_connected.ConnectID:
%   0 = not connected (only appears once)
%   N > 0 = cluster ID (appears in multiple frames)
```

**Controlled by**:
```matlab
SMF.FrameConnection.On              % Enable
SMF.FrameConnection.Method          % 'LAP-FC' (Linear Assignment)
SMF.FrameConnection.MaxSeparation   % Max distance (pixels)
SMF.FrameConnection.MaxFrameGap     % Max frame gap to bridge
SMF.FrameConnection.MinNFrameConns  % Min appearances to keep
```

**Algorithm (LAP-FC)**:
1. For consecutive frames i and i+1:
   - Build cost matrix C where C[a,b] = distance(loc_a, loc_b)
   - Add "dummy" entries for birth/death
2. Solve Linear Assignment Problem (Hungarian algorithm)
3. Accept assignments where distance < MaxSeparation
4. Bridge gaps up to MaxFrameGap frames
5. Assign unique ConnectID to each connected cluster
6. Combine localizations within cluster (weighted average)

**Example**:
```matlab
SMF.FrameConnection.On = true;
SMF.FrameConnection.Method = 'LAP-FC';
SMF.FrameConnection.MaxSeparation = 1.0;   % 1 pixel max movement
SMF.FrameConnection.MaxFrameGap = 5;       % Bridge 5-frame blinks
SMF.FrameConnection.MinNFrameConns = 2;    % Must appear ≥2 times

% Check compression ratio
unique_IDs = unique(SMD.ConnectID(SMD.ConnectID > 0));
n_connected = length(unique_IDs);
n_total = length(SMD.X);
compression = n_total / n_connected;
fprintf('Frame connection: %.1f:1 compression\n', compression);
% Typical: 2:1 to 10:1 for SMLM
```

**Benefits**:
- Improved precision (combining measurements)
- Identification of blinking behavior
- Reduced data volume
- Better photon statistics

### Stage 6B: Tracking (SPT Path)

**Input**: SMD with independent localizations
**Output**: SMD with ConnectID forming trajectories

```matlab
% What happens
SPTobj = smi.SPT(SMF);
SPTobj.SMD = SMD;
SPTobj.generateTrajectories();  % Calls genTrajFF then genTrajGC

% Two-stage tracking:
% 1. Frame-to-frame: Link nearby localizations in consecutive frames
% 2. Gap closing: Bridge blinks/missed detections
```

**Controlled by**:
```matlab
SMF.Tracking.D                 % Diffusion constant (pixels²/frame)
SMF.Tracking.K_on, K_off       % Blinking rates
SMF.Tracking.MaxDistFF         % Frame-to-frame max distance
SMF.Tracking.MaxDistGC         % Gap closing max distance
SMF.Tracking.MaxFrameGap       % Max gap to close
SMF.Tracking.MinTrackLength    % Min trajectory length
```

**Algorithm**:

**Stage 1 - Frame-to-frame tracking**:
1. For frames i and i+1:
   - Build cost matrix based on:
     - Spatial distance
     - Diffusion model
     - Blinking probabilities
2. Add costs for track birth/death
3. Solve LAP (Hungarian algorithm)
4. Create initial trajectories

**Stage 2 - Gap closing**:
1. Find trajectory starts and ends
2. Build cost matrix for possible connections:
   - Spatial distance over gap
   - Expected motion based on D
   - Gap length penalty
3. Solve LAP for optimal gap closing
4. Merge trajectories

**Example**:
```matlab
% Configure tracking parameters
SMF.Tracking.D = 0.1;                    % pixels²/frame
SMF.Tracking.K_on = 0.9;                 % High on-rate
SMF.Tracking.K_off = 0.1;                % Low off-rate (stable)
SMF.Tracking.MaxDistFF = 5;              % 5 pixels/frame max
SMF.Tracking.MaxDistGC = 10;             % Bridge larger gaps
SMF.Tracking.MaxFrameGap = 5;            % Max 5-frame blinks
SMF.Tracking.MinTrackLength = 10;        % Keep tracks ≥10 frames
SMF.Tracking.TrajwiseD = true;           % Estimate D per trajectory

% Run tracking
SPTobj = smi.SPT(SMF);
[TR, SMD, ~] = SPTobj.performFullAnalysis();

% TR is TrackingResults: array of SMD (one per trajectory)
fprintf('Found %d trajectories\n', length(TR));
```

**Output**:
- `SMD.ConnectID` - Trajectory identifier
- `TR` - Array of SMD structures (one per trajectory)
- Each TR[i] contains all localizations for trajectory i

### Stage 7: Drift Correction

**Input**: SMD with localizations (possibly connected/tracked)
**Output**: SMD with corrected positions and drift estimates

Drift correction happens in two stages for multiple datasets:

**Intra-dataset**: Correct drift within each dataset separately
**Inter-dataset**: Align all datasets to common reference

```matlab
% What happens
DC = smi_core.DriftCorrection(SMF);

% For each dataset
for i = 1:NDatasets
    SMD_i = DC.driftCorrectKNNIntra(SMD_i, i, i);
end

% After all datasets processed
SMD_all = DC.driftCorrectKNNInter(SMD_concatenated);
```

**Controlled by**:
```matlab
SMF.DriftCorrection.On          % Enable
SMF.DriftCorrection.Method      % 'DC-KNN' (K-nearest neighbors)
SMF.DriftCorrection.L_intra     % Intra-dataset smoothing (pixels)
SMF.DriftCorrection.L_inter     % Inter-dataset smoothing (pixels)
SMF.DriftCorrection.PDegree     % Polynomial degree (drift model)
```

**Algorithm (DC-KNN)**:

1. **Build reference frame**: Average positions across sliding window
2. **For each time point**:
   - Find K nearest neighbors in reference frame
   - Match to localizations at current time
   - Compute median displacement
3. **Smooth trajectory**:
   - Regularize with penalty on displacement jumps
   - Fit polynomial of degree PDegree
4. **Apply correction**: X_corrected = X_raw - DriftX

**Example**:
```matlab
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.Method = 'DC-KNN';
SMF.DriftCorrection.L_intra = 1.0;   % Moderate smoothing within dataset
SMF.DriftCorrection.L_inter = 2.0;   % More smoothing across datasets
SMF.DriftCorrection.PDegree = 1;     % Linear drift model

% After correction, examine drift
figure;
subplot(2,1,1);
plot(SMD.DriftX(:,1) * SMD.PixelSize * 1000);  % Convert to nm
xlabel('Frame'); ylabel('X Drift (nm)');
title('Estimated Drift');

subplot(2,1,2);
plot(SMD.DriftY(:,1) * SMD.PixelSize * 1000);
xlabel('Frame'); ylabel('Y Drift (nm)');

% Positions in SMD.X, SMD.Y are already corrected
% To get raw positions:
X_raw = SMD.X - SMD.DriftX(SMD.FrameNum, SMD.DatasetNum);
```

**Requirements for good drift correction**:
- Sufficient localizations per frame (>50 recommended)
- Relatively uniform spatial distribution
- Frame connection helps by providing stable reference points

### Stage 8: Output Generation

**Input**: Final SMD (or TR for SPT)
**Output**: Multiple output formats

**1. Save Results**:
```matlab
% What happens in smi.SMLM.saveResults()
SMD = obj.SMD;
SMF = obj.SMF.packageSMF();  % Convert class back to structure
save('Results.mat', 'SMD', 'SMF', '-v7.3');

% For SPT
save('Results.mat', 'TR', 'SMD', 'SMF', '-v7.3');
```

**2. Generate Super-Resolution Images**:
```matlab
% Gaussian rendering
SR_zoom = 20;  % 20× magnification
img = smi_vis.GenerateImages.gaussImage(...
    [SMD.X, SMD.Y], SMD.PixelSize, SR_zoom, [SMD.YSize, SMD.XSize]);

% Histogram rendering
img = smi_vis.GenerateImages.histImage(...
    [SMD.X, SMD.Y], SMD.PixelSize, SR_zoom, [SMD.YSize, SMD.XSize]);

% Circle rendering (size = uncertainty)
img = smi_vis.GenerateImages.circleImage(...
    [SMD.X, SMD.Y], [SMD.X_SE, SMD.Y_SE], ...
    SMD.PixelSize, SR_zoom, [SMD.YSize, SMD.XSize]);
```

**3. Generate Diagnostic Plots**:
- Photon histogram
- Precision histograms (X_SE, Y_SE, Z_SE)
- PSF sigma histogram
- Background histogram
- P-value distribution
- Localizations per frame
- Drift trajectories (X vs frame, Y vs frame)

**4. Convert to TR (for SPT)**:
```matlab
% TR is array of SMD, one per trajectory
TR = smi_core.TrackingResults.separateTraj(SMD);

% Access individual trajectories
traj_5 = TR(5);  % 5th trajectory
plot(traj_5.X, traj_5.Y);
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('Trajectory 5 (%d localizations)', length(traj_5.X)));
```

## When to Use Different Workflows

### Use SMLM Workflow When:

- Analyzing super-resolution data (PALM, STORM, PAINT)
- Emitters are mostly stationary (may blink but don't move)
- Goal: High-precision position determination
- Frame connection improves precision by averaging

```matlab
% SMLM workflow
SMF = smi_core.SingleMoleculeFitting();
SMF.FrameConnection.On = true;   % Key: Frame connection
SMF.DriftCorrection.On = true;
SMF.Tracking.Method = '';         % No tracking

SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();
```

### Use SPT Workflow When:

- Tracking moving particles
- Emitters move between frames (diffusion, active transport)
- Goal: Trajectory analysis, diffusion estimation
- Need to handle blinking and missed detections

```matlab
% SPT workflow
SMF = smi_core.SingleMoleculeFitting();
SMF.FrameConnection.On = false;  % Key: No frame connection
SMF.DriftCorrection.On = true;   % Still want drift correction
SMF.Tracking.D = 0.1;            % Configure tracking
SMF.Tracking.MaxDistFF = 5;
SMF.Tracking.MaxFrameGap = 10;

SPTobj = smi.SPT(SMF);
[TR, SMD, ~] = SPTobj.performFullAnalysis();
```

### Use Custom Workflow When:

- Need specialized processing
- Developing new methods
- Unusual data characteristics
- Want to inspect intermediate results

```matlab
% Custom modular workflow
SMF = smi_core.SingleMoleculeFitting();

% Load manually
LD = smi_core.LoadData();
[~, seq, SMF] = LD.loadRawData(SMF, 1);

% Convert to photons
photons = (seq - SMF.Data.CameraOffset) / SMF.Data.CameraGain;

% Localize
LDObj = smi_core.LocalizeData(photons, SMF);
SMD = LDObj.genLocalizations();

% Custom filtering (instead of standard thresholding)
good_idx = SMD.X_SE < 0.1 & SMD.Photons > 500;
SMD_filtered = smi_core.SingleMoleculeData.filterSMD(SMD, good_idx);

% Continue with custom processing...
```

## Complete Workflow Examples

### Example 1: Simple SMLM Analysis

```matlab
% Configure parameters
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/experiment1';
SMF.Data.FileName = {'Cell1_PAINT.h5'};
SMF.Data.PixelSize = 0.108;
SMF.Data.FrameRate = 100;
SMF.Data.CameraGain = 2.5;
SMF.Data.CameraOffset = 100;

% Box finding and fitting
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 250;
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Quality control
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

% Results saved to FileDir/Results/Cell1_PAINT_Results.mat
```

### Example 2: Multi-Dataset SPT Analysis

```matlab
% Configure for tracking
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/tracking';
SMF.Data.FileName = {'particle_tracking.h5'};
SMF.Data.DatasetList = int32(1:10);  % Analyze 10 datasets
SMF.Data.PixelSize = 0.16;
SMF.Data.FrameRate = 10;

% Detection and fitting
SMF.BoxFinding.MinPhotons = 100;  % Lower for tracking
SMF.Fitting.PSFSigma = 1.5;
SMF.Fitting.FitType = 'XYNB';

% Tracking parameters
SMF.Tracking.D = 0.05;              % Slow diffusion
SMF.Tracking.MaxDistFF = 3;         % 3 pixels per frame
SMF.Tracking.MaxDistGC = 8;         % Bridge larger gaps
SMF.Tracking.MaxFrameGap = 5;
SMF.Tracking.MinTrackLength = 15;   % Keep long tracks
SMF.Tracking.TrajwiseD = true;      # Estimate D per trajectory

% Drift correction (important for long acquisitions)
SMF.DriftCorrection.On = true;

% Run tracking
SPTobj = smi.SPT(SMF);
[TR, SMD, ~] = SPTobj.performFullAnalysis();

fprintf('Found %d trajectories\n', length(TR));

% Analyze first trajectory
traj1 = TR(1);
figure;
plot(traj1.X, traj1.Y, 'o-');
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('Trajectory 1: %d frames', length(traj1.X)));
```

### Example 3: Inspecting Data Flow

```matlab
% Step through pipeline to inspect each stage
SMF = smi_core.SingleMoleculeFitting();
% ... configure SMF ...

% 1. Load raw data
LD = smi_core.LoadData();
[~, sequence_ADU, SMF] = LD.loadRawData(SMF, 1);
fprintf('Loaded: %d x %d x %d\n', size(sequence_ADU));
figure; imagesc(sequence_ADU(:,:,1)); title('Raw ADU'); colorbar;

% 2. Convert to photons
photons = (sequence_ADU - SMF.Data.CameraOffset) / SMF.Data.CameraGain;
fprintf('Photon range: %.1f to %.1f\n', min(photons(:)), max(photons(:)));
figure; imagesc(photons(:,:,1)); title('Photons'); colorbar;

% 3. Localize
LDObj = smi_core.LocalizeData(photons, SMF);
LDObj.Verbose = 2;  % Show intermediate plots
[SMD, SMDPreThresh] = LDObj.genLocalizations();
fprintf('Found %d localizations\n', length(SMD.X));

% 4. Check thresholding
passed = sum(SMDPreThresh.ThreshFlag == 0);
total = length(SMDPreThresh.ThreshFlag);
fprintf('Passed thresholding: %d / %d (%.1f%%)\n', ...
    passed, total, 100*passed/total);

% 5. Frame connection
if SMF.FrameConnection.On
    FC = smi_core.FrameConnection(SMD, SMF);
    SMD_FC = FC.performFrameConnection();
    compression = length(SMD.X) / sum(SMD_FC.ConnectID > 0);
    fprintf('Frame connection compression: %.1f:1\n', compression);
    SMD = SMD_FC;
end

% 6. Drift correction
if SMF.DriftCorrection.On
    DC = smi_core.DriftCorrection(SMF);
    SMD = DC.driftCorrectKNNIntra(SMD, 1, 1);
    figure;
    plot(SMD.DriftX(:,1), SMD.DriftY(:,1), 'o-');
    xlabel('Drift X'); ylabel('Drift Y'); title('Drift Trajectory');
end

% 7. Generate SR image
SR_img = smi_vis.GenerateImages.gaussImage(...
    [SMD.X, SMD.Y], SMF.Data.PixelSize, 20, [SMD.YSize, SMD.XSize]);
figure; imagesc(SR_img); axis image; colormap hot;
title('Super-Resolution Image (20x)');
```

## See Also

- [Architecture Overview](architecture.md) - smite's organizational structure
- [SMF Structure](smf-structure.md) - Complete parameter reference
- [SMD Structure](smd-structure.md) - Results data format
- [SMLM Workflow](../workflows/smlm-analysis.md) - SMLM-specific details
- [SPT Workflow](../workflows/spt-tracking.md) - Tracking-specific details
- [How to Localize Molecules](../how-to/localize-molecules.md) - Fitting details
- [First Analysis Tutorial](../getting-started/first-analysis.md) - Hands-on practice
