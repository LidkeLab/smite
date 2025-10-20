---
title: "Complete SMLM Analysis Workflow"
category: "workflows"
level: "intermediate"
tags: ["smlm", "workflow", "pipeline", "super-resolution"]
prerequisites: ["../core-concepts/smf-structure.md", "../core-concepts/smd-structure.md"]
related: ["../getting-started/first-analysis.md", "../how-to/localize-molecules.md"]
summary: "Comprehensive guide to the complete SMLM analysis pipeline from raw data to super-resolution images"
estimated_time: "25 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Complete SMLM Analysis Workflow

## Purpose

This document provides a comprehensive understanding of smite's SMLM (Single Molecule Localization Microscopy) analysis pipeline. You'll learn how data flows through each processing stage, what happens at each step, how to control the process, and how to diagnose issues. This is essential for optimizing analyses and understanding your results.

## Prerequisites

- Understanding of [SMF structure](../core-concepts/smf-structure.md)
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Completion of [first analysis](../getting-started/first-analysis.md)
- Basic SMLM concepts (blinking, localization, precision)

## Overview

The SMLM workflow transforms raw camera images into super-resolution structures by:

1. **Loading** raw data from files
2. **Converting** camera ADU to photons
3. **Finding** candidate molecule locations
4. **Fitting** PSF models to estimate precise positions
5. **Thresholding** to remove poor quality localizations
6. **Connecting** the same molecule across frames
7. **Correcting** for stage drift
8. **Generating** super-resolution images and statistics

Each step is controlled by SMF parameters and produces or modifies the SMD results structure. The pipeline is designed to be both automated (via `smi.SMLM`) and modular (using individual components).

## Pipeline Architecture

### Data Flow Diagram

```
┌─────────────────┐
│ Raw Data (.h5)  │ Camera ADU values
└────────┬────────┘
         │ LoadData
         ▼
┌─────────────────┐
│  Image Stack    │ ADU (analog-to-digital units)
└────────┬────────┘
         │ DataToPhotons (gain/offset correction)
         ▼
┌─────────────────┐
│  Image Stack    │ Photons
└────────┬────────┘
         │ LocalizeData
         ├─► FindROI (box finding)
         ├─► GaussMLE (fitting)
         └─► ThresholdFits (quality filter)
         ▼
┌─────────────────┐
│  SMD (raw)      │ All localizations
└────────┬────────┘
         │ FrameConnection
         ▼
┌─────────────────┐
│  SMD (FC)       │ Connected localizations
└────────┬────────┘
         │ DriftCorrection (intra-dataset)
         ▼
┌─────────────────┐
│  SMD (DC)       │ Drift-corrected per dataset
└────────┬────────┘
         │ Concatenate datasets
         │ DriftCorrection (inter-dataset)
         ▼
┌─────────────────┐
│  SMD (final)    │ Complete results
└────────┬────────┘
         │ GenerateImages/Plots
         ▼
┌─────────────────┐
│ SR Images/Plots │ Visualizations
└─────────────────┘
```

### Control Flow

The `smi.SMLM` class orchestrates this pipeline:

```matlab
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();  % Runs complete pipeline
```

Internally, `fullAnalysis()` calls:
```
fullAnalysis()
  └─► analyzeAll()
       ├─► for each dataset in DatasetList
       │    └─► analyzeDataset()
       │         ├─► LoadData.loadRawData()
       │         ├─► DataToPhotons.convert()
       │         ├─► LocalizeData.genLocalizations()
       │         ├─► FrameConnection.connectNN()
       │         └─► DriftCorrection.correct()
       ├─► Concatenate all SMD
       └─► DriftCorrection.correct() [inter-dataset]
```

## Step-by-Step Analysis

### Step 1: Loading Raw Data

**What happens:**
- Reads .h5 or .mat file from disk
- Extracts specified datasets
- Applies ROI cropping if specified
- Returns image stack as 3D array (Y × X × Frames)

**Controlled by:**
```matlab
SMF.Data.FileDir      % Directory containing file
SMF.Data.FileName     % File name (cell array)
SMF.Data.FileType     % 'h5' or 'mat'
SMF.Data.DataVariable % Variable name for .mat files
SMF.Data.DatasetList  % Which datasets to load
SMF.Data.DataROI      % [YStart, XStart, YEnd, XEnd]
```

**Example:**
```matlab
LD = smi_core.LoadData();
[~, sequence, SMF] = LD.loadRawData(SMF, 1);  % Load dataset 1
% sequence is now Y × X × NFrames array
```

**Diagnostics:**
```matlab
% Check data dimensions
fprintf('Loaded: %d × %d × %d (Y × X × Frames)\n', size(sequence));

% View first frame
imagesc(sequence(:,:,1));
axis image; colorbar;
title('First Frame (ADU)');
```

### Step 2: Converting to Photons

**What happens:**
- Subtracts camera offset
- Divides by camera gain
- Converts ADU → photons
- Handles EMCCD vs sCMOS differences

**Controlled by:**
```matlab
SMF.Data.CameraType   % 'EMCCD' or 'SCMOS'
SMF.Data.CameraGain   % Scalar or image (ADU/photon)
SMF.Data.CameraOffset % Scalar or image (ADU)
```

**Example:**
```matlab
% For EMCCD with scalar calibration
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;   % ADU per photon
SMF.Data.CameraOffset = 100; % ADU

% Conversion happens internally:
photons = (sequence - CameraOffset) / CameraGain;
```

**For sCMOS** (pixel-wise calibration):
```matlab
% Load calibration file
load(SMF.Data.CalibrationFilePath, 'CameraGain', 'CameraOffset', 'CameraNoise');
% CameraGain, CameraOffset, CameraNoise are images (Y × X)

% Conversion:
photons = (sequence - CameraOffset) ./ CameraGain;
% Division is element-wise per pixel
```

### Step 3: Finding Molecules (Box Finding)

**What happens:**
- Scans each frame for bright spots
- Estimates local photon count
- Creates boxes around candidates exceeding threshold
- Returns box coordinates for fitting

**Controlled by:**
```matlab
SMF.BoxFinding.BoxSize      % Box size (integer, pixels, typically 7-10)
SMF.BoxFinding.MinPhotons   % Detection threshold (photons)
SMF.BoxFinding.BoxOverlap   % Allowed overlap (pixels)
```

**Algorithm:**
1. Smooth image with Gaussian filter
2. Find local maxima
3. Estimate photon count in local region
4. Keep candidates with photons > MinPhotons
5. Place BoxSize × BoxSize box centered on each candidate
6. Handle overlapping boxes based on BoxOverlap

**Example:**
```matlab
% Manual box finding for one frame
FR = smi_core.FindROI();
FR.BoxFinding = SMF.BoxFinding;
[boxes, ~] = FR.findROI(photons(:,:,1));
% boxes is N × 2 array of [Y, X] box centers

fprintf('Found %d boxes in frame 1\n', size(boxes, 1));
```

**Visualization:**
```matlab
% Show detected boxes
LDObj = smi_core.LocalizeData(photons, SMF);
LDObj.Verbose = 2;  % Shows box overlays
SMD = LDObj.genLocalizations();
```

### Step 4: Fitting PSF Models

**What happens:**
- Extracts image data from each box
- Fits 2D Gaussian PSF model using maximum likelihood
- Estimates: X, Y, photons, background (and optionally sigma, Z)
- Computes uncertainties (standard errors) via Cramér-Rao bound

**Controlled by:**
```matlab
SMF.Fitting.PSFSigma    % Initial PSF width (pixels)
SMF.Fitting.FitType     % 'XYNB', 'XYNBS', 'XYZNB', etc.
SMF.Fitting.Iterations  % Newton-Raphson iterations
```

**Fitting methods:**
- **XYNB**: Fit X, Y, photons (N), background (B) with fixed sigma
- **XYNBS**: Also fit PSF sigma
- **XYZNB**: 3D fitting using astigmatism

**Algorithm (simplified):**
1. Initialize parameters: X₀, Y₀, N₀, B₀
2. For each iteration:
   - Compute expected image from current parameters
   - Compute likelihood gradient and Hessian
   - Update parameters via Newton-Raphson
3. Compute CRLB for uncertainties
4. Return fitted parameters and standard errors

**Example:**
```matlab
% Extract one box and fit manually
box_data = photons(y:y+6, x:x+6, frame_idx);  % 7×7 box

GM = smi_core.GaussMLE();
GM.Fitting = SMF.Fitting;
[params, CRLB, LogL] = GM.gaussMLE(box_data);
% params = [X, Y, Photons, Bg] relative to box corner
% CRLB = [X_SE, Y_SE, Photons_SE, Bg_SE]
```

**GPU requirement:**
Fitting requires CUDA GPU (NVIDIA with compute capability ≥5.0) for GaussMLE operations.

### Step 5: Thresholding Localizations

**What happens:**
- Filters localizations by quality metrics
- Removes poor fits, outliers, and noise
- Sets `SMD.ThreshFlag` for rejected localizations

**Controlled by:**
```matlab
SMF.Thresholding.On              % Enable/disable
SMF.Thresholding.MaxXY_SE        % Max position uncertainty
SMF.Thresholding.MinPhotons      % Min photon count
SMF.Thresholding.MinPValue       % Min fit p-value
SMF.Thresholding.MinPSFSigma     % Min PSF width
SMF.Thresholding.MaxPSFSigma     % Max PSF width
SMF.Thresholding.AutoThreshLogL  % Auto log-likelihood threshold
```

**Filters applied:**
1. **Precision filter**: `X_SE < MaxXY_SE` and `Y_SE < MaxXY_SE`
2. **Photon filter**: `Photons > MinPhotons`
3. **Fit quality**: `PValue > MinPValue` (or auto log-likelihood)
4. **PSF range**: `MinPSFSigma < PSFSigma < MaxPSFSigma`
5. **Background**: `Bg < MaxBg`
6. **Nearest neighbor**: Reject isolated outliers

**Example:**
```matlab
% Check what was filtered
load('Results.mat', 'SMD');

% ThreshFlag meanings:
% 0  = passed
% 1  = failed precision
% 2  = failed photons
% 4  = failed p-value
% 8  = failed PSF sigma
% ... (bitwise flags can combine)

flags = SMD.ThreshFlag;
passed = sum(flags == 0);
fprintf('%d / %d (%.1f%%) passed thresholding\n', ...
    passed, length(flags), 100*passed/length(flags));

% See what failed precision filter
failed_precision = bitand(flags, 1) > 0;
fprintf('%d failed precision filter\n', sum(failed_precision));
```

### Step 6: Frame Connection

**What happens:**
- Identifies the same emitter appearing in multiple frames
- Assigns `ConnectID` to link related localizations
- Improves overall precision by combining measurements

**Controlled by:**
```matlab
SMF.FrameConnection.On              % Enable/disable
SMF.FrameConnection.Method          % 'LAP-FC' (Linear Assignment Problem)
SMF.FrameConnection.MaxSeparation   % Max distance (pixels)
SMF.FrameConnection.MaxFrameGap     % Max frame gap to bridge
SMF.FrameConnection.MinNFrameConns  % Min appearances to keep
```

**Algorithm:**
1. For consecutive frames, compute distance matrix between localizations
2. Solve linear assignment problem to find optimal pairing
3. Accept pairs with distance < MaxSeparation
4. Bridge gaps up to MaxFrameGap frames
5. Assign unique ConnectID to each cluster
6. Filter clusters with < MinNFrameConns appearances

**Example:**
```matlab
% Manual frame connection
FC = smi_core.FrameConnection();
FC.FrameConnection = SMF.FrameConnection;
SMD_connected = FC.connectNN(SMD_raw, SMF);

% Check compression ratio
unique_IDs = unique(SMD_connected.ConnectID(SMD_connected.ConnectID > 0));
compression = length(SMD_connected.X) / length(unique_IDs);
fprintf('Frame connection: %.1f localizations → 1 emitter\n', compression);
```

**Benefits:**
- Improved precision (average multiple observations)
- Reduced data size
- Better photon statistics
- Identification of blinking emitters

### Step 7: Drift Correction

**What happens:**
- Estimates stage drift from localization data
- Applies correction to align all frames
- Handles both intra-dataset and inter-dataset drift

**Controlled by:**
```matlab
SMF.DriftCorrection.On          % Enable/disable
SMF.DriftCorrection.Method      % 'DC-KNN' (K-nearest neighbors)
SMF.DriftCorrection.L_intra     % Intra-dataset smoothing (pixels)
SMF.DriftCorrection.L_inter     % Inter-dataset smoothing (pixels)
SMF.DriftCorrection.PDegree     % Polynomial degree
```

**Algorithm (DC-KNN):**
1. For each frame, find K nearest neighbors in adjacent frames
2. Compute median displacement
3. Build drift trajectory
4. Smooth with regularization parameter L
5. Fit polynomial of degree PDegree
6. Subtract drift from positions

**Two-stage correction:**
1. **Intra-dataset**: Correct drift within each dataset separately
2. **Inter-dataset**: Align all datasets to common reference

**Example:**
```matlab
% Manual drift correction
DC = smi_core.DriftCorrection();
DC.DriftCorrection = SMF.DriftCorrection;
SMD_corrected = DC.correct(SMD_raw, SMF);

% Examine drift
figure;
subplot(2,1,1);
plot(SMD_corrected.DriftX(:,1) * SMD.PixelSize * 1000);
ylabel('X Drift (nm)'); xlabel('Frame');

subplot(2,1,2);
plot(SMD_corrected.DriftY(:,1) * SMD.PixelSize * 1000);
ylabel('Y Drift (nm)'); xlabel('Frame');
```

**Drift in results:**
- `SMD.DriftX`, `SMD.DriftY`: Estimated drift per frame
- `SMD.X`, `SMD.Y`: Already corrected positions
- To get raw positions: `X_raw = SMD.X - SMD.DriftX(SMD.FrameNum, SMD.DatasetNum)`

### Step 8: Generating Visualizations

**What happens:**
- Creates super-resolution images
- Generates diagnostic histograms
- Produces drift plots and statistics
- Saves everything to Results directory

**Controlled by:**
```matlab
SMLMobj.PlotDo       % Which plots to make
SMLMobj.SRImageZoom  % SR image magnification
```

**Image types:**

1. **GaussIm**: Each localization rendered as 2D Gaussian
   ```matlab
   SR_zoom = 20;
   img = smi_vis.GenerateImages.gaussImage([SMD.X, SMD.Y], ...
       SMD.PixelSize, SR_zoom, [SMD.YSize, SMD.XSize]);
   ```

2. **HistIm**: Histogram image (bin counts)
   ```matlab
   img = smi_vis.GenerateImages.histImage([SMD.X, SMD.Y], ...
       SMD.PixelSize, SR_zoom, [SMD.YSize, SMD.XSize]);
   ```

3. **CircleIm**: Circles with size = uncertainty
   ```matlab
   SR_zoom = 25;
   img = smi_vis.GenerateImages.circleImage([SMD.X, SMD.Y], ...
       [SMD.X_SE, SMD.Y_SE], SMD.PixelSize, SR_zoom, [SMD.YSize, SMD.XSize]);
   ```

**Diagnostic plots:**
- Photon histogram
- Background histogram
- PSF sigma histogram
- X/Y/Z precision histograms
- P-value histogram
- Localizations per frame
- Drift trajectories

## Complete Workflow Examples

### Example 1: GUI-Based Analysis

```matlab
% Create SMLM object with GUI
SMLMobj = smi.SMLM();

% Use GUI to set:
% - FileDir, FileName
% - PixelSize, FrameRate
% - CameraGain, CameraOffset
% - BoxFinding parameters
% - Fitting parameters
% - Thresholding parameters

% Test on one dataset
SMLMobj.testFit();  % Verify parameters

% Run full analysis
SMLMobj.fullAnalysis();

% Results saved to FileDir/Results/
```

### Example 2: Script-Based Analysis

```matlab
% Configure SMF
SMF = smi_core.SingleMoleculeFitting();

% Data
SMF.Data.FileDir = '/data/2024-01-10';
SMF.Data.FileName = {'Cell1_PAINT.h5'};
SMF.Data.PixelSize = 0.108;   % 108 nm
SMF.Data.FrameRate = 100;     % 100 Hz
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

% Access results
SMD = SMLMobj.SMD;
fprintf('Detected %d localizations\n', length(SMD.X));
```

### Example 3: Custom Pipeline (Advanced)

```matlab
% Load data manually
LD = smi_core.LoadData();
[~, sequence, SMF] = LD.loadRawData(SMF, 1);

% Convert to photons
photons = (sequence - SMF.Data.CameraOffset) / SMF.Data.CameraGain;

% Localize
LDObj = smi_core.LocalizeData(photons, SMF);
LDObj.Verbose = 1;
SMD_raw = LDObj.genLocalizations();

% Custom filtering (instead of standard thresholding)
good = SMD_raw.X_SE < 0.1 & SMD_raw.Photons > 500 & SMD_raw.PValue > 0.05;
SMD_filtered = smi_core.SingleMoleculeData.filterSMD(SMD_raw, good);

% Frame connection
FC = smi_core.FrameConnection();
SMD_connected = FC.connectNN(SMD_filtered, SMF);

% Drift correction
DC = smi_core.DriftCorrection();
SMD_final = DC.correct(SMD_connected, SMF);

% Generate custom image
SR_img = smi_vis.GenerateImages.gaussImage([SMD_final.X, SMD_final.Y], ...
    SMF.Data.PixelSize, 30, [SMD_final.YSize, SMD_final.XSize]);

% Save
save('custom_results.mat', 'SMD_final', 'SMF');
```

## Troubleshooting

### Issue: No localizations found

**Diagnose:**
```matlab
% Check raw data
imagesc(sequence(:,:,1)); colorbar;
title('Is data visible?');

% Check after photon conversion
photons = (sequence - SMF.Data.CameraOffset) / SMF.Data.CameraGain;
imagesc(photons(:,:,1)); colorbar;
title('Photons (should be positive)');
```

**Solutions:**
- Verify CameraGain and CameraOffset are correct
- Lower BoxFinding.MinPhotons
- Check that data file has expected variable name

### Issue: Poor localization precision

**Diagnose:**
```matlab
histogram(SMD.X_SE * SMD.PixelSize * 1000);  % nm
xlabel('X Precision (nm)');
```

**Solutions:**
- Increase BoxFinding.MinPhotons (select brighter emitters)
- Check PSFSigma matches data
- Verify camera calibration is correct
- Check for focus drift (fit PSF sigma)

### Issue: Drift correction fails

**Diagnose:**
```matlab
plot(SMD.DriftX(:,1), SMD.DriftY(:,1));
xlabel('Drift X'); ylabel('Drift Y');
title('Does this look reasonable?');
```

**Solutions:**
- Ensure enough localizations per frame (>50)
- Adjust DriftCorrection.L_intra (try 0.5 or 2.0)
- Check for emitter bleaching over time
- Verify frame connection worked (need stable reference points)

### Issue: Frame connection too aggressive or too conservative

**Diagnose:**
```matlab
unique_IDs = unique(SMD.ConnectID(SMD.ConnectID > 0));
compression = length(SMD.X) / length(unique_IDs);
fprintf('Compression: %.1f:1\n', compression);
% Good range: 2:1 to 10:1 for typical SMLM
```

**Solutions:**
- Too aggressive (compression >20:1): Decrease MaxSeparation
- Too conservative (compression <2:1): Increase MaxSeparation
- Adjust MaxFrameGap based on blinking statistics

## Performance Optimization

### GPU Requirement

Verify CUDA GPU is available:
```matlab
gpuDevice  % Should show your NVIDIA GPU
```

GPU is required for fitting operations. Ensure compute capability ≥5.0.

### Parallel Processing

For batch jobs:
```matlab
parpool('local', 4);  % Start 4 workers
% Then run batch analysis
```

### Reduce Data Size

Analyze subset for testing:
```matlab
SMF.Data.DatasetList = int32(1);  % Just first dataset
SMF.Data.DataROI = [1, 1, 64, 64];  % Small ROI
```

## See Also

- [First Analysis](../getting-started/first-analysis.md) - Hands-on tutorial
- [SMF Structure](../core-concepts/smf-structure.md) - All parameters
- [SMD Structure](../core-concepts/smd-structure.md) - Results format
- [How to Load Data](../how-to/load-data.md) - Data loading details
- [How to Localize Molecules](../how-to/localize-molecules.md) - Fitting details
- MATLAB/+smi/@SMLM/README.md - Class documentation
