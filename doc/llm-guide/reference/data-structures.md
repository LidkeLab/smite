---
title: "Data Structures Reference"
category: "reference"
level: "intermediate"
tags: ["smf", "smd", "tr", "data-structures", "reference", "parameters", "results"]
prerequisites: ["../core-concepts/architecture.md"]
related: ["../core-concepts/smf-structure.md", "../core-concepts/smd-structure.md", "../core-concepts/tr-structure.md"]
summary: "Complete reference guide for all smite data structures: SMF parameters, SMD results, and TR trajectories"
estimated_time: "30 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# Data Structures Reference

## Purpose

This document serves as a comprehensive reference for all data structures in smite. Use it for quick lookup of field names, types, units, and default values. For detailed conceptual explanations and usage examples, see the related concept documents linked throughout.

## Prerequisites

- Basic understanding of [smite architecture](../core-concepts/architecture.md)
- Familiarity with MATLAB structures
- Knowledge of SMLM/SPT concepts

## Overview

smite is built around three fundamental data structures:

| Structure | Purpose | Created By | Used For |
|-----------|---------|------------|----------|
| **SMF** | Analysis parameters | User/GUI | Configuring pipeline |
| **SMD** | Localization results | Analysis workflows | All results, visualization |
| **TR** | Trajectory results | Tracking workflows | SPT analysis |

**Key relationships:**
- SMF defines what to do → SMD contains what was found
- TR is an array of SMD structures, organized by trajectory
- Same SMF + same data = same SMD (reproducibility)

## Quick Reference Tables

### SMF Sub-Structures

| Sub-Structure | Purpose | Key Fields |
|---------------|---------|------------|
| `Data` | File and camera info | FileName, PixelSize, FrameRate, CameraGain |
| `BoxFinding` | Molecule detection | BoxSize, MinPhotons, BoxOverlap |
| `Fitting` | PSF model and fitting | PSFSigma, FitType, Iterations |
| `Thresholding` | Quality filtering | MaxXY_SE, MinPhotons, MinPValue |
| `FrameConnection` | Multi-frame linking | MaxSeparation, MaxFrameGap, Method |
| `DriftCorrection` | Stage drift removal | Method, L_intra, L_inter |
| `Tracking` | Particle tracking | D, MaxDistFF, MaxFrameGap, MinTrackLength |

### SMD Field Categories

| Category | Fields | Units |
|----------|--------|-------|
| **Position** | X, Y, Z | pixels, micrometers (Z) |
| **Uncertainty** | X_SE, Y_SE, Z_SE | pixels, micrometers (Z) |
| **Photometry** | Photons, Bg | photons, photons/pixel |
| **PSF** | PSFSigma, PSFSigmaX, PSFSigmaY | pixels |
| **Temporal** | FrameNum, DatasetNum | integer indices |
| **Quality** | PValue, LogLikelihood, ThreshFlag | — |
| **Connection** | ConnectID | integer ID |
| **Drift** | DriftX, DriftY, DriftZ | pixels, micrometers (Z) |

### TR Organization

| Aspect | Description |
|--------|-------------|
| **Structure** | Array of SMD structures |
| **Indexing** | `TR(n)` = n-th trajectory |
| **Fields** | Same as SMD, all vectors |
| **ConnectID** | All points in `TR(n)` have same ID |

---

## SMF: Single Molecule Fitting

**Type:** Structure of structures
**Purpose:** Complete specification of analysis pipeline
**Created:** `SMF = smi_core.SingleMoleculeFitting()`

See [SMF Structure concept doc](../core-concepts/smf-structure.md) for detailed usage.

### SMF.Data - File and Camera Information

Specifies input data location and camera characteristics.

| Field | Type | Default | Units | Description |
|-------|------|---------|-------|-------------|
| `FileName` | cell array of char | `{}` | — | File name(s) to analyze |
| `FileDir` | char array | `''` | — | Directory containing files |
| `ResultsDir` | char array | `'FileDir/Results'` | — | Output directory |
| `AnalysisID` | char array | `''` | — | ID appended to results files |
| `FileType` | char array | auto-detect | — | 'h5', 'mat', or custom |
| `DataVariable` | char array | `'sequence'` | — | Variable name in .mat files |
| `DatasetList` | int32 array | `[]` | — | Dataset indices to analyze (empty = all) |
| `DatasetMods` | cell {incl; excl} | `{[]; []}` | — | Include/exclude dataset modifiers |
| `CameraType` | char | `'EMCCD'` | — | 'EMCCD' or 'SCMOS' |
| `CameraGain` | scalar or image | `1` | ADU/photon | Camera conversion gain |
| `CameraOffset` | scalar or image | `0` | ADU | Camera baseline offset |
| `CameraNoise` | scalar or image | `0` | electrons RMS | Read noise standard deviation |
| `CalibrationFilePath` | char | `''` | — | Path to SCMOS calibration file |
| `RegistrationFilePath` | char | `''` | — | Path to channel registration file |
| `DataROI` | 4-element vector | `[]` | pixels | [YStart, XStart, YEnd, XEnd] |
| `FrameRate` | scalar | required | Hz | Acquisition frame rate |
| `PixelSize` | scalar | required | micrometers | Camera pixel size (back-projected) |
| `SEAdjust` | scalar | `0` | pixels | Standard error inflation factor |

**Key notes:**
- `PixelSize` and `FrameRate` must be set by user
- For SCMOS cameras, use `CalibrationFilePath` instead of scalar gain/offset/noise
- `DatasetMods{1}` = datasets to include, `DatasetMods{2}` = datasets to exclude
- `DataROI` crops data before analysis: [row_start, col_start, row_end, col_end]

**Typical values:**
```matlab
SMF.Data.PixelSize = 0.108;      % 108 nm (common for 100x + 1.6x optics)
SMF.Data.FrameRate = 100;        % 100 Hz
SMF.Data.CameraGain = 2.5;       % Typical EMCCD gain
SMF.Data.CameraOffset = 100;     % Typical EMCCD offset
```

### SMF.BoxFinding - Molecule Detection

Parameters for identifying regions of interest (ROIs) containing molecules.

| Field | Type | Default | Units | Description |
|-------|------|---------|-------|-------------|
| `BoxSize` | integer | `7` | pixels | Side length of square fitting region |
| `BoxOverlap` | integer | `2` | pixels | Allowed overlap between boxes |
| `MinPhotons` | scalar | `200` | photons | Minimum signal to trigger detection |

**Guidelines:**
- **BoxSize**: Should encompass ~4-5σ of PSF. Typical: 5-10 pixels.
- **MinPhotons**: Trade-off between sensitivity (low) and false positives (high). Typical: 100-500 for SMLM, 50-200 for SPT.
- **BoxOverlap**: 0 = no overlap (conservative), 2-3 = allow nearby emitters (aggressive).

**Choosing BoxSize:**
```matlab
% Rule of thumb: BoxSize ≈ 4-5 × PSFSigma
PSFSigma = 1.3;  % pixels
BoxSize = ceil(4.5 * PSFSigma);  % ~6 pixels
```

### SMF.Fitting - PSF Model and Method

Controls how molecule positions are estimated.

| Field | Type | Default | Units | Description |
|-------|------|---------|-------|-------------|
| `PSFSigma` | scalar | `1` | pixels | Initial/fixed PSF sigma (symmetric) |
| `FitType` | char | `'XYNB'` | — | Which parameters to fit |
| `NParams` | integer | auto-set | — | Number of fit parameters |
| `Iterations` | integer | `20` | — | Newton-Raphson iterations |
| `ZFitStruct` | structure | — | — | 3D astigmatism calibration |

**FitType options:**

| FitType | Parameters | NParams | Use Case |
|---------|------------|---------|----------|
| `'XYNB'` | X, Y, Photons, Background | 4 | Standard 2D SMLM |
| `'XYNBS'` | X, Y, N, B, Sigma | 5 | Variable PSF width |
| `'XYNBSXSY'` | X, Y, N, B, SigmaX, SigmaY | 6 | Asymmetric PSF |
| `'XYZNB'` | X, Y, Z, N, B | 5 | 3D astigmatism |

**Determining PSFSigma:**
```matlab
% Theoretical: σ = 0.21 × λ / (NA × pixel_size)
% For λ=670nm, NA=1.49, pixel=108nm: σ ≈ 1.3 pixels

% Empirical: Fit some data with FitType='XYNBS', then:
sigma_empirical = median(SMD.PSFSigma);
```

**ZFitStruct fields (3D only):**
- `Ax`, `Bx`: X-sigma calibration curve parameters
- `Ay`, `By`: Y-sigma calibration curve parameters
- `Gamma`, `D`: Additional calibration parameters

### SMF.Thresholding - Quality Filtering

Removes poor quality localizations based on fit quality metrics.

| Field | Type | Default | Units | Description |
|-------|------|---------|-------|-------------|
| `On` | logical | `true` | — | Enable thresholding? |
| `MaxXY_SE` | scalar | `0.2` | pixels | Max XY localization precision |
| `MaxZ_SE` | scalar | `0.05` | micrometers | Max Z localization precision |
| `MinPValue` | scalar | `0.01` | — | Min fit p-value (goodness of fit) |
| `AutoThreshLogL` | logical | `false` | — | Auto-threshold on log-likelihood |
| `AutoThreshPrctile` | scalar | `1e-4` | — | Percentile for auto-threshold |
| `MinPSFSigma` | scalar | `0.5` | pixels | Min acceptable PSF sigma |
| `MaxPSFSigma` | scalar | `2.0` | pixels | Max acceptable PSF sigma |
| `MinPhotons` | scalar | `100` | photons | Min photons after fitting |
| `MaxBg` | scalar | `Inf` | photons/pixel | Max background level |
| `InMeanMultiplier` | scalar | `Inf` | — | Max intensity acceptance factor |
| `NNMedianMultiplier` | scalar | `3` | — | Nearest-neighbor outlier factor |
| `MinNumNeighbors` | scalar | `0` | — | Min required neighbors |

**Most important thresholds:**
```matlab
SMF.Thresholding.MaxXY_SE = 0.2;      % 20 nm for 100 nm pixels
SMF.Thresholding.MinPhotons = 100;    % Remove dim fits
SMF.Thresholding.MinPValue = 0.01;    % Remove poor fits
SMF.Thresholding.MinPSFSigma = 0.8;   % Prevent too-narrow fits
SMF.Thresholding.MaxPSFSigma = 1.8;   % Prevent too-wide fits
```

**Checking effects:**
After analysis, inspect `SMD.ThreshFlag`:
- `0` = passed all thresholds
- `>0` = failed (specific codes indicate which threshold)

### SMF.FrameConnection - Multi-Frame Linking

Links the same emitter appearing across multiple frames (for blinking SMLM).

| Field | Type | Default | Units | Description |
|-------|------|---------|-------|-------------|
| `On` | logical | `true` | — | Enable frame connection? |
| `Method` | char | `'LAP-FC'` | — | Connection algorithm |
| `MaxSeparation` | scalar | `1.0` | pixels | Max distance between frames |
| `LoS` | scalar | `0.01` | — | Level of significance for connection |
| `MaxFrameGap` | integer | `5` | frames | Max frames to bridge |
| `NSigmaDev` | scalar | `5` | — | SE multiplier for pre-clustering |
| `NNearestClusters` | integer | `2` | — | Clusters for density estimation |
| `NIterations` | integer | `1` | — | Number of iterative FC passes |
| `MinNFrameConns` | integer | `1` | — | Min appearances to retain |

**Purpose:** Improves precision by combining multiple observations of same emitter.

**Key parameters:**
- **MaxSeparation**: Should be larger than typical localization uncertainty but smaller than inter-emitter spacing
- **MaxFrameGap**: How many dark frames to bridge (typical: 5-10 for blinking fluorophores)
- **MinNFrameConns**: Require multiple appearances for stability (0 = no filtering, 2+ = stricter)

**After frame connection:**
- `SMD.ConnectID` identifies same emitter across frames
- Connected localizations can be averaged for improved precision

### SMF.DriftCorrection - Stage Drift Compensation

Corrects for stage drift during long acquisitions.

| Field | Type | Default | Units | Description |
|-------|------|---------|-------|-------------|
| `On` | logical | `true` | — | Enable drift correction? |
| `Method` | char | `'DC-KNN'` | — | Correction algorithm |
| `BFRegistration` | logical | `true` | — | Use brightfield registration? |
| `L_intra` | scalar | `1.0` | pixels | Intra-dataset smoothing threshold |
| `L_inter` | scalar | `2.0` | pixels | Inter-dataset alignment threshold |
| `PixelSizeZUnit` | scalar | `0.1` | micrometers | XY pixel size for 3D correction |
| `PDegree` | integer | `1` | — | Polynomial degree for drift model |

**Methods:**
- `'DC-KNN'`: K-nearest neighbors (default, robust)
- Other specialized methods available

**Key parameters:**
- **L_intra**: Smoothness within single dataset (lower = smoother)
- **L_inter**: Alignment strength across datasets (lower = tighter alignment)
- **PDegree**: 1 = linear drift, 2 = quadratic (accelerating drift)

**Results:**
After correction, `SMD.DriftX`, `SMD.DriftY`, `SMD.DriftZ` contain cumulative drift per frame. Positions in SMD are already corrected.

### SMF.Tracking - Single Particle Tracking

Parameters for SPT (connecting localizations into trajectories).

| Field | Type | Default | Units | Description |
|-------|------|---------|-------|-------------|
| `Method` | char | `'CostMatrix'` | — | Tracking algorithm |
| `D` | scalar | `0.01` | pixels²/frame | Diffusion constant estimate |
| `TrajwiseD` | logical | `true` | — | Estimate D per trajectory? |
| `K_on` | scalar | `0.9` | frame⁻¹ | Off→On transition rate |
| `K_off` | scalar | `0.1` | frame⁻¹ | On→Off transition rate |
| `MaxDistFF` | scalar | `5` | pixels | Max frame-to-frame distance |
| `MaxDistGC` | scalar | `10` | pixels | Max gap-closing distance |
| `MaxFrameGap` | integer | `10` | frames | Max frames to bridge |
| `MinTrackLength` | integer | `3` | frames | Min trajectory length to keep |
| `NIterMax` | integer | `5` | — | Max iterative tracking attempts |
| `NIterMaxBatch` | integer | `5` | — | Max batch iterations |
| `MaxRelativeChange` | scalar | `1e-5` | — | Convergence threshold |
| `MaxZScoreDist` | scalar | `Inf` | — | Max z-score for position jumps |
| `MaxZScorePhotons` | scalar | `Inf` | — | Max z-score for photon changes |
| `TryLowPValueLocs` | logical | `false` | — | Include low p-value localizations? |

**Converting diffusion constant:**
```matlab
% From physical units to pixels²/frame:
D_pixels = D_physical [μm²/s] / (pixel_size² [μm²] × frame_rate [Hz])

% Example: D=0.1 μm²/s, pixel=0.1 μm, rate=100 Hz
D_pixels = 0.1 / (0.01 × 100) = 0.1 pixels²/frame
```

**Key parameters:**
- **D**: Initial guess for diffusion. Used to predict likely step sizes.
- **MaxDistFF**: Maximum step between consecutive frames
- **MaxDistGC**: Maximum jump when bridging gaps (typically 2× MaxDistFF)
- **MaxFrameGap**: How many missing frames to bridge
- **MinTrackLength**: Filter short, unreliable tracks

---

## SMD: Single Molecule Data

**Type:** Structure (flat, not nested)
**Purpose:** Contains all localization results
**Created:** By `LocalizeData.genLocalizations()` and tracking workflows

See [SMD Structure concept doc](../core-concepts/smd-structure.md) for detailed usage.

### SMD Field Organization

All fields are vectors of length N (number of localizations). Think of SMD as a table where each row is a localization.

```matlab
% After analysis
load('Results.mat', 'SMD');
N_locs = length(SMD.X);  % Number of localizations
```

### Metadata Fields

Dataset-level information (scalars or matrices).

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `NDims` | integer | — | Dimensionality (2 or 3) |
| `NFrames` | integer | — | Total frames in acquisition |
| `NDatasets` | integer | — | Number of image stacks analyzed |
| `FrameRate` | scalar | Hz | Acquisition frame rate |
| `PixelSize` | scalar | micrometers | Camera pixel size (back-projected) |
| `XSize` | integer | pixels | Image width |
| `YSize` | integer | pixels | Image height |
| `ZOffset` | scalar | micrometers | Z position of focal plane |

### Position and Uncertainty Fields

Localization coordinates and precision estimates.

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `X` | vector (N×1) | pixels | X position estimates |
| `Y` | vector (N×1) | pixels | Y position estimates |
| `Z` | vector (N×1) | micrometers | Z position estimates (3D only) |
| `X_SE` | vector (N×1) | pixels | Standard error in X |
| `Y_SE` | vector (N×1) | pixels | Standard error in Y |
| `Z_SE` | vector (N×1) | micrometers | Standard error in Z (3D only) |

**Notes:**
- Origin: (1,1) is center of top-left pixel
- Positions have sub-pixel precision (real numbers)
- Already drift-corrected if `SMF.DriftCorrection.On = true`
- Standard errors from Cramér-Rao lower bound (theoretical precision)

**Typical precision:** 0.05-0.2 pixels (5-20 nm for 100 nm pixels)

### Photometric Fields

Intensity and background information.

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `Photons` | vector (N×1) | photons | Total collected photons |
| `Bg` | vector (N×1) | photons/pixel | Background per pixel |
| `Photons_SE` | vector (N×1) | photons | Standard error in photons |
| `Bg_SE` | vector (N×1) | photons/pixel | Standard error in background |

**Notes:**
- `Photons`: Integrated photons from PSF (after camera gain correction)
- `Bg`: Background per pixel in fitting region
- Signal-to-background ratio: `SBR = Photons / (Bg × PSF_area)`

### PSF Fields

Point spread function width parameters.

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `PSFSigma` | vector (N×1) | pixels | PSF width (symmetric) |
| `PSFSigmaX` | vector (N×1) | pixels | PSF width in X (asymmetric) |
| `PSFSigmaY` | vector (N×1) | pixels | PSF width in Y (asymmetric) |
| `PSFSigma_SE` | vector (N×1) | pixels | Standard error in sigma |
| `PSFSigmaX_SE` | vector (N×1) | pixels | Standard error in sigmaX |
| `PSFSigmaY_SE` | vector (N×1) | pixels | Standard error in sigmaY |

**Notes:**
- Only populated if PSF sigma was fitted (`FitType='XYNBS'` or similar)
- Otherwise, fields are empty or contain fixed value from SMF.Fitting.PSFSigma

### Temporal and Spatial Context Fields

Frame numbers and fitting box locations.

| Field | Type | Description |
|-------|------|-------------|
| `FrameNum` | vector (N×1) | Frame number (1 to NFrames) |
| `DatasetNum` | vector (N×1) | Dataset number (1 to NDatasets) |
| `XBoxCorner` | vector (N×1) | X coordinate of fitting box top-left corner |
| `YBoxCorner` | vector (N×1) | Y coordinate of fitting box top-left corner |

**Uses:**
- Temporal filtering: `early_frames = SMD.FrameNum <= 100`
- Multi-dataset separation: `dataset1 = SMD.DatasetNum == 1`
- Box corners: For reconstructing fitting regions

### Quality and Fit Metrics

Fit quality indicators.

| Field | Type | Description |
|-------|------|-------------|
| `PValue` | vector (N×1) | P-value of fit (0 to 1, higher = better) |
| `LogLikelihood` | vector (N×1) | Log-likelihood of fit (higher = better) |
| `ThreshFlag` | vector (N×1) | Threshold status (0=pass, >0=fail code) |

**PValue:** Statistical goodness-of-fit measure
- Range: 0 to 1
- Typical threshold: 0.01 to 0.05
- Low p-value = poor fit (model doesn't explain data well)

**ThreshFlag:** Indicates why localizations were filtered
- `0`: Passed all thresholds
- `>0`: Failed (specific codes identify which threshold)

### Frame Connection Fields

Multi-frame linking information.

| Field | Type | Description |
|-------|------|-------------|
| `ConnectID` | vector (N×1) | Connection identifier (integer) |
| `IndSMD` | cell array | Indices in pre-frame-connection SMD |

**ConnectID:** Links appearances of same emitter
- Same ID → same molecule across frames
- `0` or negative → not connected
- After frame connection, unique IDs identify distinct emitters

**IndSMD:** Maps frame-connected SMD back to original pre-FC SMD
- Cell array where each element contains indices
- Used for tracing localizations through processing

**Usage:**
```matlab
% Count unique emitters
unique_IDs = unique(SMD.ConnectID);
unique_IDs = unique_IDs(unique_IDs > 0);
N_emitters = length(unique_IDs);

% Average position of one emitter
emitter_1 = SMD.ConnectID == unique_IDs(1);
X_mean = mean(SMD.X(emitter_1));
Y_mean = mean(SMD.Y(emitter_1));
```

### Drift Correction Fields

Stage drift estimates per frame.

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `DriftX` | matrix (NFrames×NDatasets) | pixels | X drift per frame |
| `DriftY` | matrix (NFrames×NDatasets) | pixels | Y drift per frame |
| `DriftZ` | matrix (NFrames×NDatasets) | micrometers | Z drift per frame |

**Notes:**
- Cumulative drift from first frame (frame 1 has zero drift)
- Positions in X, Y, Z are already corrected
- Use these fields to visualize drift trajectory or manually undo correction

**Typical drift magnitude:** 10-100 nm over 1000-10000 frames

### Channel Registration Fields

Multi-color registration status.

| Field | Type | Description |
|-------|------|-------------|
| `IsTransformed` | logical | Was channel registration applied? |
| `RegError` | scalar | Registration error (pixels) |

**IsTransformed:** `true` if coordinate transformation was applied for multi-channel alignment

**RegError:** Estimated registration accuracy (typical: 0.1-0.5 pixels = 10-50 nm)

---

## TR: Tracking Results

**Type:** Array of SMD structures
**Purpose:** Organizes localizations by trajectory
**Created:** `TR = smi_core.TrackingResults.convertSMDToTR(SMD)` or by SPT workflows

See [TR Structure concept doc](../core-concepts/tr-structure.md) for detailed usage.

### TR Organization

```matlab
% TR is an array where each element is one trajectory
load('Results_SPT.mat', 'TR');
N_trajectories = length(TR);

% Each element is an SMD structure
trajectory_1 = TR(1);
positions_1 = [TR(1).X, TR(1).Y];  % All positions for trajectory 1
frames_1 = TR(1).FrameNum;         // Frames where trajectory appears
```

### TR vs SMD Comparison

| Aspect | SMD | TR |
|--------|-----|-----|
| **Type** | Single structure | Array of structures |
| **Organization** | All localizations (flat) | Grouped by trajectory |
| **Indexing** | `SMD.X(i)` = i-th localization | `TR(n).X(i)` = i-th point in n-th trajectory |
| **ConnectID** | Links locs to trajectories | All points in `TR(n)` have same ConnectID |
| **Use case** | Visualization, all-data ops | Per-trajectory analysis, SPT |

### TR Fields

Each `TR(n)` element contains **all standard SMD fields** for that trajectory:

**Position:** `TR(n).X`, `TR(n).Y`, `TR(n).Z`, `TR(n).X_SE`, `TR(n).Y_SE`, `TR(n).Z_SE`

**Photometry:** `TR(n).Photons`, `TR(n).Bg`, `TR(n).PSFSigma`

**Temporal:** `TR(n).FrameNum`, `TR(n).DatasetNum`

**Quality:** `TR(n).PValue`, `TR(n).LogLikelihood`, `TR(n).ThreshFlag`

**Connection:** `TR(n).ConnectID` (all same value for trajectory n)

**Metadata:** `TR(n).NFrames`, `TR(n).PixelSize`, `TR(n).FrameRate` (same across trajectories)

### TR Utility Methods

**Compute statistics:**
```matlab
lengths = smi_core.TrackingResults.computeTrajLengths(TR);
durations = smi_core.TrackingResults.computeTrajDurations(TR);
fidelity = smi_core.TrackingResults.computeTrajFidelity(TR);
```

**Filter trajectories:**
```matlab
TR_long = smi_core.TrackingResults.threshTrajLength(TR, MinLength);
TR_windowed = smi_core.TrackingResults.windowTR(TR, MinFrame, MaxFrame, Verbose);
TR_early = smi_core.TrackingResults.windowStartTR(TR, MaxStartFrame, Verbose);
```

**Manipulate trajectories:**
```matlab
TR_combined = smi_core.TrackingResults.catTR(TR1, TR2, CheckDims);
TR_joined = smi_core.TrackingResults.joinTraj(TR, TrajectoryIDs, Verbose);
TR_padded = smi_core.TrackingResults.padTR(TR, TRPadding);
```

**Convert between formats:**
```matlab
% TR → SMD (flatten to list)
SMD = smi_core.TrackingResults.convertTRToSMD(TR);

% SMD → TR (requires ConnectID field)
TR = smi_core.TrackingResults.convertSMDToTR(SMD);
```

---

## Common Operations

### Creating Structures

```matlab
% Create empty SMF with defaults
SMF = smi_core.SingleMoleculeFitting();

% Create empty SMD
SMD = smi_core.SingleMoleculeData.createSMD();

% Create empty TR
TR = smi_core.TrackingResults.createTR();
```

### Loading from Files

```matlab
% Load SMLM results
load('FileDir/Results/Results.mat', 'SMD', 'SMF');

% Load SPT results
load('FileDir/Results/Results_SPT.mat', 'TR', 'SMD', 'SMF');

% Load just parameters
load('my_parameters.mat', 'SMF');
```

### Saving Structures

```matlab
% Save parameters
save('SMF_parameters.mat', 'SMF');

% Save results (always include SMF for reproducibility)
save('Results.mat', 'SMD', 'SMF');

% Save tracking results
save('Results_SPT.mat', 'TR', 'SMD', 'SMF');
```

### Filtering SMD

```matlab
% By precision
good_precision = SMD.X_SE < 0.15 & SMD.Y_SE < 0.15;

% By photons
bright = SMD.Photons > 500;

% By frame range
early = SMD.FrameNum <= 100;

% By region of interest
ROI = [20, 20, 100, 100];  % [YStart, XStart, YEnd, XEnd]
in_ROI = SMD.X >= ROI(2) & SMD.X <= ROI(4) & ...
         SMD.Y >= ROI(1) & SMD.Y <= ROI(3);

% Combine filters
selected = good_precision & bright & early & in_ROI;

% Apply filter
SMD_filtered.X = SMD.X(selected);
SMD_filtered.Y = SMD.Y(selected);
% ... copy other fields ...
```

### Combining SMD Structures

```matlab
% Using built-in method
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2, SMD3);

% Manual concatenation
SMD_all.X = [SMD1.X; SMD2.X; SMD3.X];
SMD_all.Y = [SMD1.Y; SMD2.Y; SMD3.Y];
% ... etc for all fields ...
```

### Accessing Trajectory Data

```matlab
% Iterate over all trajectories
for n = 1:length(TR)
    x = TR(n).X;
    y = TR(n).Y;
    frames = TR(n).FrameNum;
    % ... analyze trajectory ...
end

% Access specific trajectory
traj_5 = TR(5);
plot(traj_5.X, traj_5.Y, 'o-');
```

### Converting Units

```matlab
% Pixels to nanometers
X_nm = SMD.X * SMD.PixelSize * 1000;
Y_nm = SMD.Y * SMD.PixelSize * 1000;

% Frames to seconds
time_s = (SMD.FrameNum - 1) / SMD.FrameRate;

% Diffusion constant: physical to pixels²/frame
D_pixels = D_um2_per_s / (SMD.PixelSize^2 * SMD.FrameRate);
```

---

## Quick Lookup: Common Field Names

### Position and Uncertainty
- `X`, `Y`, `Z` - Position estimates
- `X_SE`, `Y_SE`, `Z_SE` - Position uncertainties

### Photometry
- `Photons` - Total collected photons
- `Bg` - Background per pixel
- `PSFSigma` - PSF width (symmetric)
- `PSFSigmaX`, `PSFSigmaY` - PSF width (asymmetric)

### Temporal
- `FrameNum` - Frame index (1 to NFrames)
- `DatasetNum` - Dataset index (1 to NDatasets)
- `FrameRate` - Acquisition rate (Hz)

### Quality
- `PValue` - Fit goodness (0 to 1)
- `LogLikelihood` - Fit log-likelihood
- `ThreshFlag` - Threshold pass/fail (0=pass)

### Metadata
- `NFrames` - Total frames
- `NDatasets` - Number of datasets
- `PixelSize` - Pixel size (micrometers)
- `XSize`, `YSize` - Image dimensions (pixels)

### Connection and Tracking
- `ConnectID` - Links same emitter across frames
- `DriftX`, `DriftY`, `DriftZ` - Drift per frame

---

## Field Units Summary

| Quantity | Unit | Notes |
|----------|------|-------|
| Position (X, Y) | pixels | Origin at (1,1) = top-left pixel center |
| Position (Z) | micrometers | Relative to focal plane |
| Photons | photons | After camera gain correction |
| Background | photons/pixel | Per-pixel in fitting region |
| PSF sigma | pixels | Gaussian PSF width parameter |
| Time | frames | Index 1 to NFrames |
| Duration | seconds | Computed from frame rate |
| Precision (SE) | pixels or micrometers | Standard error (CRLB) |
| Pixel size | micrometers | Back-projected onto sample |
| Frame rate | Hz | Frames per second |
| Diffusion | pixels²/frame or μm²/s | Context-dependent |
| Drift | pixels or micrometers | Cumulative from frame 1 |

---

## Default Values Summary

### SMF Default Parameters

| Sub-Structure | Key Field | Default |
|---------------|-----------|---------|
| Data | CameraType | `'EMCCD'` |
| Data | CameraGain | `1` |
| Data | CameraOffset | `0` |
| BoxFinding | BoxSize | `7` |
| BoxFinding | MinPhotons | `200` |
| Fitting | PSFSigma | `1` |
| Fitting | FitType | `'XYNB'` |
| Fitting | Iterations | `20` |
| Thresholding | On | `true` |
| Thresholding | MaxXY_SE | `0.2` |
| Thresholding | MinPValue | `0.01` |
| FrameConnection | On | `true` |
| FrameConnection | MaxSeparation | `1.0` |
| FrameConnection | MaxFrameGap | `5` |
| DriftCorrection | On | `true` |
| DriftCorrection | Method | `'DC-KNN'` |
| Tracking | D | `0.01` |
| Tracking | MaxDistFF | `5` |
| Tracking | MinTrackLength | `3` |

---

## See Also

### Detailed Concept Documents
- [SMF Structure](../core-concepts/smf-structure.md) - Complete SMF usage guide
- [SMD Structure](../core-concepts/smd-structure.md) - Complete SMD usage guide
- [TR Structure](../core-concepts/tr-structure.md) - Complete TR usage guide
- [Architecture](../core-concepts/architecture.md) - How structures fit together

### Workflow Guides
- [SMLM Analysis](../workflows/smlm-analysis.md) - Using SMF and SMD
- [SPT Tracking](../workflows/spt-tracking.md) - Using TR
- [Frame Connection](../workflows/frame-connection.md) - Creating ConnectID

### How-To Guides
- [Load Data](../how-to/load-data.md) - Configuring SMF.Data
- [Localize Molecules](../how-to/localize-molecules.md) - Creating SMD
- [Configure Parameters](../how-to/configure-parameters.md) - Setting SMF fields

### Repository Documentation
- doc/DataStructures/SMF.md - Complete SMF field list
- doc/DataStructures/SMD.md - Complete SMD field list
- doc/DataStructures/TR.md - TR structure definition

---

## Word Count

This reference document contains approximately 2,850 words of technical content, providing comprehensive lookup tables and quick reference information for all smite data structures.
