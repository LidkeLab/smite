---
title: "Multi-Channel Registration for Two-Color SMLM"
category: "workflows"
level: "advanced"
tags: ["channel-registration", "multi-color", "fiducial", "alignment", "transform"]
prerequisites: ["../core-concepts/smf-structure.md", "../core-concepts/smd-structure.md", "smlm-analysis.md"]
related: ["drift-correction.md", "../how-to/localize-molecules.md"]
summary: "Complete guide to aligning multi-channel SMLM data using fiducial-based registration transforms for accurate two-color super-resolution imaging"
estimated_time: "30 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Multi-Channel Registration for Two-Color SMLM

## Purpose

In multi-color SMLM experiments, chromatic aberration and optical misalignment cause systematic offsets between channels, preventing accurate colocalization analysis. This document explains how to create and apply registration transforms to align multi-channel data using the `smi_core.ChannelRegistration` class. You'll learn to calibrate transforms from fiducial bead images, apply corrections to localization data, and validate alignment accuracy for quantitative two-color super-resolution microscopy.

## Prerequisites

- Understanding of [SMF structure](../core-concepts/smf-structure.md)
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Completion of [SMLM analysis workflow](smlm-analysis.md)
- Basic understanding of geometric transforms (affine, polynomial, local weighted mean)
- Multi-channel SMLM data and fiducial calibration images

## Overview

Channel registration in smite addresses spatial misalignment between imaging channels caused by:

- **Chromatic aberration**: Different wavelengths focus at different positions
- **Optical path differences**: Separate emission paths introduce distortions
- **Dichroic alignment**: Imperfect beam splitting causes shifts and rotations
- **Field curvature**: Non-uniform distortion across the field of view

The registration workflow consists of:

1. **Acquiring fiducial images**: Imaging multi-color beads simultaneously in both channels
2. **Finding localizations**: Detecting bead positions in each channel
3. **Pairing coordinates**: Matching corresponding beads between channels
4. **Computing transforms**: Fitting geometric models to align channels
5. **Validating accuracy**: Measuring registration errors
6. **Applying transforms**: Correcting experimental data
7. **Exporting results**: Saving transforms for future use

## Why Channel Registration is Essential

### The Problem

Typical registration errors without correction: 50-200 nm (0.5-2 pixels at 100 nm pixel size)

This magnitude exceeds typical SMLM localization precision (10-30 nm), making uncorrected multi-color data unsuitable for:
- Colocalization analysis
- Distance measurements between labeled structures
- Multi-color super-resolution imaging
- Quantitative FRET analysis

### The Solution

Proper registration reduces alignment errors to 5-15 nm, enabling:
- Accurate colocalization analysis
- Reliable distance measurements
- High-quality two-color overlays
- Quantitative multi-color analysis

### Visual Example

```
Before Registration:          After Registration:
  Channel 1 (green)             Both channels aligned
  Channel 2 (red)
     █  █                          █
    █  █                           █
   █  █                            █
  █  █                             █
  Misaligned                    Aligned
```

## Registration Methods and Transform Types

### Transformation Basis

The `ChannelRegistration` class supports two approaches for computing transforms:

#### 1. Coordinate-Based Registration (Default)

**Method**: Localize fiducial beads in each channel, then fit transform to paired coordinates.

**When to use**:
- High SNR fiducial images
- Discrete, well-separated beads
- Best accuracy for sparse features
- Standard multi-color SMLM applications

**Process**:
1. Run `LocalizeData` on fiducial images
2. Pair localizations between channels
3. Fit geometric transform to coordinate pairs

**Advantages**:
- Sub-pixel accuracy
- Robust to noise
- Handles sparse fiducials

#### 2. Image-Based Registration

**Method**: Register fiducial images directly using intensity-based methods.

**When to use**:
- Low SNR or dense features
- Continuous structures (not point-like)
- Quick registration without localization step

**Process**:
1. Use MATLAB's `imregtform` directly on images
2. Optimize intensity-based metric

**Limitations**:
- Less accurate than coordinate-based
- Requires more computational resources
- Limited to simpler transform types

### Transform Types

Different geometric models accommodate various types of distortion:

#### Global Transforms

**Non-reflective similarity** (`'nonreflectivesimilarity'`):
- Translation + rotation + isotropic scaling
- 4 degrees of freedom
- Simple chromatic shifts
- Example: uniform offset + small rotation

**Similarity** (`'similarity'`):
- Translation + rotation + isotropic scaling + reflection
- 4-5 degrees of freedom
- Handles mirror configurations

**Affine** (`'affine'`):
- Translation + rotation + scaling + shear
- 6 degrees of freedom
- Moderate field distortions
- Most common for standard microscopy

**Projective** (`'projective'`):
- Perspective distortion
- 8 degrees of freedom
- Rarely needed for SMLM

#### Polynomial Transforms

**Polynomial** (`'polynomial'`):
- 2nd, 3rd, or 4th order polynomial
- Degree 2: 12 DOF, Degree 3: 20 DOF, Degree 4: 30 DOF
- Moderate non-linear field curvature
- Good for optical aberrations

**Recommended**: Polynomial degree 2 for most applications

#### Local Transforms

**Local Weighted Mean (LWM)** (`'lwm'`, **Default**):
- Locally varying transform
- Uses K nearest neighbor points
- Handles arbitrary non-linear distortions
- Best accuracy for complex aberrations

**Parameters**:
- `NNeighborPoints`: Number of neighbors (default 12, minimum 6)
- Higher values = smoother transform
- Lower values = more flexible, higher risk of overfitting

**Recommended**: LWM with 12 neighbor points for typical SMLM systems

**Piecewise Linear (PWL)** (`'pwl'`):
- Triangle-based tessellation
- Another local approach
- Similar to LWM but different algorithm

### Choosing Transform Type

**Decision guide**:

```
Is distortion uniform across FOV?
├─ YES → Try 'affine' first
│         If insufficient, try 'polynomial' degree 2
│
└─ NO → Use 'lwm' (default)
         Adjust NNeighborPoints based on spatial frequency of distortion
```

**Field size considerations**:
- Small FOV (< 50 μm): 'affine' often sufficient
- Medium FOV (50-100 μm): 'polynomial' or 'lwm'
- Large FOV (> 100 μm): 'lwm' with careful validation

## Complete Registration Workflow

### Step 1: Acquire Fiducial Images

**Requirements**:
- Multi-color fluorescent beads (e.g., TetraSpeck, 100-200 nm diameter)
- Same imaging conditions as experiments (magnification, filters, camera settings)
- Sufficient bead density: 10-100 beads per FOV
- Good SNR: > 10:1 signal to background

**Imaging protocol**:
1. Prepare fiducial sample (beads on coverslip)
2. Image in both channels with identical exposure
3. Use same optical path as experimental data
4. Save as separate files per channel OR single file with split-view

**File formats**:
- **Option A**: Two files, one per channel (`Fiducials_Ch1.h5`, `Fiducials_Ch2.h5`)
- **Option B**: One file with side-by-side channels (`Fiducials_SplitView.h5`)

### Step 2: Create ChannelRegistration Object

```matlab
% Create SMF for fiducial localization
SMF = smi_core.SingleMoleculeFitting();

% Configure for fiducial detection
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 500;  % Higher for bright beads
SMF.Fitting.FitType = 'XYNBS';    % Fit position, photons, bg, sigma
SMF.Fitting.PSFSigma = 1.3;       % Typical PSF width

% Set fiducial data path
SMF.Data.FileDir = '/path/to/fiducials';
SMF.Data.FileName = {'Fiducials_Ch1.h5', 'Fiducials_Ch2.h5'};

% Create ChannelRegistration object
CR = smi_core.ChannelRegistration([], [], SMF);
CR.Verbose = 1;  % Enable progress updates
```

### Step 3: Configure Registration Parameters

```matlab
% Choose transformation basis (default: 'coordinates')
CR.TransformationBasis = 'coordinates';  % or 'images'

% Choose transform type (default: 'lwm')
CR.TransformationType = 'lwm';  % Options: 'affine', 'polynomial', 'lwm', etc.

% Set LWM parameters (if using 'lwm')
CR.NNeighborPoints = 12;  % Number of neighbors (6-20 typical)

% Set polynomial degree (if using 'polynomial')
CR.PolynomialDegree = 2;  % Options: 2, 3, or 4

% Set pairing threshold
CR.SeparationThreshold = 2.0;  % Max distance for pairing (pixels)

% Enable manual cull to remove bad pairs
CR.ManualCull = true;  % Interactive quality control

% Auto-scale fiducial intensities
CR.AutoscaleFiducials = true;  % Simplifies processing
```

**Key parameters explained**:

- **`SeparationThreshold`**: Maximum distance (pixels) to pair beads between channels. Set based on expected misalignment (typically 1-5 pixels).
- **`ManualCull`**: Enables interactive GUI to remove incorrect pairings before computing transform.
- **`AutoscaleFiducials`**: Automatically normalizes image intensities, avoiding need for precise gain/offset correction.

### Step 4: Configure ROI Splitting (If Needed)

If channels are split side-by-side in one image, define the split format:

#### Option A: Horizontal Split (Default)

```matlab
% Channels arranged horizontally: [Ch1 | Ch2]
CR.SplitFormat = [1, 2];  % Default

% This splits the image in half vertically
% Left half = Channel 1 (reference)
% Right half = Channel 2 (moving)
```

#### Option B: Vertical Split

```matlab
% Channels arranged vertically:
% [Ch1]
% [Ch2]
CR.SplitFormat = [1; 2];

% This splits the image in half horizontally
% Top half = Channel 1 (reference)
% Bottom half = Channel 2 (moving)
```

#### Option C: Quadrant Split

```matlab
% Four quadrants:
% [Ch1 | Ch3]
% [Ch2 | Ch4]
CR.SplitFormat = [1, 3; 2, 4];

% Useful for 4-channel systems or beam splitter configurations
```

#### Option D: Manual ROI Definition

```matlab
% Manually specify exact ROI coordinates
CR.SplitFormat = [];  % Disable automatic splitting

% Define ROIs: [YStart, XStart, YEnd, XEnd, ZStart, ZPeriod]
CR.FiducialROI = [1, 1, 256, 256, 1, 1;      % Channel 1 (reference)
                  1, 257, 256, 512, 1, 1];   % Channel 2 (moving)

% First row is always the reference channel
% All other rows are transformed to match the reference
```

### Step 5: Compute Registration Transform

```matlab
% Find the transform
RegistrationTransform = CR.findTransform();

% The transform is stored in CR.RegistrationTransform
% Coordinates used for fitting are in CR.Coordinates
```

**What happens internally**:

1. **Load fiducials**: Images loaded from files or split according to ROI settings
2. **Rescale images**: Apply gain/offset correction or auto-scaling
3. **Localize beads** (if `TransformationBasis = 'coordinates'`):
   - Run `LocalizeData` on each channel
   - Find bead positions with sub-pixel accuracy
4. **Pair coordinates**:
   - Match beads between channels using nearest neighbor
   - Apply `SeparationThreshold` constraint
   - Remove unpaired or ambiguous matches
5. **Manual cull** (if `ManualCull = true`):
   - Display interactive GUI
   - User clicks to remove bad pairs
6. **Fit transform**:
   - Use MATLAB's `fitgeotrans` (coordinate-based) or `imregtform` (image-based)
   - Compute optimal transformation parameters
7. **Store result**:
   - Save transform object in `CR.RegistrationTransform{2}`
   - Save paired coordinates in `CR.Coordinates{2}`

**Manual culling interface**:

When `ManualCull = true`, a GUI appears showing:
- Overlay of both channels with paired beads marked
- Click on bad pairs to remove them
- Press 'Enter' when satisfied
- Good practice: Remove outliers, edge beads, and mismatched pairs

### Step 6: Validate Registration Accuracy

#### Compute Registration Errors

```matlab
% Estimate registration error on calibration data
MovingCoords = CR.Coordinates{2}(:, :, 2);  % Channel 2 (moving)
FixedCoords = CR.Coordinates{2}(:, :, 1);   % Channel 1 (reference)

% Standard registration error
SquaredError = CR.estimateRegistrationError(...
    CR.RegistrationTransform{2}, MovingCoords, FixedCoords);
RMSE = sqrt(mean(SquaredError));

fprintf('Registration RMSE: %.3f pixels (%.1f nm)\n', ...
    RMSE, RMSE * SMF.Data.PixelSize * 1000);

% Leave-one-out error (more robust estimate)
RegParams = {CR.NNeighborPoints};  % For LWM
SquaredErrorLOO = CR.estimateRegErrorLOO(...
    CR.TransformationType, RegParams, MovingCoords, FixedCoords);
RMSE_LOO = sqrt(mean(SquaredErrorLOO));

fprintf('Leave-one-out RMSE: %.3f pixels (%.1f nm)\n', ...
    RMSE_LOO, RMSE_LOO * SMF.Data.PixelSize * 1000);
```

**Target accuracy**:
- Excellent: RMSE < 0.1 pixels (< 10 nm at 100 nm pixel size)
- Good: RMSE < 0.2 pixels (< 20 nm)
- Acceptable: RMSE < 0.5 pixels (< 50 nm)
- Poor: RMSE > 0.5 pixels (re-evaluate transform type or fiducials)

**Leave-one-out error**: More conservative estimate that predicts performance on new data by iteratively excluding each point during fitting.

#### Visualize Registration Results

```matlab
% Show before/after comparison
PlotFigure = figure('Position', [100, 100, 1200, 500]);

CR.visualizeRegistrationResults(PlotFigure, ...
    CR.RegistrationTransform{2}, ...
    MovingCoords, FixedCoords, ...
    CR.FiducialImages(:, :, 2), CR.FiducialImages(:, :, 1));

% This creates side-by-side plots:
% Left: Before registration (green + magenta overlay)
% Right: After registration (aligned overlay)
% Titles show RMSE before and after
```

**What to look for**:
- Magenta and green markers should overlap after registration
- RMSE should decrease significantly (typically 5-10x reduction)
- Residual errors should be spatially uniform (no systematic patterns)

#### Visualize Error Distribution

```matlab
% Show spatial distribution of registration errors
figure('Position', [100, 100, 800, 600]);

CR.visualizeRegistrationError(gca, ...
    CR.RegistrationTransform{2}, ...
    MovingCoords, FixedCoords, ...
    [256, 256], 10);  % FOV size, grid spacing

% This shows:
% - Color map of registration error magnitude across FOV
% - Quiver plot showing error vectors
% - Identifies regions with higher residual errors
```

**Diagnosis**:
- **Uniform low error**: Excellent registration
- **Edge artifacts**: Insufficient fiducials near edges, consider tighter FOV
- **Systematic patterns**: Wrong transform type, try higher order polynomial or LWM
- **Random spikes**: Bad fiducial localizations, re-run with stricter thresholds

#### Visualize Transform Distortion

```matlab
% Show how transform warps a regular grid
figure('Position', [100, 100, 800, 600]);

CR.visualizeCoordTransform(gcf, ...
    CR.RegistrationTransform{2}, ...
    [256, 256], 10);  % FOV size, grid spacing

% This shows:
% - Regular grid before and after transformation
% - Magnitude and direction of local distortions
% - Useful for understanding optical aberrations
```

### Step 7: Export Transform

```matlab
% Save transform to file for later use
FilePath = CR.exportTransform(SMF.Data.FileDir);

fprintf('Transform saved to: %s\n', FilePath);
```

**Saved file contents** (`RegistrationTransform_*.mat`):
- `RegistrationTransform`: Cell array of transform objects
- `Coordinates`: Paired coordinates used for fitting
- `FiducialROI`: ROI definitions
- `TransformationType`, `TransformationBasis`: Method parameters
- `NNeighborPoints`, `PolynomialDegree`: Transform-specific parameters
- `SeparationThreshold`: Pairing threshold
- `RegistrationError`, `RegistrationErrorLOO`: Accuracy metrics
- `SMF`: Complete SMF structure used for registration

**Naming convention**: `RegistrationTransform_<FiducialFileName>.mat`

**Reuse**: Load this file to apply the same transform to other datasets acquired with identical optical setup.

### Step 8: Apply Transform to Experimental Data

#### Apply to SMD Structures

```matlab
% Load experimental data (already localized)
load('ExperimentData_Ch2.mat', 'SMD_Ch2');  % Moving channel

% Load the registration transform
load('RegistrationTransform_Fiducials_Ch1.mat', 'RegistrationTransform');

% Apply transform to align Ch2 to Ch1
SMD_Ch2_Aligned = smi_core.ChannelRegistration.transformSMD(...
    RegistrationTransform{2}, SMD_Ch2);

% Now SMD_Ch2_Aligned coordinates match Ch1 reference frame
% Save aligned data
save('ExperimentData_Ch2_Aligned.mat', 'SMD_Ch2_Aligned');
```

**What changes**:
- `SMD.X`, `SMD.Y`: Transformed to reference channel coordinates
- `SMD.XRegCorrection`, `SMD.YRegCorrection`: Correction vectors applied
- `SMD.IsTransformed`: Set to `true` (prevents re-transformation)

**Important**: Transformation is applied to the `SMD` structure, NOT the raw images. Localization must be done before registration.

#### Apply to TR Structures (Tracking Data)

```matlab
% For single particle tracking data
load('TrackingData_Ch2.mat', 'TR_Ch2');

% Apply transform
TR_Ch2_Aligned = smi_core.ChannelRegistration.transformTR(...
    RegistrationTransform{2}, TR_Ch2);

% Trajectories now aligned to Ch1 reference
save('TrackingData_Ch2_Aligned.mat', 'TR_Ch2_Aligned');
```

#### Apply to Raw Coordinates

```matlab
% Transform arbitrary coordinate arrays
MovingCoords = [X_Ch2, Y_Ch2];  % N × 2 array

% Apply transform
AlignedCoords = smi_core.ChannelRegistration.transformCoords(...
    RegistrationTransform{2}, MovingCoords);

X_Ch2_Aligned = AlignedCoords(:, 1);
Y_Ch2_Aligned = AlignedCoords(:, 2);
```

#### Apply to Images (Advanced)

```matlab
% Transform entire images (less common, mostly for visualization)
MovingImage = imread('Image_Ch2.tif');

% Apply transform
AlignedImage = smi_core.ChannelRegistration.transformImages(...
    RegistrationTransform{2}, MovingImage);

% Save aligned image
imwrite(uint16(AlignedImage), 'Image_Ch2_Aligned.tif');
```

**Note**: Image transformation is mainly for visualization. Localizations should be aligned, not raw images.

## Complete Example Workflows

### Example 1: Two-File Configuration

```matlab
% Two separate fiducial files, one per channel

% Step 1: Configure SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/2024-03-15/calibration';
SMF.Data.FileName = {'Fiducials_Ch1.h5', 'Fiducials_Ch2.h5'};
SMF.Data.PixelSize = 0.100;  % 100 nm pixels
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 800;
SMF.Fitting.FitType = 'XYNBS';
SMF.Fitting.PSFSigma = 1.3;

% Step 2: Create registration object
CR = smi_core.ChannelRegistration([], [], SMF);
CR.TransformationType = 'lwm';
CR.NNeighborPoints = 12;
CR.SeparationThreshold = 2.5;
CR.ManualCull = true;
CR.Verbose = 1;

% Step 3: Compute transform
RegistrationTransform = CR.findTransform();

% Step 4: Validate
MovingCoords = CR.Coordinates{2}(:, :, 2);
FixedCoords = CR.Coordinates{2}(:, :, 1);
SquaredError = CR.estimateRegistrationError(...
    CR.RegistrationTransform{2}, MovingCoords, FixedCoords);
fprintf('RMSE: %.3f pixels\n', sqrt(mean(SquaredError)));

% Step 5: Visualize
figure;
CR.visualizeRegistrationResults(gcf, ...
    CR.RegistrationTransform{2}, MovingCoords, FixedCoords, ...
    CR.FiducialImages(:, :, 2), CR.FiducialImages(:, :, 1));

% Step 6: Export
FilePath = CR.exportTransform();

% Step 7: Apply to experimental data
load('Experiment_Ch2_Results.mat', 'SMD');
SMD_Aligned = smi_core.ChannelRegistration.transformSMD(...
    RegistrationTransform{2}, SMD);
save('Experiment_Ch2_Aligned.mat', 'SMD_Aligned');
```

### Example 2: Split-View Configuration

```matlab
% Single file with side-by-side channels

% Step 1: Configure SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/2024-03-15/calibration';
SMF.Data.FileName = {'Fiducials_SplitView.h5'};
SMF.Data.PixelSize = 0.108;  % 108 nm pixels
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 600;
SMF.Fitting.FitType = 'XYNB';  % Fixed sigma
SMF.Fitting.PSFSigma = 1.4;

% Step 2: Create registration object
CR = smi_core.ChannelRegistration([], [], SMF);
CR.TransformationType = 'polynomial';
CR.PolynomialDegree = 2;
CR.SeparationThreshold = 3.0;
CR.ManualCull = true;

% Step 3: Configure horizontal split
CR.SplitFormat = [1, 2];  % Left = Ch1 (reference), Right = Ch2 (moving)

% Step 4: Compute transform
RegistrationTransform = CR.findTransform();

% Step 5: Validate
MovingCoords = CR.Coordinates{2}(:, :, 2);
FixedCoords = CR.Coordinates{2}(:, :, 1);
SquaredError = CR.estimateRegistrationError(...
    CR.RegistrationTransform{2}, MovingCoords, FixedCoords);
RMSE = sqrt(mean(SquaredError));
fprintf('Registration accuracy: %.1f nm\n', RMSE * SMF.Data.PixelSize * 1000);

% Step 6: Export
FilePath = CR.exportTransform();

% Step 7: Load and apply to multiple datasets
for ii = 1:5
    filename = sprintf('Cell%d_Ch2.mat', ii);
    load(filename, 'SMD');
    SMD_Aligned = smi_core.ChannelRegistration.transformSMD(...
        RegistrationTransform{2}, SMD);
    save(sprintf('Cell%d_Ch2_Aligned.mat', ii), 'SMD_Aligned');
end
```

### Example 3: Image-Based Registration

```matlab
% Use image-based registration (no localization step)

% Step 1: Configure SMF (minimal settings)
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/calibration';
SMF.Data.FileName = {'Fiducials_Ch1.h5', 'Fiducials_Ch2.h5'};
SMF.Data.PixelSize = 0.100;

% Step 2: Create registration object
CR = smi_core.ChannelRegistration([], [], SMF);
CR.TransformationBasis = 'images';  % Use images directly
CR.TransformationType = 'affine';   % Simpler transform for image method
CR.ManualCull = false;  % Not applicable for image-based
CR.AutoscaleFiducials = true;
CR.Verbose = 1;

% Step 3: Compute transform
RegistrationTransform = CR.findTransform();

% Step 4: Export
FilePath = CR.exportTransform();

% Step 5: Apply to data
load('Experiment_Ch2.mat', 'SMD');
SMD_Aligned = smi_core.ChannelRegistration.transformSMD(...
    RegistrationTransform{2}, SMD);
save('Experiment_Ch2_Aligned.mat', 'SMD_Aligned');
```

### Example 4: Batch Processing Multiple Sessions

```matlab
% Register multiple experimental sessions with same optical setup

% Step 1: Create transform once from fiducials
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/calibration';
SMF.Data.FileName = {'Fiducials_Ch1.h5', 'Fiducials_Ch2.h5'};
SMF.Data.PixelSize = 0.100;
SMF.BoxFinding.MinPhotons = 700;
SMF.Fitting.FitType = 'XYNBS';
SMF.Fitting.PSFSigma = 1.3;

CR = smi_core.ChannelRegistration([], [], SMF);
CR.TransformationType = 'lwm';
CR.NNeighborPoints = 12;
CR.findTransform();

% Export transform
TransformPath = CR.exportTransform();

% Step 2: Apply to all experimental sessions
SessionDirs = {'/data/2024-03-15', '/data/2024-03-16', '/data/2024-03-17'};

for session_idx = 1:length(SessionDirs)
    SessionDir = SessionDirs{session_idx};
    fprintf('Processing session: %s\n', SessionDir);

    % Load transform (same for all sessions)
    load(TransformPath, 'RegistrationTransform');

    % Find all Ch2 results in this session
    ResultFiles = dir(fullfile(SessionDir, 'Results', '*_Ch2.mat'));

    for file_idx = 1:length(ResultFiles)
        % Load Ch2 data
        FilePath = fullfile(ResultFiles(file_idx).folder, ...
            ResultFiles(file_idx).name);
        load(FilePath, 'SMD');

        % Apply transform
        SMD_Aligned = smi_core.ChannelRegistration.transformSMD(...
            RegistrationTransform{2}, SMD);

        % Save with '_Aligned' suffix
        [~, name, ext] = fileparts(FilePath);
        SavePath = fullfile(ResultFiles(file_idx).folder, ...
            [name, '_Aligned', ext]);
        save(SavePath, 'SMD_Aligned');

        fprintf('  Aligned: %s\n', ResultFiles(file_idx).name);
    end
end

fprintf('Batch registration complete!\n');
```

## Troubleshooting

### Issue: Poor registration accuracy (high RMSE)

**Diagnose**:
```matlab
% Check number of paired beads
NPairs = size(CR.Coordinates{2}, 1);
fprintf('Number of paired beads: %d\n', NPairs);

% Need at least:
% - Affine: 3 pairs minimum, 10+ recommended
% - Polynomial degree 2: 6 pairs minimum, 20+ recommended
% - LWM: 10+ pairs minimum, 50+ recommended

% Visualize paired coordinates
figure;
plot(CR.Coordinates{2}(:, 1, 1), CR.Coordinates{2}(:, 2, 1), 'go', ...
    'MarkerSize', 10, 'DisplayName', 'Channel 1');
hold on;
plot(CR.Coordinates{2}(:, 1, 2), CR.Coordinates{2}(:, 2, 2), 'r+', ...
    'MarkerSize', 10, 'DisplayName', 'Channel 2');
legend; title('Paired bead coordinates');
```

**Solutions**:
1. **Insufficient fiducials**: Image more beads or use larger FOV
2. **Wrong transform type**: Try higher-order polynomial or LWM
3. **Bad localizations**: Increase `BoxFinding.MinPhotons`, use `ManualCull = true`
4. **Pairing errors**: Decrease `SeparationThreshold`, enable `ManualCull`
5. **Poor fiducial quality**: Use brighter beads, improve focus, increase exposure

### Issue: Transform fails to compute

**Error**: "Not enough points to compute transform"

**Solution**:
```matlab
% Check if enough points were paired
NPairs = sum(~isnan(CR.Coordinates{2}(:, 1, 1)));
fprintf('Valid pairs: %d\n', NPairs);

% Minimum required points:
% - 'similarity': 2
% - 'affine': 3
% - 'projective': 4
% - 'polynomial' degree 2: 6
% - 'lwm': NNeighborPoints + 1

% If insufficient, relax pairing threshold
CR.SeparationThreshold = 5.0;  % Increase from default
CR.findTransform();
```

### Issue: Transform applied incorrectly

**Symptom**: Channels more misaligned after transformation

**Diagnose**:
```matlab
% Check if SMD was already transformed
if isfield(SMD, 'IsTransformed') && SMD.IsTransformed
    warning('SMD already transformed!');
end

% Check transform indices
load('RegistrationTransform_*.mat', 'RegistrationTransform');
fprintf('Number of transforms: %d\n', length(RegistrationTransform));

% RegistrationTransform{1} is always empty (reference channel)
% RegistrationTransform{2} transforms Ch2 → Ch1
% RegistrationTransform{3} transforms Ch3 → Ch1, etc.
```

**Solutions**:
1. **Already transformed**: Use original untransformed SMD
2. **Wrong transform index**: Use `RegistrationTransform{2}` for Channel 2
3. **Wrong direction**: Ensure moving channel gets transformed, not reference
4. **Reference mismatch**: Verify same reference channel for calibration and data

### Issue: Different results each time

**Symptom**: Registration RMSE varies significantly between runs

**Cause**: Localization randomness, manual culling differences

**Solution**:
```matlab
% Disable manual culling for reproducibility
CR.ManualCull = false;

% Use stricter automated pairing
CR.SeparationThreshold = 1.5;  % Tighter threshold

% Use higher photon threshold for consistent localizations
SMF.BoxFinding.MinPhotons = 1000;  % Higher threshold

% Re-run
RegistrationTransform = CR.findTransform();

% Save settings and transform together
save('RegisterConfig.mat', 'CR', 'SMF', 'RegistrationTransform');
```

### Issue: Edge artifacts in registration

**Symptom**: Higher errors near FOV edges

**Diagnose**:
```matlab
% Check spatial distribution of fiducials
figure;
scatter(CR.Coordinates{2}(:, 1, 1), CR.Coordinates{2}(:, 2, 1), ...
    50, 'filled');
axis equal; title('Fiducial distribution');
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Look for sparse coverage near edges
```

**Solutions**:
1. **Use smaller FOV**: Crop edges with sparse coverage
2. **Image more fiducials**: Increase bead density
3. **Use lower-order transform**: Polynomial extrapolates poorly, try 'affine'
4. **Validate on interior region**: Exclude edge localizations from analysis

### Issue: Systematic colocalization bias

**Symptom**: All colocalizations show consistent offset in one direction

**Diagnose**:
```matlab
% Compute mean offset between channels
load('Experiment_Ch1.mat', 'SMD_Ch1');
load('Experiment_Ch2_Aligned.mat', 'SMD_Ch2_Aligned');

% Match nearby localizations (simple NN)
D = pdist2([SMD_Ch1.X, SMD_Ch1.Y], ...
           [SMD_Ch2_Aligned.X, SMD_Ch2_Aligned.Y]);
[MinDist, NN_idx] = min(D, [], 2);

% Compute offset for nearby pairs (< 100 nm)
NearbyMask = MinDist < 0.1 / SMF.Data.PixelSize;  % 100 nm
OffsetX = mean(SMD_Ch2_Aligned.X(NN_idx(NearbyMask)) - ...
              SMD_Ch1.X(NearbyMask));
OffsetY = mean(SMD_Ch2_Aligned.Y(NN_idx(NearbyMask)) - ...
              SMD_Ch1.Y(NearbyMask));

fprintf('Residual offset: (%.1f, %.1f) nm\n', ...
    OffsetX * SMF.Data.PixelSize * 1000, ...
    OffsetY * SMF.Data.PixelSize * 1000);
```

**Solutions**:
1. **Re-calibrate with fresh fiducials**: Old calibration may not match current alignment
2. **Check for stage drift**: Fiducials and experiment may have different drift
3. **Verify same optical setup**: Filter changes, focus shifts invalidate calibration
4. **Add constant offset correction**: Apply residual offset manually if systematic

## Best Practices

### Fiducial Imaging

1. **Bead selection**: 100-200 nm TetraSpeck beads (Thermo Fisher)
2. **Density**: 10-100 beads per FOV (sparse enough to resolve individually)
3. **SNR**: > 10:1 signal-to-background ratio
4. **Coverage**: Uniform distribution across FOV, including edges
5. **Mounting**: Same coverslip and mounting medium as experiments
6. **Imaging**: Same microscope settings, exposure time, and focus as experiments

### Transform Selection

1. **Start simple**: Try 'affine' first, then 'polynomial', then 'lwm'
2. **Small FOV**: 'affine' often sufficient
3. **Large FOV or complex optics**: Use 'lwm' with 10-15 neighbor points
4. **Avoid overfitting**: More complex transforms require more fiducials
5. **Validate on independent data**: Test transform on separate fiducial images

### Workflow Organization

1. **One calibration per optical configuration**: Different magnifications, filters, or beam paths need separate calibrations
2. **Regular recalibration**: Weekly or after microscope adjustments
3. **Document calibration conditions**: Note date, settings, focus, fiducial lot
4. **Version control transforms**: Save with descriptive names including date
5. **Validate before experiments**: Quick test with fresh fiducials before each session

### Quality Control

1. **Check RMSE**: Should be < 0.2 pixels (< 20 nm at 100 nm pixel size)
2. **Inspect visualization**: Verify alignment in `visualizeRegistrationResults`
3. **Leave-one-out validation**: More robust than standard RMSE
4. **Spatial uniformity**: Errors should be uniform across FOV
5. **Test on control samples**: Known colocalized structures as positive control

## Registration File Format

The exported `.mat` file contains:

```matlab
% Load registration file
load('RegistrationTransform_Fiducials_Ch1.mat');

% Key variables:
% RegistrationTransform{2}  - Transform object (empty for {1})
% Coordinates{2}            - Paired coords: (:,:,1) = fixed, (:,:,2) = moving
% TransformationType        - 'lwm', 'affine', 'polynomial', etc.
% TransformationBasis       - 'coordinates' or 'images'
% NNeighborPoints           - For LWM transform
% PolynomialDegree          - For polynomial transform
% SeparationThreshold       - Pairing threshold (pixels)
% FiducialROI               - ROI definitions
% SplitFormat               - Split configuration
% SMF                       - Full SMF structure
% RegistrationError         - RMSE (pixels)
% RegistrationErrorLOO      - Leave-one-out RMSE (pixels)
```

**Usage**:
```matlab
% Apply saved transform
load('RegistrationTransform_Fiducials_Ch1.mat', 'RegistrationTransform');
SMD_Aligned = smi_core.ChannelRegistration.transformSMD(...
    RegistrationTransform{2}, SMD_Moving);
```

## Integration with SMLM Workflow

Channel registration typically occurs after localization and before colocalization analysis:

```
1. Acquire multi-channel data
2. Localize each channel independently (smi.SMLM)
3. Create/load registration transform (smi_core.ChannelRegistration)
4. Apply transform to non-reference channels
5. Perform colocalization analysis (smi_cluster.PairCorrelation, etc.)
6. Generate multi-color overlays (smi_vis.GenerateImages)
```

**Complete integrated workflow**:

```matlab
% Calibration phase (once per optical setup)
SMF_Fid = smi_core.SingleMoleculeFitting();
SMF_Fid.Data.FileDir = '/data/calibration';
SMF_Fid.Data.FileName = {'Fiducials_Ch1.h5', 'Fiducials_Ch2.h5'};
% ... configure SMF_Fid ...

CR = smi_core.ChannelRegistration([], [], SMF_Fid);
% ... configure CR ...
CR.findTransform();
TransformPath = CR.exportTransform();

% Experimental phase (for each dataset)
% 1. Localize Ch1 (reference)
SMF_Ch1 = smi_core.SingleMoleculeFitting();
SMF_Ch1.Data.FileName = {'Experiment_Ch1.h5'};
% ... configure SMF_Ch1 ...
SMLMobj_Ch1 = smi.SMLM(SMF_Ch1);
SMLMobj_Ch1.fullAnalysis();
SMD_Ch1 = SMLMobj_Ch1.SMD;

% 2. Localize Ch2 (moving)
SMF_Ch2 = smi_core.SingleMoleculeFitting();
SMF_Ch2.Data.FileName = {'Experiment_Ch2.h5'};
% ... configure SMF_Ch2 ...
SMLMobj_Ch2 = smi.SMLM(SMF_Ch2);
SMLMobj_Ch2.fullAnalysis();
SMD_Ch2 = SMLMobj_Ch2.SMD;

% 3. Register Ch2 to Ch1
load(TransformPath, 'RegistrationTransform');
SMD_Ch2_Aligned = smi_core.ChannelRegistration.transformSMD(...
    RegistrationTransform{2}, SMD_Ch2);

% 4. Perform colocalization analysis
% (example: pair correlation)
PC = smi_cluster.PairCorrelation();
PC.analyzeCorrelation(SMD_Ch1, SMD_Ch2_Aligned);

% 5. Generate two-color overlay
SR_zoom = 20;
Img_Ch1 = smi_vis.GenerateImages.gaussImage(...
    [SMD_Ch1.X, SMD_Ch1.Y], ...
    SMD_Ch1.PixelSize, SR_zoom, [SMD_Ch1.YSize, SMD_Ch1.XSize]);
Img_Ch2 = smi_vis.GenerateImages.gaussImage(...
    [SMD_Ch2_Aligned.X, SMD_Ch2_Aligned.Y], ...
    SMD_Ch2_Aligned.PixelSize, SR_zoom, [SMD_Ch2_Aligned.YSize, SMD_Ch2_Aligned.XSize]);

% Create RGB overlay
RGB = zeros([size(Img_Ch1), 3]);
RGB(:,:,2) = Img_Ch1 / max(Img_Ch1(:));  % Green = Ch1
RGB(:,:,1) = Img_Ch2 / max(Img_Ch2(:));  % Red = Ch2
figure; imshow(RGB); title('Two-color overlay (aligned)');
```

## See Also

- [SMLM Analysis Workflow](smlm-analysis.md) - Localization pipeline
- [Drift Correction](drift-correction.md) - Temporal alignment
- [SMF Structure](../core-concepts/smf-structure.md) - Parameter configuration
- [SMD Structure](../core-concepts/smd-structure.md) - Results format
- MATLAB/+smi_core/@ChannelRegistration/README.md - Class documentation
