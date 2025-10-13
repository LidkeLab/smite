---
title: "Drift Correction in SMLM"
category: "workflows"
level: "intermediate"
tags: ["drift-correction", "knn", "brightfield", "registration", "alignment"]
prerequisites: ["../core-concepts/smf-structure.md", "../core-concepts/smd-structure.md", "smlm-analysis.md"]
related: ["../how-to/localize-molecules.md", "../examples/advanced-processing.md"]
summary: "Comprehensive guide to correcting stage drift during long SMLM acquisitions using fiducial-free and brightfield-based methods"
estimated_time: "20 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Drift Correction in SMLM

## Purpose

Stage drift during long SMLM acquisitions can introduce significant errors in localization positions, degrading super-resolution image quality. This document explains why drift correction is essential, describes the available methods in smite, and shows how to configure and apply drift correction to achieve accurate, aligned localizations across entire datasets.

## Prerequisites

- Understanding of [SMF structure](../core-concepts/smf-structure.md)
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Completion of [SMLM analysis workflow](smlm-analysis.md)
- Basic understanding of image registration concepts

## Overview

Drift correction in smite addresses two related problems:

1. **Intra-dataset drift**: Gradual stage drift within a single acquisition
2. **Inter-dataset drift**: Alignment between multiple acquisitions of the same region

The toolbox provides two main approaches:

- **DC-KNN**: Fiducial-free drift correction using K-nearest neighbors
- **Brightfield registration**: Drift correction using brightfield/focus images

Both methods produce drift trajectories stored in `SMD.DriftX` and `SMD.DriftY`, and correct positions are stored in `SMD.X` and `SMD.Y`.

## Why Drift Correction is Needed

### The Problem

During SMLM acquisitions lasting minutes to hours, the microscope stage drifts due to:

- Thermal expansion/contraction
- Mechanical relaxation
- Environmental vibrations
- Temperature gradients

Typical drift rates: 5-50 nm/minute (0.05-0.5 pixels/minute at 100 nm pixel size)

Over a 30-minute acquisition, this produces 150-1500 nm of drift, comparable to or exceeding the super-resolution precision.

### The Effect

**Without drift correction:**
- Super-resolution images appear blurred
- Structures artificially elongated in drift direction
- Multiple acquisitions cannot be aligned
- Resolution degraded from ~20 nm to >100 nm

**With drift correction:**
- Sharp super-resolution images
- True structural features preserved
- Multiple datasets properly aligned
- Full resolution recovered

### Visual Example

```
Before Drift Correction:     After Drift Correction:
   Blurred, elongated           Sharp, accurate
        ████                         ██
       ████                          ██
      ████                           ██
     ████                            ██
```

## Drift Correction Methods

### Method 1: DC-KNN (Fiducial-Free)

The DC-KNN method estimates drift directly from localization data without requiring fiducial markers.

#### How It Works

**Basic principle**: Molecules that appear in multiple frames should maintain constant relative positions. Systematic shifts indicate stage drift.

**Algorithm (simplified)**:

1. For each localization, find K nearest neighbors in adjacent frames
2. Compute median displacement vectors
3. Build smooth drift trajectory using polynomial fitting
4. Apply corrections to all positions

**Two-stage process**:

**Intra-dataset correction**:
- Estimates smooth drift within each dataset
- Fits polynomial of degree PDegree to drift vs. time
- Handles gradual, continuous drift

**Inter-dataset correction**:
- Aligns multiple datasets to common reference
- Estimates constant offset between datasets
- Handles discrete jumps between acquisitions

#### Mathematical Details

The cost function minimizes the sum of nearest neighbor distances:

```
minimize: Σ min(||x_i(t) - drift(t) - NN(x_i(t))||, L)
```

where:
- `x_i(t)` = position of localization i at time t
- `drift(t)` = drift model (polynomial for intra, constant for inter)
- `NN(x_i(t))` = nearest neighbor position
- `L` = threshold to limit influence of distant neighbors

The threshold `L` prevents spurious long-distance matches from dominating the optimization, critical for sparse datasets.

#### Advantages

- No fiducial markers required
- Works with any sample
- Estimates smooth drift trajectories
- Handles both 2D and 3D data
- Robust to emitter blinking

#### Limitations

- Requires sufficient localization density (>50 per frame)
- Assumes stable emitter population
- Can fail with severe drift (>5 pixels between frames)
- Sensitive to non-uniform bleaching

#### When to Use

- Dense PAINT or STORM data
- No fiducials available
- Continuous acquisitions
- Smooth, gradual drift expected

### Method 2: Brightfield Registration

This method uses brightfield or phase contrast images taken before/after each dataset to estimate drift.

#### How It Works

1. Acquire brightfield image before dataset starts
2. Acquire brightfield image after dataset completes
3. Register these images to reference using cross-correlation
4. Interpolate drift linearly between pre/post images
5. Apply to localizations based on frame number

**Intra-dataset**: Linear interpolation between pre and post images
**Inter-dataset**: Registration of all pre-images to common reference

#### Advantages

- Independent of localization density
- Works with sparse samples
- Fast computation
- Reliable for large drifts
- No assumptions about emitter stability

#### Limitations

- Requires brightfield/focus images in data file
- Assumes linear drift within datasets
- Sample must have visible features in brightfield
- Cannot detect nonlinear drift

#### When to Use

- Sparse localization data
- Large drift expected
- Brightfield images available
- Quick processing needed

## Configuring Drift Correction in SMF

### Basic Configuration

```matlab
% Enable drift correction
SMF.DriftCorrection.On = true;

% Choose method
SMF.DriftCorrection.Method = 'DC-KNN';  % or not used for BF

% Brightfield registration flag
SMF.DriftCorrection.BFRegistration = true;  % Use BF images if available

% Intra-dataset smoothing threshold
SMF.DriftCorrection.L_intra = 1;  % pixels

% Inter-dataset smoothing threshold
SMF.DriftCorrection.L_inter = 2;  % pixels

% Polynomial degree for drift model
SMF.DriftCorrection.PDegree = 1;  % 1 = linear drift

% Pixel size for 3D (only needed for Z drift)
SMF.DriftCorrection.PixelSizeZUnit = 0.1;  % micrometers
```

### Parameter Guide

#### `On` (logical)
- `true`: Enable drift correction
- `false`: Disable drift correction
- Default: `true`

#### `Method` (string)
- `'DC-KNN'`: K-nearest neighbors method
- Default: `'DC-KNN'`

#### `BFRegistration` (logical)
- `true`: Use brightfield images if available
- `false`: Use DC-KNN on localizations
- Default: `true`
- Note: If true but no BF images found, falls back to DC-KNN

#### `L_intra` (scalar, pixels)
- Threshold for intra-dataset nearest neighbor distances
- Larger values = more smoothing, less detail
- Smaller values = follow rapid changes, more noise sensitive
- Typical range: 0.5 - 2.0 pixels
- Default: 1 pixel

**How to choose**:
- High density data (>200 locs/frame): 0.5 - 1.0
- Medium density (50-200 locs/frame): 1.0 - 1.5
- Low density (<50 locs/frame): 1.5 - 2.0

#### `L_inter` (scalar, pixels)
- Threshold for inter-dataset nearest neighbor distances
- Should be larger than L_intra (handles discrete jumps)
- Typical range: 1.0 - 3.0 pixels
- Default: 2 pixels

**How to choose**:
- Small jumps expected: 1.0 - 2.0
- Large jumps or shifts: 2.0 - 3.0

#### `PDegree` (integer)
- Polynomial degree for intra-dataset drift model
- 0 = constant (no intra-dataset correction)
- 1 = linear drift (most common)
- 2 = quadratic drift (rarely needed)
- Default: 1

**How to choose**:
- Smooth, steady drift: 1
- Accelerating drift: 2
- Very short acquisitions: 0

#### `PixelSizeZUnit` (scalar, micrometers)
- XY pixel size in same units as Z
- Only needed for 3D drift correction
- Used to convert Z coordinates to pixels
- Default: 0.1 micrometers

## Applying Drift Correction

### Automatic Correction (via SMLM Workflow)

Drift correction happens automatically when enabled:

```matlab
% Configure SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'cell1.h5'};

% Enable drift correction
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.L_intra = 1.0;
SMF.DriftCorrection.L_inter = 2.0;

% Run full analysis
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();

% Results already drift-corrected
SMD = SMLMobj.SMD;
```

### Manual Correction (Advanced)

For more control, apply drift correction manually:

```matlab
% After localization and frame connection
SMD_raw = ...; % From LocalizeData and FrameConnection

% Create drift correction object
DC = smi_core.DriftCorrection(SMF, SMD_raw);

% Apply DC-KNN method
[SMD_corrected, Statistics] = DC.driftCorrectKNN(SMD_raw);

% Examine statistics
fprintf('Intra-dataset correction time: %.2f sec\n', ...
    Statistics.Intra_elapsedTime);
fprintf('Inter-dataset correction time: %.2f sec\n', ...
    Statistics.Inter_elapsedTime);
```

### Separated Intra/Inter Correction

For batch processing of multiple datasets:

```matlab
% Initialize
DC = smi_core.DriftCorrection(SMF, SMD);
SMDIntra = [];

% Process each dataset separately (intra-dataset only)
for i = 1:NDatasets
    SMD_i = ...; % Extract dataset i
    [SMD_i_corrected, Stats_i] = DC.driftCorrectKNNIntra(SMD_i, i, i);
    SMDIntra = smi_core.SingleMoleculeData.catSMD(SMDIntra, SMD_i_corrected, false);
end

% Align all datasets (inter-dataset)
[SMD_final, StatsInter] = DC.driftCorrectKNNInter(SMDIntra);
```

### Brightfield-Based Correction

When brightfield images are stored in the data file:

```matlab
% Brightfield images should be in 'FocusImages' group in HDF5
% Format: BFStruct(n).Data.PreSeqImages, PostSeqImages

% Apply brightfield correction
[SMD_corrected, BFStruct] = smi_core.DriftCorrection.driftCorrectBF(...
    SMD, SMF, [], [], []);

% First dataset pre-image used as reference by default
% Can specify custom reference image as 3rd argument
```

## Understanding Drift Correction Results

### Drift Fields in SMD

After drift correction, SMD contains:

```matlab
SMD.X           % Drift-corrected X positions
SMD.Y           % Drift-corrected Y positions
SMD.Z           % Drift-corrected Z positions (if 3D)

SMD.DriftX      % X drift estimate (NFrames × NDatasets)
SMD.DriftY      % Y drift estimate (NFrames × NDatasets)
SMD.DriftZ      % Z drift estimate (NFrames × NDatasets, if 3D)
```

### Sign Convention

Drift corrections follow this relationship:

```matlab
% For each localization:
i = SMD.FrameNum(k);
j = SMD.DatasetNum(k);

% Corrected = Drifted - Drift
SMD.X(k) = X_drifted(k) - SMD.DriftX(i, j);
SMD.Y(k) = Y_drifted(k) - SMD.DriftY(i, j);

% To recover drifted positions:
X_drifted(k) = SMD.X(k) + SMD.DriftX(i, j);
```

### Interpreting Drift Values

```matlab
% Convert to physical units (nanometers)
DriftX_nm = SMD.DriftX * SMD.PixelSize * 1000;
DriftY_nm = SMD.DriftY * SMD.PixelSize * 1000;

% Total drift magnitude over acquisition
total_drift_x = range(DriftX_nm(:, 1));
total_drift_y = range(DriftY_nm(:, 1));
fprintf('Total drift: X = %.1f nm, Y = %.1f nm\n', ...
    total_drift_x, total_drift_y);

% Drift rate (nm per frame)
drift_rate_x = mean(diff(DriftX_nm(:, 1)));
drift_rate_y = mean(diff(DriftY_nm(:, 1)));
fprintf('Drift rate: X = %.2f nm/frame, Y = %.2f nm/frame\n', ...
    drift_rate_x, drift_rate_y);
```

## Visualizing Drift Correction

### Plotting Drift Trajectories

```matlab
% Create drift correction object
DC = smi_core.DriftCorrection(SMF, SMD);

% Generate drift plot
DC_fig = DC.plotDriftCorrection(SMD, 'A');
% 'A' = absolute drift
% 'R' = relative drift between datasets
% '1' = initial values only

% Customize plot
figure(DC_fig);
title('Stage Drift During Acquisition');
xlabel('Drift X (nm)');
ylabel('Drift Y (nm)');
```

The plot shows:
- Color gradient (blue to red) = time progression
- Tapering line width = direction of time
- Separate trajectories per dataset

### Drift vs. Time

```matlab
% Extract drift for first dataset
frames = 1:SMD.NFrames;
drift_x = SMD.DriftX(:, 1) * SMD.PixelSize * 1000;  % nm
drift_y = SMD.DriftY(:, 1) * SMD.PixelSize * 1000;  % nm

% Plot drift components
figure;
subplot(2,1,1);
plot(frames, drift_x, 'b-', 'LineWidth', 1.5);
xlabel('Frame Number');
ylabel('X Drift (nm)');
title('X Drift vs. Time');
grid on;

subplot(2,1,2);
plot(frames, drift_y, 'r-', 'LineWidth', 1.5);
xlabel('Frame Number');
ylabel('Y Drift (nm)');
title('Y Drift vs. Time');
grid on;
```

### Before/After Comparison

```matlab
% Generate images before and after correction
SRImageZoom = 20;

% Before: Use raw positions (add drift back)
SMD_before = SMD;
for k = 1:length(SMD.X)
    i = SMD.FrameNum(k);
    j = SMD.DatasetNum(k);
    SMD_before.X(k) = SMD.X(k) + SMD.DriftX(i, j);
    SMD_before.Y(k) = SMD.Y(k) + SMD.DriftY(i, j);
end

img_before = smi_vis.GenerateImages.gaussImage(...
    [SMD_before.X, SMD_before.Y], SMD.PixelSize, SRImageZoom, ...
    [SMD.YSize, SMD.XSize]);

% After: Use corrected positions
img_after = smi_vis.GenerateImages.gaussImage(...
    [SMD.X, SMD.Y], SMD.PixelSize, SRImageZoom, ...
    [SMD.YSize, SMD.XSize]);

% Display side by side
figure;
subplot(1,2,1);
imagesc(img_before); axis image; colormap(gray);
title('Before Drift Correction');

subplot(1,2,2);
imagesc(img_after); axis image; colormap(gray);
title('After Drift Correction');
```

### Cumulative Drift Plot

```matlab
% Plot cumulative drift magnitude
DC_fig = smi_core.DriftCorrection.plotCumDrift(SMD, 'DriftX');
figure(DC_fig);
title('Cumulative X Drift');

DC_fig = smi_core.DriftCorrection.plotCumDrift(SMD, 'DriftY');
figure(DC_fig);
title('Cumulative Y Drift');
```

### Parametric Drift Plot

```matlab
% XY parametric plot showing drift path
DC_fig = smi_core.DriftCorrection.plotXYDriftParametric(SMD);
figure(DC_fig);
title('Drift Trajectory (Parametric)');
```

## Troubleshooting

### Issue: Drift correction fails or produces poor results

**Diagnose**:

```matlab
% Check localization density
locs_per_frame = histcounts(SMD.FrameNum, 1:SMD.NFrames+1);
mean_density = mean(locs_per_frame);
fprintf('Average localizations per frame: %.1f\n', mean_density);

% Plot density over time
figure;
plot(locs_per_frame);
xlabel('Frame'); ylabel('Localizations');
title('Localization Density');
```

**Solutions**:

- **Low density (<50/frame)**:
  - Increase `L_intra` to 1.5-2.0 for more smoothing
  - Consider using brightfield registration instead
  - Try reducing `PDegree` to 0 (constant drift)

- **Decreasing density (bleaching)**:
  - Split into smaller datasets with reorganization
  - Use brightfield registration
  - Increase `L_intra` to be more robust

- **Uneven density**:
  - Frame connection can help stabilize
  - Increase smoothing parameters

### Issue: Drift correction is too smooth (under-corrected)

**Diagnose**:

```matlab
% Check if drift model matches data
DC_fig = DC.plotDriftCorrection(SMD, 'A');

% Look for:
% - Smooth polynomial not capturing rapid changes
% - Obvious drift remaining in images
```

**Solutions**:

- Decrease `L_intra` to 0.5-0.7 (less smoothing, more responsive)
- Increase `PDegree` to 2 if drift is nonlinear
- Split acquisition into more datasets for better inter-dataset alignment

### Issue: Drift correction is too noisy (over-fitted)

**Diagnose**:

```matlab
% Examine drift trajectory
plot(SMD.DriftX(:, 1));
% Should be smooth, not jagged
```

**Solutions**:

- Increase `L_intra` to 1.5-2.5 (more smoothing)
- Check if frame connection was applied (helps stability)
- Reduce `PDegree` if set too high

### Issue: Inter-dataset alignment fails

**Diagnose**:

```matlab
% Visualize datasets separately
for i = 1:SMD.NDatasets
    mask = SMD.DatasetNum == i;
    figure;
    plot(SMD.X(mask), SMD.Y(mask), '.');
    title(sprintf('Dataset %d', i));
    axis equal;
end

% Look for misalignment between datasets
```

**Solutions**:

- Increase `L_inter` to 3.0 for more tolerance
- Use brightfield registration for independent alignment
- Check that datasets have overlapping FOV
- Verify sufficient common features between datasets

### Issue: Brightfield registration not working

**Diagnose**:

```matlab
% Check if brightfield images exist
try
    BFStruct = smi_core.LoadData.readH5File(...
        fullfile(SMF.Data.FileDir, SMF.Data.FileName{1}), 'FocusImages');
    fprintf('Found %d brightfield image sets\n', length(BFStruct));
catch
    fprintf('ERROR: No brightfield images in file\n');
end
```

**Solutions**:

- Verify 'FocusImages' group exists in HDF5 file
- Check that pre/post sequence images were acquired
- Ensure brightfield images have visible features
- Fall back to DC-KNN: `SMF.DriftCorrection.BFRegistration = false;`

### Issue: 3D drift correction fails

**Diagnose**:

```matlab
% Check Z coordinate units
fprintf('Z range: %.3f to %.3f\n', min(SMD.Z), max(SMD.Z));
fprintf('XY range: %.1f to %.1f pixels\n', min(SMD.X), max(SMD.X));
fprintf('PixelSizeZUnit: %.3f um\n', SMF.DriftCorrection.PixelSizeZUnit);
```

**Solutions**:

- Verify `PixelSizeZUnit` matches pixel size (typically 0.1 um)
- Check that Z coordinates are in micrometers
- Ensure 3D localization quality is sufficient
- Z drift correction requires higher precision than XY

## Advanced Topics

### Reorganizing Datasets

For optimal drift correction, you can reorganize frames into different dataset divisions:

```matlab
% Original: 2 datasets × 1000 frames = 2000 total frames
% Reorganize to: 10 datasets × 200 frames

SMF = smi_core.SingleMoleculeFitting();

% Specify new organization
SMF.DriftCorrection.NDatasets = 10;
% or equivalently:
% SMF.DriftCorrection.NFrames = 200;

% This internally reorganizes but preserves original numbering
DC = smi_core.DriftCorrection(SMF, SMD);
[SMD_corrected, Statistics] = DC.driftCorrectKNN(SMD);

% SMD_corrected.Internal contains reorganized numbering
% SMD_corrected still reports original NDatasets/NFrames
```

**When to use**:
- Very long acquisitions (>500 frames per dataset)
- Non-uniform drift within datasets
- Severe bleaching over time
- Better inter-dataset alignment needed

### Registration Via Drift Correction

Align two differently labeled datasets of the same structure:

```matlab
% Load two datasets (e.g., different color channels)
SMD1 = ...; % First label
SMD2 = ...; % Second label

% Compute optimal alignment
[delta12, Statistics] = smi_core.DriftCorrection.regViaDC(SMD1, SMD2);

% delta12 contains [deltaX, deltaY] shift to align SMD2 to SMD1
fprintf('Shift: X = %.3f pixels, Y = %.3f pixels\n', delta12(1), delta12(2));

% Apply shift
SMD2.X = SMD2.X + delta12(1);
SMD2.Y = SMD2.Y + delta12(2);
```

### Changing Inter-Dataset Reference

After drift correction, change which dataset is the reference:

```matlab
% Change reference to dataset 3
SMD_reref = smi_core.DriftCorrection.changeInterRef(SMD, 3);

% All coordinates now relative to dataset 3 as origin
```

### Custom Drift Models

For specialized needs, you can modify the drift model:

```matlab
% Access internal drift minimization function
DC = smi_core.DriftCorrection(SMF, SMD);

% Custom cost function for drift fitting
% See DriftCorrection.m, minD() method
% Modify optimization parameters:
DC.TolFun_intra = 1e-3;  % Tighter tolerance
DC.TolX_intra = 1e-5;    % Higher precision
```

## Best Practices

### Data Acquisition

1. **Include brightfield images**: Acquire before/after each dataset for robust registration
2. **Minimize drift sources**: Temperature stabilization, vibration isolation
3. **Optimize acquisition time**: Balance between statistics and drift accumulation
4. **Use fiducials when possible**: Gold nanoparticles provide ground truth

### Parameter Selection

1. **Start with defaults**: Work well for typical data
2. **Adjust based on density**: Higher density → lower L_intra
3. **Monitor drift plots**: Verify smoothness and realism
4. **Compare methods**: Try both DC-KNN and brightfield, use better result

### Processing Strategy

1. **Apply frame connection first**: Improves stability for DC-KNN
2. **Check localization quality**: Poor localizations corrupt drift estimation
3. **Visualize before/after**: Confirm improvement
4. **Save drift data**: Keep SMD.DriftX/Y for later analysis

### Quality Control

1. **Inspect drift magnitude**: Should be <2-3 pixels total
2. **Check drift smoothness**: Should be continuous, not jumpy
3. **Verify correction**: Images should look sharper
4. **Compare datasets**: Should align within localization precision

## Performance Considerations

### Computation Time

DC-KNN scales with:
- Number of localizations: O(N log N) for neighbor search
- Number of frames: O(F) for polynomial fitting
- Number of iterations: Typically 10-50 fminsearch iterations

Typical times (on modern CPU):
- 100,000 localizations, 1000 frames: ~5-10 seconds
- 1,000,000 localizations, 5000 frames: ~30-60 seconds

Brightfield registration is faster:
- Independent of localization count
- Dominated by image cross-correlation
- Typical time: <1 second per dataset

### Memory Usage

DC-KNN requires:
- Neighbor search structure: ~50 MB per 100,000 localizations
- Temporary arrays: ~20 MB
- Total: Modest, not limiting for typical datasets

### Optimization Tips

1. **Use frame connection**: Reduces data size, improves stability
2. **Filter low-quality localizations**: Speeds processing, improves results
3. **Split very large datasets**: Process in chunks if memory limited
4. **Parallel processing**: Process multiple datasets independently

## Citation

If you use drift correction in publications, please cite:

Michael J. Wester, David J. Schodt, Hanieh Mazloom-Farsibaf, Mohamadreza Fazel, Sandeep Pallikkuth and Keith A. Lidke, "Robust, fiducial-free drift correction for super-resolution imaging", Scientific Reports, Volume 11, Article 23672, December 8, 2021, 1-14, [https://www.nature.com/articles/s41598-021-02850-7](https://www.nature.com/articles/s41598-021-02850-7) (DOI: 10.1038/s41598-021-02850-7).

## See Also

- [SMLM Analysis Workflow](smlm-analysis.md) - Complete pipeline including drift correction
- [SMF Structure](../core-concepts/smf-structure.md) - All drift correction parameters
- [SMD Structure](../core-concepts/smd-structure.md) - Drift fields in results
- MATLAB/+smi_core/@DriftCorrection/README.md - Technical documentation
- MATLAB/+smi_core/@DriftCorrection/driftCorrectKNN.m - Algorithm implementation
