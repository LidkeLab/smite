---
title: "SMD Structure: Localization Results"
category: "core-concepts"
level: "intermediate"
tags: ["smd", "results", "data-structures", "localizations"]
prerequisites: ["architecture.md", "smf-structure.md"]
related: ["../workflows/smlm-analysis.md", "../how-to/localize-molecules.md"]
summary: "Complete reference to the Single Molecule Data (SMD) structure that contains all localization results"
estimated_time: "15 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# SMD Structure: Localization Results

## Purpose

The SMD (Single Molecule Data) structure is smite's standard format for storing localization results. Every analysis produces an SMD containing positions, uncertainties, photon counts, and metadata for all detected molecules. Understanding SMD is essential for analyzing results, creating visualizations, and performing downstream analysis like clustering or tracking.

## Prerequisites

- Understanding of [smite architecture](architecture.md)
- Familiarity with [SMF structure](smf-structure.md)
- Basic SMLM/SPT concepts

## Overview

SMD is a structure (not an array) where each field is a vector containing values for all localizations. Think of it as a table where:
- Each **row** is a localization
- Each **column** is a property (X, Y, Photons, etc.)
- **Number of rows** = `length(SMD.X)`

Key characteristics:
- **Flat structure**: All fields are vectors of the same length
- **Standard format**: Same fields across all smite analyses
- **Self-contained**: Includes metadata about acquisition
- **Easy to filter**: Use logical indexing on any field

```matlab
% After analysis
load('Results.mat', 'SMD');

% How many localizations?
N = length(SMD.X);

% Plot positions
plot(SMD.X, SMD.Y, '.');
```

## Creating and Accessing SMD

### From Analysis

SMD is typically created by analysis workflows:

```matlab
% From SMLM analysis
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();
SMD = SMLMobj.SMD;  % SMD is now populated

% From LocalizeData
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();

% From tracking
SPT = smi.SPT(SMF);
SPT.performFullAnalysis();
load(resultsFile, 'SMD', 'TR');  % Both SMD and TR available
```

### Loading from Files

Results are saved as .mat files:

```matlab
% Load from Results.mat
load('FileDir/Results/Results.mat', 'SMD', 'SMF');

% Now access fields
positions = [SMD.X, SMD.Y];
photons = SMD.Photons;
```

### Creating Empty SMD

For custom workflows:

```matlab
% Create empty SMD with correct structure
SMD = smi_core.SingleMoleculeData.createSMD();

% Manually populate
SMD.X = [10.5; 20.3; 15.7];  % 3 localizations
SMD.Y = [12.1; 18.9; 25.4];
SMD.Photons = [1000; 850; 1200];
% ... etc.
```

## SMD Fields Reference

### Metadata Fields

Information about the dataset:

| Field | Type | Description |
|-------|------|-------------|
| `NDims` | integer | Dimensionality (2 or 3) |
| `NFrames` | integer | Total frames in raw data |
| `NDatasets` | integer | Number of image stacks |
| `FrameRate` | scalar | Acquisition frame rate (Hz) |
| `PixelSize` | scalar | Camera pixel size (micrometers) |
| `XSize` | integer | Image width (pixels) |
| `YSize` | integer | Image height (pixels) |
| `ZOffset` | scalar | Z position of focal plane |

**Example:**

```matlab
% Check acquisition parameters
fprintf('Acquired %d frames at %.1f Hz\n', SMD.NFrames, SMD.FrameRate);
fprintf('Image size: %d × %d pixels\n', SMD.XSize, SMD.YSize);
fprintf('Pixel size: %.1f nm\n', SMD.PixelSize * 1000);
```

### Position Fields

Localization coordinates and uncertainties:

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `X` | vector (N×1) | pixels | X position estimates |
| `Y` | vector (N×1) | pixels | Y position estimates |
| `Z` | vector (N×1) | micrometers | Z position estimates (3D only) |
| `X_SE` | vector (N×1) | pixels | Standard error in X |
| `Y_SE` | vector (N×1) | pixels | Standard error in Y |
| `Z_SE` | vector (N×1) | micrometers | Standard error in Z (3D only) |

**Position coordinates:**
- Origin: (1, 1) is center of top-left pixel
- Sub-pixel precision: Positions are real numbers, not integers
- Already drift-corrected if `SMF.DriftCorrection.On = true`

**Standard errors (uncertainties):**
- Estimated from Cramér-Rao lower bound (CRLB)
- Smaller = more precise localization
- Typical values: 0.05-0.2 pixels (5-20 nm for 100 nm pixels)

**Example:**

```matlab
% Convert to nanometers
X_nm = SMD.X * SMD.PixelSize * 1000;
Y_nm = SMD.Y * SMD.PixelSize * 1000;
X_SE_nm = SMD.X_SE * SMD.PixelSize * 1000;

% Plot with error bars (sample 1000 points)
idx = randsample(length(SMD.X), min(1000, length(SMD.X)));
errorbar(X_nm(idx), Y_nm(idx), X_SE_nm(idx), 'horizontal', '.');

% Check precision distribution
median_precision_nm = median(SMD.X_SE) * SMD.PixelSize * 1000;
fprintf('Median precision: %.1f nm\n', median_precision_nm);
```

### Photometric Fields

Intensity and background information:

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `Photons` | vector (N×1) | photons | Total collected photons |
| `Bg` | vector (N×1) | photons/pixel | Background level |
| `Photons_SE` | vector (N×1) | photons | Standard error in photons |
| `Bg_SE` | vector (N×1) | photons/pixel | Standard error in background |

**Photons**: Total photons detected from emitter (after camera conversion)

**Background**: Photons per pixel from background (autofluorescence, scatter, etc.)

**Example:**

```matlab
% Photon count statistics
fprintf('Photon counts: %.0f ± %.0f (mean ± std)\n', ...
    mean(SMD.Photons), std(SMD.Photons));

% Signal-to-background ratio
SBR = SMD.Photons ./ (SMD.Bg * pi * SMF.Fitting.PSFSigma^2);
fprintf('Median SBR: %.1f\n', median(SBR));

% Plot photon histogram
histogram(SMD.Photons, 50);
xlabel('Photons'); ylabel('Count');
title('Photon Distribution');
```

### PSF Fields

Point spread function parameters:

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `PSFSigma` | vector (N×1) | pixels | PSF width (symmetric) |
| `PSFSigmaX` | vector (N×1) | pixels | PSF width in X (asymmetric) |
| `PSFSigmaY` | vector (N×1) | pixels | PSF width in Y (asymmetric) |
| `PSFSigma_SE` | vector (N×1) | pixels | Standard error in PSF sigma |
| `PSFSigmaX_SE` | vector (N×1) | pixels | Standard error in X sigma |
| `PSFSigmaY_SE` | vector (N×1) | pixels | Standard error in Y sigma |

**Note**: Only populated if PSF sigma was a fit parameter (`FitType='XYNBS'` or similar). Otherwise, these fields are empty or contain the fixed value from `SMF.Fitting.PSFSigma`.

**Example:**

```matlab
% Check if PSF was fit
if ~isempty(SMD.PSFSigma) && std(SMD.PSFSigma) > 0.01
    % PSF sigma was fitted
    histogram(SMD.PSFSigma, 30);
    xlabel('PSF Sigma (pixels)');
    title('Fitted PSF Width Distribution');

    % Check for aberrations (asymmetric PSF)
    if isfield(SMD, 'PSFSigmaX') && ~isempty(SMD.PSFSigmaX)
        scatter(SMD.PSFSigmaX, SMD.PSFSigmaY, '.');
        xlabel('Sigma X'); ylabel('Sigma Y');
        axis equal;
    end
end
```

### Frame and Dataset Fields

Temporal and organizational information:

| Field | Type | Description |
|-------|------|-------------|
| `FrameNum` | vector (N×1) | Frame number (1 to NFrames) |
| `DatasetNum` | vector (N×1) | Dataset number (1 to NDatasets) |
| `XBoxCorner` | vector (N×1) | X coordinate of fitting box corner |
| `YBoxCorner` | vector (N×1) | Y coordinate of fitting box corner |

**FrameNum**: Which frame the localization came from

**DatasetNum**: Which dataset (if analyzing multiple .h5 files)

**Box corners**: Top-left corner of the fitting box

**Example:**

```matlab
% Localizations per frame
frames = 1:SMD.NFrames;
counts = histcounts(SMD.FrameNum, [frames, SMD.NFrames+1]);
plot(frames, counts);
xlabel('Frame'); ylabel('Localizations');
title('Localizations per Frame');

% Filter by dataset
if SMD.NDatasets > 1
    dataset1_idx = SMD.DatasetNum == 1;
    SMD_dataset1.X = SMD.X(dataset1_idx);
    SMD_dataset1.Y = SMD.Y(dataset1_idx);
end

% Temporal filtering
early_frames = SMD.FrameNum <= 100;
late_frames = SMD.FrameNum > 100;
figure; subplot(1,2,1);
plot(SMD.X(early_frames), SMD.Y(early_frames), '.');
title('Frames 1-100');
subplot(1,2,2);
plot(SMD.X(late_frames), SMD.Y(late_frames), '.');
title('Frames 101+');
```

### Quality and Fit Fields

Fit quality and filtering information:

| Field | Type | Description |
|-------|------|-------------|
| `PValue` | vector (N×1) | P-value of fit (goodness of fit) |
| `LogLikelihood` | vector (N×1) | Log-likelihood of fit |
| `ThreshFlag` | vector (N×1) | Threshold flag (0=pass, >0=fail) |

**PValue**: Statistical measure of fit quality
- Range: 0 to 1
- Higher = better fit
- Threshold typically 0.01-0.05

**LogLikelihood**: Log probability of data given model
- Higher = better fit
- Used for automatic thresholding

**ThreshFlag**: Indicates which threshold failed
- 0: Passed all thresholds
- >0: Failed (specific codes indicate which threshold)

**Example:**

```matlab
% Fit quality assessment
fprintf('P-values: %.3f (median), %.3f (min)\n', ...
    median(SMD.PValue), min(SMD.PValue));

% Identify poor fits
poor_fits = SMD.PValue < 0.01;
fprintf('%.1f%% of fits below p=0.01\n', 100*mean(poor_fits));

% Check thresholding
if isfield(SMD, 'ThreshFlag')
    passed = sum(SMD.ThreshFlag == 0);
    total = length(SMD.ThreshFlag);
    fprintf('%d / %d (%.1f%%) passed thresholds\n', ...
        passed, total, 100*passed/total);
end
```

### Frame Connection Fields

Linking localizations across frames:

| Field | Type | Description |
|-------|------|-------------|
| `ConnectID` | vector (N×1) | Connection identifier |
| `IndSMD` | cell array | Indices in original SMD pre-frame-connection |

**ConnectID**: Integer identifying same emitter across frames
- Same ID = same molecule
- Different ID = different molecules
- 0 or negative = not connected

**IndSMD**: Maps back to pre-frame-connection SMD

**Example:**

```matlab
% How many unique emitters?
unique_emitters = unique(SMD.ConnectID);
unique_emitters = unique_emitters(unique_emitters > 0);
fprintf('Detected %d unique emitters\n', length(unique_emitters));

% Average appearances per emitter
appearances = histcounts(SMD.ConnectID, [unique_emitters; max(unique_emitters)+1]);
fprintf('Average appearances: %.1f frames\n', mean(appearances));

% Find emitters appearing many times
long_blinkers = unique_emitters(appearances > 10);
fprintf('%d emitters appeared >10 times\n', length(long_blinkers));

% Combine localizations from same emitter for precision improvement
emitter_id = long_blinkers(1);  % Pick one
same_emitter = SMD.ConnectID == emitter_id;
X_mean = mean(SMD.X(same_emitter));
Y_mean = mean(SMD.Y(same_emitter));
X_SE_combined = std(SMD.X(same_emitter)) / sqrt(sum(same_emitter));
fprintf('Emitter %d: (%.2f, %.2f) ± %.3f pixels\n', ...
    emitter_id, X_mean, Y_mean, X_SE_combined);
```

### Drift Correction Fields

Stage drift estimates:

| Field | Type | Units | Description |
|-------|------|-------|-------------|
| `DriftX` | matrix (NFrames×NDatasets) | pixels | X drift per frame |
| `DriftY` | matrix (NFrames×NDatasets) | pixels | Y drift per frame |
| `DriftZ` | matrix (NFrames×NDatasets) | micrometers | Z drift per frame |

**Note**: Drift values are cumulative offsets from the first frame. Positions in `SMD.X`, `SMD.Y` are already corrected if `SMF.DriftCorrection.On = true`.

**Example:**

```matlab
% Plot drift trajectory
if isfield(SMD, 'DriftX') && ~isempty(SMD.DriftX)
    figure;
    plot(SMD.DriftX(:,1), SMD.DriftY(:,1), '-');
    xlabel('Drift X (pixels)'); ylabel('Drift Y (pixels)');
    title('Stage Drift Trajectory');
    axis equal;

    % Drift magnitude over time
    figure;
    drift_mag = sqrt(SMD.DriftX(:,1).^2 + SMD.DriftY(:,1).^2);
    plot(drift_mag * SMD.PixelSize * 1000);  % Convert to nm
    xlabel('Frame'); ylabel('Cumulative Drift (nm)');
    title('Drift Magnitude');
end
```

### Channel Registration Fields

Multi-color registration information:

| Field | Type | Description |
|-------|------|-------------|
| `IsTransformed` | logical | Was channel registration applied? |
| `RegError` | scalar | Registration error (pixels) |

**Example:**

```matlab
if isfield(SMD, 'IsTransformed') && SMD.IsTransformed
    fprintf('Channel registration applied\n');
    fprintf('Registration error: %.3f pixels (%.1f nm)\n', ...
        SMD.RegError, SMD.RegError * SMD.PixelSize * 1000);
end
```

## Common Operations

### Filtering SMD

Create filtered SMD based on criteria:

```matlab
% Filter by precision
good_precision = SMD.X_SE < 0.15 & SMD.Y_SE < 0.15;
SMD_filtered.X = SMD.X(good_precision);
SMD_filtered.Y = SMD.Y(good_precision);
SMD_filtered.Photons = SMD.Photons(good_precision);
% ... copy other fields ...

% Filter by photons
bright = SMD.Photons > 500;

% Filter by region of interest
ROI = [20, 20, 100, 100];  % [YStart, XStart, YEnd, XEnd]
in_ROI = SMD.X >= ROI(2) & SMD.X <= ROI(4) & ...
         SMD.Y >= ROI(1) & SMD.Y <= ROI(3);

% Combine filters
selected = good_precision & bright & in_ROI;
```

### Combining SMD Structures

Merge multiple SMD structures:

```matlab
% Using built-in method
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2, SMD3);

% Manual concatenation
SMD_all.X = [SMD1.X; SMD2.X; SMD3.X];
SMD_all.Y = [SMD1.Y; SMD2.Y; SMD3.Y];
% ... etc for all fields ...
```

### Exporting SMD

Convert to other formats:

```matlab
% To table (for easier manipulation)
T = table(SMD.X, SMD.Y, SMD.Photons, SMD.X_SE, SMD.Y_SE, SMD.FrameNum, ...
    'VariableNames', {'X', 'Y', 'Photons', 'X_SE', 'Y_SE', 'Frame'});

% To CSV
writetable(T, 'localizations.csv');

% To array for plotting
positions = [SMD.X, SMD.Y];
scatter(positions(:,1), positions(:,2), 5, SMD.Photons, 'filled');
colorbar; xlabel('X (pixels)'); ylabel('Y (pixels)');
```

### Computing Derived Quantities

Calculate additional metrics:

```matlab
% Localization precision (combined X and Y)
precision = sqrt(SMD.X_SE.^2 + SMD.Y_SE.^2);

% Signal-to-noise ratio
PSF_area = pi * SMF.Fitting.PSFSigma^2;
SNR = SMD.Photons ./ sqrt(SMD.Photons + SMD.Bg * PSF_area);

% Localization density (localizations per μm²)
area_um2 = SMD.XSize * SMD.YSize * SMD.PixelSize^2;
density = length(SMD.X) / area_um2;
fprintf('Density: %.1f localizations/μm²\n', density);

% Time in seconds
time_s = (SMD.FrameNum - 1) / SMD.FrameRate;
```

## Relationship to TR (Tracking Results)

For SPT analyses, results are organized as **TR** (Tracking Results), which is an array of SMD structures:

```matlab
% TR is an array where each element is one trajectory
load('SPT_Results.mat', 'TR');

% Number of trajectories
N_traj = length(TR);

% Each TR(i) is an SMD for trajectory i
traj_1_X = TR(1).X;  % X positions for first trajectory
traj_1_Y = TR(1).Y;  % Y positions for first trajectory

% Plot all trajectories
figure; hold on;
for i = 1:length(TR)
    plot(TR(i).X, TR(i).Y, '-');
end
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('%d Trajectories', length(TR)));
```

See tracking documentation for more details on TR.

## Complete Example

Comprehensive SMD analysis:

```matlab
% Load results
load('Results.mat', 'SMD', 'SMF');

% Basic statistics
fprintf('=== Dataset Summary ===\n');
fprintf('Localizations: %d\n', length(SMD.X));
fprintf('Frames: %d\n', SMD.NFrames);
fprintf('Datasets: %d\n', SMD.NDatasets);
fprintf('Image size: %d × %d pixels\n', SMD.XSize, SMD.YSize);

fprintf('\n=== Localization Quality ===\n');
fprintf('X precision: %.1f ± %.1f nm (mean ± std)\n', ...
    mean(SMD.X_SE)*SMD.PixelSize*1000, std(SMD.X_SE)*SMD.PixelSize*1000);
fprintf('Photons: %.0f ± %.0f\n', mean(SMD.Photons), std(SMD.Photons));
fprintf('Background: %.1f ± %.1f photons/pixel\n', mean(SMD.Bg), std(SMD.Bg));

% Frame connection statistics
if isfield(SMD, 'ConnectID')
    unique_IDs = unique(SMD.ConnectID);
    unique_IDs = unique_IDs(unique_IDs > 0);
    fprintf('\n=== Frame Connection ===\n');
    fprintf('Unique emitters: %d\n', length(unique_IDs));
    fprintf('Compression ratio: %.1f:1\n', length(SMD.X) / length(unique_IDs));
end

% Drift correction
if isfield(SMD, 'DriftX') && ~isempty(SMD.DriftX)
    total_drift = sqrt(SMD.DriftX(end,1)^2 + SMD.DriftY(end,1)^2);
    fprintf('\n=== Drift Correction ===\n');
    fprintf('Total drift: %.1f nm\n', total_drift * SMD.PixelSize * 1000);
end

% Create super-resolution image
SR_zoom = 20;  % 20× zoom
SR_image = smi_vis.GenerateImages.gaussImage([SMD.X, SMD.Y], ...
    SMD.PixelSize, SR_zoom, [SMD.YSize, SMD.XSize]);
figure; imagesc(SR_image); axis image; colormap hot;
title('Super-Resolution Image');
```

## See Also

- [SMF Structure](smf-structure.md) - Analysis parameters
- [Architecture Overview](architecture.md) - How SMD fits into smite
- [SMLM Workflow](../workflows/smlm-analysis.md) - Generating SMD
- [How to Localize Molecules](../how-to/localize-molecules.md) - Creating SMD
- doc/DataStructures/SMD.md (in repository) - Complete field reference
- doc/DataStructures/TR.md - Tracking results format
