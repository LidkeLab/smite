---
title: "ROI Analysis in SMLM"
category: "workflows"
level: "intermediate"
tags: ["roi", "region-of-interest", "spatial-analysis", "subregion", "segmentation"]
prerequisites: ["../core-concepts/smd-structure.md", "smlm-analysis.md"]
related: ["../how-to/visualize-results.md", "../examples/spatial-analysis.md"]
summary: "Comprehensive guide to selecting, extracting, and analyzing Regions of Interest (ROIs) for spatially-restricted SMLM analysis"
estimated_time: "20 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# ROI Analysis in SMLM

## Purpose

Region of Interest (ROI) analysis allows you to spatially restrict your SMLM analysis to specific cellular structures, organelles, or features. This document explains how to define ROIs using automatic or manual methods, extract localization data within ROIs, perform ROI-based analysis, and combine results from multiple ROIs. ROI analysis is essential for quantifying spatial heterogeneity, comparing different cellular regions, and reducing computational complexity for large datasets.

## Prerequisites

- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Completion of [SMLM analysis workflow](smlm-analysis.md)
- Basic understanding of MATLAB plotting and image display
- Familiarity with coordinate systems in microscopy

## Overview

smite provides two complementary approaches to ROI analysis:

1. **Automatic ROI finding** (`smi_core.FindROI`): GPU-accelerated detection of candidate molecule locations during localization
2. **Manual ROI selection** (`smi_helpers.ROITools`): Interactive selection of spatial regions from localization or image data

Both approaches integrate seamlessly with the SMLM workflow, allowing you to:
- Reduce computational load by analyzing only relevant regions
- Quantify spatial heterogeneity across different cellular compartments
- Compare multiple structures within the same dataset
- Combine results from separately processed ROIs

## ROI Methods Overview

### Automatic ROI Finding: smi_core.FindROI

**Purpose**: Rapidly identify candidate emitter locations in raw image stacks using GPU-accelerated filtering.

**How it works**:
1. Apply difference-of-Gaussians (DoG) filter to raw images
2. Find local maxima above photon threshold
3. Extract small boxes around each candidate location
4. Pass boxes to fitting algorithm for precise localization

**When to use**:
- During standard SMLM localization workflow
- High-density data with many emitters per frame
- Need fast, automated box finding
- GPU available for acceleration

**Key features**:
- GPU-accelerated using CUDA
- Recursive Gaussian filtering (Young-van Vliet method)
- Adjustable box size and overlap
- Minimum photon threshold filtering

### Manual ROI Selection: smi_helpers.ROITools

**Purpose**: Interactively select spatial regions from localization data or super-resolution images.

**How it works**:
1. Display localizations as point cloud or Gaussian image
2. Use mouse to draw rectangular ROIs
3. Extract coordinates within each ROI
4. Save ROI definitions for later use

**When to use**:
- Analyzing specific cellular structures
- Comparing different regions within cells
- Processing dense regions separately
- Multi-label spatial analysis

**Key features**:
- Interactive GUI-based selection
- Point cloud or Gaussian image display
- Fixed-size or adjustable rectangular regions
- Multi-label ROI compatibility
- Coordinate transformations for multi-channel data

## Automatic ROI Finding Workflow

### Basic Configuration

Automatic ROI finding is configured through `SMF.BoxFinding` parameters:

```matlab
% Create SMF and configure box finding
SMF = smi_core.SingleMoleculeFitting();

% Box finding parameters
SMF.BoxFinding.BoxSize = 7;           % Box size in pixels (odd number)
SMF.BoxFinding.BoxOverlap = 2;        % Allowed overlap between boxes (pixels)
SMF.BoxFinding.MinPhotons = 200;      % Minimum photons to be candidate

% Fitting parameters (affects ROI detection threshold)
SMF.Fitting.PSFSigma = [1.3, 1.3];    % Expected PSF width (pixels)
```

### Parameter Guide

#### BoxSize (integer, pixels)
- Linear dimension of extracted boxes
- Must be odd number (7, 9, 11, etc.)
- Larger boxes: capture more background, slower fitting
- Smaller boxes: faster but may clip PSF tails
- Typical range: 5-11 pixels
- **Recommendation**: Start with 7 pixels

#### BoxOverlap (integer, pixels)
- Maximum allowed overlap between adjacent boxes
- Prevents duplicate detections of same emitter
- Should be less than BoxSize/2
- Larger overlap: more redundant detections
- Smaller overlap: may miss closely-spaced emitters
- **Recommendation**: 1-3 pixels

#### MinPhotons (scalar, photons)
- Minimum integrated photons for candidate detection
- Higher threshold: fewer false positives, may miss dim emitters
- Lower threshold: more candidates, more computational load
- Depends on signal-to-noise ratio
- **Recommendation**: 100-300 photons for typical STORM data

### Using FindROI Directly

For more control or debugging, use `FindROI` independently:

```matlab
% Load raw data
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'data.h5'};

LD = smi_core.LoadData();
[~, RawData] = LD.loadRawData(SMF, 1);

% Convert to photons
DTP = smi_core.DataToPhotons(SMF);
DTP.RawData = RawData;
ImageStack = DTP.convertData();

% Create FindROI object and find boxes
FR = smi_core.FindROI(SMF, ImageStack);
FR.Verbose = 3;  % Show visualization
[ROIStack, SMD] = FR.findROI();

% ROIStack: BoxSize x BoxSize x NBoxes array of image boxes
% SMD: Contains XBoxCorner, YBoxCorner, FrameNum for each box
```

### Visualizing Detected ROIs

```matlab
% Plot boxes on a specific frame
Frame = 10;
FR.plotBox(SMD, ImageStack, Frame, SMF.BoxFinding.BoxSize);
title(sprintf('Detected Boxes in Frame %d', Frame));

% Create video showing boxes on all frames (requires Verbose >= 3)
FR.ResultsDir = fullfile(SMF.Data.FileDir, 'Results');
[ROIStack, SMD] = FR.findROI();
% Video saved as 'RawDataFromSMD.mp4' in ResultsDir
```

### Extracting ROIs from Existing SMD

After localization, you can reload raw data and extract the original ROI boxes:

```matlab
% Load existing SMD from analysis
SMD = load(fullfile(SMF.Data.FileDir, 'Results', 'SMD.mat'));
SMD = SMD.SMD;

% Extract ROI boxes from raw data
[ROIStack, CameraGain, CameraOffset, CameraReadNoise] = ...
    smi_core.FindROI.extractROIs(SMD, SMF, true);

% ROIStack(y, x, n) = nth ROI box
% Corrected for gain/offset if third argument is true

% Access individual boxes
BoxIndex = 100;
figure;
imagesc(ROIStack(:, :, BoxIndex));
axis image; colormap(gray);
title(sprintf('ROI %d (Frame %d)', BoxIndex, SMD.FrameNum(BoxIndex)));
```

## Manual ROI Selection Workflow

### Basic ROI Selection from SMD

The simplest way to select ROIs from localization data:

```matlab
% Load SMD from analysis
SMD = load(fullfile(SMF.Data.FileDir, 'Results', 'SMD.mat'));
SMD = SMD.SMD;

% Create ROITools object and configure
RT = smi_helpers.ROITools();
RT.ROI_sizes = [3000, 3000];  % Fixed ROI size in nm
RT.Pixel2nm = SMD.PixelSize * 1000;  % Pixel size in nm

% Select ROIs interactively
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Select ROIs for Analysis');

% RoI is a cell array with one element per ROI:
% RoI{i}.ROI = [xmin, xmax, ymin, ymax] in nm
% RoI{i}.X{1} = X coordinates in ROI (nm)
% RoI{i}.Y{1} = Y coordinates in ROI (nm)
% RoI{i}.X_SE{1} = X localization errors (nm)
% RoI{i}.Y_SE{1} = Y localization errors (nm)
% RoI{i}.SMD{1} = SMD structure for ROI
```

### Interactive Controls

When the ROI selection window appears:

**Mouse controls**:
- **Left click**: Place fixed-size rectangular ROI (size set by `ROI_sizes`)
- **Right click**: Draw adjustable rectangular ROI (drag to desired size)
- **Middle click**: Ignored (no action)

**Keyboard controls**:
- **Backspace** or **Delete**: Remove the most recently added ROI
- **Any other key**: Finish selection and close window

**Visual feedback**:
- ROIs outlined in red
- ROI numbers displayed at center of each box
- All localizations shown as colored points

### Display Options

#### Point Cloud Display (default)

```matlab
RT = smi_helpers.ROITools();
RT.GaussIm = false;        % Use point display
RT.OriginLLvsUL = true;    % Lower-left origin
RT.Msize = 7;              % Marker size

[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Point Cloud Display');
```

#### Gaussian Image Display

For single-label data, display as rendered super-resolution image:

```matlab
RT = smi_helpers.ROITools();
RT.GaussIm = true;          % Use Gaussian rendering
RT.OriginLLvsUL = false;    % Upper-left origin (consistent with images)
RT.SRzoom = 20;             % Super-resolution zoom factor

[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Gaussian Image Display');
```

**Note**: `GaussIm = true` currently only works for single-label data. Use point display for multi-label ROI selection.

### Multi-Label ROI Selection

Select ROIs that contain localizations from multiple labels (e.g., different color channels):

```matlab
% Load SMD for two labels
SMD1 = load('Cell_01_Label_01_Results/SMD.mat');
SMD1 = SMD1.SMD;
SMD2 = load('Cell_01_Label_02_Results/SMD.mat');
SMD2 = SMD2.SMD;

% Configure ROITools for multi-label
RT = smi_helpers.ROITools();
RT.Color = ['g', 'm'];      % Colors for each label
RT.Order = [1, 2];          % Plotting order
RT.ROI_sizes = [2000, 2000];

% Select ROIs from both labels
[n_ROIs, RoI, XYsize] = RT.getROI({SMD1, SMD2}, 'Dual-Color ROIs');

% Access multi-label data
for i = 1:n_ROIs
    fprintf('ROI %d:\n', i);
    fprintf('  Label 1: %d localizations\n', length(RoI{i}.X{1}));
    fprintf('  Label 2: %d localizations\n', length(RoI{i}.X{2}));

    % Extract SMD for each label
    SMD1_ROI = RoI{i}.SMD{1};
    SMD2_ROI = RoI{i}.SMD{2};
end
```

**Important**: ROIs are only accepted if all labels have at least one localization within the region. Empty ROIs are automatically rejected.

### Channel Registration and Coordinate Transforms

For multi-color data with chromatic aberration correction:

```matlab
RT = smi_helpers.ROITools();

% Define transformation from label 2 to label 1 coordinates
% Transform{1} is identity (label 1 is reference)
RT.Transform{1} = [];
% Transform{2} is affine transformation matrix
RT.Transform{2} = ChannelRegTransform;  % From channel registration

[n_ROIs, RoI, XYsize] = RT.getROI({SMD1, SMD2}, 'Registered Channels');
```

### Using Coordinate Masks

Restrict ROI selection to specific regions using a binary mask:

```matlab
RT = smi_helpers.ROITools();

% Create or load binary mask (logical array)
% mask(y, x) = true for valid regions
mask = SMD.Y > 5 & SMD.Y < 15 & SMD.X > 10 & SMD.X < 20;

% Apply mask
RT.Mask{1} = mask;  % One mask per label

[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Masked ROIs');
% Only localizations where mask == true will be displayed/selected
```

## Extracting Localizations from ROIs

### From RoI Structure

The `RoI` cell array from `ROITools.getROI()` contains complete information:

```matlab
% Get ROIs
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Analysis ROIs');

% Extract localization data for specific ROI
roi_index = 3;
X_nm = RoI{roi_index}.X{1};        % X coordinates (nm)
Y_nm = RoI{roi_index}.Y{1};        % Y coordinates (nm)
X_SE_nm = RoI{roi_index}.X_SE{1};  % X errors (nm)
Y_SE_nm = RoI{roi_index}.Y_SE{1};  % Y errors (nm)

% Get full SMD structure for this ROI
SMD_ROI = RoI{roi_index}.SMD{1};

% SMD_ROI contains all fields (Photons, Bg, FrameNum, etc.)
fprintf('ROI %d: %d localizations\n', roi_index, length(SMD_ROI.X));
fprintf('  Mean photons: %.1f\n', mean(SMD_ROI.Photons));
fprintf('  Frames: %d to %d\n', min(SMD_ROI.FrameNum), max(SMD_ROI.FrameNum));
```

### Directly from SMD Using isolateSubROI

Extract ROI directly from SMD without interactive selection:

```matlab
% Define ROI boundaries in pixels
% ROI = [YStart, XStart, YEnd, XEnd]
ROI_pixels = [50, 100, 150, 200];  % Y: 50-150, X: 100-200

% Extract localizations within ROI
SMD_ROI = smi_core.SingleMoleculeData.isolateSubROI(SMD, ROI_pixels);

% SMD_ROI contains only localizations within the specified region
fprintf('Original: %d localizations\n', length(SMD.X));
fprintf('In ROI: %d localizations\n', length(SMD_ROI.X));
```

**Coordinate convention**:
- Pixel edges: left edge of pixel n has coordinate n-1, right edge has coordinate n
- ROI boundaries are inclusive
- Example: ROI = [1, 1, 10, 10] includes pixels (1,1) through (10,10)

### Using Boolean Masks

Extract arbitrary subsets using logical indexing:

```matlab
% Create boolean mask for desired localizations
% Example: High-photon, central region localizations
mask = (SMD.Photons > 500) & ...
       (SMD.X > 50) & (SMD.X < 200) & ...
       (SMD.Y > 50) & (SMD.Y < 200);

% Extract subset
SMD_subset = smi_core.SingleMoleculeData.isolateSubSMD(SMD, mask);

% Equivalent to:
SMD_subset.X = SMD.X(mask);
SMD_subset.Y = SMD.Y(mask);
SMD_subset.Photons = SMD.Photons(mask);
% ... (all vector fields indexed by mask)
```

## ROI-Based Analysis Workflows

### Processing Dense Regions Separately

For very dense datasets, process ROIs independently to reduce memory and computational load:

```matlab
% Step 1: Initial localization to get overview
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'dense_data.h5'};

% Use permissive thresholds for overview
SMF.Thresholding.MinPhotons = 100;
SMLMobj = smi.SMLM(SMF);
SMLMobj.genLocalizations();
SMD_overview = SMLMobj.SMD;

% Step 2: Select ROIs from overview
RT = smi_helpers.ROITools();
RT.Pixel2nm = SMD_overview.PixelSize * 1000;
RT.ROI_sizes = [5000, 5000];  % 5 micron squares
[n_ROIs, RoI, XYsize] = RT.getROI(SMD_overview, 'Select Dense Regions');

% Save ROI definitions
save('Analysis_ROIs.mat', 'n_ROIs', 'RoI', 'XYsize', ...
    'Pixel2nm', 'ResultsFile', '-v7.3');

% Step 3: Process each ROI independently with strict thresholds
for i = 1:n_ROIs
    fprintf('Processing ROI %d/%d...\n', i, n_ROIs);

    % Configure for this ROI
    SMF_ROI = SMF;  % Copy parameters
    SMF_ROI.Thresholding.MinPhotons = 300;  % Stricter threshold

    % Define spatial restriction (implementation depends on data loading)
    % For now, process full dataset then extract
    SMLMobj_ROI = smi.SMLM(SMF_ROI);
    SMLMobj_ROI.genLocalizations();

    % Extract ROI
    ROI_px = RoI{i}.ROI / (SMD_overview.PixelSize * 1000);
    ROI_px = ROI_px([3, 1, 4, 2]);  % Convert to [YStart, XStart, YEnd, XEnd]
    SMD_ROI = smi_core.SingleMoleculeData.isolateSubROI(...
        SMLMobj_ROI.SMD, ROI_px);

    % Save ROI result
    save(sprintf('SMD_ROI_%02d.mat', i), 'SMD_ROI', '-v7.3');
end
```

### Comparing Statistics Across ROIs

Quantify spatial heterogeneity by comparing ROI properties:

```matlab
% Get ROIs
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Compartment Analysis');

% Initialize results
roi_stats = struct();

for i = 1:n_ROIs
    SMD_ROI = RoI{i}.SMD{1};

    % Compute statistics
    roi_stats(i).n_locs = length(SMD_ROI.X);
    roi_stats(i).mean_photons = mean(SMD_ROI.Photons);
    roi_stats(i).std_photons = std(SMD_ROI.Photons);
    roi_stats(i).mean_bg = mean(SMD_ROI.Bg);
    roi_stats(i).density = roi_stats(i).n_locs / ...
        (diff(RoI{i}.ROI(1:2)) * diff(RoI{i}.ROI(3:4)) / 1e6);  % per um^2

    % Localization precision
    roi_stats(i).mean_precision_x = mean(SMD_ROI.X_SE);
    roi_stats(i).mean_precision_y = mean(SMD_ROI.Y_SE);
end

% Display comparison
fprintf('\nROI Statistics Summary:\n');
fprintf('%-10s %10s %12s %10s %12s\n', ...
    'ROI', 'N_locs', 'Photons', 'Bg', 'Density');
for i = 1:n_ROIs
    fprintf('%-10d %10d %12.1f %10.1f %12.2f\n', ...
        i, roi_stats(i).n_locs, roi_stats(i).mean_photons, ...
        roi_stats(i).mean_bg, roi_stats(i).density);
end

% Statistical comparison
photons_per_roi = [roi_stats.mean_photons];
[h, p] = ttest2(photons_per_roi(1:2), photons_per_roi(3:4));
fprintf('\nPhoton comparison p-value: %.4f\n', p);
```

### Spatial Clustering Within ROIs

Perform clustering analysis restricted to ROI:

```matlab
% Select ROI
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Cluster Analysis Region');
SMD_ROI = RoI{1}.SMD{1};

% DBSCAN clustering
Eps = 0.020;  % 20 nm search radius in pixels
MinPts = 5;   % Minimum points per cluster

% Run clustering
ClustSMD = smi_cluster.Clustering();
ClustSMD.E = Eps;
ClustSMD.minPts = MinPts;
ClustSMD.ClusterMethod = 'DBSCAN_Daszykowski';

[SMD_clustered, ClusterSMD] = ClustSMD.cluster(SMD_ROI);

% Analyze clusters
n_clusters = max(SMD_clustered.ClusterID);
fprintf('Found %d clusters in ROI\n', n_clusters);

% Plot clustered data
figure;
gscatter(SMD_clustered.X, SMD_clustered.Y, SMD_clustered.ClusterID);
axis equal;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('Clusters in ROI (n=%d)', n_clusters));
```

### Drift Correction Per ROI

For datasets with spatially-varying drift, correct each ROI independently:

```matlab
% Select multiple ROIs
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Drift Correction Regions');

% Configure drift correction
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.L_intra = 1.0;

for i = 1:n_ROIs
    SMD_ROI = RoI{i}.SMD{1};

    % Apply drift correction to this ROI
    DC = smi_core.DriftCorrection(SMF, SMD_ROI);
    [SMD_ROI_corrected, Stats] = DC.driftCorrectKNN(SMD_ROI);

    % Save corrected ROI
    RoI{i}.SMD_corrected{1} = SMD_ROI_corrected;

    fprintf('ROI %d drift: X = %.2f nm, Y = %.2f nm\n', i, ...
        range(SMD_ROI_corrected.DriftX(:,1)) * SMD.PixelSize * 1000, ...
        range(SMD_ROI_corrected.DriftY(:,1)) * SMD.PixelSize * 1000);
end
```

## Combining ROI Results

### Concatenating ROI-Processed Data

After separately processing ROIs, combine into single dataset:

```matlab
% Initialize combined SMD
SMD_combined = smi_core.SingleMoleculeData.createSMD();

% Load and concatenate individual ROIs
n_ROIs = 5;
for i = 1:n_ROIs
    % Load ROI result
    data = load(sprintf('SMD_ROI_%02d.mat', i));
    SMD_ROI = data.SMD_ROI;

    % Add ROI identifier (optional)
    SMD_ROI.ROI_ID = i * ones(size(SMD_ROI.X));

    % Concatenate
    SMD_combined = smi_core.SingleMoleculeData.catSMD(...
        SMD_combined, SMD_ROI, false);
end

fprintf('Combined %d ROIs: %d total localizations\n', ...
    n_ROIs, length(SMD_combined.X));

% Verify no overlap
if isfield(SMD_combined, 'ROI_ID')
    for i = 1:n_ROIs
        n_in_roi = sum(SMD_combined.ROI_ID == i);
        fprintf('  ROI %d: %d localizations\n', i, n_in_roi);
    end
end
```

### Combining BaGoL ROI Results

For BaGoL (Bayesian Grouping of Localizations) analysis of ROIs:

```matlab
% This function combines separately-processed BaGoL ROIs
% See MATLAB/+smi_cluster/@ClusterInterface/combineBaGoLROIs.m

% Paths to results
pathnameR = '/path/to/SR_Results';   % Original SR ROI definitions
filesR = {'Cell_01_Label_01_Results_ROIs.mat'};

pathnameB = '/path/to/BaGoL_Results';  % BaGoL processed ROIs
filesB = {
    'MAPN_Cell_01_Label_01_Results_ROI_01.mat',
    'MAPN_Cell_01_Label_01_Results_ROI_02.mat',
    'MAPN_Cell_01_Label_02_Results_ROI_01.mat',
    'MAPN_Cell_01_Label_02_Results_ROI_02.mat'
};

% Combine ROIs
smi_cluster.ClusterInterface.combineBaGoLROIs(...
    pathnameR, filesR, pathnameB, filesB, true, false, false);

% Output: pathnameB/Analysis/Cell_01_Results_BaGoL_ROIs.mat
```

**Parameters**:
- `pathnameR`: Path to original SR ROI files
- `filesR`: Cell array of `_Results_ROIs.mat` files
- `pathnameB`: Path to BaGoL ROI files
- `filesB`: Cell array of BaGoL MAPN files
- `MAPNfile`: true if using MAPN files, false for BaGoL_Results files
- `keep_numbering`: Retain ROI numbers even if some ROIs missing
- `GaussIm`: true if original ROIs selected from Gaussian image

## Visualization Techniques

### Displaying ROI Boundaries

```matlab
% Generate super-resolution image
SRImageZoom = 20;
img = smi_vis.GenerateImages.gaussImage(...
    [SMD.X, SMD.Y], SMD.PixelSize, SRImageZoom, [SMD.YSize, SMD.XSize]);

% Display with ROI overlays
figure;
imagesc(img); axis image; colormap(gray);
hold on;

% Overlay ROI rectangles
for i = 1:n_ROIs
    roi_nm = RoI{i}.ROI;  % [xmin, xmax, ymin, ymax] in nm
    roi_px = roi_nm / (SMD.PixelSize * 1000) * SRImageZoom;

    rectangle('Position', [roi_px(1), roi_px(3), ...
        roi_px(2)-roi_px(1), roi_px(4)-roi_px(3)], ...
        'EdgeColor', 'r', 'LineWidth', 2);

    % Label ROI
    text(mean(roi_px([1,2])), mean(roi_px([3,4])), sprintf('%d', i), ...
        'Color', 'yellow', 'FontSize', 14, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end
hold off;
title('ROI Locations on SR Image');
xlabel('X (SR pixels)'); ylabel('Y (SR pixels)');
```

### Side-by-Side ROI Comparison

```matlab
% Generate images for each ROI
n_cols = 3;
n_rows = ceil(n_ROIs / n_cols);

figure('Position', [100, 100, 1200, 400*n_rows]);
for i = 1:n_ROIs
    subplot(n_rows, n_cols, i);

    SMD_ROI = RoI{i}.SMD{1};

    % Generate ROI image
    % Adjust coordinates to ROI-local frame
    X_local = SMD_ROI.X - min(SMD_ROI.X);
    Y_local = SMD_ROI.Y - min(SMD_ROI.Y);

    img_roi = smi_vis.GenerateImages.gaussImage(...
        [X_local, Y_local], SMD.PixelSize, 20, ...
        [max(Y_local), max(X_local)]);

    imagesc(img_roi); axis image; colormap(gray);
    title(sprintf('ROI %d (N=%d)', i, length(SMD_ROI.X)));
    xlabel('X (SR pixels)'); ylabel('Y (SR pixels)');
end
```

### Density Heatmaps Per ROI

```matlab
% Compute and display localization density for each ROI
figure('Position', [100, 100, 1200, 400]);

for i = 1:min(3, n_ROIs)  % Show first 3 ROIs
    subplot(1, 3, i);

    SMD_ROI = RoI{i}.SMD{1};

    % Compute density image
    DensityIm = smi_core.SingleMoleculeData.computeDensityImage(...
        SMD_ROI, 0.5, 10);  % SigmaSmooth=0.5 pixels, Mag=10x

    imagesc(DensityIm); axis image;
    colormap(hot); colorbar;
    title(sprintf('ROI %d Density', i));
    xlabel('X'); ylabel('Y');
end
```

### Multi-Label ROI Overlays

```matlab
% Display two-color ROIs with different colors
figure;

% Plot both labels
colors = ['g', 'r'];
for j = 1:2  % Two labels
    for i = 1:n_ROIs
        SMD_ROI = RoI{i}.SMD{j};
        plot(SMD_ROI.X, SMD_ROI.Y, [colors(j), '.'], 'MarkerSize', 5);
        hold on;
    end
end

% Overlay ROI boundaries
for i = 1:n_ROIs
    roi_px = RoI{i}.ROI / (SMD.PixelSize * 1000);
    roi_px = roi_px([1, 2, 2, 1, 1; 3, 3, 4, 4, 3]);
    plot(roi_px(1,:), roi_px(2,:), 'w-', 'LineWidth', 2);
end

axis equal;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Dual-Color ROIs');
legend('Label 1', 'Label 2');
hold off;
```

## Troubleshooting

### Issue: FindROI detects too many false positives

**Diagnose**:
```matlab
% Check detected boxes on sample frame
FR.Verbose = 3;
FR.PlotBoxFrame = 10;  % Frame to visualize
[ROIStack, SMD] = FR.findROI();
```

**Solutions**:
- Increase `SMF.BoxFinding.MinPhotons` (e.g., 200 → 400)
- Adjust `SMF.Fitting.PSFSigma` to match actual PSF
- Check camera calibration (gain, offset, readnoise)
- Verify data is properly converted to photons

### Issue: FindROI misses dim emitters

**Diagnose**:
```matlab
% Check photon distribution
figure;
histogram(SMD.Photons, 50);
xlabel('Photons'); ylabel('Count');
title('Photon Distribution');
```

**Solutions**:
- Decrease `SMF.BoxFinding.MinPhotons` threshold
- Increase `SMF.BoxFinding.BoxOverlap` to allow closer detections
- Verify PSFSigma is not too large (inflates threshold)

### Issue: ROI selection displays empty or incorrect image

**Diagnose**:
```matlab
% Check SMD coordinate range
fprintf('X range: %.2f to %.2f pixels\n', min(SMD.X), max(SMD.X));
fprintf('Y range: %.2f to %.2f pixels\n', min(SMD.Y), max(SMD.Y));
fprintf('N localizations: %d\n', length(SMD.X));

% Verify pixel size
fprintf('PixelSize: %.4f um\n', SMD.PixelSize);
```

**Solutions**:
- Ensure SMD has valid coordinates (not all zeros/NaN)
- Check `RT.Pixel2nm` matches data (typically 100 nm)
- For Gaussian images, verify `RT.SRzoom` is reasonable (10-30)
- Set `RT.OriginLLvsUL` appropriately for data type

### Issue: isolateSubROI returns empty SMD

**Diagnose**:
```matlab
% Check ROI bounds relative to data
ROI_px = [50, 100, 150, 200];
fprintf('ROI: Y=[%d,%d], X=[%d,%d]\n', ROI_px([1,3]), ROI_px([2,4]));
fprintf('Data: Y=[%.1f,%.1f], X=[%.1f,%.1f]\n', ...
    min(SMD.Y), max(SMD.Y), min(SMD.X), max(SMD.X));

% Count localizations in ROI
in_roi = (SMD.Y >= ROI_px(1)) & (SMD.Y <= ROI_px(3)) & ...
         (SMD.X >= ROI_px(2)) & (SMD.X <= ROI_px(4));
fprintf('Localizations in ROI: %d\n', sum(in_roi));
```

**Solutions**:
- Verify ROI coordinates overlap with data range
- Check coordinate convention: ROI is [YStart, XStart, YEnd, XEnd]
- Ensure units match (pixels, not nm)
- ROI bounds are inclusive

### Issue: Multi-label ROI selection rejects all ROIs

**Diagnose**:
```matlab
% Check if labels overlap spatially
figure;
subplot(1,2,1);
plot(SMD1.X, SMD1.Y, 'g.', 'MarkerSize', 5);
title('Label 1'); axis equal;

subplot(1,2,2);
plot(SMD2.X, SMD2.Y, 'r.', 'MarkerSize', 5);
title('Label 2'); axis equal;
```

**Solutions**:
- Verify both labels are properly aligned (channel registration)
- Apply coordinate transform: `RT.Transform{2} = RegistrationMatrix;`
- Check that labels actually colocalize in selected regions
- ROI must contain at least one localization from EACH label

### Issue: Memory errors when processing large ROIs

**Solutions**:
- Split large ROIs into smaller sub-regions
- Process ROIs sequentially, clearing variables between iterations:
  ```matlab
  for i = 1:n_ROIs
      SMD_ROI = RoI{i}.SMD{1};
      % ... process ...
      save(sprintf('ROI_%02d.mat', i), 'SMD_ROI');
      clear SMD_ROI;  % Free memory
  end
  ```
- Reduce super-resolution zoom factor
- Use sparse representations when possible

## Best Practices

### ROI Design

1. **Size appropriately**: Balance detail vs. statistics
   - Too small: Poor statistics, boundary effects
   - Too large: Loss of spatial specificity
   - Typical: 1-5 μm for cellular structures

2. **Avoid boundaries**: Keep ROIs away from image edges
   - Edge artifacts in localization
   - Incomplete PSFs near borders

3. **Consider symmetry**: Use matched ROI sizes for comparisons
   - Same size enables fair statistical comparison
   - Fixed-size ROIs with left-click in ROITools

4. **Plan for overlap**: Decide if ROIs should overlap or tile
   - Non-overlapping: Independent measurements
   - Overlapping: Smooth spatial analysis

### Processing Strategy

1. **Overview first**: Generate low-quality overview before ROI selection
   - Fast initial localization with permissive thresholds
   - Select ROIs from overview
   - Re-process with strict thresholds

2. **Save ROI definitions**: Keep ROI boundaries for reproducibility
   ```matlab
   save('Analysis_ROIs.mat', 'n_ROIs', 'RoI', 'XYsize', ...
       'Pixel2nm', 'SMF', '-v7.3');
   ```

3. **Document parameters**: Record all processing parameters per ROI
   - Save SMF alongside results
   - Include date, software version
   - Note any manual interventions

4. **Validate extraction**: Verify ROI data makes sense
   - Check localization counts
   - Visualize extracted regions
   - Compare statistics to full dataset

### Analysis Workflow

1. **Process independently**: Treat each ROI as separate experiment
   - Independent thresholding
   - Independent drift correction
   - Independent clustering

2. **Control for artifacts**: Check for systematic biases
   - Edge effects
   - Position-dependent photon counts
   - Drift variations across field

3. **Statistical rigor**: Use appropriate tests for ROI comparisons
   - Account for multiple comparisons (Bonferroni, FDR)
   - Check assumptions (normality, equal variance)
   - Report effect sizes, not just p-values

## Performance Considerations

### FindROI Performance

**GPU acceleration**: FindROI uses CUDA for filtering
- Speed: ~10-100 ms per frame (256×256)
- Requires: NVIDIA GPU with compute capability ≥5.0
- Memory: Automatically chunks data to fit GPU memory

**Optimization tips**:
- Use `IsSCMOS = 0` if not using sCMOS camera (faster)
- Reduce `BoxSize` if possible (fewer pixels per ROI)
- Increase `MinPhotons` to reduce candidate count

### ROITools Performance

**Display speed**:
- Point cloud: Fast (direct plotting)
- Gaussian image: Slower (rendering step)
  - Rendering time scales with localization count
  - Typical: 1-5 seconds for 100,000 localizations

**Memory usage**:
- Proportional to number of localizations
- Typical: ~1 KB per localization
- 1 million localizations ≈ 1 GB memory

## See Also

- [SMLM Analysis Workflow](smlm-analysis.md) - Complete pipeline overview
- [SMD Structure](../core-concepts/smd-structure.md) - Understanding localization data
- [Visualization Guide](../how-to/visualize-results.md) - Plotting and image generation
- MATLAB/+smi_core/@FindROI/README.md - FindROI technical details
- MATLAB/+smi_helpers/@ROITools/README.md - ROITools technical details
- MATLAB/+smi_core/@SingleMoleculeData/isolateSubROI.m - ROI extraction implementation
