---
title: "Working with ROIs"
category: "how-to"
level: "beginner"
tags: ["roi", "region-of-interest", "roi-tools", "spatial-selection", "data-extraction"]
prerequisites: ["../core-concepts/smd-structure.md", "../getting-started/first-analysis.md"]
related: ["../workflows/roi-analysis.md", "visualize-results.md"]
summary: "Practical guide to creating, modifying, saving, and loading Regions of Interest (ROIs) for spatial data selection"
estimated_time: "12 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Working with ROIs

## Purpose

Regions of Interest (ROIs) let you select specific spatial areas from your localization data for focused analysis. This guide covers the practical operations of creating ROIs, understanding the ROI data structure, converting between formats, extracting localization data, and visualizing ROI boundaries. For complete workflow-level ROI analysis including clustering and statistics, see the [ROI Analysis workflow](../workflows/roi-analysis.md).

## Prerequisites

- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Completed at least one [SMLM analysis](../getting-started/first-analysis.md)
- Basic MATLAB plotting knowledge

## Overview

smite uses the `smi_helpers.ROITools` class for interactive ROI operations. The basic workflow is:

1. Create an ROITools object and configure display options
2. Interactively select ROIs from localizations or images
3. Extract coordinates and data from each ROI
4. Save ROI definitions for later use
5. Visualize ROI boundaries on images

ROIs are stored in a cell array structure (`RoI`) that contains coordinates, uncertainties, and SMD subsets for each selected region.

## The ROI Data Structure

### RoI Cell Array Format

When you create ROIs using `ROITools.getROI()`, the output is a cell array where each element represents one ROI:

```matlab
% Create ROIs
RT = smi_helpers.ROITools();
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'My ROIs');

% RoI structure for single label:
% RoI{i}.ROI         [xmin, xmax, ymin, ymax] in nm
% RoI{i}.X{1}        X coordinates in ROI (nm) [Nx1 vector]
% RoI{i}.Y{1}        Y coordinates in ROI (nm) [Nx1 vector]
% RoI{i}.X_SE{1}     X localization errors (nm) [Nx1 vector]
% RoI{i}.Y_SE{1}     Y localization errors (nm) [Nx1 vector]
% RoI{i}.SMD{1}      Full SMD structure for this ROI
```

### Understanding the Structure

The `RoI` cell array uses nested cells to support multi-label data:

```matlab
% Access first ROI
roi_index = 1;

% Boundary box (nm)
bounds = RoI{roi_index}.ROI;  % [xmin, xmax, ymin, ymax]
fprintf('ROI bounds: X=[%.1f, %.1f] Y=[%.1f, %.1f] nm\n', ...
    bounds(1), bounds(2), bounds(3), bounds(4));

% Coordinates for label 1
X_coords = RoI{roi_index}.X{1};      % X positions (nm)
Y_coords = RoI{roi_index}.Y{1};      % Y positions (nm)
X_errors = RoI{roi_index}.X_SE{1};   % X uncertainties (nm)
Y_errors = RoI{roi_index}.Y_SE{1};   % Y uncertainties (nm)

fprintf('ROI contains %d localizations\n', length(X_coords));
```

### Multi-Label ROI Structure

For multi-label data (e.g., dual-color imaging), the coordinate cells contain data for each label:

```matlab
% Create dual-label ROIs
[n_ROIs, RoI, XYsize] = RT.getROI({SMD1, SMD2}, 'Dual-Color ROIs');

% Access both labels in first ROI
roi_index = 1;

% Label 1 data
X_label1 = RoI{roi_index}.X{1};
Y_label1 = RoI{roi_index}.Y{1};
SMD_label1 = RoI{roi_index}.SMD{1};

% Label 2 data
X_label2 = RoI{roi_index}.X{2};
Y_label2 = RoI{roi_index}.Y{2};
SMD_label2 = RoI{roi_index}.SMD{2};

fprintf('Label 1: %d locs, Label 2: %d locs\n', ...
    length(X_label1), length(X_label2));
```

## Creating ROIs Interactively

### Basic ROI Selection

The simplest way to create ROIs from localization data:

```matlab
% Load your localization results
SMD = load('Results/SMD.mat');
SMD = SMD.SMD;

% Create ROITools object
RT = smi_helpers.ROITools();

% Set ROI size (nm)
RT.ROI_sizes = [2000, 2000];  % 2 x 2 micron squares

% Set pixel conversion (from SMD or manually)
RT.Pixel2nm = SMD.PixelSize * 1000;  % Convert micrometers to nm

% Interactive selection
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Select Analysis Regions');

fprintf('Created %d ROIs\n', n_ROIs);
```

### Interactive Controls

When the ROI selection window opens:

**Mouse Actions**:
- **Left click**: Place fixed-size rectangular ROI centered at cursor
  - Size determined by `RT.ROI_sizes`
  - Quick and consistent for comparing regions
- **Right click and drag**: Draw adjustable rectangular ROI
  - Click at one corner, drag to opposite corner
  - Use when ROIs need different sizes
- **Middle click**: Ignored (no action)

**Keyboard Actions**:
- **Backspace** or **Delete**: Remove most recently added ROI
  - Outline disappears immediately
  - Can delete multiple times to remove several ROIs
- **Any other key** (Enter, Space, etc.): Finish selection
  - Closes window and returns ROI data
  - ROIs are not automatically saved

### Configuring Display Options

Choose how localizations are displayed during selection:

```matlab
RT = smi_helpers.ROITools();

% Point cloud display (default, fastest)
RT.GaussIm = false;
RT.OriginLLvsUL = true;   % Lower-left origin
RT.Msize = 7;              % Marker size for points
RT.Color = ['g', 'm', 'c'];  % Colors for multiple labels

% Gaussian image display (single label only)
RT.GaussIm = true;
RT.OriginLLvsUL = false;   % Upper-left origin (image convention)
RT.SRzoom = 20;            % Super-resolution zoom factor

[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Gaussian Image ROIs');
```

**Display recommendations**:
- Use point cloud for quick selection and multi-label data
- Use Gaussian image for better visualization of structures
- Gaussian rendering takes longer for large datasets
- `GaussIm = true` currently only supports single labels

### Creating ROIs from Different Sources

ROITools can load data from multiple sources:

```matlab
RT = smi_helpers.ROITools();

% From SMD structure
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'From SMD');

% From results file
[n_ROIs, RoI, XYsize] = RT.getROI('Results/SMD.mat', 'From File');

% From coordinate arrays (N x 2 matrix in nm)
coords = [X_nm, Y_nm];  % X, Y coordinates
[n_ROIs, RoI, XYsize] = RT.getROI(coords, 'From Coordinates');

% From coordinates with errors (N x 4 matrix in nm)
coords_with_errors = [X_nm, Y_nm, X_SE_nm, Y_SE_nm];
[n_ROIs, RoI, XYsize] = RT.getROI(coords_with_errors, 'With Errors');

% From BaGoL results
[n_ROIs, RoI, XYsize] = RT.getROI('MAPN_data.mat', 'BaGoL ROIs');
```

## Extracting Data from ROIs

### Accessing Coordinates

Extract position and uncertainty data from ROI structures:

```matlab
% Select a specific ROI
roi_num = 3;

% Get coordinates (nm)
X_nm = RoI{roi_num}.X{1};
Y_nm = RoI{roi_num}.Y{1};

% Get localization precision (nm)
X_SE_nm = RoI{roi_num}.X_SE{1};
Y_SE_nm = RoI{roi_num}.Y_SE{1};

% ROI boundaries
roi_box = RoI{roi_num}.ROI;  % [xmin, xmax, ymin, ymax] (nm)

% Compute ROI properties
n_locs = length(X_nm);
roi_width = roi_box(2) - roi_box(1);   % nm
roi_height = roi_box(4) - roi_box(3);  % nm
roi_area = roi_width * roi_height / 1e6;  % square microns

fprintf('ROI %d: %d localizations in %.2f um^2\n', ...
    roi_num, n_locs, roi_area);
```

### Extracting Full SMD Data

The ROI structure includes complete SMD with all fields:

```matlab
% Get SMD for specific ROI
SMD_roi = RoI{roi_num}.SMD{1};

% All SMD fields available
fprintf('Photons: %.1f +/- %.1f\n', ...
    mean(SMD_roi.Photons), std(SMD_roi.Photons));
fprintf('Background: %.1f +/- %.1f\n', ...
    mean(SMD_roi.Bg), std(SMD_roi.Bg));
fprintf('Frame range: %d to %d\n', ...
    min(SMD_roi.FrameNum), max(SMD_roi.FrameNum));

% Localization precision
fprintf('Precision X: %.2f nm\n', mean(SMD_roi.X_SE) * SMD.PixelSize * 1000);
fprintf('Precision Y: %.2f nm\n', mean(SMD_roi.Y_SE) * SMD.PixelSize * 1000);

% Can use SMD_roi for further analysis
% (clustering, drift correction, frame connection, etc.)
```

### Processing All ROIs

Loop through ROIs to extract or analyze all regions:

```matlab
% Initialize storage
roi_summary = struct();

for i = 1:n_ROIs
    % Extract SMD
    SMD_roi = RoI{i}.SMD{1};

    % Compute statistics
    roi_summary(i).n_locs = length(SMD_roi.X);
    roi_summary(i).mean_photons = mean(SMD_roi.Photons);
    roi_summary(i).mean_bg = mean(SMD_roi.Bg);
    roi_summary(i).bounds = RoI{i}.ROI;

    % Calculate density
    area = diff(RoI{i}.ROI(1:2)) * diff(RoI{i}.ROI(3:4)) / 1e6;  % um^2
    roi_summary(i).density = roi_summary(i).n_locs / area;
end

% Display summary table
fprintf('%-8s %10s %12s %10s %12s\n', ...
    'ROI', 'N_locs', 'Photons', 'Bg', 'Density');
for i = 1:n_ROIs
    fprintf('%-8d %10d %12.1f %10.1f %12.2f\n', ...
        i, roi_summary(i).n_locs, roi_summary(i).mean_photons, ...
        roi_summary(i).mean_bg, roi_summary(i).density);
end
```

## Converting ROI Formats

### Converting ROI Bounds to Pixel Coordinates

The `RoI` structure stores boundaries in nm. Convert to pixels for use with `isolateSubROI`:

```matlab
% ROI bounds from RoI structure (nm)
roi_nm = RoI{1}.ROI;  % [xmin, xmax, ymin, ymax]

% Convert to pixels
pixel_size_nm = SMD.PixelSize * 1000;  % micrometers to nm
roi_px = roi_nm / pixel_size_nm;

% Reorder for isolateSubROI: [YStart, XStart, YEnd, XEnd]
roi_px_ordered = roi_px([3, 1, 4, 2]);

% Extract using isolateSubROI
SMD_extracted = smi_core.SingleMoleculeData.isolateSubROI(SMD, roi_px_ordered);
```

### Creating ROI Structure Programmatically

Build RoI cell array manually without interactive selection:

```matlab
% Define ROI boundaries (nm)
roi_bounds = [
    1000, 3000, 1000, 3000;   % ROI 1
    4000, 6000, 2000, 4000;   % ROI 2
    5000, 7000, 5000, 7000    % ROI 3
];

n_ROIs = size(roi_bounds, 1);
RoI = cell(1, n_ROIs);

for i = 1:n_ROIs
    % Store boundary
    RoI{i}.ROI = roi_bounds(i, :);  % [xmin, xmax, ymin, ymax]

    % Find localizations in ROI
    in_roi = (SMD.X * SMD.PixelSize * 1000 >= roi_bounds(i, 1)) & ...
             (SMD.X * SMD.PixelSize * 1000 <= roi_bounds(i, 2)) & ...
             (SMD.Y * SMD.PixelSize * 1000 >= roi_bounds(i, 3)) & ...
             (SMD.Y * SMD.PixelSize * 1000 <= roi_bounds(i, 4));

    % Extract coordinates (nm)
    RoI{i}.X{1} = SMD.X(in_roi) * SMD.PixelSize * 1000;
    RoI{i}.Y{1} = SMD.Y(in_roi) * SMD.PixelSize * 1000;
    RoI{i}.X_SE{1} = SMD.X_SE(in_roi) * SMD.PixelSize * 1000;
    RoI{i}.Y_SE{1} = SMD.Y_SE(in_roi) * SMD.PixelSize * 1000;

    % Extract full SMD
    roi_px = roi_bounds(i, :) / (SMD.PixelSize * 1000);
    roi_px = roi_px([3, 1, 4, 2]);
    RoI{i}.SMD{1} = smi_core.SingleMoleculeData.isolateSubROI(SMD, roi_px);
end

fprintf('Created %d programmatic ROIs\n', n_ROIs);
```

### Converting Between Coordinate Systems

Handle coordinate origin conversions for different data types:

```matlab
% Lower-left origin (typical for SMD) to upper-left origin (images)
X_ll = SMD.X * SMD.PixelSize * 1000;  % nm, lower-left
Y_ll = SMD.Y * SMD.PixelSize * 1000;  % nm, lower-left

image_height_nm = SMD.YSize * SMD.PixelSize * 1000;
Y_ul = image_height_nm - Y_ll;  % Convert to upper-left

% Upper-left back to lower-left
Y_ll_back = image_height_nm - Y_ul;

% Apply to ROI bounds
roi_ll = [xmin, xmax, ymin, ymax];  % Lower-left origin
roi_ul = [xmin, xmax, image_height_nm - ymax, image_height_nm - ymin];
```

## Saving and Loading ROIs

### Saving ROI Definitions

Save ROI structures for reproducibility and later use:

```matlab
% After creating ROIs
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Analysis ROIs');

% Save complete ROI data
Pixel2nm = SMD.PixelSize * 1000;
OriginLLvsUL = RT.OriginLLvsUL;
ResultsFile = fullfile(SMF.Data.FileDir, 'Results', 'SMD.mat');

save('MyAnalysis_ROIs.mat', 'n_ROIs', 'RoI', 'XYsize', ...
    'Pixel2nm', 'OriginLLvsUL', 'ResultsFile', '-v7.3');

fprintf('Saved %d ROIs to MyAnalysis_ROIs.mat\n', n_ROIs);
```

### Loading Saved ROIs

Reload ROI definitions from file:

```matlab
% Load ROI data
roi_data = load('MyAnalysis_ROIs.mat');

n_ROIs = roi_data.n_ROIs;
RoI = roi_data.RoI;
XYsize = roi_data.XYsize;
Pixel2nm = roi_data.Pixel2nm;

fprintf('Loaded %d ROIs\n', n_ROIs);

% Display ROI information
for i = 1:n_ROIs
    bounds = RoI{i}.ROI;
    n_locs = length(RoI{i}.X{1});
    fprintf('ROI %d: [%.0f,%.0f,%.0f,%.0f] nm, %d locs\n', ...
        i, bounds(1), bounds(2), bounds(3), bounds(4), n_locs);
end
```

### Saving ROI Bounds Only

Save just the spatial boundaries for lightweight storage:

```matlab
% Extract boundaries
roi_bounds = zeros(n_ROIs, 4);
for i = 1:n_ROIs
    roi_bounds(i, :) = RoI{i}.ROI;  % [xmin, xmax, ymin, ymax]
end

% Save to text file
fid = fopen('ROI_bounds.txt', 'w');
fprintf(fid, '# ROI_num xmin xmax ymin ymax (nm)\n');
for i = 1:n_ROIs
    fprintf(fid, '%d %.3f %.3f %.3f %.3f\n', i, roi_bounds(i, :));
end
fclose(fid);

% Save to MAT file
save('ROI_bounds.mat', 'roi_bounds', 'Pixel2nm');
```

### Loading and Applying ROI Bounds

Reload boundaries and apply to new data:

```matlab
% Load saved bounds
bounds_data = load('ROI_bounds.mat');
roi_bounds = bounds_data.roi_bounds;  % [n_ROIs x 4]

% Load new SMD data
SMD_new = load('NewData_Results/SMD.mat');
SMD_new = SMD_new.SMD;

% Apply each ROI boundary
n_ROIs = size(roi_bounds, 1);
for i = 1:n_ROIs
    % Convert to pixels and extract
    roi_px = roi_bounds(i, :) / (SMD_new.PixelSize * 1000);
    roi_px = roi_px([3, 1, 4, 2]);  % Reorder

    SMD_roi = smi_core.SingleMoleculeData.isolateSubROI(SMD_new, roi_px);

    % Save or process
    save(sprintf('NewData_ROI_%02d.mat', i), 'SMD_roi');
end
```

## Visualizing ROI Boundaries

### Overlaying ROIs on Images

Display ROI boundaries on super-resolution images:

```matlab
% Generate super-resolution image
SRzoom = 20;
img = smi_vis.GenerateImages.gaussImage(...
    [SMD.X, SMD.Y], SMD.PixelSize, SRzoom, [SMD.YSize, SMD.XSize]);

% Display image
figure('Position', [100, 100, 800, 800]);
imagesc(img);
axis image;
colormap(gray);
hold on;

% Overlay ROI rectangles
for i = 1:n_ROIs
    % Get bounds in nm
    roi_nm = RoI{i}.ROI;  % [xmin, xmax, ymin, ymax]

    % Convert to SR image pixels
    roi_sr = roi_nm / (SMD.PixelSize * 1000) * SRzoom;

    % Draw rectangle
    width = roi_sr(2) - roi_sr(1);
    height = roi_sr(4) - roi_sr(3);
    rectangle('Position', [roi_sr(1), roi_sr(3), width, height], ...
        'EdgeColor', 'red', 'LineWidth', 2);

    % Add label
    text(mean(roi_sr([1, 2])), mean(roi_sr([3, 4])), sprintf('%d', i), ...
        'Color', 'yellow', 'FontSize', 12, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center');
end

hold off;
title('ROI Locations');
xlabel('X (SR pixels)');
ylabel('Y (SR pixels)');
```

### Plotting ROI Boundaries on Point Cloud

Overlay ROI boxes on scatter plot of localizations:

```matlab
figure('Position', [100, 100, 800, 800]);

% Plot all localizations
plot(SMD.X, SMD.Y, 'k.', 'MarkerSize', 2);
hold on;
axis equal;

% Overlay ROIs
for i = 1:n_ROIs
    % Get bounds in nm and convert to pixels
    roi_nm = RoI{i}.ROI;
    roi_px = roi_nm / (SMD.PixelSize * 1000);

    % Define rectangle corners
    x_rect = [roi_px(1), roi_px(2), roi_px(2), roi_px(1), roi_px(1)];
    y_rect = [roi_px(3), roi_px(3), roi_px(4), roi_px(4), roi_px(3)];

    % Draw boundary
    plot(x_rect, y_rect, 'r-', 'LineWidth', 2);

    % Label ROI
    text(mean(roi_px([1, 2])), mean(roi_px([3, 4])), sprintf('%d', i), ...
        'Color', 'red', 'FontSize', 14, 'FontWeight', 'bold', ...
        'HorizontalAlignment', 'center', 'BackgroundColor', 'white');
end

hold off;
xlabel('X (pixels)');
ylabel('Y (pixels)');
title(sprintf('Localizations with %d ROIs', n_ROIs));
```

### Highlighting ROI Contents

Show which localizations belong to each ROI:

```matlab
figure('Position', [100, 100, 1200, 400]);

% Plot 1: All data with ROI boundaries
subplot(1, 3, 1);
plot(SMD.X, SMD.Y, 'k.', 'MarkerSize', 2);
hold on;
for i = 1:min(3, n_ROIs)
    roi_px = RoI{i}.ROI / (SMD.PixelSize * 1000);
    rectangle('Position', [roi_px(1), roi_px(3), ...
        roi_px(2)-roi_px(1), roi_px(4)-roi_px(3)], ...
        'EdgeColor', 'r', 'LineWidth', 2);
end
hold off;
axis equal;
title('All Data');

% Plot 2 & 3: Individual ROI contents
for i = 1:min(2, n_ROIs)
    subplot(1, 3, i+1);

    % Get ROI data
    X_roi = RoI{i}.X{1} / (SMD.PixelSize * 1000);  % Convert to pixels
    Y_roi = RoI{i}.Y{1} / (SMD.PixelSize * 1000);

    % Plot
    plot(X_roi, Y_roi, 'b.', 'MarkerSize', 5);
    axis equal;
    title(sprintf('ROI %d (%d locs)', i, length(X_roi)));
    xlabel('X (pixels)');
    ylabel('Y (pixels)');
end
```

## Practical Examples

### Example 1: Select and Analyze Dense Regions

```matlab
% Load data
SMD = load('Results/SMD.mat');
SMD = SMD.SMD;

% Configure ROITools
RT = smi_helpers.ROITools();
RT.ROI_sizes = [3000, 3000];  % 3 micron squares
RT.Pixel2nm = SMD.PixelSize * 1000;

% Select regions of interest
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Select Dense Regions');

% Analyze each ROI
for i = 1:n_ROIs
    SMD_roi = RoI{i}.SMD{1};

    % Compute local density
    area_um2 = prod(diff(reshape(RoI{i}.ROI, 2, 2), 1, 1)) / 1e6;
    density = length(SMD_roi.X) / area_um2;

    fprintf('ROI %d: %.1f locs/um^2\n', i, density);
end

% Save results
save('DenseRegions_ROIs.mat', 'n_ROIs', 'RoI', 'XYsize', '-v7.3');
```

### Example 2: Extract Specific Cellular Structures

```matlab
% Load and display overview
SMD = load('Cell_01_Results/SMD.mat');
SMD = SMD.SMD;

RT = smi_helpers.ROITools();
RT.ROI_sizes = [2000, 2000];
RT.Pixel2nm = SMD.PixelSize * 1000;
RT.GaussIm = true;   % Use Gaussian rendering
RT.OriginLLvsUL = false;
RT.SRzoom = 20;

% Select structures of interest
[n_ROIs, RoI, XYsize] = RT.getROI(SMD, 'Select Structures');

% Process each structure separately
for i = 1:n_ROIs
    SMD_structure = RoI{i}.SMD{1};

    % Perform clustering on this structure
    Clust = smi_cluster.Clustering();
    Clust.E = 0.020;  % 20 nm
    Clust.minPts = 5;
    [SMD_clust, ~] = Clust.cluster(SMD_structure);

    % Save result
    save(sprintf('Structure_%02d_clustered.mat', i), 'SMD_clust');
end
```

### Example 3: Reuse ROI Definitions

```matlab
% Define ROIs once from first dataset
SMD1 = load('Condition_A_Cell_01/SMD.mat');
SMD1 = SMD1.SMD;

RT = smi_helpers.ROITools();
RT.ROI_sizes = [4000, 4000];
RT.Pixel2nm = SMD1.PixelSize * 1000;

[n_ROIs, RoI, XYsize] = RT.getROI(SMD1, 'Define ROIs');

% Save ROI boundaries
roi_bounds = zeros(n_ROIs, 4);
for i = 1:n_ROIs
    roi_bounds(i, :) = RoI{i}.ROI;
end
save('Analysis_ROI_bounds.mat', 'roi_bounds', 'RT');

% Apply to multiple datasets
datasets = {
    'Condition_A_Cell_02/SMD.mat',
    'Condition_A_Cell_03/SMD.mat',
    'Condition_B_Cell_01/SMD.mat'
};

for d = 1:length(datasets)
    SMD = load(datasets{d});
    SMD = SMD.SMD;

    for i = 1:n_ROIs
        roi_px = roi_bounds(i, :) / (SMD.PixelSize * 1000);
        roi_px = roi_px([3, 1, 4, 2]);

        SMD_roi = smi_core.SingleMoleculeData.isolateSubROI(SMD, roi_px);

        [~, name, ~] = fileparts(datasets{d});
        save(sprintf('%s_ROI_%02d.mat', name, i), 'SMD_roi');
    end
end
```

## Troubleshooting

### Issue: ROI selection window shows no points

**Check**: Verify SMD contains valid data
```matlab
fprintf('N localizations: %d\n', length(SMD.X));
fprintf('X range: %.2f to %.2f\n', min(SMD.X), max(SMD.X));
fprintf('Y range: %.2f to %.2f\n', min(SMD.Y), max(SMD.Y));
```

**Solutions**:
- Ensure SMD has non-empty X, Y fields
- Check Pixel2nm matches your data (typically 100-108 nm)
- Verify coordinate units are correct

### Issue: ROI contains no localizations

**Check**: Compare ROI bounds to data range
```matlab
roi_nm = RoI{1}.ROI;
fprintf('ROI: X=[%.0f,%.0f] Y=[%.0f,%.0f] nm\n', roi_nm);
fprintf('Data: X=[%.0f,%.0f] Y=[%.0f,%.0f] nm\n', ...
    min(SMD.X)*RT.Pixel2nm, max(SMD.X)*RT.Pixel2nm, ...
    min(SMD.Y)*RT.Pixel2nm, max(SMD.Y)*RT.Pixel2nm);
```

**Solutions**:
- ROI may be outside data bounds
- Check that Pixel2nm is set correctly
- Verify coordinate origin (OriginLLvsUL) matches data

### Issue: Can't delete ROIs during selection

**Solution**:
- Ensure figure window has focus
- Press Backspace or Delete keys (not mouse button)
- Each press removes the most recent ROI

### Issue: Gaussian image display fails

**Solutions**:
- `GaussIm = true` only works for single-label data
- Use point display (`GaussIm = false`) for multi-label
- Reduce SRzoom if rendering is too slow
- Check that SMD structure is provided to getROI

## Best Practices

### ROI Size Selection

- **Too small**: Poor statistics, sensitive to boundaries
- **Too large**: Loss of spatial specificity
- **Recommendation**: 1-5 microns for cellular features
- Use fixed sizes (`RT.ROI_sizes`) for fair comparisons

### Saving ROI Data

Always save:
- ROI boundaries (`RoI{i}.ROI`)
- Pixel conversion factor (`Pixel2nm`)
- Origin convention (`OriginLLvsUL`)
- Source file path (`ResultsFile`)
- Processing date and parameters

### ROI Selection Strategy

1. Generate overview image first (low quality, fast)
2. Identify regions of interest
3. Select ROIs using consistent size
4. Verify ROI placement visually
5. Save ROI definitions immediately
6. Apply strict thresholds during analysis

### Multi-Label ROI Selection

- Ensure channels are registered before selection
- Use `RT.Transform{}` for chromatic aberration correction
- ROI must contain localizations from ALL labels
- Empty labels cause ROI rejection

## Performance Tips

- Point cloud display is faster than Gaussian images
- Reduce marker size (`RT.Msize`) for dense data
- Lower `RT.SRzoom` for faster Gaussian rendering
- Process ROIs sequentially for large datasets
- Clear variables between ROI processing iterations

## See Also

- [ROI Analysis Workflow](../workflows/roi-analysis.md) - Complete ROI analysis pipeline
- [SMD Structure](../core-concepts/smd-structure.md) - Understanding localization data
- [Visualize Results](visualize-results.md) - Plotting and image generation
- `MATLAB/+smi_helpers/@ROITools/ROITools.m` - ROITools class documentation
- `MATLAB/examples/simpleROIcluster.m` - ROI clustering example
