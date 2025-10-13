---
title: "Image Coordinate System in smite"
category: "core-concepts"
level: "beginner"
tags: ["coordinates", "pixels", "positions", "conventions", "transformations"]
prerequisites: ["architecture.md"]
related: ["smd-structure.md", "../how-to/localize-molecules.md"]
summary: "Understanding smite's coordinate system, pixel indexing, sub-pixel precision, and coordinate transformations"
estimated_time: "12 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Image Coordinate System in smite

## Purpose

Understanding coordinate systems is fundamental to working with imaging data in smite. This guide explains how smite represents positions in images, how to work with pixel coordinates and sub-pixel localization precision, and how to perform coordinate transformations. Mastering these concepts is essential for accurate data analysis, visualization, and interpretation of super-resolution results.

## Prerequisites

- Basic understanding of [smite architecture](architecture.md)
- Familiarity with MATLAB array indexing
- Basic concepts of digital imaging

## Overview

smite follows MATLAB's image coordinate conventions with specific rules for pixel positioning and sub-pixel precision. Key concepts:

- **Pixel-centered coordinates**: Coordinate (1,1) is at the center of the top-left pixel, not the corner
- **Column-major ordering**: Arrays are indexed as image(row, column) or image(Y, X)
- **Sub-pixel precision**: Localization positions are real numbers (e.g., 10.347 pixels), not integers
- **Units**: Positions in SMD are in pixels; conversion to physical units requires PixelSize

Understanding these conventions prevents off-by-one errors, ensures correct image rendering, and enables accurate spatial analysis.

## Coordinate System Basics

### Pixel-Centered Coordinates

In smite (and MATLAB), coordinates refer to pixel **centers**, not edges or corners:

```
        1.0       2.0       3.0       4.0       5.0
         |         |         |         |         |
    +----+----+----+----+----+----+----+----+----+----+
    |         |         |         |         |         |
1.0 |  (1,1)  |  (1,2)  |  (1,3)  |  (1,4)  |  (1,5)  |
    |         |         |         |         |         |
    +----+----+----+----+----+----+----+----+----+----+
    |         |         |         |         |         |
2.0 |  (2,1)  |  (2,2)  |  (2,3)  |  (2,4)  |  (2,5)  |
    |         |         |         |         |         |
    +----+----+----+----+----+----+----+----+----+----+
    |         |         |         |         |         |
3.0 |  (3,1)  |  (3,2)  |  (3,3)  |  (3,4)  |  (3,5)  |
    |         |         |         |         |         |
    +----+----+----+----+----+----+----+----+----+----+

Legend:
  + = pixel corner (at half-integer coordinates: 0.5, 1.5, 2.5, etc.)
  (Y,X) = pixel center coordinates
  | and - = pixel boundaries
```

**Key points:**
- Coordinate (1,1) is at the **center** of the top-left pixel
- Pixel corners are at half-integer positions (0.5, 1.5, 2.5, etc.)
- The image extends from 0.5 to (size+0.5) along each dimension

**Example:**
```matlab
% Create a 5×5 image
image = rand(5, 5);

% The pixel at coordinate (1,1) contains image(1,1)
pixel_value = image(1, 1);

% The physical extent is from 0.5 to 5.5 in both dimensions
x_range = [0.5, 5.5];
y_range = [0.5, 5.5];

% Display with correct axis limits
imagesc(x_range, y_range, image);
axis equal tight;
```

### Array Indexing: Row-Column Order

MATLAB arrays use **row-major indexing**: the first index is the row (Y-axis), the second is the column (X-axis):

```matlab
% Array indexing: array(row, column)
%                      (Y,    X)

% Create image
image = magic(5);

% Access pixel at row 2, column 3
% This is Y=2, X=3 in coordinate terms
value = image(2, 3);

% IMPORTANT: Y comes first, X comes second
Y_index = 2;
X_index = 3;
value = image(Y_index, X_index);  % Correct order
```

**Axes directions:**
```
        X increasing →
    +------------------------+
Y   | (1,1)  (1,2)  (1,3)   |
↓   | (2,1)  (2,2)  (2,3)   |
    | (3,1)  (3,2)  (3,3)   |
    +------------------------+
```

### SMD Position Fields

In SMD structures, X and Y are separate vectors:

```matlab
% SMD uses separate X and Y fields
% X: horizontal position (column)
% Y: vertical position (row)

% Example: Create SMD with 3 localizations
SMD.X = [10.5; 20.3; 15.7];  % X positions (columns)
SMD.Y = [12.1; 18.9; 25.4];  % Y positions (rows)

% Plot positions
figure;
plot(SMD.X, SMD.Y, 'r.', 'MarkerSize', 20);
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Localization Positions');
axis equal;

% Note: plot() uses (X, Y) order naturally
% This matches SMD.X and SMD.Y
```

## Sub-Pixel Precision

### Localization Coordinates are Real Numbers

Unlike pixel indices (which are integers), localization positions are **real-valued** with sub-pixel precision:

```matlab
% Pixel indices are integers
pixel_row = 10;      % Integer
pixel_col = 15;      % Integer
pixel_value = image(pixel_row, pixel_col);

% Localization positions are real numbers
loc_Y = 10.347;      % Real number with sub-pixel precision
loc_X = 15.892;      % Real number with sub-pixel precision

% This molecule is 0.347 pixels down from the center of pixel (10, 15)
% and 0.892 pixels right from that center
```

**Physical interpretation:**

A localization at (X=15.892, Y=10.347) means:
- The molecule is centered at X=15.892 pixels from the left edge
- And Y=10.347 pixels from the top edge
- With precision typically 0.05-0.2 pixels (5-20 nm for 100 nm pixels)

```matlab
% Typical localization precision
SMD.X_SE = 0.08;  % Standard error ~8% of pixel
SMD.Y_SE = 0.07;  % ~7% of pixel

% Convert to nanometers
pixel_size_um = 0.1;  % 100 nm pixels
precision_nm = SMD.X_SE * pixel_size_um * 1000;  % ~8 nm
fprintf('Localization precision: %.1f nm\n', precision_nm);
```

### Interpolation and Sub-Pixel Operations

When working with sub-pixel coordinates, use interpolation:

```matlab
% Extract value at sub-pixel position
X_loc = 15.892;
Y_loc = 10.347;

% Use interp2 for sub-pixel interpolation
value_interp = interp2(image, X_loc, Y_loc, 'linear');

% Or for many positions
X_locs = SMD.X;
Y_locs = SMD.Y;
values = interp2(image, X_locs, Y_locs, 'linear');
```

### Box Corners and ROIs

When fitting molecules, smite extracts small boxes (ROIs) around each detected spot. The box corner positions are stored:

```matlab
% Box fitting parameters
BoxSize = 11;  % 11×11 pixel box

% Molecule at sub-pixel position
X_center = 25.4;
Y_center = 18.7;

% Box corner calculation (top-left corner)
% Center the box on the localization
XBoxCorner = round(X_center - BoxSize/2);  % = 20
YBoxCorner = round(Y_center - BoxSize/2);  % = 13

% Box extends from corner to corner+BoxSize
X_box_range = XBoxCorner : (XBoxCorner + BoxSize - 1);  % [20...30]
Y_box_range = YBoxCorner : (YBoxCorner + BoxSize - 1);  % [13...23]

% Extract box from image
box = image(Y_box_range, X_box_range);

% Position relative to box corner
X_in_box = X_center - XBoxCorner;  % = 5.4 pixels from left edge
Y_in_box = Y_center - YBoxCorner;  % = 5.7 pixels from top edge
```

## Coordinate Transformations

### Pixel to Physical Units

Convert pixel coordinates to physical units (micrometers or nanometers):

```matlab
% Load results
load('Results.mat', 'SMD');

% Pixel coordinates
X_pixels = SMD.X;  % pixels
Y_pixels = SMD.Y;  % pixels

% Physical units
pixel_size = SMD.PixelSize;  % micrometers (e.g., 0.1 um = 100 nm)

% Convert to micrometers
X_um = X_pixels * pixel_size;
Y_um = Y_pixels * pixel_size;

% Convert to nanometers
X_nm = X_pixels * pixel_size * 1000;
Y_nm = Y_pixels * pixel_size * 1000;

% Example
fprintf('Position: (%.2f, %.2f) pixels\n', X_pixels(1), Y_pixels(1));
fprintf('         = (%.3f, %.3f) μm\n', X_um(1), Y_um(1));
fprintf('         = (%.1f, %.1f) nm\n', X_nm(1), Y_nm(1));
```

### Physical to Pixel Units

Convert physical measurements back to pixels:

```matlab
% Define ROI in micrometers
ROI_um = [5.0, 5.0, 10.0, 10.0];  % [Y_min, X_min, Y_max, X_max] in μm

% Convert to pixels
pixel_size = SMD.PixelSize;  % μm
ROI_pixels = ROI_um / pixel_size;

% Select localizations in ROI
in_ROI = (SMD.X >= ROI_pixels(2)) & (SMD.X <= ROI_pixels(4)) & ...
         (SMD.Y >= ROI_pixels(1)) & (SMD.Y <= ROI_pixels(3));

fprintf('%d localizations in ROI\n', sum(in_ROI));
```

### Drift Correction

Drift correction adjusts all positions by frame-dependent offsets:

```matlab
% Drift fields in SMD
% DriftX, DriftY: NFrames × NDatasets matrices

% Positions in SMD.X and SMD.Y are ALREADY drift-corrected
% if SMF.DriftCorrection.On = true

% To see uncorrected positions
if isfield(SMD, 'DriftX') && ~isempty(SMD.DriftX)
    % Get drift for each localization
    for i = 1:length(SMD.X)
        frame = SMD.FrameNum(i);
        dataset = SMD.DatasetNum(i);
        drift_x = SMD.DriftX(frame, dataset);
        drift_y = SMD.DriftY(frame, dataset);

        % Uncorrected position
        X_uncorrected = SMD.X(i) - drift_x;
        Y_uncorrected = SMD.Y(i) - drift_y;
    end

    % Visualize drift
    figure;
    plot(SMD.DriftX(:,1), SMD.DriftY(:,1), '-o');
    xlabel('Drift X (pixels)'); ylabel('Drift Y (pixels)');
    title('Stage Drift Over Time');
    axis equal;
end
```

### Channel Registration

For multi-color imaging, coordinates may be transformed between channels:

```matlab
% Channel registration transforms coordinates from one channel to another
% using a spatial transformation (affine, polynomial, etc.)

% After registration, SMD.IsTransformed = true
if isfield(SMD, 'IsTransformed') && SMD.IsTransformed
    fprintf('Positions are in registered coordinate system\n');
    fprintf('Registration error: %.3f pixels (%.1f nm)\n', ...
        SMD.RegError, SMD.RegError * SMD.PixelSize * 1000);
end

% For manual transformation (not typical)
% Use MATLAB's spatial transformation functions
tform = affine2d([1 0 0; 0 1 0; 5 10 1]);  % Example: shift by (5,10)
[X_new, Y_new] = transformPointsForward(tform, SMD.X, SMD.Y);
```

### Image Coordinate to Display Coordinate

When displaying images with imagesc, coordinates may need adjustment:

```matlab
% Create super-resolution image
SR_zoom = 20;  % 20× magnification
SR_size = [SMD.YSize, SMD.XSize] * SR_zoom;

% Positions in SR coordinates
X_SR = SMD.X * SR_zoom;
Y_SR = SMD.Y * SR_zoom;

% Create image
SR_image = zeros(SR_size);
% ... populate SR_image ...

% Display with correct extent
x_extent = [0.5, SMD.XSize + 0.5];
y_extent = [0.5, SMD.YSize + 0.5];

figure;
imagesc(x_extent, y_extent, SR_image);
hold on;
plot(SMD.X, SMD.Y, 'r.', 'MarkerSize', 10);
xlabel('X (pixels)'); ylabel('Y (pixels)');
axis equal tight;
colormap hot;
```

## Common Pitfalls and Solutions

### Pitfall 1: Swapping X and Y

**Problem:** Confusing row/column with X/Y

```matlab
% WRONG: Swapped X and Y
wrong_value = image(SMD.X(i), SMD.Y(i));  % X should be column!

% CORRECT: Y is row, X is column
correct_value = interp2(image, SMD.X(i), SMD.Y(i));
```

### Pitfall 2: Integer vs. Sub-Pixel Positions

**Problem:** Rounding localization positions inappropriately

```matlab
% WRONG: Loses sub-pixel precision
X_rounded = round(SMD.X);
Y_rounded = round(SMD.Y);

% CORRECT: Keep sub-pixel precision
X_exact = SMD.X;  % Real-valued positions
Y_exact = SMD.Y;
```

### Pitfall 3: Incorrect Image Extent

**Problem:** Plotting without accounting for pixel centers

```matlab
% WRONG: Assumes pixels go from 1 to N
imagesc(image);  % Default extent is [1, N]

% CORRECT: Account for pixel centers at 0.5 offset
imagesc([0.5, SMD.XSize+0.5], [0.5, SMD.YSize+0.5], image);
```

### Pitfall 4: Unit Confusion

**Problem:** Mixing pixels and physical units

```matlab
% WRONG: Mixing units
ROI_nm = [100, 100, 500, 500];  % In nanometers
in_ROI = (SMD.X >= ROI_nm(2));  % SMD.X is in pixels!

% CORRECT: Convert to same units
ROI_pixels = ROI_nm / (SMD.PixelSize * 1000);  % Convert nm to pixels
in_ROI = (SMD.X >= ROI_pixels(2));
```

## Practical Examples

### Example 1: Extracting Sub-Pixel Intensities

```matlab
% Load data
load('Results.mat', 'SMD');
image = imread('raw_frame.tif');

% Get sub-pixel intensities at localization positions
intensities = zeros(length(SMD.X), 1);
for i = 1:length(SMD.X)
    intensities(i) = interp2(double(image), SMD.X(i), SMD.Y(i), 'cubic');
end

% Plot intensity distribution
histogram(intensities, 50);
xlabel('Interpolated Intensity'); ylabel('Count');
title('Sub-Pixel Intensities at Localizations');
```

### Example 2: Creating a Precision Map

```matlab
% Create image showing localization precision across field of view
SR_zoom = 10;
SR_size = [SMD.YSize, SMD.XSize] * SR_zoom;

% Initialize precision map
precision_map = zeros(SR_size);
count_map = zeros(SR_size);

% Accumulate precision values
for i = 1:length(SMD.X)
    X_SR = round(SMD.X(i) * SR_zoom);
    Y_SR = round(SMD.Y(i) * SR_zoom);

    if X_SR >= 1 && X_SR <= SR_size(2) && Y_SR >= 1 && Y_SR <= SR_size(1)
        combined_SE = sqrt(SMD.X_SE(i)^2 + SMD.Y_SE(i)^2);
        precision_map(Y_SR, X_SR) = precision_map(Y_SR, X_SR) + combined_SE;
        count_map(Y_SR, X_SR) = count_map(Y_SR, X_SR) + 1;
    end
end

% Average precision per pixel
precision_map(count_map > 0) = precision_map(count_map > 0) ./ count_map(count_map > 0);

% Convert to nm
precision_map_nm = precision_map * SMD.PixelSize * 1000;

% Display
figure;
imagesc(precision_map_nm);
colorbar;
xlabel('X (SR pixels)'); ylabel('Y (SR pixels)');
title('Localization Precision Map (nm)');
axis equal tight;
```

### Example 3: ROI Selection with Physical Units

```matlab
% Define ROI in physical units (micrometers)
ROI_center_um = [5.0, 5.0];  % [X, Y] center in μm
ROI_size_um = 2.0;            % 2 μm × 2 μm

% Convert to pixel coordinates
pixel_size = SMD.PixelSize;
ROI_center_px = ROI_center_um / pixel_size;
ROI_half_size_px = (ROI_size_um / 2) / pixel_size;

% Define bounds
X_min = ROI_center_px(1) - ROI_half_size_px;
X_max = ROI_center_px(1) + ROI_half_size_px;
Y_min = ROI_center_px(2) - ROI_half_size_px;
Y_max = ROI_center_px(2) + ROI_half_size_px;

% Filter localizations
in_ROI = (SMD.X >= X_min) & (SMD.X <= X_max) & ...
         (SMD.Y >= Y_min) & (SMD.Y <= Y_max);

% Extract ROI data
SMD_ROI.X = SMD.X(in_ROI);
SMD_ROI.Y = SMD.Y(in_ROI);
SMD_ROI.Photons = SMD.Photons(in_ROI);
SMD_ROI.PixelSize = SMD.PixelSize;

fprintf('Selected %d / %d localizations in ROI\n', sum(in_ROI), length(SMD.X));

% Visualize
figure;
plot(SMD.X, SMD.Y, 'k.', 'MarkerSize', 2); hold on;
plot(SMD_ROI.X, SMD_ROI.Y, 'r.', 'MarkerSize', 10);
rectangle('Position', [X_min, Y_min, 2*ROI_half_size_px, 2*ROI_half_size_px], ...
    'EdgeColor', 'b', 'LineWidth', 2);
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('ROI Selection');
axis equal;
legend('All localizations', 'ROI localizations', 'ROI boundary');
```

### Example 4: Coordinate System Verification

```matlab
% Verify coordinate system understanding
fprintf('=== Coordinate System Verification ===\n\n');

% Image properties
fprintf('Image size: %d × %d pixels\n', SMD.YSize, SMD.XSize);
fprintf('Pixel size: %.1f nm\n', SMD.PixelSize * 1000);

% Coordinate ranges
X_min = min(SMD.X);
X_max = max(SMD.X);
Y_min = min(SMD.Y);
Y_max = max(SMD.Y);

fprintf('\nLocalization coordinate ranges:\n');
fprintf('  X: %.2f to %.2f pixels\n', X_min, X_max);
fprintf('  Y: %.2f to %.2f pixels\n', Y_min, Y_max);

% Check if within image bounds
in_bounds = (SMD.X >= 0.5) & (SMD.X <= SMD.XSize+0.5) & ...
            (SMD.Y >= 0.5) & (SMD.Y <= SMD.YSize+0.5);
fprintf('\n%.1f%% of localizations within image bounds\n', 100*mean(in_bounds));

% Sub-pixel precision check
X_fractional = abs(SMD.X - round(SMD.X));
has_subpixel = mean(X_fractional > 0.01);
fprintf('\n%.1f%% of positions have sub-pixel precision\n', 100*has_subpixel);

% Precision statistics
median_precision_nm = median(sqrt(SMD.X_SE.^2 + SMD.Y_SE.^2)) * SMD.PixelSize * 1000;
fprintf('Median localization precision: %.1f nm\n', median_precision_nm);
```

## Key Takeaways

1. **Pixel centers**: Coordinate (1,1) is at the center of the top-left pixel, not a corner

2. **Row-column order**: Arrays use image(Y, X) or image(row, column) indexing

3. **Sub-pixel precision**: Localization coordinates are real numbers, not integers

4. **SMD conventions**: SMD.X is horizontal (column), SMD.Y is vertical (row)

5. **Unit conversions**: Multiply pixel coordinates by PixelSize to get physical units

6. **Image extent**: Images span from 0.5 to (size+0.5) along each dimension

7. **Always verify**: Check coordinate ranges and units to catch errors early

## See Also

- [SMD Structure](smd-structure.md) - Position fields and metadata
- [Architecture Overview](architecture.md) - Overall coordinate system context
- [How to Localize Molecules](../how-to/localize-molecules.md) - Generating sub-pixel positions
- doc/DataStructures/SMD.md (in repository) - Complete SMD field reference
