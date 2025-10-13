---
title: "How to Visualize SMLM/SPT Results"
category: "how-to"
level: "beginner"
tags: ["visualization", "images", "plots", "export", "publication"]
prerequisites: ["localize-molecules.md", "../core-concepts/smd-structure.md"]
related: ["../workflows/smlm-analysis.md", "../workflows/spt-tracking.md", "../examples/basic-localization.md"]
summary: "Comprehensive guide to creating publication-quality visualizations from SMLM and SPT data"
estimated_time: "20 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# How to Visualize SMLM/SPT Results

## Purpose

Visualization transforms raw localization data into interpretable images and plots. This guide covers all visualization methods in smite, from basic scatter plots to publication-ready super-resolution images. You'll learn how to render localizations, apply colormaps, add scale bars, visualize drift, plot trajectories, and export high-quality figures.

## Prerequisites

- SMD or TR structure with localization data
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Completed localization or tracking analysis
- Basic MATLAB plotting knowledge

## Overview

smite provides visualization tools for:

1. **Super-resolution images**: Gaussian rendering, histograms, circles
2. **Drift visualization**: Color-coded temporal information
3. **Trajectory plots**: SPT motion visualization
4. **Quality assessment**: Distribution plots and statistics
5. **Multi-channel overlays**: Color-combined images
6. **Publication export**: Scale bars, colormaps, high-resolution output

All methods are in the `smi_vis.GenerateImages` class and use SMD/TR structures as input.

## Basic Visualization: Scatter Plots

### Simple Localization Scatter Plot

The simplest visualization is plotting X and Y coordinates:

```matlab
% Load results
load('Results.mat', 'SMD');

% Basic scatter plot
figure;
plot(SMD.X, SMD.Y, '.', 'MarkerSize', 2);
axis equal; axis tight;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('%d Localizations', length(SMD.X)));
```

### Color by Frame Number

Show temporal information:

```matlab
figure;
scatter(SMD.X, SMD.Y, 5, SMD.FrameNum, 'filled');
axis equal; axis tight;
colormap jet; colorbar;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Localizations Colored by Frame');
```

### Color by Precision

Visualize localization quality:

```matlab
% Compute precision in nm
precision_nm = SMD.X_SE * SMD.PixelSize * 1000;

figure;
scatter(SMD.X, SMD.Y, 10, precision_nm, 'filled');
axis equal; axis tight;
colormap hot; colorbar;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Localizations Colored by Precision (nm)');
caxis([0, 30]);  % Set color scale
```

## Super-Resolution Images

Super-resolution images render localizations at higher resolution than raw camera pixels.

### Gaussian Rendering

Creates smooth, publication-quality images by rendering each localization as a Gaussian blob:

```matlab
% Load results
load('Results.mat', 'SMD');

% Gaussian rendering
SR_zoom = 20;         % 20× zoom over raw pixel size
scalebar_um = 10;     % 10 μm scale bar

% Generate image
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, SR_zoom, scalebar_um);

% Display
figure;
imshow(GaussIm);
title(sprintf('Gaussian SR Image (%d× zoom)', SR_zoom));
```

**How it works:**
- Each localization rendered as 2D Gaussian
- Width = localization precision (X_SE, Y_SE)
- Intensity = uniform (all equal weight)
- Images scaled to [0, 1] and colored with hot colormap

**Choosing zoom factor:**

```matlab
% Rule of thumb: effective pixel size 5-10 nm
effective_pixel_nm = SMD.PixelSize * 1000 / SR_zoom;

% For 108 nm pixels:
SR_zoom = 20;  % → 5.4 nm effective pixels
SR_zoom = 10;  % → 10.8 nm effective pixels

fprintf('Effective SR pixel size: %.1f nm\n', effective_pixel_nm);
```

**Without scale bar:**

```matlab
% Set scalebar length to 0 to disable
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, SR_zoom, 0);
```

### Histogram Images

Creates pixelated images by binning localizations:

```matlab
% Histogram rendering
SR_zoom = 10;
colormap_name = 'hot';

% Generate image (both grayscale and RGB)
[HistIm, RgbHistIm] = smi_vis.GenerateImages.histogramImage(...
    SMD, SR_zoom, colormap_name);

% Display
figure;
subplot(1,2,1);
imshow(HistIm, []); colormap hot; colorbar;
title('Histogram Image (Grayscale)');

subplot(1,2,2);
imshow(RgbHistIm);
title('Histogram Image (RGB)');
```

**How it works:**
- Localizations binned into pixels
- Each pixel counts number of localizations
- Fast, good for dense data

**Available colormaps:**

```matlab
% Common colormaps
'hot'     % Black → red → yellow → white (default)
'jet'     % Blue → cyan → yellow → red
'parula'  % MATLAB default blue → yellow
'gray'    % Grayscale
'viridis' % Perceptually uniform

% Example with different colormaps
figure;
for i = 1:4
    cmaps = {'hot', 'jet', 'parula', 'gray'};
    [~, RgbIm] = smi_vis.GenerateImages.histogramImage(SMD, 10, cmaps{i});
    subplot(2,2,i);
    imshow(RgbIm);
    title(cmaps{i});
end
```

### Comparison: Gaussian vs Histogram

```matlab
% Generate both
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 0);
[~, HistIm] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');

% Display side-by-side
figure('Position', [100, 100, 1200, 500]);

subplot(1,2,1);
imshow(GaussIm);
title('Gaussian Rendering (Smooth)');

subplot(1,2,2);
imshow(HistIm);
title('Histogram Rendering (Pixelated)');
```

**When to use which:**
- **Gaussian**: Publication figures, smooth appearance, precision-weighted
- **Histogram**: Quick preview, dense data, exact counts

## Circle Images

Visualizes each localization as a circle with radius proportional to precision:

```matlab
% Load results
load('Results.mat', 'SMD');

% Generate circle image
SR_zoom = 20;
ColorMap = [1, 0, 0];  % Red circles (RGB normalized to 1)

[CircleImage, CircleImageRGB] = smi_vis.GenerateImages.circleImage(...
    SMD, ColorMap, SR_zoom);

% Display
figure;
subplot(1,2,1);
imshow(CircleImage);
title('Circle Image (Binary)');

subplot(1,2,2);
imshow(CircleImageRGB);
title('Circle Image (Color)');
```

**Understanding circle sizes:**

Circles have radius = mean(X_SE, Y_SE), representing localization uncertainty:

```matlab
% Visualize precision
figure;
imshow(CircleImageRGB);
title('Circle sizes = Localization Precision');

% Larger circles = higher uncertainty
% Smaller circles = better precision
```

**Custom circle colors:**

```matlab
% Color by photon count
Photons_normalized = (SMD.Photons - min(SMD.Photons)) / ...
    (max(SMD.Photons) - min(SMD.Photons));
ColorMap = [Photons_normalized, zeros(size(Photons_normalized)), ...
    1-Photons_normalized];  % Red (bright) to blue (dim)

[~, CircleImageRGB] = smi_vis.GenerateImages.circleImage(...
    SMD, ColorMap, 20);

figure; imshow(CircleImageRGB);
title('Circles Colored by Photon Count');
```

**Adjusting circle sizes:**

```matlab
% Scale factor makes circles larger/smaller
SEScaleFactor = 3;  % 3× larger circles

[~, CircleImageRGB] = smi_vis.GenerateImages.circleImage(...
    SMD, [1, 0, 0], 20, [], SEScaleFactor);

figure; imshow(CircleImageRGB);
title('Circles Scaled 3×');
```

**Auto-sizing circles:**

```matlab
% Automatically determine zoom based on desired circle size
MinPixelsPerCircle = 16;  % Minimum circle size in pixels
SR_zoom = [];             % Let function determine zoom

[CircleImage, CircleImageRGB, ActualZoom] = ...
    smi_vis.GenerateImages.circleImage(SMD, [1, 0, 0], SR_zoom, ...
    MinPixelsPerCircle, 1);

fprintf('Automatically selected zoom: %d\n', ActualZoom);
```

## Drift Visualization

Visualizes temporal drift by color-coding localizations:

```matlab
% Load results
load('Results.mat', 'SMD');

% Generate drift image
SR_zoom = 10;
[DriftIm, DriftImRGB] = smi_vis.GenerateImages.driftImage(SMD, SR_zoom);

% Display
figure;
imshow(DriftImRGB);
title('Drift Image (Color = Time)');
colormap jet; colorbar;

% Early frames = blue
% Middle frames = green/yellow
% Late frames = red
```

**How it works:**
- Uses jet colormap to encode time
- Blue → cyan → green → yellow → red across frames
- Visualizes stage drift or sample movement

**Circle drift images:**

Combines circles with temporal colors:

```matlab
% Circle drift image
SR_zoom = 20;
CircleDriftIm = smi_vis.GenerateImages.circleDriftImage(SMD, SR_zoom);

figure;
imshow(CircleDriftIm);
title('Circle Drift Image');

% Circle sizes = precision
% Circle colors = time
```

**Interpreting drift:**
- Rainbow streaks → uncorrected drift
- Single color region → stable
- Color shifts → temporal changes

## Trajectory Visualization (SPT)

For tracking results (TR structures):

### Basic Trajectory Plot

```matlab
% Load tracking results
load('TrackingResults.mat', 'TR');

% Plot all trajectories
figure;
hold on;
for i = 1:length(TR)
    plot(TR(i).X, TR(i).Y, '-');
end
axis equal; axis tight;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('%d Trajectories', length(TR)));
```

### Color by Trajectory

```matlab
% Each trajectory a different color
figure; hold on;
colors = lines(length(TR));  % Distinct colors

for i = 1:length(TR)
    plot(TR(i).X, TR(i).Y, '-', 'Color', colors(i,:), 'LineWidth', 1.5);
end

axis equal; axis tight;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Trajectories (Color = ID)');
```

### 3D Trajectory Plot (Time as Z-axis)

```matlab
% Plot X, Y, Time
figure; hold on;

for i = 1:length(TR)
    plot3(TR(i).X, TR(i).Y, TR(i).FrameNum, '-', 'LineWidth', 2);
end

xlabel('X (pixels)'); ylabel('Y (pixels)'); zlabel('Frame');
title('Trajectories in Space-Time');
grid on; rotate3d on;
```

### Trajectory Movies

Using `smi_vis.GenerateMovies`:

```matlab
% Create movie showing trajectories building over time
GM = smi_vis.GenerateMovies();

% Load raw data and tracking results
load('RawData.mat', 'sequence');
load('TrackingResults.mat', 'TR');

% Configure movie parameters
GM.SMD = TR;              % Tracking results
GM.Sequence = sequence;   % Raw data
GM.FrameRate = 10;        % 10 fps output
GM.TrailLength = 20;      % Show 20-frame trail

% Generate and save movie
GM.generateMovie();
GM.saveMovie('trajectories.mp4');
```

## Scale Bars and Annotations

### Adding Scale Bars

```matlab
% Load and render image
load('Results.mat', 'SMD');
SR_zoom = 20;
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, SR_zoom, 0);  % No auto scale bar

% Add custom scale bar
PixelSize = SMD.PixelSize / SR_zoom;  % Effective pixel size
ScalebarLength_um = 5;                % 5 μm bar
Location = 'bottomright';

ImageWithBar = smi_vis.GenerateImages.scalebar(...
    GaussIm, PixelSize, ScalebarLength_um, Location);

figure; imshow(ImageWithBar);
title('Image with 5 μm Scale Bar');
```

**Scale bar locations:**

```matlab
% Available locations
'bottomright'  % Default
'bottomleft'
'topright'
'topleft'

% Example: all four corners
figure;
for i = 1:4
    locs = {'bottomright', 'bottomleft', 'topright', 'topleft'};
    subplot(2,2,i);
    ImWithBar = smi_vis.GenerateImages.scalebar(GaussIm, PixelSize, 5, locs{i});
    imshow(ImWithBar);
    title(locs{i});
end
```

### Custom Annotations

```matlab
% Start with image
figure; imshow(ImageWithBar);

% Add text annotation
text(50, 50, 'Region of Interest', ...
    'Color', 'white', 'FontSize', 14, 'FontWeight', 'bold');

% Add arrow pointing to feature
annotation('arrow', [0.3, 0.4], [0.7, 0.6], ...
    'Color', 'yellow', 'LineWidth', 2);

% Add rectangle ROI
rectangle('Position', [100, 100, 200, 200], ...
    'EdgeColor', 'cyan', 'LineWidth', 2);
```

## Multi-Channel Overlays

Combine multiple channels with distinct colors:

```matlab
% Load two-channel data
load('Channel1_Results.mat', 'SMD');
SMD1 = SMD;
load('Channel2_Results.mat', 'SMD');
SMD2 = SMD;

% Generate images for each channel
SR_zoom = 20;
Im1 = smi_vis.GenerateImages.gaussianImage(SMD1, SR_zoom, 0);
Im2 = smi_vis.GenerateImages.gaussianImage(SMD2, SR_zoom, 0);

% Convert to grayscale if needed (remove colormap)
Im1_gray = mat2gray(Im1);
Im2_gray = mat2gray(Im2);

% Stack images
ImageStack = cat(3, Im1_gray, Im2_gray);

% Create overlay (green and magenta)
[OverlayIm, ColorOrder] = smi_vis.GenerateImages.overlayNImages(ImageStack);

figure; imshow(OverlayIm);
title(sprintf('Two-Channel Overlay (%s)', ColorOrder));
% ColorOrder = 'GM' means channel 1=green, channel 2=magenta
```

**Three-channel overlay:**

```matlab
% Load three channels
load('Channel1_Results.mat', 'SMD'); SMD1 = SMD;
load('Channel2_Results.mat', 'SMD'); SMD2 = SMD;
load('Channel3_Results.mat', 'SMD'); SMD3 = SMD;

% Generate images
Im1 = smi_vis.GenerateImages.gaussianImage(SMD1, 20, 0);
Im2 = smi_vis.GenerateImages.gaussianImage(SMD2, 20, 0);
Im3 = smi_vis.GenerateImages.gaussianImage(SMD3, 20, 0);

% Stack and overlay (green, magenta, cyan)
ImageStack = cat(3, mat2gray(Im1), mat2gray(Im2), mat2gray(Im3));
[OverlayIm, ColorOrder] = smi_vis.GenerateImages.overlayNImages(ImageStack);

figure; imshow(OverlayIm);
title(sprintf('Three-Channel Overlay (%s)', ColorOrder));
% ColorOrder = 'GMC' means: 1=green, 2=magenta, 3=cyan
```

**Color schemes:**
- 2 channels: Green + Magenta
- 3 channels: Green + Magenta + Cyan
- 4 channels: Green + Magenta + Cyan + Yellow

## Custom Colormaps

### Applying Custom Colors

```matlab
% Generate grayscale histogram
[HistIm, ~] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');

% Apply custom colormap
CustomMap = [linspace(0, 1, 256)', zeros(256, 1), linspace(1, 0, 256)'];  % Blue to red

RgbCustom = smi_vis.GenerateImages.colorImage(HistIm, CustomMap);

figure; imshow(RgbCustom);
title('Custom Blue-to-Red Colormap');
```

### Colormap Examples

```matlab
% Perceptually uniform (recommended for quantitative data)
viridis_map = viridis(256);
[~, Im_viridis] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'parula');

% High contrast (good for presentations)
hot_map = hot(256);
[~, Im_hot] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');

% Diverging (for showing deviations)
N = 256;
RedBlue = [linspace(1, 0, N/2), linspace(0, 1, N/2); ...
           linspace(0, 0, N/2), linspace(0, 0, N/2); ...
           linspace(0, 1, N/2), linspace(1, 0, N/2)]';
```

### Manual Color Scaling

```matlab
% Generate grayscale image
[HistIm, ~] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');

% Manual min/max scaling
MinMax = [0, 10];  % Only show 0-10 localizations per pixel
ColorMap = hot(256);

RgbScaled = smi_vis.GenerateImages.colorImage(HistIm, ColorMap, MinMax);

figure; imshow(RgbScaled);
title('Manual Color Scaling [0, 10]');
```

## Quality Assessment Plots

### Precision Distribution

```matlab
% Load results
load('Results.mat', 'SMD');

% Convert to nm
precision_nm = SMD.X_SE * SMD.PixelSize * 1000;

% Histogram
figure;
histogram(precision_nm, 50, 'FaceColor', [0.3, 0.6, 0.9]);
xlabel('Localization Precision (nm)');
ylabel('Count');
title('Precision Distribution');
grid on;

% Add median line
median_prec = median(precision_nm);
xline(median_prec, 'r--', 'LineWidth', 2);
text(median_prec+1, max(ylim)*0.9, ...
    sprintf('Median: %.1f nm', median_prec), ...
    'FontSize', 12, 'Color', 'red');
```

### Photon Distribution

```matlab
figure;
histogram(SMD.Photons, 50, 'FaceColor', [0.9, 0.4, 0.3]);
xlabel('Detected Photons');
ylabel('Count');
title('Photon Distribution');
grid on;

% Add statistics
mean_photons = mean(SMD.Photons);
xline(mean_photons, 'k--', 'LineWidth', 2);
text(mean_photons+100, max(ylim)*0.9, ...
    sprintf('Mean: %.0f', mean_photons), 'FontSize', 12);
```

### Localization Density Map

```matlab
% 2D histogram of localization density
figure;
histogram2(SMD.X, SMD.Y, 100, 'DisplayStyle', 'tile');
colormap hot; colorbar;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Localization Density Map');
axis equal; axis tight;
```

### Multi-Panel Quality Figure

```matlab
figure('Position', [100, 100, 1400, 800]);

% Panel 1: Localizations
subplot(2,3,1);
plot(SMD.X, SMD.Y, '.', 'MarkerSize', 1);
axis equal; axis tight;
title(sprintf('%d Localizations', length(SMD.X)));

% Panel 2: Precision
subplot(2,3,2);
precision_nm = SMD.X_SE * SMD.PixelSize * 1000;
histogram(precision_nm, 50);
xlabel('Precision (nm)');
title(sprintf('Median: %.1f nm', median(precision_nm)));

% Panel 3: Photons
subplot(2,3,3);
histogram(SMD.Photons, 50);
xlabel('Photons');
title(sprintf('Mean: %.0f', mean(SMD.Photons)));

% Panel 4: Background
subplot(2,3,4);
histogram(SMD.Bg, 50);
xlabel('Background (photons/pixel)');
title(sprintf('Mean: %.1f', mean(SMD.Bg)));

% Panel 5: Frame distribution
subplot(2,3,5);
histogram(SMD.FrameNum, 50);
xlabel('Frame Number');
ylabel('Localizations');
title('Localizations per Frame');

% Panel 6: Precision vs Photons
subplot(2,3,6);
scatter(SMD.Photons, precision_nm, 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Photons'); ylabel('Precision (nm)');
title('Precision vs Photons');
grid on;
```

## Export for Publication

### High-Resolution PNG

```matlab
% Generate image
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 10);

% Create figure with no borders
fig = figure('Color', 'white');
imshow(GaussIm);
axis off;

% Export high-resolution PNG
exportgraphics(fig, 'SR_Image.png', 'Resolution', 300);  % 300 DPI
```

### TIFF Export

```matlab
% Generate RGB image
[~, RgbIm] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');

% Save as TIFF
imwrite(RgbIm, 'SR_Image.tif');
```

### Vector Graphics (EPS/PDF)

```matlab
% Create figure
fig = figure('Color', 'white');
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 10);
imshow(GaussIm);

% Export as vector graphics
print(fig, 'SR_Image.eps', '-depsc', '-r300');  % EPS
print(fig, 'SR_Image.pdf', '-dpdf', '-r300');   % PDF
```

### Multi-Panel Publication Figure

```matlab
% Create figure with specific size (in inches)
fig = figure('Units', 'inches', 'Position', [0, 0, 7, 5], 'Color', 'white');

% Panel A: Gaussian rendering
subplot(2,3,1);
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 5);
imshow(GaussIm);
title('A', 'FontSize', 14, 'FontWeight', 'bold');

% Panel B: Histogram rendering
subplot(2,3,2);
[~, HistIm] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');
imshow(HistIm);
title('B', 'FontSize', 14, 'FontWeight', 'bold');

% Panel C: Circle image
subplot(2,3,3);
[~, CircleIm] = smi_vis.GenerateImages.circleImage(SMD, [1,0,0], 20);
imshow(CircleIm);
title('C', 'FontSize', 14, 'FontWeight', 'bold');

% Panel D: Precision distribution
subplot(2,3,4);
precision_nm = SMD.X_SE * SMD.PixelSize * 1000;
histogram(precision_nm, 30);
xlabel('Precision (nm)', 'FontSize', 10);
ylabel('Count', 'FontSize', 10);
title('D', 'FontSize', 14, 'FontWeight', 'bold');
box off;

% Panel E: Photon distribution
subplot(2,3,5);
histogram(SMD.Photons, 30);
xlabel('Photons', 'FontSize', 10);
ylabel('Count', 'FontSize', 10);
title('E', 'FontSize', 14, 'FontWeight', 'bold');
box off;

% Panel F: Localizations per frame
subplot(2,3,6);
histogram(SMD.FrameNum, 50);
xlabel('Frame', 'FontSize', 10);
ylabel('Localizations', 'FontSize', 10);
title('F', 'FontSize', 14, 'FontWeight', 'bold');
box off;

% Export
exportgraphics(fig, 'Figure1.png', 'Resolution', 300);
exportgraphics(fig, 'Figure1.pdf', 'ContentType', 'vector');
```

## Complete Visualization Workflow

```matlab
%% Complete Visualization Example
% Demonstrates all major visualization types

% Load data
load('Results.mat', 'SMD', 'SMF');

fprintf('=== Generating Visualizations ===\n');
fprintf('Total localizations: %d\n', length(SMD.X));

%% 1. Super-Resolution Images
fprintf('\n1. Creating SR images...\n');

% Gaussian rendering (publication quality)
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 10);
figure('Name', 'Gaussian SR'); imshow(GaussIm);

% Histogram rendering (fast preview)
[HistIm, HistRGB] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');
figure('Name', 'Histogram SR'); imshow(HistRGB);

% Circle image (show precision)
[~, CircleIm] = smi_vis.GenerateImages.circleImage(SMD, [1,0,0], 20);
figure('Name', 'Circle Image'); imshow(CircleIm);

%% 2. Drift Visualization
fprintf('2. Creating drift image...\n');

[DriftIm, DriftRGB] = smi_vis.GenerateImages.driftImage(SMD, 20);
figure('Name', 'Drift Image'); imshow(DriftRGB);
title('Drift Visualization (Color = Time)');

%% 3. Quality Plots
fprintf('3. Creating quality assessment plots...\n');

figure('Name', 'Quality Assessment', 'Position', [100, 100, 1200, 800]);

% Precision
subplot(2,3,1);
precision_nm = SMD.X_SE * SMD.PixelSize * 1000;
histogram(precision_nm, 50);
xlabel('Precision (nm)'); ylabel('Count');
title(sprintf('Median: %.1f nm', median(precision_nm)));

% Photons
subplot(2,3,2);
histogram(SMD.Photons, 50);
xlabel('Photons'); ylabel('Count');
title(sprintf('Mean: %.0f', mean(SMD.Photons)));

% Background
subplot(2,3,3);
histogram(SMD.Bg, 50);
xlabel('Background (photons/pixel)'); ylabel('Count');
title(sprintf('Mean: %.1f', mean(SMD.Bg)));

% Spatial distribution
subplot(2,3,4);
plot(SMD.X, SMD.Y, '.', 'MarkerSize', 1);
axis equal; axis tight;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Spatial Distribution');

% Frame distribution
subplot(2,3,5);
histogram(SMD.FrameNum, 50);
xlabel('Frame'); ylabel('Localizations');
title('Temporal Distribution');

% Precision vs Photons
subplot(2,3,6);
scatter(SMD.Photons, precision_nm, 5, 'filled', 'MarkerFaceAlpha', 0.2);
xlabel('Photons'); ylabel('Precision (nm)');
title('Precision vs Photons');
grid on;

%% 4. Custom Scale Bars
fprintf('4. Adding scale bars...\n');

% Different locations
figure('Position', [100, 100, 1000, 800]);
locations = {'bottomright', 'bottomleft', 'topright', 'topleft'};
for i = 1:4
    subplot(2,2,i);
    ImWithBar = smi_vis.GenerateImages.scalebar(GaussIm, ...
        SMD.PixelSize/20, 5, locations{i});
    imshow(ImWithBar);
    title(locations{i});
end

%% 5. Export Publication Figure
fprintf('5. Exporting publication figure...\n');

fig = figure('Color', 'white', 'Units', 'inches', 'Position', [0, 0, 7, 5]);

subplot(1,2,1);
imshow(GaussIm);
title('A. Gaussian Rendering', 'FontSize', 12, 'FontWeight', 'bold');

subplot(1,2,2);
precision_nm = SMD.X_SE * SMD.PixelSize * 1000;
histogram(precision_nm, 50, 'FaceColor', [0.3, 0.6, 0.9]);
xlabel('Localization Precision (nm)', 'FontSize', 10);
ylabel('Count', 'FontSize', 10);
title('B. Precision Distribution', 'FontSize', 12, 'FontWeight', 'bold');
box off; grid on;

% Export
exportgraphics(fig, 'Publication_Figure.png', 'Resolution', 300);
fprintf('Saved: Publication_Figure.png\n');

fprintf('\nVisualization complete!\n');
```

## Tips and Best Practices

### Choosing Visualization Methods

**For presentations:**
- Use Gaussian rendering (smooth, professional)
- High contrast colormaps (hot, jet)
- Scale bars in all images
- Large fonts

**For publications:**
- Gaussian or circle rendering
- Perceptually uniform colormaps (viridis, parula)
- Vector graphics when possible (PDF, EPS)
- Multi-panel figures with labels

**For analysis:**
- Histogram images (fast)
- Scatter plots with color coding
- Distribution histograms
- Drift images to check quality

### Optimizing Image Quality

**Resolution guidelines:**

```matlab
% Rule of thumb: 5-10 nm effective pixel size
desired_pixel_nm = 5;
SR_zoom = round((SMD.PixelSize * 1000) / desired_pixel_nm);

fprintf('Recommended zoom: %d×\n', SR_zoom);
```

**Memory considerations:**

Large zoom factors create huge images:

```matlab
% Image size calculation
ImageSize_pixels = [SMD.YSize, SMD.XSize] * SR_zoom;
Memory_MB = prod(ImageSize_pixels) * 3 / (1024^2);  % RGB

fprintf('Image size: %d × %d pixels\n', ImageSize_pixels);
fprintf('Memory: %.1f MB\n', Memory_MB);

% If too large, reduce zoom
if Memory_MB > 1000  % 1 GB
    warning('Image very large, consider reducing zoom');
end
```

### Common Issues

**Issue: Images too dark**

```matlab
% Adjust color scaling manually
[HistIm, ~] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');
MinMax = [0, prctile(HistIm(:), 99)];  % Use 99th percentile
RgbScaled = smi_vis.GenerateImages.colorImage(HistIm, hot(256), MinMax);
```

**Issue: Circles too small**

```matlab
% Increase scale factor
SEScaleFactor = 5;  % Make circles 5× larger
[~, CircleIm] = smi_vis.GenerateImages.circleImage(SMD, [1,0,0], 20, [], SEScaleFactor);
```

**Issue: Scale bar wrong size**

```matlab
% Verify pixel size
fprintf('Pixel size: %.3f μm\n', SMD.PixelSize);
fprintf('Effective SR pixel: %.1f nm\n', SMD.PixelSize*1000/SR_zoom);

% Adjust scale bar
ScalebarLength_um = 1;  % Try smaller bar
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, ScalebarLength_um);
```

## See Also

- [How to Localize Molecules](localize-molecules.md) - Generate data for visualization
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete analysis pipeline
- [SPT Workflow](../workflows/spt-tracking.md) - Trajectory visualization
- [Basic Localization Example](../examples/basic-localization.md) - Working example with visualization
- MATLAB/+smi_vis/@GenerateImages/ - Source code for all methods
- MATLAB/+smi_vis/@GenerateMovies/ - Movie generation tools
