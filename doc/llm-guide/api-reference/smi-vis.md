---
title: "API Reference: +smi_vis Namespace"
category: "api-reference"
level: "beginner"
tags: ["api", "smi_vis", "visualization", "images", "plots", "movies", "super-resolution"]
prerequisites: ["../core-concepts/smd-structure.md", "../core-concepts/tr-structure.md", "../how-to/localize-molecules.md"]
related: ["../how-to/visualize-results.md", "../workflows/smlm-analysis.md", "../workflows/spt-tracking.md"]
summary: "Complete API reference for the +smi_vis namespace covering image generation, movie creation, and visualization tools for SMLM and SPT data"
estimated_time: "25 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# API Reference: +smi_vis Namespace

## Purpose

The `+smi_vis` namespace provides comprehensive visualization tools for single molecule localization and tracking data. This reference documents all methods for generating publication-quality super-resolution images, diagnostic plots, multi-channel overlays, trajectory movies, and interactive data exploration. These tools transform raw localization coordinates into interpretable visual representations suitable for analysis, presentations, and publications.

## Prerequisites

- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Understanding of [TR structure](../core-concepts/tr-structure.md) for trajectory visualization
- Completed [localization analysis](../how-to/localize-molecules.md)
- Basic MATLAB image processing knowledge

## Overview

The `+smi_vis` namespace contains three main classes:

- **`GenerateImages`**: Static methods for creating super-resolution images, histograms, and overlays from localization data
- **`GenerateMovies`**: Class for creating trajectory movies with raw data underlays
- **`InspectResults`**: Interactive GUI for exploring localization results

All visualization methods accept SMD or TR structures as primary inputs and produce images, figures, or movies suitable for analysis and publication.

---

## smi_vis.GenerateImages

### Description

`GenerateImages` is a collection of static methods for generating various types of super-resolution images and visualizations. All methods operate on SMD structures (or compatible subsets) and produce images at specified magnifications. This is the primary class for creating publication figures.

### Class Definition

```matlab
classdef GenerateImages
    methods(Static)
        % Image generation methods
    end
end
```

All methods are static and called without instantiation:
```matlab
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 10);
```

---

### Core Image Generation Methods

#### gaussianImage

**Purpose:** Creates smooth super-resolution image by rendering each localization as a 2D Gaussian blob.

**Signature:**
```matlab
GaussIm = gaussianImage(SMD, SRImageZoom, ScalebarLength)
```

**Parameters:**
- `SMD` (struct, required): SMD structure with required fields:
  - `X`, `Y`: Localization coordinates (pixels)
  - `X_SE`, `Y_SE`: Standard errors (pixels) - used as Gaussian widths
  - `XSize`, `YSize`: Raw image dimensions (pixels)
  - `PixelSize`: Pixel size (micrometers)
  - `Photons`: Photon counts (optional, set to 1 for uniform intensity)
  - `Bg`: Background (optional, set to 0)
  - `FrameNum`: Frame numbers (optional)
- `SRImageZoom` (numeric, optional): Magnification factor relative to raw pixels (Default: 10)
- `ScalebarLength` (numeric, optional): Scale bar length in micrometers (Default: 10, set to 0 to disable)

**Returns:**
- `GaussIm` (array): RGB super-resolution image scaled to [0, 1] with hot colormap

**Algorithm:**
1. Scales coordinates and standard errors by `SRImageZoom`
2. Renders each localization as 2D Gaussian with width = standard error
3. All localizations have uniform photon count (equal weighting)
4. Scales final image to [0, 1]
5. Applies hot colormap (black → red → yellow → white)
6. Adds scale bar if requested

**Example:**
```matlab
% Load results
load('Results.mat', 'SMD');

% Generate Gaussian SR image with 20× zoom
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 10);

% Display
figure; imshow(GaussIm);
title('Gaussian Super-Resolution Image (20× zoom, 10 μm scale bar)');

% Without scale bar
GaussIm_NoBar = smi_vis.GenerateImages.gaussianImage(SMD, 20, 0);
```

**Best Practices:**
- Choose zoom to achieve 5-10 nm effective pixel size:
  ```matlab
  effective_pixel_nm = SMD.PixelSize * 1000 / SRImageZoom;
  % For 108 nm pixels: zoom=20 → 5.4 nm, zoom=10 → 10.8 nm
  ```
- Higher zoom = smoother appearance but larger memory usage
- Standard errors determine blur amount (better precision = sharper)
- Use for publication figures and final presentations

**See Also:** `histogramImage`, `grayscaleImage`

---

#### histogramImage

**Purpose:** Creates pixelated super-resolution image by binning localizations into pixels.

**Signature:**
```matlab
[HistIm, RgbHistIm] = histogramImage(SMD, SRImageZoom, ColorMap)
```

**Parameters:**
- `SMD` (struct, required): SMD structure with fields `X`, `Y`, `XSize`, `YSize`
- `SRImageZoom` (numeric, optional): Magnification factor (Default: 10)
- `ColorMap` (string, optional): MATLAB colormap name (Default: 'hot')
  - Options: 'hot', 'jet', 'parula', 'gray', 'viridis', 'turbo', etc.

**Returns:**
- `HistIm` (array): Grayscale histogram image (counts per pixel)
- `RgbHistIm` (array): RGB image with applied colormap

**Algorithm:**
1. Bins localization coordinates into super-resolution pixels
2. Each pixel value = count of localizations
3. Grayscale output: raw counts
4. RGB output: counts mapped to specified colormap

**Example:**
```matlab
% Generate histogram with hot colormap
[HistGray, HistRGB] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');

% Display side-by-side
figure;
subplot(1,2,1);
imshow(HistGray, []); colormap hot; colorbar;
title('Histogram (Grayscale)');

subplot(1,2,2);
imshow(HistRGB);
title('Histogram (Hot Colormap)');

% Compare different colormaps
figure;
cmaps = {'hot', 'jet', 'parula', 'gray'};
for ii = 1:4
    [~, Im] = smi_vis.GenerateImages.histogramImage(SMD, 20, cmaps{ii});
    subplot(2,2,ii);
    imshow(Im);
    title(cmaps{ii});
end
```

**Comparison with Gaussian:**
- **Histogram**: Fast, exact counts, pixelated appearance
- **Gaussian**: Smooth, publication quality, precision-weighted

**Best Practices:**
- Use for quick previews and dense data
- Good for quantitative analysis (exact counts)
- Faster than Gaussian rendering
- Pixelation visible at lower zoom factors

---

#### circleImage

**Purpose:** Visualizes each localization as a circle with radius proportional to localization precision.

**Signature:**
```matlab
[CircleImage, CircleImageRGB, SRImageZoom] = circleImage(SMD, ColorMap, ...
    SRImageZoom, MinPixelsPerCircle, SEScaleFactor)
```

**Parameters:**
- `SMD` (struct, required): SMD structure with fields:
  - `X`, `Y`: Coordinates (pixels)
  - `X_SE`, `Y_SE`: Standard errors (pixels) - determine circle radii
  - `XSize`, `YSize`: Image dimensions (pixels)
- `ColorMap` (array, optional): Color specification (Default: [1, 0, 0] = red)
  - Single color: N×3 array with one row (RGB values 0-1)
  - Per-localization colors: N×3 array where N = number of localizations
- `SRImageZoom` (numeric, optional): Magnification factor (Default: 20)
  - Set to `[]` to auto-calculate based on `MinPixelsPerCircle`
- `MinPixelsPerCircle` (numeric, optional): Minimum circle size in pixels (Default: 16)
  - Only used if `SRImageZoom` is empty
- `SEScaleFactor` (numeric, optional): Multiplicative scaling of circle radii (Default: 1)

**Returns:**
- `CircleImage` (array): Binary image (0 or 1) with circle outlines
- `CircleImageRGB` (array): RGB color image with colored circles
- `SRImageZoom` (numeric): Actual zoom used (may differ if requested zoom too large)

**Algorithm:**
1. Circle radius = mean(X_SE, Y_SE) × SEScaleFactor × SRImageZoom
2. Each circle drawn with ~8π×radius points for smoothness
3. Binary image: all circles white (1) on black (0)
4. RGB image: circles colored according to `ColorMap`

**Example:**
```matlab
% Basic red circles
[CircBinary, CircRGB] = smi_vis.GenerateImages.circleImage(SMD, [1,0,0], 20);
figure; imshow(CircRGB);

% Color by photon count
Photons_norm = (SMD.Photons - min(SMD.Photons)) / ...
    (max(SMD.Photons) - min(SMD.Photons));
ColorMap = [Photons_norm, zeros(size(Photons_norm)), 1-Photons_norm];
[~, CircRGB] = smi_vis.GenerateImages.circleImage(SMD, ColorMap, 20);
figure; imshow(CircRGB);
title('Circles Colored by Photon Count (Red=bright, Blue=dim)');

% Larger circles for visibility
[~, CircRGB] = smi_vis.GenerateImages.circleImage(SMD, [1,0,0], 20, [], 3);
figure; imshow(CircRGB);
title('Circles Scaled 3× Larger');

% Auto-calculate zoom for minimum circle size
[~, CircRGB, ActualZoom] = smi_vis.GenerateImages.circleImage(...
    SMD, [1,0,0], [], 16, 1);
fprintf('Auto-selected zoom: %d\n', ActualZoom);
```

**Interpreting Circle Sizes:**
- Larger circles = higher uncertainty (worse precision)
- Smaller circles = better precision
- Circle radius directly represents localization standard error

**Best Practices:**
- Use to visualize localization quality spatially
- Combine with color coding for additional information
- Increase `SEScaleFactor` if circles too small to see
- Use `MinPixelsPerCircle` for automatic sizing

**See Also:** `circleDriftImage`

---

#### driftImage

**Purpose:** Creates super-resolution image color-coded by time to visualize drift.

**Signature:**
```matlab
[DriftIm, DriftImRGB] = driftImage(SMD, SRImageZoom)
```

**Parameters:**
- `SMD` (struct, required): SMD structure with fields:
  - `X`, `Y`: Coordinates (pixels)
  - `FrameNum`: Frame numbers (required for temporal coloring)
  - `X_SE`, `Y_SE`: Standard errors (pixels)
  - `XSize`, `YSize`: Image dimensions (pixels)
- `SRImageZoom` (numeric, optional): Magnification factor (Default: 10)

**Returns:**
- `DriftIm` (array): Grayscale drift image
- `DriftImRGB` (array): RGB drift image with jet colormap (blue→cyan→green→yellow→red)

**Color Encoding:**
- **Blue**: Early frames
- **Cyan**: Early-middle frames
- **Green**: Middle frames
- **Yellow**: Middle-late frames
- **Red**: Late frames

**Example:**
```matlab
% Generate drift image
[DriftGray, DriftRGB] = smi_vis.GenerateImages.driftImage(SMD, 20);

figure; imshow(DriftRGB);
title('Drift Visualization (Color = Time)');
colormap jet; colorbar;

% Interpretation
% - Single color region: No drift or stable structure
% - Rainbow streaks: Uncorrected drift
% - Color gradients: Temporal evolution or movement
```

**Interpreting Results:**
- **Rainbow streaks**: Indicates uncorrected stage drift
- **Uniform color**: Good drift correction or stable structure
- **Color patches**: Different structures imaged at different times
- **Radial color patterns**: Possible drift correction artifacts

**Use Cases:**
- Verify drift correction quality
- Identify temporal artifacts
- Visualize dynamic processes
- Quality control for long acquisitions

**See Also:** `circleDriftImage`, `smi_core.DriftCorrection.plotDriftCorrection`

---

#### circleDriftImage

**Purpose:** Combines circle rendering with temporal color coding.

**Signature:**
```matlab
[CircleDriftImage, SRImageZoom] = circleDriftImage(SMD, SRImageZoom, ...
    MinPixelsPerCircle, SEScaleFactor)
```

**Parameters:**
- `SMD` (struct, required): SMD structure with `X`, `Y`, `X_SE`, `Y_SE`, `FrameNum`
- `SRImageZoom` (numeric, optional): Magnification factor (Default: 20)
- `MinPixelsPerCircle` (numeric, optional): Minimum circle size (Default: 16)
- `SEScaleFactor` (numeric, optional): Circle size scaling (Default: 1)

**Returns:**
- `CircleDriftImage` (array): RGB image with time-colored circles
- `SRImageZoom` (numeric): Actual zoom factor used

**Example:**
```matlab
% Generate circle drift image
CircDriftIm = smi_vis.GenerateImages.circleDriftImage(SMD, 20);

figure; imshow(CircDriftIm);
title('Circle Drift Image: Size=Precision, Color=Time');
colormap jet; colorbar;
```

**Interpretation:**
- Circle size = localization precision
- Circle color = frame number
- Combines quality assessment with drift visualization

---

### Color and Overlay Methods

#### colorImage

**Purpose:** Applies colormap to grayscale image with optional intensity scaling.

**Signature:**
```matlab
RGBimage = colorImage(Image, ColorMap, MinMax)
```

**Parameters:**
- `Image` (array, required): Grayscale image
- `ColorMap` (array, optional): N×3 colormap array (Default: hot(256))
- `MinMax` (array, optional): [Min, Max] intensity range for scaling (Default: auto from image)

**Returns:**
- `RGBimage` (array): RGB image with applied colormap

**Example:**
```matlab
% Generate histogram
[HistIm, ~] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'gray');

% Apply custom colormap
CustomMap = [linspace(0,1,256)', zeros(256,1), linspace(1,0,256)'];
RgbCustom = smi_vis.GenerateImages.colorImage(HistIm, CustomMap);

figure; imshow(RgbCustom);
title('Custom Blue-to-Red Colormap');

% Manual intensity scaling
MinMax = [0, 10];  % Only 0-10 localizations per pixel
RgbScaled = smi_vis.GenerateImages.colorImage(HistIm, hot(256), MinMax);

figure; imshow(RgbScaled);
title('Manual Color Scaling [0, 10]');
```

**Use Cases:**
- Apply consistent colormaps across multiple images
- Custom color schemes for specific needs
- Manual intensity scaling for better contrast

---

#### overlayNImages

**Purpose:** Combines 2-4 grayscale images into multi-color overlay.

**Signature:**
```matlab
[OverlayImage, ColorOrderTag] = overlayNImages(ImageStack)
```

**Parameters:**
- `ImageStack` (array, required): M×N×K array where K = 2, 3, or 4 images

**Returns:**
- `OverlayImage` (array): M×N×3 RGB overlay image
- `ColorOrderTag` (string): Color assignment code

**Color Schemes:**
- **2 images**: Green (G) + Magenta (M)
- **3 images**: Green (G) + Magenta (M) + Cyan (C)
- **4 images**: Green (G) + Magenta (M) + Cyan (C) + Yellow (Y)

**Example:**
```matlab
% Load two-channel data
load('Channel1_Results.mat', 'SMD'); SMD1 = SMD;
load('Channel2_Results.mat', 'SMD'); SMD2 = SMD;

% Generate SR images
Im1 = smi_vis.GenerateImages.gaussianImage(SMD1, 20, 0);
Im2 = smi_vis.GenerateImages.gaussianImage(SMD2, 20, 0);

% Convert to grayscale
Im1_gray = mat2gray(Im1);
Im2_gray = mat2gray(Im2);

% Create overlay
ImageStack = cat(3, Im1_gray, Im2_gray);
[OverlayIm, ColorOrder] = smi_vis.GenerateImages.overlayNImages(ImageStack);

figure; imshow(OverlayIm);
title(sprintf('Two-Channel Overlay (%s)', ColorOrder));
% ColorOrder = 'GM' means: channel 1 = green, channel 2 = magenta

% Three-channel overlay
load('Channel3_Results.mat', 'SMD'); SMD3 = SMD;
Im3 = smi_vis.GenerateImages.gaussianImage(SMD3, 20, 0);
ImageStack3 = cat(3, mat2gray(Im1), mat2gray(Im2), mat2gray(Im3));
[Overlay3, Colors3] = smi_vis.GenerateImages.overlayNImages(ImageStack3);

figure; imshow(Overlay3);
title(sprintf('Three-Channel Overlay (%s)', Colors3));
% ColorOrder = 'GMC' means: 1=green, 2=magenta, 3=cyan
```

**Best Practices:**
- Images should be normalized to [0, 1] or similar range
- All images must be same size
- Each color channel scaled independently
- White = all channels overlapping

**See Also:** `rgbImage`

---

#### rgbImage

**Purpose:** Combines separate R, G, B channel images into RGB image.

**Signature:**
```matlab
RGBimage = rgbImage(R, G, B)
```

**Parameters:**
- `R` (array, required): Red channel (grayscale)
- `G` (array, required): Green channel (grayscale)
- `B` (array, required): Blue channel (grayscale)

**Returns:**
- `RGBimage` (array): M×N×3 RGB image

**Example:**
```matlab
% Create custom RGB from three channels
R = HistIm1;  % Red channel
G = HistIm2;  % Green channel
B = HistIm3;  % Blue channel

RGBIm = smi_vis.GenerateImages.rgbImage(R, G, B);
figure; imshow(RGBIm);
```

---

### Annotation Methods

#### scalebar

**Purpose:** Adds scale bar to image.

**Signature:**
```matlab
[ImageOut, Image] = scalebar(Image, PixelSize, Length, Location)
```

**Parameters:**
- `Image` (array, required): Input image (grayscale or RGB)
- `PixelSize` (numeric, required): Physical size of one pixel (micrometers)
- `Length` (numeric, required): Desired scale bar length (micrometers)
- `Location` (string, optional): Scale bar position (Default: 'bottomright')
  - Options: 'bottomright', 'bottomleft', 'topright', 'topleft'

**Returns:**
- `ImageOut` (array): Image with added scale bar
- `Image` (array): Original image (unchanged)

**Example:**
```matlab
% Add scale bar to Gaussian image
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 0);
PixelSize_SR = SMD.PixelSize / 20;  % Effective SR pixel size

ImWithBar = smi_vis.GenerateImages.scalebar(GaussIm, PixelSize_SR, 5, ...
    'bottomright');

figure; imshow(ImWithBar);
title('Image with 5 μm Scale Bar');

% Different positions
figure;
locations = {'bottomright', 'bottomleft', 'topright', 'topleft'};
for ii = 1:4
    subplot(2,2,ii);
    ImWithBar = smi_vis.GenerateImages.scalebar(GaussIm, PixelSize_SR, ...
        5, locations{ii});
    imshow(ImWithBar);
    title(locations{ii});
end
```

**Scale Bar Appearance:**
- White horizontal line
- Positioned near specified corner
- Length accurate to within one pixel
- Works with grayscale and RGB images

---

### Diagnostic and Helper Methods

#### plotHistogram

**Purpose:** Creates histogram plot of a field from SMD structure.

**Signature:**
```matlab
FigHandle = plotHistogram(Vector_in, Hist_Name)
```

**Parameters:**
- `Vector_in` (array, required): Data vector to histogram
- `Hist_Name` (string, required): Histogram title

**Returns:**
- `FigHandle` (handle): Figure handle

**Example:**
```matlab
% Histogram of photon counts
FigHandle = smi_vis.GenerateImages.plotHistogram(SMD.Photons, 'Photon Count');

% Histogram of precision
precision_nm = SMD.X_SE * SMD.PixelSize * 1000;
FigHandle = smi_vis.GenerateImages.plotHistogram(precision_nm, ...
    'Localization Precision (nm)');
```

---

#### blobColorOverlay

**Purpose:** Creates color overlay of fitted emitters (green) on raw data (red).

**Signature:**
```matlab
OverlayImage = blobColorOverlay(Sequence, SMD)
```

**Parameters:**
- `Sequence` (array, required): Raw image stack (Y×X×Frames)
- `SMD` (struct, required): SMD structure with localization coordinates

**Returns:**
- `OverlayImage` (array): RGB overlay image

**Example:**
```matlab
% Load raw data and results
load('RawData.mat', 'sequence');
load('Results.mat', 'SMD');

% Create overlay
OverlayIm = smi_vis.GenerateImages.blobColorOverlay(sequence, SMD);

figure; imshow(OverlayIm);
title('Raw Data (red) + Localizations (green)');
```

**Use Cases:**
- Verify localization accuracy
- Check for missed detections
- Quality control visualization

---

#### showim / dispIm

**Purpose:** Display images with flexible controls.

**Signature:**
```matlab
showim(Image)
dispIm()  % Opens GUI
```

**Parameters:**
- `Image` (array, optional): Image to display

**Description:**
- `showim`: Simple display function
- `dispIm`: Opens interactive GUI for exploring multiple images

**Example:**
```matlab
% Simple display
smi_vis.GenerateImages.showim(GaussIm);

% Interactive GUI
smi_vis.GenerateImages.dispIm();
```

---

### Complete Workflow Example

```matlab
%% Complete Visualization Workflow
% Load results
load('SMLM_Results.mat', 'SMD', 'SMF');

fprintf('=== Generating All Visualization Types ===\n');
fprintf('Total localizations: %d\n', length(SMD.X));

%% 1. Super-Resolution Images
% Gaussian rendering (publication quality)
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 10);
figure('Name', 'Gaussian SR');
imshow(GaussIm);
saveas(gcf, 'GaussianImage.png');

% Histogram rendering (fast preview)
[HistGray, HistRGB] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');
figure('Name', 'Histogram SR');
imshow(HistRGB);
saveas(gcf, 'HistogramImage.png');

% Circle image (precision visualization)
[~, CircleIm] = smi_vis.GenerateImages.circleImage(SMD, [1,0,0], 20);
figure('Name', 'Circle Image');
imshow(CircleIm);
saveas(gcf, 'CircleImage.png');

%% 2. Drift Visualization
[~, DriftRGB] = smi_vis.GenerateImages.driftImage(SMD, 20);
figure('Name', 'Drift Image');
imshow(DriftRGB);
colormap jet; colorbar;
title('Temporal Drift (Blue=early, Red=late)');
saveas(gcf, 'DriftImage.png');

%% 3. Multi-panel Publication Figure
fig = figure('Color', 'white', 'Units', 'inches', ...
    'Position', [0, 0, 7, 5]);

% Panel A: Gaussian rendering
subplot(2,3,1);
imshow(GaussIm);
title('A. Gaussian Rendering', 'FontSize', 12, 'FontWeight', 'bold');

% Panel B: Histogram rendering
subplot(2,3,2);
imshow(HistRGB);
title('B. Histogram Rendering', 'FontSize', 12, 'FontWeight', 'bold');

% Panel C: Circle image
subplot(2,3,3);
imshow(CircleIm);
title('C. Precision Circles', 'FontSize', 12, 'FontWeight', 'bold');

% Panel D: Photon distribution
subplot(2,3,4);
histogram(SMD.Photons, 50, 'FaceColor', [0.3, 0.6, 0.9]);
xlabel('Photons', 'FontSize', 10);
ylabel('Count', 'FontSize', 10);
title('D. Photon Distribution', 'FontSize', 12, 'FontWeight', 'bold');
box off;

% Panel E: Precision distribution
subplot(2,3,5);
precision_nm = SMD.X_SE * SMD.PixelSize * 1000;
histogram(precision_nm, 50, 'FaceColor', [0.9, 0.4, 0.3]);
xlabel('Precision (nm)', 'FontSize', 10);
ylabel('Count', 'FontSize', 10);
title('E. Precision Distribution', 'FontSize', 12, 'FontWeight', 'bold');
box off;

% Panel F: Drift trajectory
subplot(2,3,6);
if isfield(SMD, 'DriftX') && ~isempty(SMD.DriftX)
    plot(SMD.DriftX(:,1), 'b-', 'LineWidth', 1.5); hold on;
    plot(SMD.DriftY(:,1), 'r-', 'LineWidth', 1.5);
    xlabel('Frame', 'FontSize', 10);
    ylabel('Drift (pixels)', 'FontSize', 10);
    legend('X', 'Y', 'Location', 'best');
    title('F. Drift Correction', 'FontSize', 12, 'FontWeight', 'bold');
    box off; grid on;
else
    text(0.5, 0.5, 'No drift data', 'HorizontalAlignment', 'center');
    axis off;
end

% Export
exportgraphics(fig, 'Publication_Figure.png', 'Resolution', 300);
exportgraphics(fig, 'Publication_Figure.pdf', 'ContentType', 'vector');

fprintf('\nVisualization complete! Figures saved.\n');
```

---

## smi_vis.GenerateMovies

### Description

`GenerateMovies` creates movies of trajectory data overlaid on raw microscopy images. Designed primarily for single particle tracking (SPT) visualization, it shows particle motion over time with optional trails, markers, and annotations. Movies can be played interactively or saved to video files.

### Class Definition

```matlab
classdef GenerateMovies < handle
```

**Key Concept:** This is a handle class - properties are modified in place. Create an instance, configure properties, then call `generateMovie()` or `gui()`.

---

### Constructor

```matlab
GM = smi_vis.GenerateMovies()
GM = smi_vis.GenerateMovies(MovieParams)
```

**Parameters:**
- `MovieParams` (struct, optional): Parameter structure (see `prepDefaults()`)

**Returns:**
- `GM` (object): GenerateMovies instance

**Example:**
```matlab
% Create with defaults
GM = smi_vis.GenerateMovies();

% Create with custom parameters
Params = smi_vis.GenerateMovies.prepDefaults();
Params.FrameRate = 20;  % 20 fps
GM = smi_vis.GenerateMovies(Params);
```

---

### Properties

#### Data Properties

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `RawData` | array | [] | Raw data (Y×X×3×Frames) or (Y×X×Frames) |
| `SMD` | struct | empty | SMD structure with points to mark |
| `TR` | struct array | empty | Tracking Results (array of SMD, one per trajectory) |
| `SMF` | struct | empty | SMF structure (for pixel size, framerate) |

#### Parameter Properties

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `Params` | struct | see prepDefaults() | All movie parameters in one structure |

**Key `Params` Fields** (see `prepDefaults()` for complete list):
- `UnitFlag` (boolean): false for pixels/frames, true for μm/seconds (Default: false)
- `MaxTrajLength` (numeric): Max trajectory points shown (Default: inf)
- `XPixels`, `YPixels` (array): ROI for movie ([Min, Max] pixels) (Default: full image)
- `ZFrames` (array): Frame range ([MinFrame, MaxFrame]) (Default: all frames)
- `FrameRate` (numeric): Playback framerate (fps) (Default: 10)
- `PlotMarker` (string): Marker for trajectory points ('o', 'x', 'none') (Default: 'none')
- `SMDMarker` (string): Marker for SMD points (Default: 'none')
- `AddTimeStamp` (boolean): Add frame/time stamp (Default: false)
- `AutoCrop` (boolean): Auto-crop to trajectory extents (Default: false)
- `Resolution` (numeric): Output movie resolution (dpi) (Default: 0 = screen resolution)

#### Protected Properties

| Property | Type | Description |
|----------|------|-------------|
| `ScaledData` | array | Preprocessed/cropped raw data |
| `DataIsPrepped` | boolean | Flag indicating data is ready |
| `VideoObject` | object | MATLAB VideoWriter for saving |

---

### Methods

#### Primary Workflow Methods

##### generateMovie

**Purpose:** Generates and displays the movie.

**Signature:**
```matlab
GM.generateMovie()
```

**Workflow:**
1. Prepares raw data (scaling, cropping)
2. Initializes figure and axes
3. Loops through frames, creating each frame with trajectories
4. Displays interactively or writes to video file (if `VideoObject` set)

**Example:**
```matlab
% Load data
load('TrackingResults.mat', 'TR', 'SMF');
load('RawData.mat', 'sequence');

% Configure movie
GM = smi_vis.GenerateMovies();
GM.TR = TR;
GM.RawData = sequence;
GM.SMF = SMF;

% Set parameters
GM.Params.FrameRate = 15;  % 15 fps
GM.Params.MaxTrajLength = 20;  % Show 20-frame trails
GM.Params.AddTimeStamp = true;

% Generate and display
GM.generateMovie();
```

---

##### saveMovie

**Purpose:** Saves movie to video file.

**Signature:**
```matlab
GM.saveMovie(SavePath)
```

**Parameters:**
- `SavePath` (string, required): Full path to output video file
  - Extension determines format: `.mp4`, `.avi`, `.mj2`

**Example:**
```matlab
% Generate and save movie
GM = smi_vis.GenerateMovies();
GM.TR = TR;
GM.RawData = sequence;
GM.SMF = SMF;

% Configure
GM.Params.FrameRate = 10;
GM.Params.Resolution = 150;  % 150 dpi

% Save
GM.saveMovie('trajectories.mp4');
% Movie saved and also displayed during generation
```

**Supported Formats:**
- `.mp4`: MPEG-4 (recommended, widely compatible)
- `.avi`: Audio Video Interleave
- `.mj2`: Motion JPEG 2000

---

##### gui

**Purpose:** Opens interactive GUI for movie configuration and generation.

**Signature:**
```matlab
GM.gui()
```

**Features:**
- Load data files (TR, SMD, raw data)
- Configure all parameters via controls
- Preview movie interactively
- Save movie to file
- Adjust crop regions
- Set trajectory colors

**Example:**
```matlab
GM = smi_vis.GenerateMovies();
GM.gui();
% Use GUI to:
% - Load data
% - Set parameters
% - Generate movie
% - Save to file
```

---

#### Data Preparation Methods

##### prepRawData

**Purpose:** Prepares raw data (scaling, cropping) for movie generation.

**Signature:**
```matlab
GM.prepRawData()
```

**Called automatically by:** `generateMovie()`

**Actions:**
- Applies percentile clipping to raw data
- Crops to specified ROI (if set)
- Applies `AutoCrop` (if enabled)
- Converts to RGB if needed
- Updates `GM.ScaledData`

---

#### Static Methods

##### playMovie

**Purpose:** Static method to play movie from prepared data.

**Signature:**
```matlab
smi_vis.GenerateMovies.playMovie(PlotAxes, TR, ScaledData, Params, ...
    SMF, SMD, VideoObject)
```

**Parameters:**
- `PlotAxes` (axes handle): Axes to plot in
- `TR` (struct): Tracking Results
- `ScaledData` (array): Prepared raw data
- `Params` (struct): Movie parameters
- `SMF` (struct): SMF structure
- `SMD` (struct): Optional SMD for markers
- `VideoObject` (object): Optional VideoWriter

**Use Case:**
- Advanced users who want full control over movie generation
- Integration into custom workflows
- Typically called internally by `generateMovie()`

---

##### prepDefaults

**Purpose:** Creates default `Params` structure.

**Signature:**
```matlab
Params = smi_vis.GenerateMovies.prepDefaults()
```

**Returns:**
- `Params` (struct): Complete parameter structure with defaults

**Example:**
```matlab
% Get defaults
Params = smi_vis.GenerateMovies.prepDefaults();

% Modify as needed
Params.FrameRate = 20;
Params.MaxTrajLength = 15;
Params.AddTimeStamp = true;
Params.AutoCrop = true;

% Create movie object with custom params
GM = smi_vis.GenerateMovies(Params);
```

**Key Parameters in Structure:**

**Display Settings:**
- `UnitFlag`: Physical (true) vs camera (false) units
- `MaxTrajLength`: Max trajectory points shown (frames)
- `XPixels`, `YPixels`, `ZFrames`: Crop regions

**Raw Data Processing:**
- `MinScaleIntensity`: Minimum scaling intensity
- `PercentileCeiling`: Upper percentile clip (%)
- `PercentileFloor`: Lower percentile clip (%)

**Trajectory Appearance:**
- `PlotMarker`: Marker type for trajectory points
- `TrajColor`: Color for each trajectory (N×3 array)

**SMD Markers:**
- `SMDMarker`: Marker for SMD points
- `SMDColor`: Color for SMD markers

**Output Settings:**
- `FrameRate`: Playback rate (fps)
- `Resolution`: Output resolution (dpi)
- `AddTimeStamp`: Show frame/time stamp

**Auto-Cropping:**
- `AutoCrop`: Enable auto-cropping to trajectories
- `MinXYRange`: Minimum XY range (pixels)
- `NPadPixels`: Padding pixels around data
- `NPadFrames`: Padding frames around trajectories

---

##### saveRawDataMovie

**Purpose:** Save raw data only (no trajectories) as movie.

**Signature:**
```matlab
smi_vis.GenerateMovies.saveRawDataMovie(RawData, FilePath, Params, FrameRate)
```

**Parameters:**
- `RawData` (array): Raw image stack
- `FilePath` (string): Output file path
- `Params` (struct): Movie parameters
- `FrameRate` (numeric): Framerate (fps)

**Example:**
```matlab
% Save raw data only
Params = smi_vis.GenerateMovies.prepDefaults();
smi_vis.GenerateMovies.saveRawDataMovie(sequence, 'raw_movie.mp4', ...
    Params, 10);
```

---

### Complete Movie Workflow Example

```matlab
%% Complete Movie Generation Workflow
% For single particle tracking results

% Load data
load('SPT_Results.mat', 'TR', 'SMF');
load('RawData.mat', 'sequence');

fprintf('=== Generating Trajectory Movie ===\n');
fprintf('Number of trajectories: %d\n', length(TR));

%% 1. Basic Movie
% Create movie object
GM = smi_vis.GenerateMovies();

% Set data
GM.TR = TR;
GM.RawData = sequence;
GM.SMF = SMF;

% Set basic parameters
GM.Params.FrameRate = 10;          % 10 fps playback
GM.Params.MaxTrajLength = 20;       % 20-frame trails
GM.Params.AddTimeStamp = true;      % Show time

% Generate and display
GM.generateMovie();

%% 2. High-Resolution Movie for Publication
GM2 = smi_vis.GenerateMovies();
GM2.TR = TR;
GM2.RawData = sequence;
GM2.SMF = SMF;

% High quality settings
GM2.Params.FrameRate = 15;          % Smooth playback
GM2.Params.Resolution = 300;        % 300 dpi
GM2.Params.MaxTrajLength = 30;      % Longer trails
GM2.Params.AddTimeStamp = true;
GM2.Params.UnitFlag = true;         % Physical units (μm, s)

% Custom trajectory colors (rainbow)
NTraj = length(TR);
GM2.Params.TrajColor = jet(NTraj);

% Save
GM2.saveMovie('trajectories_publication.mp4');
fprintf('Saved: trajectories_publication.mp4\n');

%% 3. Cropped Movie (Region of Interest)
GM3 = smi_vis.GenerateMovies();
GM3.TR = TR;
GM3.RawData = sequence;
GM3.SMF = SMF;

% Manual crop
GM3.Params.XPixels = [50, 150];     % X range
GM3.Params.YPixels = [50, 150];     % Y range
GM3.Params.ZFrames = [100, 200];    % Frames 100-200

% Generate
GM3.generateMovie();

%% 4. Auto-Cropped Movie
GM4 = smi_vis.GenerateMovies();
GM4.TR = TR;
GM4.RawData = sequence;
GM4.SMF = SMF;

% Auto-crop to trajectory extents
GM4.Params.AutoCrop = true;
GM4.Params.MinXYRange = 30;         % Min 30 pixel range
GM4.Params.NPadPixels = 10;         % 10 pixel padding
GM4.Params.NPadFrames = 5;          % 5 frame padding

% Generate
GM4.generateMovie();

%% 5. Movie with Markers
GM5 = smi_vis.GenerateMovies();
GM5.TR = TR;
GM5.RawData = sequence;
GM5.SMF = SMF;

% Add markers at trajectory points
GM5.Params.PlotMarker = 'o';        % Circle markers
GM5.Params.MaxTrajLength = 10;      % Shorter trails

% Also mark all localizations
load('SPT_Results.mat', 'SMD');
GM5.SMD = SMD;
GM5.Params.SMDMarker = 'x';         % X markers
GM5.Params.SMDColor = [1, 0, 0];    % Red

% Generate
GM5.generateMovie();

%% 6. GUI-Based Movie
% For interactive exploration
GM_gui = smi_vis.GenerateMovies();
GM_gui.gui();
% Use GUI to load data, configure, and save

fprintf('\nMovie generation complete!\n');
```

---

## Best Practices

### Choosing Visualization Methods

**For Publications:**
- Use `gaussianImage()` for smooth, professional appearance
- Export as PDF (vector graphics) when possible
- Add scale bars to all images
- Use perceptually uniform colormaps (parula, viridis)
- High resolution (300 dpi minimum)

**For Presentations:**
- `gaussianImage()` or `histogramImage()` with high contrast colormaps (hot, jet)
- Larger fonts and scale bars
- Movies for dynamic data
- Multi-panel figures with labels

**For Analysis:**
- `histogramImage()` for fast previews
- `driftImage()` for quality control
- `circleImage()` for precision assessment
- Distribution histograms for quantitative analysis

### Image Resolution Guidelines

**Choosing zoom factor:**
```matlab
% Rule of thumb: 5-10 nm effective pixel size
desired_pixel_nm = 5;  % Target 5 nm pixels
SRZoom = round((SMD.PixelSize * 1000) / desired_pixel_nm);

fprintf('Recommended zoom: %d×\n', SRZoom);
fprintf('Effective pixel: %.1f nm\n', SMD.PixelSize*1000/SRZoom);
```

**Memory considerations:**
```matlab
% Estimate memory usage
ImageSize = [SMD.YSize, SMD.XSize] * SRZoom;
Memory_MB = prod(ImageSize) * 3 / (1024^2);  % RGB

fprintf('Image size: %d × %d pixels\n', ImageSize);
fprintf('Memory: %.1f MB\n', Memory_MB);

if Memory_MB > 1000
    warning('Image very large (%.0f MB). Consider reducing zoom.', Memory_MB);
end
```

### Movie Generation Best Practices

**Framerate selection:**
- Playback: 10-15 fps (smooth, easy to follow)
- Publication: 15-30 fps (professional appearance)
- Analysis: 5-10 fps (easier to see details)

**Trail length:**
- Short (5-10 frames): High density data, multiple trajectories
- Medium (15-25 frames): Standard tracking
- Long (30+ frames): Sparse data, long-range motion

**Resolution:**
- Screen preview: Default (0)
- Presentation: 150 dpi
- Publication: 300 dpi

**File formats:**
- MP4: Best compression, widely compatible (recommended)
- AVI: Larger files, universal support
- MJ2: High quality, less common

---

## Common Patterns

### Pattern 1: Complete Visualization Suite

```matlab
% Generate all standard visualizations
SMD = load('Results.mat', 'SMD');

% SR images
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 10);
[~, HistIm] = smi_vis.GenerateImages.histogramImage(SMD, 20, 'hot');
[~, CircleIm] = smi_vis.GenerateImages.circleImage(SMD, [1,0,0], 20);

% Drift
[~, DriftIm] = smi_vis.GenerateImages.driftImage(SMD, 20);

% Save all
imwrite(GaussIm, 'Gaussian.png');
imwrite(HistIm, 'Histogram.png');
imwrite(CircleIm, 'Circles.png');
imwrite(DriftIm, 'Drift.png');
```

### Pattern 2: Multi-Channel Analysis

```matlab
% Load channels
load('Ch1.mat', 'SMD'); SMD1 = SMD;
load('Ch2.mat', 'SMD'); SMD2 = SMD;

% Generate images
Im1 = smi_vis.GenerateImages.gaussianImage(SMD1, 20, 0);
Im2 = smi_vis.GenerateImages.gaussianImage(SMD2, 20, 0);

% Overlay
Stack = cat(3, mat2gray(Im1), mat2gray(Im2));
[Overlay, Colors] = smi_vis.GenerateImages.overlayNImages(Stack);

imwrite(Overlay, 'TwoChannel_Overlay.png');
```

### Pattern 3: Publication Figure with Statistics

```matlab
% Multi-panel figure
fig = figure('Units', 'inches', 'Position', [0,0,7,5], 'Color', 'white');

% SR image
subplot(2,2,1);
GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, 20, 5);
imshow(GaussIm);
title('A', 'FontSize', 14, 'FontWeight', 'bold');

% Precision
subplot(2,2,2);
prec = SMD.X_SE * SMD.PixelSize * 1000;
histogram(prec, 50);
xlabel('Precision (nm)'); ylabel('Count');
title('B', 'FontSize', 14, 'FontWeight', 'bold');

% Photons
subplot(2,2,3);
histogram(SMD.Photons, 50);
xlabel('Photons'); ylabel('Count');
title('C', 'FontSize', 14, 'FontWeight', 'bold');

% Drift
subplot(2,2,4);
[~, DriftIm] = smi_vis.GenerateImages.driftImage(SMD, 20);
imshow(DriftIm);
title('D', 'FontSize', 14, 'FontWeight', 'bold');

% Export
exportgraphics(fig, 'Figure1.pdf', 'ContentType', 'vector');
exportgraphics(fig, 'Figure1.png', 'Resolution', 300);
```

---

## See Also

### Related Documentation
- [How to Visualize Results](../how-to/visualize-results.md) - Comprehensive visualization guide
- [SMD Structure](../core-concepts/smd-structure.md) - Data structure reference
- [TR Structure](../core-concepts/tr-structure.md) - Trajectory structure reference
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete analysis pipeline
- [SPT Workflow](../workflows/spt-tracking.md) - Tracking workflow with movies

### Related Classes
- `smi.SMLM` - SMLM analysis (generates SMD)
- `smi.SPT` - SPT analysis (generates TR)
- `smi_core.DriftCorrection` - Drift visualization methods
- `smi_sim.GaussBlobs` - Underlying Gaussian rendering

---

## Summary

The `+smi_vis` namespace provides comprehensive visualization capabilities:

**GenerateImages** offers static methods for:
- Super-resolution image rendering (Gaussian, histogram, circle)
- Drift visualization (temporal color coding)
- Multi-channel overlays (2-4 channels)
- Custom colormaps and annotations
- Scale bars and quality assessment plots

**GenerateMovies** provides:
- Trajectory movie generation with raw data underlays
- Interactive GUI for movie configuration
- Flexible parameter control
- Multiple output formats
- Auto-cropping and custom ROIs

These tools transform raw localization data into publication-ready visualizations, enabling both qualitative assessment and quantitative analysis of single molecule experiments. All methods are designed for ease of use while providing extensive customization options for advanced users.
