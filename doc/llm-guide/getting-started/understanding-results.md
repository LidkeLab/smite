---
title: "Understanding Your Results"
category: "getting-started"
level: "beginner"
tags: ["results", "interpretation", "quality-control", "smd", "tr", "analysis"]
prerequisites: ["first-analysis.md", "../core-concepts/smd-structure.md"]
related: ["../how-to/visualize-results.md", "../how-to/threshold-results.md"]
summary: "Comprehensive guide for beginners to interpret localization and tracking results, assess quality, and identify common problems"
estimated_time: "20 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# Understanding Your Results

## Purpose

After running SMLM or SPT analysis, you have a wealth of data stored in SMD (Single Molecule Data) or TR (Tracking Results) structures. But what do all these numbers mean? How do you know if your results are good? This guide teaches you how to read and interpret your results, assess localization quality, identify common problems, and determine when your analysis succeeded or needs adjustment.

## Prerequisites

- Completed at least one analysis (see [first-analysis.md](first-analysis.md))
- Basic understanding of [SMD structure](../core-concepts/smd-structure.md)
- Results file (typically Results.mat) to examine

## Overview

Understanding results involves four key activities:

1. **Reading the data** - What each SMD field represents
2. **Assessing quality** - Interpreting precision, photons, and p-values
3. **Visual inspection** - What good and bad results look like
4. **Problem diagnosis** - Recognizing and fixing common issues

By the end of this guide, you'll be able to confidently evaluate whether your analysis produced reliable, publication-quality data.

## Loading and Inspecting Your Results

### Load Your Results File

After running an SMLM or SPT analysis, results are saved in your data directory:

```matlab
% Typical location: FileDir/Results/Results.mat
ResultsFile = 'C:/path/to/your/data/Results/Results.mat';
load(ResultsFile, 'SMD', 'SMF');

fprintf('=== Results Summary ===\n');
fprintf('Total localizations: %d\n', length(SMD.X));
fprintf('Frames analyzed: %d\n', SMD.NFrames);
fprintf('Image size: %d x %d pixels\n', SMD.XSize, SMD.YSize);
fprintf('Pixel size: %.1f nm\n', SMD.PixelSize * 1000);
```

**What you see:**
- SMD contains all localization data
- SMF contains all analysis parameters (what settings you used)
- Number of localizations can range from hundreds to millions depending on your data

### Quick Visual Check

Before diving into numbers, look at your localizations spatially:

```matlab
% Simple scatter plot
figure;
plot(SMD.X, SMD.Y, '.', 'MarkerSize', 2);
axis equal; axis tight;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('%d Localizations', length(SMD.X)));
```

**What to look for:**
- Do you see structure that makes biological sense?
- Are localizations distributed across the field of view or clustered in specific regions?
- Are there obvious artifacts (lines, edges, hot pixels)?

## Understanding SMD Fields

SMD contains many fields, but these are the most important for assessing quality:

### Position Fields

**SMD.X, SMD.Y** (pixels)
- The fitted position of each molecule
- Sub-pixel precision: values like 45.327 pixels are normal
- Origin is (1, 1) at the center of the top-left pixel

```matlab
% Check position range
fprintf('X range: %.1f to %.1f pixels\n', min(SMD.X), max(SMD.X));
fprintf('Y range: %.1f to %.1f pixels\n', min(SMD.Y), max(SMD.Y));

% Convert to physical units (nanometers)
X_nm = SMD.X * SMD.PixelSize * 1000;
Y_nm = SMD.Y * SMD.PixelSize * 1000;
```

**SMD.X_SE, SMD.Y_SE** (pixels)
- Standard error (uncertainty) of the position
- Lower values = more precise localization
- Derived from Cramer-Rao Lower Bound (CRLB)

```matlab
% Check precision distribution
precision_nm = SMD.X_SE * SMD.PixelSize * 1000;
fprintf('Precision (X): %.1f ± %.1f nm (mean ± std)\n', ...
    mean(precision_nm), std(precision_nm));
fprintf('Median precision: %.1f nm\n', median(precision_nm));
```

**Good precision values:**
- Excellent: 5-10 nm
- Good: 10-20 nm
- Acceptable: 20-30 nm
- Poor: >30 nm

### Photometric Fields

**SMD.Photons**
- Total photons detected from each emitter
- More photons = better signal, more precise localization
- Typical range: 100-5000 photons

```matlab
% Check photon statistics
fprintf('Photons: %.0f ± %.0f (mean ± std)\n', ...
    mean(SMD.Photons), std(SMD.Photons));
fprintf('Photon range: %.0f to %.0f\n', ...
    min(SMD.Photons), max(SMD.Photons));

% Visualize distribution
figure;
histogram(SMD.Photons, 50);
xlabel('Detected Photons');
ylabel('Count');
title('Photon Distribution');
```

**Interpreting photon counts:**
- >1000 photons: Excellent signal
- 500-1000 photons: Good signal
- 200-500 photons: Moderate signal
- 100-200 photons: Marginal
- <100 photons: Likely noise or very dim emitters

**SMD.Bg** (photons per pixel)
- Background level at each localization
- Includes autofluorescence, scattered light, camera noise
- Lower is better

```matlab
% Background statistics
fprintf('Background: %.1f ± %.1f photons/pixel\n', ...
    mean(SMD.Bg), std(SMD.Bg));
```

**Typical background:**
- Low background: 5-15 photons/pixel
- Moderate: 15-30 photons/pixel
- High: 30-50 photons/pixel
- Very high: >50 photons/pixel (may indicate problems)

### Quality Metrics

**SMD.PValue**
- P-value from chi-squared goodness-of-fit test
- Range: 0 to 1
- Higher = better fit to Gaussian PSF model

```matlab
% P-value distribution
fprintf('P-value: %.3f (median), %.3f (mean)\n', ...
    median(SMD.PValue), mean(SMD.PValue));
fprintf('P-values < 0.01: %.1f%%\n', ...
    100 * sum(SMD.PValue < 0.01) / length(SMD.PValue));

% Visualize
figure;
histogram(SMD.PValue, 50);
xlabel('P-value');
ylabel('Count');
title('Fit Quality Distribution');
xline(0.01, 'r--', 'LineWidth', 2, 'Label', 'Typical threshold');
```

**Interpreting p-values:**
- P > 0.1: Excellent fit
- P = 0.01-0.1: Good fit
- P < 0.01: Questionable fit (may be overlapping molecules, wrong PSF model)

**Warning:** A few low p-values are normal. If >20% have p-value < 0.01, you may have systematic problems.

### Temporal Fields

**SMD.FrameNum**
- Which frame each localization came from
- Range: 1 to SMD.NFrames

```matlab
% Localizations per frame
figure;
histogram(SMD.FrameNum, 50);
xlabel('Frame Number');
ylabel('Localizations');
title('Temporal Distribution');

% Average density
avg_per_frame = length(SMD.X) / SMD.NFrames;
fprintf('Average localizations per frame: %.1f\n', avg_per_frame);
```

**Interpreting temporal distribution:**
- Should be relatively uniform unless emitters are photobleaching
- Sharp drops indicate photobleaching or focus drift
- Peaks may indicate blinking events or motion artifacts

### Frame Connection Fields (if applicable)

**SMD.ConnectID**
- Links localizations from the same molecule across frames
- Same ID = same molecule
- Only present if frame connection was performed

```matlab
if isfield(SMD, 'ConnectID')
    unique_emitters = unique(SMD.ConnectID);
    unique_emitters = unique_emitters(unique_emitters > 0);

    fprintf('Unique molecules detected: %d\n', length(unique_emitters));
    fprintf('Average appearances per molecule: %.1f frames\n', ...
        length(SMD.X) / length(unique_emitters));

    % Blinking statistics
    appearances = histcounts(SMD.ConnectID, ...
        [unique_emitters; max(unique_emitters)+1]);
    fprintf('Molecules appearing >5 times: %d\n', ...
        sum(appearances > 5));
end
```

## Assessing Localization Quality

### The Quality Triangle: Photons, Precision, and Background

Localization quality depends on three interconnected factors:

```matlab
% Calculate signal-to-background ratio
PSF_area = pi * SMF.Fitting.PSFSigma^2;  % Approximate PSF area
SBR = SMD.Photons ./ (SMD.Bg * PSF_area);

% Relationship analysis
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
scatter(SMD.Photons, precision_nm, 5, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Photons'); ylabel('Precision (nm)');
title('Photons vs Precision');
grid on;
% Expected: inverse relationship (more photons = lower precision value = better)

subplot(1,3,2);
scatter(SMD.Bg, precision_nm, 5, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Background (photons/pixel)'); ylabel('Precision (nm)');
title('Background vs Precision');
grid on;
% Expected: positive relationship (higher background = worse precision)

subplot(1,3,3);
scatter(SBR, precision_nm, 5, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Signal-to-Background Ratio'); ylabel('Precision (nm)');
title('SBR vs Precision');
grid on;
% Expected: inverse relationship (higher SBR = better precision)

fprintf('Median SBR: %.1f\n', median(SBR));
```

**Good data shows:**
- More photons correlate with better (lower) precision
- Higher background correlates with worse (higher) precision
- SBR >10 produces good localizations

### Multi-Panel Quality Report

Create a comprehensive quality assessment:

```matlab
figure('Position', [100, 100, 1400, 900]);

% Panel 1: Spatial distribution
subplot(2,3,1);
plot(SMD.X, SMD.Y, '.', 'MarkerSize', 1);
axis equal; axis tight;
title(sprintf('%d Localizations', length(SMD.X)));
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Panel 2: Precision distribution
subplot(2,3,2);
histogram(precision_nm, 50, 'FaceColor', [0.3, 0.6, 0.9]);
xlabel('Precision (nm)');
ylabel('Count');
title(sprintf('Precision: %.1f nm (median)', median(precision_nm)));
xline(median(precision_nm), 'r--', 'LineWidth', 2);
xlim([0, 50]);  % Focus on reasonable range

% Panel 3: Photon distribution
subplot(2,3,3);
histogram(SMD.Photons, 50, 'FaceColor', [0.9, 0.4, 0.3]);
xlabel('Photons');
ylabel('Count');
title(sprintf('Photons: %.0f (median)', median(SMD.Photons)));
xline(median(SMD.Photons), 'r--', 'LineWidth', 2);

% Panel 4: Background distribution
subplot(2,3,4);
histogram(SMD.Bg, 50, 'FaceColor', [0.3, 0.8, 0.4]);
xlabel('Background (photons/pixel)');
ylabel('Count');
title(sprintf('Background: %.1f (median)', median(SMD.Bg)));
xline(median(SMD.Bg), 'r--', 'LineWidth', 2);

% Panel 5: P-value distribution
subplot(2,3,5);
histogram(SMD.PValue, 50, 'FaceColor', [0.8, 0.3, 0.8]);
xlabel('P-value');
ylabel('Count');
title(sprintf('P-value: %.3f (median)', median(SMD.PValue)));
xline(0.01, 'r--', 'LineWidth', 2, 'Label', 'Threshold');

% Panel 6: Localizations per frame
subplot(2,3,6);
counts = histcounts(SMD.FrameNum, 1:SMD.NFrames+1);
plot(1:SMD.NFrames, counts, 'LineWidth', 1.5);
xlabel('Frame Number');
ylabel('Localizations');
title(sprintf('Avg: %.1f per frame', mean(counts)));
grid on;
```

### Quality Scores

Create an overall quality assessment:

```matlab
fprintf('\n=== QUALITY ASSESSMENT ===\n\n');

% Precision score
median_prec = median(precision_nm);
if median_prec < 15
    prec_score = 'Excellent';
elseif median_prec < 25
    prec_score = 'Good';
elseif median_prec < 35
    prec_score = 'Fair';
else
    prec_score = 'Poor';
end
fprintf('Precision: %.1f nm - %s\n', median_prec, prec_score);

% Photon score
median_photons = median(SMD.Photons);
if median_photons > 1000
    photon_score = 'Excellent';
elseif median_photons > 500
    photon_score = 'Good';
elseif median_photons > 200
    photon_score = 'Fair';
else
    photon_score = 'Poor';
end
fprintf('Photons: %.0f - %s\n', median_photons, photon_score);

% Background score
median_bg = median(SMD.Bg);
if median_bg < 15
    bg_score = 'Excellent (low)';
elseif median_bg < 30
    bg_score = 'Good';
elseif median_bg < 50
    bg_score = 'Fair';
else
    bg_score = 'Poor (high)';
end
fprintf('Background: %.1f photons/pixel - %s\n', median_bg, bg_score);

% P-value score
low_pval_frac = sum(SMD.PValue < 0.01) / length(SMD.PValue);
if low_pval_frac < 0.05
    pval_score = 'Excellent';
elseif low_pval_frac < 0.15
    pval_score = 'Good';
elseif low_pval_frac < 0.30
    pval_score = 'Fair';
else
    pval_score = 'Poor';
end
fprintf('Fit quality: %.1f%% low p-values - %s\n', ...
    low_pval_frac * 100, pval_score);

% Localization density
area_um2 = SMD.XSize * SMD.YSize * SMD.PixelSize^2;
density = length(SMD.X) / area_um2;
fprintf('\nLocalization density: %.1f per square micron\n', density);

% Overall assessment
fprintf('\n=== OVERALL ASSESSMENT ===\n');
if strcmp(prec_score, 'Excellent') && strcmp(photon_score, 'Excellent')
    fprintf('Status: Publication quality\n');
elseif ~strcmp(prec_score, 'Poor') && ~strcmp(photon_score, 'Poor')
    fprintf('Status: Good quality, suitable for analysis\n');
else
    fprintf('Status: Quality issues detected - review settings\n');
end
```

## Visual Quality Assessment

Numbers tell part of the story, but visualizing your results reveals the full picture.

### Super-Resolution Images

Generate a super-resolution image to see if structure is resolved:

```matlab
% Create Gaussian-rendered SR image
SR_zoom = 20;  % 20x magnification
scalebar_um = 5;  % 5 micron scale bar

GaussIm = smi_vis.GenerateImages.gaussianImage(SMD, SR_zoom, scalebar_um);

figure;
imshow(GaussIm);
title('Super-Resolution Image');
```

**What to look for:**

**Good results show:**
- Clear, well-defined structures
- Sharp features at expected biological scales
- Uniform intensity across field
- Smooth, continuous structures

**Problems to watch for:**
- Blurry, indistinct features (poor precision or drift)
- Patchy, discontinuous structures (missing localizations)
- Obvious grid patterns (stage drift artifacts)
- Bright spots or lines (hot pixels, dust, artifacts)

### Drift Visualization

Check for uncorrected drift using temporal color coding:

```matlab
% Create drift-colored image
[DriftIm, DriftImRGB] = smi_vis.GenerateImages.driftImage(SMD, 20);

figure;
imshow(DriftImRGB);
title('Drift Image (Color = Time)');
% Blue = early frames, Red = late frames
```

**Interpreting drift images:**
- **Single color structure:** Good - no visible drift or drift well-corrected
- **Rainbow streaks:** Bad - significant uncorrected drift
- **Color separation:** Stage moved during acquisition

If drift is present, enable or adjust drift correction:

```matlab
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.Method = 'DC-KNN';  % or 'RCC', 'GCC', 'RDC'
```

### Precision Mapping

Visualize precision spatially to identify regions with poor fits:

```matlab
% Color localizations by precision
figure;
scatter(SMD.X, SMD.Y, 10, precision_nm, 'filled');
axis equal; axis tight;
colormap hot; colorbar;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Precision Map (color = precision in nm)');
caxis([0, 40]);  % Set color scale
```

**What good data looks like:**
- Uniform colors across the field (consistent quality)
- Blue/purple colors (low precision values = good)

**Warning signs:**
- Red regions (high precision = poor quality)
- Systematic patterns (edge effects, illumination problems)
- Clusters of poor precision (potential artifacts)

## Understanding Tracking Results (TR)

If you performed SPT analysis, you'll have TR (Tracking Results) instead of or in addition to SMD.

### Basic TR Inspection

```matlab
% Load tracking results
load('Results.mat', 'TR', 'SMD', 'SMF');

% Basic statistics
N_trajectories = length(TR);
fprintf('=== Tracking Results ===\n');
fprintf('Total trajectories: %d\n', N_trajectories);

% Trajectory lengths
lengths = arrayfun(@(x) length(x.X), TR);
fprintf('Trajectory length: %.1f ± %.1f frames (mean ± std)\n', ...
    mean(lengths), std(lengths));
fprintf('Shortest trajectory: %d frames\n', min(lengths));
fprintf('Longest trajectory: %d frames\n', max(lengths));

% Duration statistics
if isfield(TR, 'FrameRate') && TR(1).FrameRate > 0
    durations = lengths / TR(1).FrameRate;
    fprintf('Trajectory duration: %.2f ± %.2f seconds\n', ...
        mean(durations), std(durations));
end
```

### Visualizing Trajectories

```matlab
% Plot all trajectories
figure;
hold on;
for i = 1:min(length(TR), 100)  % Limit to first 100 for clarity
    plot(TR(i).X, TR(i).Y, '-', 'LineWidth', 1.5);
end
axis equal; axis tight;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('%d Trajectories', min(length(TR), 100)));

% Trajectory length distribution
figure;
histogram(lengths, 50);
xlabel('Trajectory Length (frames)');
ylabel('Count');
title('Trajectory Length Distribution');
```

**Good tracking shows:**
- Reasonable trajectory lengths (>10 frames is good)
- Smooth, continuous paths
- Biologically plausible motion

**Problems to watch for:**
- Most trajectories only 2-3 frames (poor tracking)
- Extremely long trajectories (may be multiple molecules incorrectly linked)
- Jumps or discontinuities (linking errors)

## Common Problems and Solutions

### Problem 1: Very Poor Precision (>40 nm median)

**Symptoms:**
- Precision values consistently above 30-40 nm
- Blurry super-resolution images
- Low photon counts

**Causes and solutions:**

```matlab
% Check photon counts
if median(SMD.Photons) < 200
    fprintf('Issue: Low photon counts (%.0f median)\n', median(SMD.Photons));
    fprintf('Solutions:\n');
    fprintf('  - Increase laser power\n');
    fprintf('  - Increase camera exposure time\n');
    fprintf('  - Use brighter fluorophores\n');
end

% Check background
if median(SMD.Bg) > 30
    fprintf('Issue: High background (%.1f median)\n', median(SMD.Bg));
    fprintf('Solutions:\n');
    fprintf('  - Reduce autofluorescence (better buffer, washing)\n');
    fprintf('  - Improve illumination uniformity\n');
    fprintf('  - Reduce ambient light contamination\n');
end
```

### Problem 2: Many Low P-Values

**Symptoms:**
- >20% of localizations have p-value < 0.01
- Fitting struggles with data

**Diagnosis:**

```matlab
low_pval_frac = sum(SMD.PValue < 0.01) / length(SMD.PValue);
fprintf('Low p-value fraction: %.1f%%\n', low_pval_frac * 100);

if low_pval_frac > 0.2
    fprintf('Likely causes:\n');
    fprintf('  1. Overlapping molecules (too high density)\n');
    fprintf('  2. Wrong PSF model\n');
    fprintf('  3. Aberrations or non-Gaussian PSFs\n');

    % Check for overlaps by examining spatial density
    if density > 5  % Localizations per square micron
        fprintf('  -> High density detected: %.1f/um^2\n', density);
        fprintf('     Consider reducing emitter density\n');
    end
end
```

**Solutions:**
- Reduce emitter density (lower dye concentration, faster photobleaching)
- Try different PSF model (astigmatism for 3D, asymmetric Gaussian)
- Enable thresholding to remove poor fits automatically

### Problem 3: Drift Artifacts

**Symptoms:**
- Rainbow streaks in drift visualization
- Structures appear blurred or doubled
- Drift plots show large movements

**Check drift correction:**

```matlab
if isfield(SMD, 'DriftX') && ~isempty(SMD.DriftX)
    total_drift = sqrt(SMD.DriftX(end,1)^2 + SMD.DriftY(end,1)^2);
    drift_nm = total_drift * SMD.PixelSize * 1000;

    fprintf('Total drift: %.1f nm\n', drift_nm);

    if drift_nm > 100
        fprintf('Significant drift detected!\n');
        fprintf('Solutions:\n');
        fprintf('  - Enable drift correction: SMF.DriftCorrection.On = true\n');
        fprintf('  - Try different method: DC-KNN, RCC, or GCC\n');
        fprintf('  - Check microscope stability\n');
    end

    % Visualize drift trajectory
    figure;
    plot(SMD.DriftX(:,1) * SMD.PixelSize * 1000, ...
         SMD.DriftY(:,1) * SMD.PixelSize * 1000, '-o');
    xlabel('Drift X (nm)'); ylabel('Drift Y (nm)');
    title('Stage Drift Trajectory');
    axis equal; grid on;
end
```

### Problem 4: Too Few Localizations

**Symptoms:**
- Sparse super-resolution image
- Much fewer localizations than expected
- Empty regions

**Diagnosis:**

```matlab
fprintf('Localizations per frame: %.1f\n', ...
    length(SMD.X) / SMD.NFrames);

if length(SMD.X) / SMD.NFrames < 10
    fprintf('Very sparse localization density\n');
    fprintf('Check:\n');

    % Check if thresholding removed too many
    if isfield(SMD, 'ThreshFlag')
        rejected = sum(SMD.ThreshFlag > 0);
        if rejected > 0.5 * length(SMD.ThreshFlag)
            fprintf('  -> Thresholding removed %.1f%% of fits\n', ...
                100 * rejected / length(SMD.ThreshFlag));
            fprintf('     Solution: Relax threshold parameters\n');
        end
    end

    % Check detection sensitivity
    fprintf('  -> Detection threshold (MinPhotons): %.0f\n', ...
        SMF.BoxFinding.MinPhotons);
    fprintf('     Try lowering to 100-150 if currently higher\n');
end
```

### Problem 5: Hot Pixels or Artifacts

**Symptoms:**
- Bright spots in same location across many frames
- Lines or patterns in super-resolution image
- Outliers in position distribution

**Detection:**

```matlab
% Find potential hot pixels (same position repeatedly)
[unique_pos, ~, idx] = unique([round(SMD.X), round(SMD.Y)], 'rows');
counts = histcounts(idx, 1:max(idx)+1);
hot_pixel_threshold = SMD.NFrames * 0.1;  % >10% of frames

if any(counts > hot_pixel_threshold)
    fprintf('Warning: Potential hot pixels detected\n');
    fprintf('  %d positions appear in >%.0f frames\n', ...
        sum(counts > hot_pixel_threshold), hot_pixel_threshold);
    fprintf('Solution: Apply median filter or bad pixel mask\n');
end
```

## Complete Quality Check Script

Here's a comprehensive script to run on any results:

```matlab
%% Complete Results Quality Check
% Run this after loading SMD and SMF

fprintf('========================================\n');
fprintf('         SMITE RESULTS QUALITY CHECK    \n');
fprintf('========================================\n\n');

%% 1. Basic Information
fprintf('--- DATASET INFORMATION ---\n');
fprintf('Localizations: %d\n', length(SMD.X));
fprintf('Frames: %d\n', SMD.NFrames);
fprintf('Image size: %d x %d pixels\n', SMD.XSize, SMD.YSize);
fprintf('Pixel size: %.1f nm\n', SMD.PixelSize * 1000);
fprintf('Frame rate: %.1f Hz\n', SMD.FrameRate);
area_um2 = SMD.XSize * SMD.YSize * SMD.PixelSize^2;
density = length(SMD.X) / area_um2;
fprintf('Density: %.1f localizations/um^2\n', density);
fprintf('Average per frame: %.1f\n\n', length(SMD.X) / SMD.NFrames);

%% 2. Precision Assessment
fprintf('--- PRECISION ASSESSMENT ---\n');
precision_nm = SMD.X_SE * SMD.PixelSize * 1000;
fprintf('Median: %.1f nm\n', median(precision_nm));
fprintf('Mean ± Std: %.1f ± %.1f nm\n', mean(precision_nm), std(precision_nm));
fprintf('Range: %.1f to %.1f nm\n', min(precision_nm), max(precision_nm));
if median(precision_nm) < 15
    fprintf('Quality: EXCELLENT\n\n');
elseif median(precision_nm) < 25
    fprintf('Quality: GOOD\n\n');
elseif median(precision_nm) < 35
    fprintf('Quality: FAIR\n\n');
else
    fprintf('Quality: POOR - Consider improving SNR\n\n');
end

%% 3. Photon Statistics
fprintf('--- PHOTON STATISTICS ---\n');
fprintf('Median: %.0f photons\n', median(SMD.Photons));
fprintf('Mean ± Std: %.0f ± %.0f photons\n', ...
    mean(SMD.Photons), std(SMD.Photons));
if median(SMD.Photons) > 500
    fprintf('Quality: GOOD signal strength\n\n');
elseif median(SMD.Photons) > 200
    fprintf('Quality: ADEQUATE signal strength\n\n');
else
    fprintf('Quality: LOW signal - consider brighter conditions\n\n');
end

%% 4. Background Analysis
fprintf('--- BACKGROUND ANALYSIS ---\n');
fprintf('Median: %.1f photons/pixel\n', median(SMD.Bg));
fprintf('Mean ± Std: %.1f ± %.1f photons/pixel\n', ...
    mean(SMD.Bg), std(SMD.Bg));
if median(SMD.Bg) < 20
    fprintf('Quality: LOW background (excellent)\n\n');
elseif median(SMD.Bg) < 40
    fprintf('Quality: MODERATE background\n\n');
else
    fprintf('Quality: HIGH background - may affect precision\n\n');
end

%% 5. Fit Quality
fprintf('--- FIT QUALITY ---\n');
fprintf('Median p-value: %.3f\n', median(SMD.PValue));
low_pval = sum(SMD.PValue < 0.01) / length(SMD.PValue) * 100;
fprintf('Low p-values (<0.01): %.1f%%\n', low_pval);
if low_pval < 10
    fprintf('Quality: GOOD fits\n\n');
elseif low_pval < 25
    fprintf('Quality: ACCEPTABLE fits\n\n');
else
    fprintf('Quality: MANY POOR FITS - check for overlaps\n\n');
end

%% 6. Drift Check
if isfield(SMD, 'DriftX') && ~isempty(SMD.DriftX)
    fprintf('--- DRIFT ANALYSIS ---\n');
    total_drift = sqrt(SMD.DriftX(end,1)^2 + SMD.DriftY(end,1)^2);
    drift_nm = total_drift * SMD.PixelSize * 1000;
    fprintf('Total drift: %.1f nm\n', drift_nm);
    if drift_nm < 50
        fprintf('Drift: MINIMAL\n\n');
    elseif drift_nm < 150
        fprintf('Drift: MODERATE (corrected)\n\n');
    else
        fprintf('Drift: SIGNIFICANT - verify correction worked\n\n');
    end
end

%% 7. Overall Recommendation
fprintf('========================================\n');
fprintf('--- OVERALL ASSESSMENT ---\n');
if median(precision_nm) < 20 && median(SMD.Photons) > 300 && low_pval < 15
    fprintf('Status: PUBLICATION QUALITY\n');
    fprintf('These results are suitable for publication.\n');
elseif median(precision_nm) < 30 && median(SMD.Photons) > 200
    fprintf('Status: GOOD QUALITY\n');
    fprintf('These results are suitable for analysis.\n');
else
    fprintf('Status: NEEDS IMPROVEMENT\n');
    fprintf('Consider optimizing acquisition or analysis parameters.\n');
end
fprintf('========================================\n');
```

## Next Steps

After understanding your results:

1. **Good results:** Generate publication figures with [visualize-results.md](../how-to/visualize-results.md)
2. **Quality issues:** Adjust parameters and reanalyze
3. **Need filtering:** Apply quality thresholds with [threshold-results.md](../how-to/threshold-results.md)
4. **Ready for analysis:** Explore clustering, tracking, or advanced analysis

## See Also

- [How to Visualize Results](../how-to/visualize-results.md) - Create publication-quality figures
- [How to Apply Quality Filters](../how-to/threshold-results.md) - Remove poor quality localizations
- [SMD Structure Reference](../core-concepts/smd-structure.md) - Complete field documentation
- [SMLM Analysis Workflow](../workflows/smlm-analysis.md) - Full pipeline details
- [SPT Tracking Workflow](../workflows/spt-tracking.md) - Trajectory analysis
