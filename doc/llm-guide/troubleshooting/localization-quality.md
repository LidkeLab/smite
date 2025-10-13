---
title: "Diagnosing and Fixing Poor Localization Quality"
category: "troubleshooting"
level: "intermediate"
tags: ["troubleshooting", "quality", "precision", "detection", "psf", "drift", "fitting", "parameters"]
prerequisites: ["../getting-started/understanding-results.md", "../how-to/tune-parameters.md", "../how-to/localize-molecules.md"]
related: ["../how-to/threshold-results.md", "../workflows/smlm-analysis.md", "../core-concepts/smf-structure.md"]
summary: "Comprehensive troubleshooting guide for diagnosing and fixing common localization quality issues including poor precision, low detection rates, false positives, bad fit quality, wrong PSF models, and drift artifacts"
estimated_time: "30 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# Diagnosing and Fixing Poor Localization Quality

## Purpose

When localization results don't meet expectations, systematic diagnosis and correction are essential. This guide helps you identify the root causes of poor localization quality and provides concrete solutions for six major problem categories: poor precision, low detection rate, excessive false positives, poor fit quality, wrong PSF model, and drift artifacts. Each problem section covers symptoms, diagnosis methods, root causes, and step-by-step solutions.

## Prerequisites

- Understanding of [result interpretation](../getting-started/understanding-results.md)
- Familiarity with [parameter tuning](../how-to/tune-parameters.md)
- Basic knowledge of the [localization process](../how-to/localize-molecules.md)
- Results file (SMD/SMF) showing quality issues

## Overview

Quality issues typically manifest in specific patterns:

1. **Poor Precision** (>30 nm): Localization uncertainties too large for meaningful analysis
2. **Low Detection Rate**: Missing emitters that should be visible in your data
3. **False Positives**: Detecting noise or artifacts as real molecules
4. **Bad P-values**: Poor fit quality suggesting model mismatch or overlapping emitters
5. **Wrong PSF Model**: Systematic biases from incorrect point spread function assumptions
6. **Drift Artifacts**: Blurred or doubled structures from uncorrected stage movement

Understanding which problem affects your data is the first step to solving it. Many issues are interrelated, so fixing one often improves others.

## Problem 1: Poor Precision (>30 nm)

### Symptoms

Poor localization precision manifests as:

```matlab
% Load your results
load('Results.mat', 'SMD', 'SMF');
precision_nm = SMD.X_SE * SMF.Data.PixelSize * 1000;

fprintf('Median precision: %.1f nm\n', median(precision_nm));
fprintf('Mean precision: %.1f nm\n', mean(precision_nm));
fprintf('Precision range: %.1f to %.1f nm\n', min(precision_nm), max(precision_nm));
```

**Key indicators:**
- Median precision >30 nm (target: <20 nm for high-quality SMLM)
- Wide precision distribution (large standard deviation)
- Super-resolution images appear blurry or lacking fine detail
- Difficulty resolving structures separated by <100 nm

**Visual check:**
```matlab
% Precision histogram
figure;
histogram(precision_nm, 50);
xlabel('Localization Precision (nm)');
ylabel('Count');
title('Precision Distribution');
xline(20, 'r--', 'LineWidth', 2, 'Label', 'Target');
xline(30, 'r-', 'LineWidth', 2, 'Label', 'Poor threshold');

% Spatial precision map
figure;
scatter(SMD.X, SMD.Y, 10, precision_nm, 'filled');
colormap hot; colorbar;
title('Precision Map (hot colors = poor precision)');
axis equal; axis tight;
caxis([0, 50]);
```

### Diagnosis

Precision depends on three fundamental factors described by the Cramer-Rao Lower Bound:

**1. Check photon counts** - More photons enable more precise localization:

```matlab
fprintf('Median photons: %.0f\n', median(SMD.Photons));
fprintf('Mean photons: %.0f\n', mean(SMD.Photons));
fprintf('Photon range: %.0f to %.0f\n', min(SMD.Photons), max(SMD.Photons));

% Analyze photon-precision relationship
figure;
scatter(SMD.Photons, precision_nm, 5, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Detected Photons'); ylabel('Precision (nm)');
title('Photons vs Precision (should show inverse relationship)');
grid on;
```

**Problem indicators:**
- Median photons <200: Very low, expect poor precision
- Median photons 200-500: Marginal, precision will suffer
- No clear inverse relationship: Suggests other factors dominate

**2. Check background levels** - High background adds noise:

```matlab
fprintf('Median background: %.1f photons/pixel\n', median(SMD.Bg));
fprintf('Mean background: %.1f photons/pixel\n', mean(SMD.Bg));

% Calculate signal-to-background ratio
PSF_area = pi * SMF.Fitting.PSFSigma^2;
SBR = SMD.Photons ./ (SMD.Bg * PSF_area);
fprintf('Median SBR: %.1f\n', median(SBR));

% Target SBR >5 for decent precision, >10 for excellent
if median(SBR) < 5
    fprintf('WARNING: Low signal-to-background ratio\n');
end

% Background-precision relationship
figure;
scatter(SMD.Bg, precision_nm, 5, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Background (photons/pixel)'); ylabel('Precision (nm)');
title('Background vs Precision (should show positive correlation)');
grid on;
```

**Problem indicators:**
- Median background >30-40 photons/pixel: High, will degrade precision
- SBR <5: Too low for good localization
- Spatially varying background: Indicates illumination or sample issues

**3. Check PSF parameters** - Wrong PSF model degrades precision estimates:

```matlab
% If PSFSigma was fitted (FitType includes 'S')
if isfield(SMD, 'PSFSigma')
    fprintf('Fitted PSFSigma: %.3f +/- %.3f pixels (mean +/- std)\n', ...
        mean(SMD.PSFSigma), std(SMD.PSFSigma));
    fprintf('SMF PSFSigma setting: %.3f pixels\n', SMF.Fitting.PSFSigma);

    % Check for mismatch
    sigma_bias = mean(SMD.PSFSigma) - SMF.Fitting.PSFSigma;
    if abs(sigma_bias) > 0.15
        fprintf('WARNING: PSFSigma mismatch of %.3f pixels detected\n', sigma_bias);
    end
end

% Pixel size check - affects precision in physical units
fprintf('Pixel size: %.1f nm\n', SMF.Data.PixelSize * 1000);
% Typical: 100-160 nm for high magnification (100x-150x objectives)
```

### Solutions

**Solution 1: Increase signal (photon counts)**

If median photons <500:

```matlab
% Option A: Increase laser power (if not already saturated)
% Check your microscope settings
fprintf('Current exposure time: %.0f ms\n', 1000 / SMF.Data.FrameRate);
% Increase laser power by 20-50%% or increase exposure time

% Option B: Use brighter fluorophores
% Consider switching to brighter dyes (e.g., Alexa647 for dSTORM)
% For DNA-PAINT: Increase imager strand concentration

% Option C: Improve collection efficiency
% Check objective NA, dichroic/filter transmission
% Clean optical path components
```

**Test the effect:**
```matlab
% Theoretical precision improvement from increasing photons
current_photons = median(SMD.Photons);
target_photons = 1000;  % Aim for this
improvement_factor = sqrt(target_photons / current_photons);

fprintf('Expected precision improvement: %.1fx\n', improvement_factor);
fprintf('Current: %.1f nm -> Predicted: %.1f nm\n', ...
    median(precision_nm), median(precision_nm) / improvement_factor);
```

**Solution 2: Reduce background**

If median background >30 photons/pixel or SBR <5:

```matlab
% Option A: Improve sample preparation
% - Better washing to remove unbound fluorophores
% - Reduce autofluorescence (optimize buffer, use quenchers)
% - Minimize coverslip fluorescence (use high-quality glass)

% Option B: Optimize imaging conditions
% - Reduce ambient light leakage
% - Check for scattered excitation light
% - Use proper emission filters
% - Minimize exposure time while maintaining signal

% Option C: Computational background correction
% Enable background subtraction in SMF
SMF.Fitting.FitType = 'XYNB';  % Fit background per localization (default)

% Or pre-subtract background before localization
if false  % Set to true to test
    % This modifies raw data, use carefully
    for i = 1:size(rawData, 3)
        frame = rawData(:,:,i);
        bg_estimate = median(frame(:));  % Simple median background
        rawData(:,:,i) = max(frame - bg_estimate, 0);
    end
end
```

**Solution 3: Verify and correct PSFSigma**

Critical parameter that affects all estimates:

```matlab
% Step 1: Empirically determine correct PSFSigma
SMF_test = SMF;
SMF_test.Fitting.FitType = 'XYNBS';  % Fit sigma as free parameter
SMF_test.Fitting.PSFSigma = 1.5;     % Initial guess only

% Re-run localization on subset
LD = smi_core.LocalizeData();
[rawData, ~, SMF_test] = LD.loadRawData(SMF_test, 1);
subset = rawData(:,:,1:min(100, size(rawData,3)));
LD = smi_core.LocalizeData(subset, SMF_test);
SMD_test = LD.genLocalizations();

% Extract fitted sigma
fitted_sigma = median(SMD_test.PSFSigma);
fprintf('Optimal PSFSigma: %.3f pixels\n', fitted_sigma);

% Step 2: Visualize distribution to check for aberrations
figure;
histogram(SMD_test.PSFSigma, 50);
xlabel('Fitted PSF Sigma (pixels)');
ylabel('Count');
title(sprintf('PSF Sigma Distribution (median: %.3f)', fitted_sigma));
xline(fitted_sigma, 'r--', 'LineWidth', 2);

% Check for spatial variation (aberrations)
figure;
scatter(SMD_test.X, SMD_test.Y, 10, SMD_test.PSFSigma, 'filled');
colormap hot; colorbar;
title('Spatial PSF Sigma Variation');
axis equal; axis tight;

% Step 3: Update SMF with correct value
SMF.Fitting.PSFSigma = fitted_sigma;
SMF.Fitting.FitType = 'XYNB';  % Fix sigma for speed

% Set reasonable thresholds
SMF.Thresholding.MinPSFSigma = fitted_sigma * 0.7;
SMF.Thresholding.MaxPSFSigma = fitted_sigma * 1.3;

% Step 4: Re-run full analysis
fprintf('Reanalyzing with corrected PSFSigma...\n');
LD = smi_core.LocalizeData(rawData, SMF);
SMD_corrected = LD.genLocalizations();

% Check improvement
precision_nm_new = SMD_corrected.X_SE * SMF.Data.PixelSize * 1000;
fprintf('Old median precision: %.1f nm\n', median(precision_nm));
fprintf('New median precision: %.1f nm\n', median(precision_nm_new));
fprintf('Improvement: %.1f%%\n', ...
    100 * (median(precision_nm) - median(precision_nm_new)) / median(precision_nm));
```

**Solution 4: Optimize box size**

BoxSize affects photon collection and background inclusion:

```matlab
% Current box size
fprintf('Current BoxSize: %d pixels\n', SMF.BoxFinding.BoxSize);

% Rule of thumb: BoxSize should be 4-5x PSFSigma
recommended_size = ceil(4.5 * SMF.Fitting.PSFSigma);
fprintf('Recommended BoxSize: %d pixels\n', recommended_size);

if abs(SMF.BoxFinding.BoxSize - recommended_size) > 1
    % Test different box sizes
    box_sizes = [recommended_size-2, recommended_size, recommended_size+2];

    for i = 1:length(box_sizes)
        SMF_test = SMF;
        SMF_test.BoxFinding.BoxSize = box_sizes(i);

        LD = smi_core.LocalizeData(subset, SMF_test);
        SMD_test = LD.genLocalizations();

        test_precision = median(SMD_test.X_SE) * SMF.Data.PixelSize * 1000;
        fprintf('BoxSize=%d: %.1f nm precision, %.0f photons\n', ...
            box_sizes(i), test_precision, median(SMD_test.Photons));
    end

    % Use size that minimizes precision
    SMF.BoxFinding.BoxSize = recommended_size;
end
```

### Expected Outcomes

After implementing solutions:
- Median precision <20 nm (excellent) or 20-30 nm (good)
- Precision distribution narrower (smaller std)
- Clear inverse relationship between photons and precision
- Super-resolution images show finer structural detail
- Improved SBR (>10 for best results)

## Problem 2: Low Detection Rate

### Symptoms

Missing expected localizations:

```matlab
% Check detection statistics
avg_per_frame = length(SMD.X) / SMD.NFrames;
fprintf('Localizations per frame: %.1f\n', avg_per_frame);
fprintf('Total localizations: %d\n', length(SMD.X));

% Temporal distribution
figure;
histogram(SMD.FrameNum, 50);
xlabel('Frame Number'); ylabel('Localizations');
title('Temporal Detection Distribution');

% If ground truth available (from simulation)
if exist('SMD_true', 'var')
    detection_rate = 100 * length(SMD.X) / length(SMD_true.X);
    fprintf('Detection rate: %.1f%%\n', detection_rate);
    if detection_rate < 85
        fprintf('WARNING: Low detection rate (<85%%)\n');
    end
end
```

**Key indicators:**
- Sparse super-resolution images with obvious gaps
- Fewer localizations than expected from sample labeling density
- Structures appear incomplete or fragmented
- Localizations per frame much lower than experimental design suggests
- Temporal drops suggesting missed emitters

**Visual check:**
```matlab
% Compare to raw data
LD = smi_core.LoadData();
[rawData, ~, ~] = LD.loadRawData(SMF, 1);
typical_frame = 50;

figure;
subplot(1,2,1);
imagesc(rawData(:,:,typical_frame));
colormap gray; axis image;
title('Raw Frame');
colorbar;

subplot(1,2,2);
imagesc(rawData(:,:,typical_frame));
hold on;
same_frame = (SMD.FrameNum == typical_frame);
plot(SMD.X(same_frame), SMD.Y(same_frame), 'r+', 'MarkerSize', 10, 'LineWidth', 1);
axis image;
title(sprintf('Frame %d: %d detections', typical_frame, sum(same_frame)));
% Are obvious bright spots being missed?
```

### Diagnosis

**1. Check detection threshold** - Too high MinPhotons rejects dim emitters:

```matlab
fprintf('Current MinPhotons threshold: %.0f\n', SMF.BoxFinding.MinPhotons);
fprintf('Detected photons: %.0f (median), %.0f (mean)\n', ...
    median(SMD.Photons), mean(SMD.Photons));
fprintf('Photon range: %.0f to %.0f\n', min(SMD.Photons), max(SMD.Photons));

% The threshold should be well below median detected photons
threshold_ratio = SMF.BoxFinding.MinPhotons / median(SMD.Photons);
fprintf('Threshold / Median ratio: %.2f\n', threshold_ratio);

if threshold_ratio > 0.7
    fprintf('WARNING: Detection threshold too high (>70%% of median)\n');
end
```

**2. Check post-fit filtering** - Thresholding may reject too many:

```matlab
if SMF.Thresholding.On
    fprintf('\nThresholding is ON:\n');
    fprintf('  MinPhotons: %.0f\n', SMF.Thresholding.MinPhotons);
    fprintf('  MaxXY_SE: %.3f pixels (%.1f nm)\n', ...
        SMF.Thresholding.MaxXY_SE, SMF.Thresholding.MaxXY_SE * SMF.Data.PixelSize * 1000);
    fprintf('  MinPValue: %.3f\n', SMF.Thresholding.MinPValue);

    % Check rejection rates
    if isfield(SMD, 'ThreshFlag')
        rejected = sum(SMD.ThreshFlag > 0);
        total = length(SMD.ThreshFlag);
        rejection_rate = 100 * rejected / total;
        fprintf('  Rejection rate: %.1f%% (%d/%d)\n', rejection_rate, rejected, total);

        if rejection_rate > 40
            fprintf('WARNING: >40%% of fits rejected by thresholding\n');
        end
    end
end
```

**3. Check camera calibration** - Wrong gain/offset causes photon miscalculation:

```matlab
fprintf('\nCamera calibration:\n');
fprintf('  CameraGain: %.3f\n', SMF.Data.CameraGain);
fprintf('  CameraOffset: %.0f\n', SMF.Data.CameraOffset);
fprintf('  CameraType: %s\n', SMF.Data.CameraType);

% Check if raw data values make sense
sample_frame = rawData(:,:,typical_frame);
fprintf('  Raw data range: %.0f to %.0f ADU\n', min(sample_frame(:)), max(sample_frame(:)));
fprintf('  Raw data median: %.0f ADU\n', median(sample_frame(:)));

% After offset/gain correction, photon values should be reasonable
% Background should be ~5-30 photons/pixel, peaks ~100-5000 photons
```

**4. Check box size** - Too small misses extended PSFs:

```matlab
fprintf('\nBoxFinding parameters:\n');
fprintf('  BoxSize: %d pixels\n', SMF.BoxFinding.BoxSize);
fprintf('  PSFSigma: %.3f pixels\n', SMF.Fitting.PSFSigma);

ratio = SMF.BoxFinding.BoxSize / SMF.Fitting.PSFSigma;
fprintf('  BoxSize/PSFSigma ratio: %.1f\n', ratio);

if ratio < 4.0
    fprintf('WARNING: BoxSize may be too small (ratio <4.0)\n');
end
```

### Solutions

**Solution 1: Lower detection threshold**

If MinPhotons is too high:

```matlab
% Current threshold
current_threshold = SMF.BoxFinding.MinPhotons;

% Recommended: 50-70%% of expected photon count
% For first pass, try 60%% of current median
if ~isempty(SMD.Photons)
    recommended_threshold = 0.6 * median(SMD.Photons);
else
    % If no localizations yet, start with conservative value
    recommended_threshold = 150;
end

fprintf('Recommended MinPhotons: %.0f (current: %.0f)\n', ...
    recommended_threshold, current_threshold);

% Test on subset
SMF_test = SMF;
SMF_test.BoxFinding.MinPhotons = recommended_threshold;

LD = smi_core.LocalizeData(subset, SMF_test);
SMD_test = LD.genLocalizations();

fprintf('Before: %d detections\n', length(SMD.X));
fprintf('After: %d detections (%.1f%% increase)\n', ...
    length(SMD_test.X), 100 * (length(SMD_test.X) - length(SMD.X)) / length(SMD.X));

% If improvement is good, update SMF
if length(SMD_test.X) > length(SMD.X)
    SMF.BoxFinding.MinPhotons = recommended_threshold;
    fprintf('Updated MinPhotons to %.0f\n', recommended_threshold);
end
```

**Solution 2: Relax post-fit thresholds**

If thresholding rejects too many:

```matlab
if SMF.Thresholding.On && rejection_rate > 30
    fprintf('Relaxing threshold parameters...\n');

    % Option A: Disable thresholding temporarily to see maximum detections
    SMF_test = SMF;
    SMF_test.Thresholding.On = false;

    LD = smi_core.LocalizeData(subset, SMF_test);
    SMD_test = LD.genLocalizations();

    fprintf('Without thresholding: %d localizations\n', length(SMD_test.X));
    fprintf('With thresholding: %d localizations\n', length(SMD.X));

    % Option B: Adjust thresholds to reject only worst ~10-20%%
    precision_threshold_90th = prctile(SMD_test.X_SE, 90);
    photon_threshold_10th = prctile(SMD_test.Photons, 10);

    SMF.Thresholding.MaxXY_SE = precision_threshold_90th;
    SMF.Thresholding.MinPhotons = photon_threshold_10th;

    fprintf('Updated thresholds:\n');
    fprintf('  MaxXY_SE: %.3f pixels (90th percentile)\n', SMF.Thresholding.MaxXY_SE);
    fprintf('  MinPhotons: %.0f (10th percentile)\n', SMF.Thresholding.MinPhotons);
end
```

**Solution 3: Correct camera calibration**

If camera parameters are wrong:

```matlab
% Empirically determine camera gain and offset
% Method: Compare noise in raw data to expected Poisson statistics

% Take a uniform region (e.g., background)
background_region = rawData(20:40, 20:40, 1:min(100, size(rawData,3)));
background_values = background_region(:);

measured_mean = mean(background_values);
measured_var = var(background_values);

fprintf('\nBackground statistics:\n');
fprintf('  Mean: %.1f ADU\n', measured_mean);
fprintf('  Variance: %.1f ADU^2\n', measured_var);
fprintf('  Ratio (Var/Mean): %.2f\n', measured_var / measured_mean);

% For Poisson process in photons: Var = Mean
% With gain G and offset O: ADU = G*photons + O
% Var(ADU) = G^2 * Var(photons) = G^2 * photons = G * (ADU - O)
% Therefore: G = Var(ADU) / (Mean(ADU) - O)

% If offset is known (from dark frames):
% estimated_gain = measured_var / (measured_mean - SMF.Data.CameraOffset);

% If offset unknown, this is tricky - consult camera documentation
% Typical values:
%   EMCCD: Gain = 1-50, Offset = 0-2000
%   sCMOS: Gain = 0.3-2, Offset = 100-500

fprintf('\nCurrent calibration:\n');
fprintf('  Gain: %.3f\n', SMF.Data.CameraGain);
fprintf('  Offset: %.0f\n', SMF.Data.CameraOffset);
fprintf('\nIf these seem wrong, measure empirically or consult camera manual.\n');
```

**Solution 4: Increase box size**

If BoxSize is too small:

```matlab
% Calculate optimal box size
optimal_box = ceil(4.5 * SMF.Fitting.PSFSigma);
fprintf('Optimal BoxSize: %d (current: %d)\n', optimal_box, SMF.BoxFinding.BoxSize);

if SMF.BoxFinding.BoxSize < optimal_box
    SMF_test = SMF;
    SMF_test.BoxFinding.BoxSize = optimal_box;

    LD = smi_core.LocalizeData(subset, SMF_test);
    SMD_test = LD.genLocalizations();

    fprintf('BoxSize %d: %d detections, %.0f median photons\n', ...
        SMF.BoxFinding.BoxSize, length(SMD.X), median(SMD.Photons));
    fprintf('BoxSize %d: %d detections, %.0f median photons\n', ...
        optimal_box, length(SMD_test.X), median(SMD_test.Photons));

    if length(SMD_test.X) > length(SMD.X) * 1.1  % >10%% improvement
        SMF.BoxFinding.BoxSize = optimal_box;
        fprintf('Updated BoxSize to %d\n', optimal_box);
    end
end
```

### Expected Outcomes

After implementing solutions:
- Detection rate >90% (if ground truth available)
- More complete structures in super-resolution images
- Uniform detection across frames (no systematic gaps)
- Localizations per frame matches experimental design
- Biological structures appear continuous rather than fragmented

## Problem 3: Many False Positives

### Symptoms

Detecting noise or artifacts as molecules:

```matlab
% Check for excess detections
fprintf('Total localizations: %d\n', length(SMD.X));
fprintf('Localizations per frame: %.1f\n', length(SMD.X) / SMD.NFrames);

% With ground truth
if exist('SMD_true', 'var')
    false_positive_count = length(SMD.X) - length(SMD_true.X);
    false_positive_rate = 100 * false_positive_count / length(SMD.X);
    fprintf('Estimated false positives: %d (%.1f%% of detections)\n', ...
        false_positive_count, false_positive_rate);

    if false_positive_rate > 10
        fprintf('WARNING: High false positive rate (>10%%)\n');
    end
end
```

**Key indicators:**
- Many localizations with very low photon counts
- Random scattered localizations not forming biological structures
- "Sparkly" appearance in super-resolution image (noise pepper)
- Localizations in regions that should be empty
- Wide distribution of photon counts extending to very low values

**Visual check:**
```matlab
% Look at spatial distribution
figure;
subplot(1,2,1);
plot(SMD.X, SMD.Y, '.', 'MarkerSize', 1);
axis equal; axis tight;
title('All Localizations');

% Color by photon count
subplot(1,2,2);
scatter(SMD.X, SMD.Y, 5, SMD.Photons, 'filled');
colormap hot; colorbar;
title('Photon Count (low values = likely false positives)');
axis equal; axis tight;
caxis([0, 500]);

% Check photon distribution for tail at low values
figure;
histogram(SMD.Photons, 50);
xlabel('Photons'); ylabel('Count');
title('Photon Distribution (check for low-photon tail)');
xline(SMF.BoxFinding.MinPhotons, 'r--', 'LineWidth', 2, 'Label', 'Detection threshold');
```

### Diagnosis

**1. Check detection stringency** - MinPhotons too low admits noise:

```matlab
fprintf('Detection threshold: %.0f photons\n', SMF.BoxFinding.MinPhotons);
fprintf('Detected photon range: %.0f to %.0f\n', min(SMD.Photons), max(SMD.Photons));

% Fraction near threshold suggests false positives
near_threshold = sum(SMD.Photons < SMF.BoxFinding.MinPhotons * 1.2);
fprintf('Within 20%% of threshold: %d (%.1f%%)\n', ...
    near_threshold, 100 * near_threshold / length(SMD.Photons));

if 100 * near_threshold / length(SMD.Photons) > 30
    fprintf('WARNING: Many detections near threshold (likely noise)\n');
end
```

**2. Check background level** - High background increases false positives:

```matlab
fprintf('Background: %.1f +/- %.1f photons/pixel (mean +/- std)\n', ...
    mean(SMD.Bg), std(SMD.Bg));

% High background generates noise peaks
if mean(SMD.Bg) > 40
    fprintf('WARNING: High background level (>40 photons/pixel)\n');
    fprintf('  This increases probability of false positive detections\n');
end
```

**3. Check p-values** - Poor fits suggest noise detections:

```matlab
fprintf('P-value distribution:\n');
fprintf('  Median: %.3f\n', median(SMD.PValue));
fprintf('  Mean: %.3f\n', mean(SMD.PValue));
fprintf('  < 0.01: %.1f%%\n', 100 * sum(SMD.PValue < 0.01) / length(SMD.PValue));

% Many very low p-values (e.g., <0.001) suggest fitting noise
very_low_pval = sum(SMD.PValue < 0.001) / length(SMD.PValue);
if very_low_pval > 0.1
    fprintf('WARNING: %.1f%% have p-value <0.001 (likely fitting noise)\n', ...
        100 * very_low_pval);
end
```

**4. Check spatial randomness** - False positives lack structure:

```matlab
% Pair correlation analysis (if available)
% False positives show no spatial correlation

% Simple check: Compare to ground truth if available
if exist('SMD_true', 'var')
    % Find unmatched localizations (potential false positives)
    max_dist = 2;  % pixels
    matched_idx = false(length(SMD.X), 1);

    for i = 1:length(SMD_true.X)
        dist = sqrt((SMD.X - SMD_true.X(i)).^2 + (SMD.Y - SMD_true.Y(i)).^2);
        matched_idx = matched_idx | (dist < max_dist & SMD.FrameNum == SMD_true.FrameNum(i));
    end

    false_positives = ~matched_idx;
    fprintf('Potential false positives: %d (%.1f%%)\n', ...
        sum(false_positives), 100 * sum(false_positives) / length(SMD.X));

    % Characteristics of false positives
    fprintf('False positive characteristics:\n');
    fprintf('  Photons: %.0f (vs %.0f for true)\n', ...
        mean(SMD.Photons(false_positives)), mean(SMD.Photons(~false_positives)));
    fprintf('  Background: %.1f (vs %.1f for true)\n', ...
        mean(SMD.Bg(false_positives)), mean(SMD.Bg(~false_positives)));
end
```

### Solutions

**Solution 1: Increase detection threshold**

Raise MinPhotons to exclude noise:

```matlab
% Analyze photon distribution to find separation between noise and signal
figure;
histogram(SMD.Photons, 100);
xlabel('Photons'); ylabel('Count');
title('Find valley between noise and signal peaks');
set(gca, 'YScale', 'log');  % Log scale helps see distribution

% Strategy: Set threshold above noise peak
% Typically: Identify lowest mode, set threshold 2-3x higher

photon_percentiles = prctile(SMD.Photons, [10, 20, 30, 40, 50]);
fprintf('Photon percentiles: 10th=%.0f, 20th=%.0f, 30th=%.0f, 40th=%.0f, 50th=%.0f\n', ...
    photon_percentiles);

% Conservative approach: Set threshold at 30-40th percentile
% This removes lowest 30-40%% which are likely noise
new_threshold = photon_percentiles(3);  % 30th percentile
fprintf('Suggested MinPhotons: %.0f (current: %.0f)\n', ...
    new_threshold, SMF.BoxFinding.MinPhotons);

if new_threshold > SMF.BoxFinding.MinPhotons
    SMF_test = SMF;
    SMF_test.BoxFinding.MinPhotons = new_threshold;

    LD = smi_core.LocalizeData(subset, SMF_test);
    SMD_test = LD.genLocalizations();

    fprintf('Before: %d localizations\n', length(SMD.X));
    fprintf('After: %d localizations (%.1f%% reduction)\n', ...
        length(SMD_test.X), 100 * (length(SMD.X) - length(SMD_test.X)) / length(SMD.X));

    % Verify false positives reduced without losing true signal
    % Check that super-resolution image still shows expected structures

    SMF.BoxFinding.MinPhotons = new_threshold;
end
```

**Solution 2: Enable or strengthen thresholding**

Post-fit quality filtering removes noise:

```matlab
% Enable thresholding if not already on
if ~SMF.Thresholding.On
    fprintf('Enabling thresholding...\n');
    SMF.Thresholding.On = true;
end

% Set aggressive thresholds to remove noise
% Target: Remove worst ~20-30%% by multiple criteria

% Photon threshold (post-fit, can be more stringent than detection)
SMF.Thresholding.MinPhotons = prctile(SMD.Photons, 25);  % Remove lowest 25%%

% Precision threshold
SMF.Thresholding.MaxXY_SE = prctile(SMD.X_SE, 75);  % Remove worst 25%%

% P-value threshold (standard)
SMF.Thresholding.MinPValue = 0.01;

% PSF sigma bounds (if fitted)
if isfield(SMD, 'PSFSigma')
    median_sigma = median(SMD.PSFSigma);
    SMF.Thresholding.MinPSFSigma = median_sigma * 0.6;
    SMF.Thresholding.MaxPSFSigma = median_sigma * 1.4;
end

% Alternative: Use log-likelihood auto-thresholding
SMF.Thresholding.AutoThreshLogL = true;
% This adaptively removes outliers based on fit quality

fprintf('Updated thresholding parameters:\n');
fprintf('  MinPhotons: %.0f\n', SMF.Thresholding.MinPhotons);
fprintf('  MaxXY_SE: %.3f pixels\n', SMF.Thresholding.MaxXY_SE);
fprintf('  MinPValue: %.3f\n', SMF.Thresholding.MinPValue);
fprintf('  AutoThreshLogL: %d\n', SMF.Thresholding.AutoThreshLogL);

% Test on subset
LD = smi_core.LocalizeData(subset, SMF);
SMD_test = LD.genLocalizations();

if isfield(SMD_test, 'ThreshFlag')
    rejected = sum(SMD_test.ThreshFlag > 0);
    fprintf('Thresholding rejected: %d (%.1f%%)\n', ...
        rejected, 100 * rejected / length(SMD_test.ThreshFlag));
end
```

**Solution 3: Reduce background**

Lower background reduces noise detections:

```matlab
% If background is high (>30-40 photons/pixel), consider:

% Option A: Sample preparation improvements
fprintf('Current median background: %.1f photons/pixel\n', median(SMD.Bg));
fprintf('Consider:\n');
fprintf('  - Better washing protocol\n');
fprintf('  - Oxygen scavenging system (for dSTORM)\n');
fprintf('  - Reduce unbound fluorophore concentration\n');
fprintf('  - Use background suppression techniques\n');

% Option B: Spatially-varying background subtraction
% For structured background (e.g., autofluorescence patterns)
if false  % Enable for testing
    % Create background model from median of all frames
    background_model = median(rawData, 3);

    % Subtract from each frame
    rawData_corrected = rawData;
    for i = 1:size(rawData, 3)
        rawData_corrected(:,:,i) = max(rawData(:,:,i) - background_model, 0);
    end

    % Reanalyze
    LD = smi_core.LocalizeData(rawData_corrected, SMF);
    SMD_test = LD.genLocalizations();

    fprintf('After background subtraction:\n');
    fprintf('  Median background: %.1f photons/pixel\n', median(SMD_test.Bg));
    fprintf('  Localizations: %d (was %d)\n', length(SMD_test.X), length(SMD.X));
end
```

**Solution 4: Verify detections are molecule-like**

Use shape and size filtering:

```matlab
% P-value filtering already helps (rejects non-Gaussian shapes)
% But can add additional checks:

% Check PSF sigma (if fitted)
if isfield(SMD, 'PSFSigma')
    % Noise peaks often have wrong PSF size
    figure;
    subplot(2,1,1);
    histogram(SMD.PSFSigma, 50);
    xlabel('Fitted PSF Sigma (pixels)');
    title('PSF Sigma Distribution');

    subplot(2,1,2);
    scatter(SMD.Photons, SMD.PSFSigma, 5, 'filled', 'MarkerFaceAlpha', 0.3);
    xlabel('Photons'); ylabel('PSF Sigma');
    title('Photons vs PSF Sigma (noise often has extreme values)');

    % Filter outliers
    median_sigma = median(SMD.PSFSigma);
    sigma_std = std(SMD.PSFSigma);
    valid = abs(SMD.PSFSigma - median_sigma) < 2 * sigma_std;
    fprintf('PSF sigma filtering removes: %d (%.1f%%)\n', ...
        sum(~valid), 100 * sum(~valid) / length(valid));
end
```

### Expected Outcomes

After implementing solutions:
- False positive rate <5% (if ground truth available)
- Super-resolution image shows biological structures without noise speckles
- Photon distribution has clear peak, no extended low-photon tail
- P-value distribution more uniform (fewer very low values)
- Spatial distribution shows structure rather than randomness

## Problem 4: Bad P-Values

### Symptoms

Poor fit quality indicated by p-values:

```matlab
% P-value assessment
fprintf('P-value statistics:\n');
fprintf('  Median: %.3f\n', median(SMD.PValue));
fprintf('  Mean: %.3f\n', mean(SMD.PValue));
fprintf('  < 0.01: %.1f%%\n', 100 * sum(SMD.PValue < 0.01) / length(SMD.PValue));
fprintf('  < 0.001: %.1f%%\n', 100 * sum(SMD.PValue < 0.001) / length(SMD.PValue));

% P-value distribution
figure;
histogram(SMD.PValue, 50);
xlabel('P-value'); ylabel('Count');
title('Fit Quality (P-value) Distribution');
xline(0.01, 'r--', 'LineWidth', 2, 'Label', 'Typical threshold');

% For good fits, p-values should be roughly uniform
% Excess low p-values indicates systematic problems
```

**Key indicators:**
- More than 20% of localizations with p-value <0.01
- Many very low p-values (<0.001)
- P-value distribution strongly skewed toward zero
- Lower median p-value (<0.1)

**Visual check:**
```matlab
% Spatial distribution of poor fits
low_pval = SMD.PValue < 0.01;

figure;
subplot(1,2,1);
plot(SMD.X, SMD.Y, '.', 'MarkerSize', 1, 'Color', [0.7 0.7 0.7]);
hold on;
plot(SMD.X(low_pval), SMD.Y(low_pval), 'r.', 'MarkerSize', 3);
axis equal; axis tight;
title('Spatial Distribution (red = low p-value)');
legend('Good fits', 'Poor fits (p<0.01)');

% Correlation with other parameters
subplot(1,2,2);
scatter(SMD.Photons, SMD.PValue, 5, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Photons'); ylabel('P-value');
title('Photons vs P-value');
set(gca, 'YScale', 'log');
```

### Diagnosis

P-value is goodness-of-fit for PSF model. Low values indicate mismatch.

**1. Check for overlapping molecules** - Most common cause:

```matlab
% High density causes overlaps
area_um2 = SMD.XSize * SMD.YSize * SMF.Data.PixelSize^2;
density = length(SMD.X) / area_um2;
fprintf('Localization density: %.1f per square micron\n', density);

if density > 5
    fprintf('WARNING: High density (>5/um^2) increases overlap probability\n');
end

% Check localizations per frame
avg_per_frame = length(SMD.X) / SMD.NFrames;
fprintf('Average localizations per frame: %.1f\n', avg_per_frame);

% Estimate overlap probability (crude)
% For random distribution, probability of overlap within PSF area:
PSF_area_pixels = pi * (2 * SMF.Fitting.PSFSigma)^2;
overlap_prob = 1 - exp(-avg_per_frame * PSF_area_pixels / (SMD.XSize * SMD.YSize));
fprintf('Estimated overlap probability: %.1f%%\n', 100 * overlap_prob);

if overlap_prob > 0.1
    fprintf('WARNING: >10%% overlap probability expected\n');
end
```

**2. Check PSFSigma mismatch** - Wrong model parameter:

```matlab
% If sigma was fitted, compare to assumed value
if isfield(SMD, 'PSFSigma')
    fprintf('Assumed PSFSigma: %.3f pixels\n', SMF.Fitting.PSFSigma);
    fprintf('Fitted PSFSigma: %.3f +/- %.3f pixels (mean +/- std)\n', ...
        mean(SMD.PSFSigma), std(SMD.PSFSigma));

    % Check correlation between sigma mismatch and p-value
    if SMF.Fitting.PSFSigma > 0
        sigma_error = abs(SMD.PSFSigma - SMF.Fitting.PSFSigma);

        figure;
        scatter(sigma_error, SMD.PValue, 5, 'filled', 'MarkerFaceAlpha', 0.3);
        xlabel('|Fitted - Assumed| PSF Sigma (pixels)');
        ylabel('P-value');
        title('PSF Sigma Error vs Fit Quality');
        set(gca, 'YScale', 'log');

        % Strong negative correlation indicates sigma mismatch problem
    end
end
```

**3. Check for non-Gaussian PSF** - Aberrations or astigmatism:

```matlab
% Non-Gaussian PSFs give poor fits
% Check if problem is spatially varying (aberrations)

figure;
scatter(SMD.X, SMD.Y, 10, log10(SMD.PValue), 'filled');
colormap jet; colorbar;
title('Spatial P-value Map (log10)');
axis equal; axis tight;

% Systematic patterns suggest:
%   - Center-edge variation: Spherical aberration
%   - Directional patterns: Astigmatism
%   - Random: Overlapping molecules

% Check for asymmetry (if PSFSigmaX/Y available)
if isfield(SMD, 'PSFSigmaX') && isfield(SMD, 'PSFSigmaY')
    asymmetry = abs(SMD.PSFSigmaX - SMD.PSFSigmaY) ./ ...
                ((SMD.PSFSigmaX + SMD.PSFSigmaY) / 2);
    fprintf('PSF asymmetry: %.1f%% +/- %.1f%% (mean +/- std)\n', ...
        100 * mean(asymmetry), 100 * std(asymmetry));

    if mean(asymmetry) > 0.15
        fprintf('WARNING: Significant PSF asymmetry detected\n');
    end
end
```

**4. Check residuals** - Pattern in residuals indicates model problem:

```matlab
% Unfortunately SMD doesn't store residuals directly
% But can refit a few localizations to check

% This requires access to raw data and is more involved
% Key indicator: Structured residuals (not random noise)
```

### Solutions

**Solution 1: Reduce emitter density (for overlaps)**

If overlap is the issue:

```matlab
% This requires new data acquisition
fprintf('To reduce overlaps:\n');
fprintf('  DNA-PAINT: Lower imager strand concentration\n');
fprintf('  dSTORM: Faster photobleaching (higher laser power initially)\n');
fprintf('  PALM: Use sparser labeling or faster photoconversion\n');
fprintf('\nTarget density: <3 localizations per square micron\n');
fprintf('Current density: %.1f per square micron\n', density);

% For existing data, can try multi-emitter fitting (advanced)
% Or accept higher p-value threshold for overlapping regions
```

**Solution 2: Correct PSFSigma**

Critical for good fits:

```matlab
% Determine correct PSFSigma empirically
fprintf('Determining optimal PSFSigma...\n');

SMF_test = SMF;
SMF_test.Fitting.FitType = 'XYNBS';  % Fit sigma
SMF_test.Fitting.PSFSigma = 1.5;     % Initial guess

LD = smi_core.LocalizeData(subset, SMF_test);
SMD_test = LD.genLocalizations();

% Find optimal value
optimal_sigma = median(SMD_test.PSFSigma);
fprintf('Optimal PSFSigma: %.3f pixels (was %.3f)\n', ...
    optimal_sigma, SMF.Fitting.PSFSigma);

% Check p-value improvement
fprintf('Median p-value with fitted sigma: %.3f\n', median(SMD_test.PValue));

% Update SMF
SMF.Fitting.PSFSigma = optimal_sigma;
SMF.Fitting.FitType = 'XYNB';  % Fix sigma at optimal value

% Reanalyze
LD = smi_core.LocalizeData(subset, SMF);
SMD_corrected = LD.genLocalizations();

fprintf('Median p-value with fixed optimal sigma: %.3f\n', median(SMD_corrected.PValue));
fprintf('Localizations with p<0.01: %.1f%% (was %.1f%%)\n', ...
    100 * sum(SMD_corrected.PValue < 0.01) / length(SMD_corrected.PValue), ...
    100 * sum(SMD.PValue < 0.01) / length(SMD.PValue));
```

**Solution 3: Use asymmetric PSF model**

For astigmatic or aberrated PSFs:

```matlab
% Try fitting separate X and Y sigmas
SMF_test = SMF;
SMF_test.Fitting.FitType = 'XYNBSXSY';  % Fit independent sigmaX, sigmaY

LD = smi_core.LocalizeData(subset, SMF_test);
SMD_test = LD.genLocalizations();

fprintf('Asymmetric PSF fitting:\n');
fprintf('  Median p-value: %.3f (was %.3f)\n', ...
    median(SMD_test.PValue), median(SMD.PValue));
fprintf('  SigmaX: %.3f +/- %.3f\n', mean(SMD_test.PSFSigmaX), std(SMD_test.PSFSigmaX));
fprintf('  SigmaY: %.3f +/- %.3f\n', mean(SMD_test.PSFSigmaY), std(SMD_test.PSFSigmaY));

% If p-values improve significantly, use this model
improvement = (median(SMD.PValue) - median(SMD_test.PValue)) / median(SMD.PValue);
if improvement > 0.2  % >20%% improvement
    fprintf('Asymmetric PSF model improves fits by %.1f%%\n', 100 * improvement);
    SMF.Fitting.FitType = 'XYNBSXSY';
end
```

**Solution 4: Accept lower p-values with appropriate filtering**

If overlaps are unavoidable:

```matlab
% For high-density data, many low p-values are expected
% Use appropriate threshold rather than default 0.01

% Analyze p-value distribution to find reasonable cutoff
p_percentiles = prctile(SMD.PValue, [5, 10, 20, 30]);
fprintf('P-value percentiles: 5th=%.3f, 10th=%.3f, 20th=%.3f, 30th=%.3f\n', ...
    p_percentiles);

% Set threshold to exclude worst ~10%%
new_threshold = p_percentiles(2);  % 10th percentile
fprintf('Suggested MinPValue threshold: %.4f (current: %.4f)\n', ...
    new_threshold, SMF.Thresholding.MinPValue);

SMF.Thresholding.MinPValue = max(new_threshold, 0.001);  % Don't go below 0.001

% Also enable log-likelihood auto-thresholding
SMF.Thresholding.AutoThreshLogL = true;
fprintf('Enabled AutoThreshLogL for adaptive quality filtering\n');
```

### Expected Outcomes

After implementing solutions:
- Fewer than 10-15% of localizations with p-value <0.01
- Median p-value >0.1
- More uniform p-value distribution
- No strong spatial patterns in p-value map
- Better agreement between fitted and assumed PSF parameters

## Problem 5: Wrong PSF Model

### Symptoms

Incorrect PSF assumptions cause systematic biases:

```matlab
% Check for PSF-related issues
fprintf('Current PSF model:\n');
fprintf('  PSFSigma: %.3f pixels\n', SMF.Fitting.PSFSigma);
fprintf('  FitType: %s\n', SMF.Fitting.FitType);

% If sigma was fitted, check distribution
if contains(SMF.Fitting.FitType, 'S')
    fprintf('  Fitted sigma range: %.3f to %.3f pixels\n', ...
        min(SMD.PSFSigma), max(SMD.PSFSigma));
    fprintf('  Fitted sigma std: %.3f pixels\n', std(SMD.PSFSigma));

    % High variability suggests wrong model
    if std(SMD.PSFSigma) / mean(SMD.PSFSigma) > 0.2
        fprintf('WARNING: High PSF sigma variability (>20%%)\n');
    end
end
```

**Key indicators:**
- Biased photon estimates (systematic over/under-estimation)
- Poor p-values even at low density
- Spatial variation in fitted PSF parameters
- Precision doesn't improve with more photons as expected
- Different results in different regions of field of view

**Visual checks:**
```matlab
% Spatial PSF variation
if isfield(SMD, 'PSFSigma')
    figure;
    subplot(1,2,1);
    scatter(SMD.X, SMD.Y, 10, SMD.PSFSigma, 'filled');
    colormap hot; colorbar;
    title('Spatial PSF Sigma Variation');
    axis equal; axis tight;

    % Radial dependence (aberrations)
    center_x = SMD.XSize / 2;
    center_y = SMD.YSize / 2;
    radius = sqrt((SMD.X - center_x).^2 + (SMD.Y - center_y).^2);

    subplot(1,2,2);
    scatter(radius, SMD.PSFSigma, 5, 'filled', 'MarkerFaceAlpha', 0.3);
    xlabel('Distance from Center (pixels)');
    ylabel('PSF Sigma (pixels)');
    title('Radial PSF Variation (flat = good)');
    grid on;
end

% Asymmetry check
if isfield(SMD, 'PSFSigmaX') && isfield(SMD, 'PSFSigmaY')
    figure;
    plot(SMD.PSFSigmaX, SMD.PSFSigmaY, '.', 'MarkerSize', 2);
    hold on;
    plot([0.8, 2], [0.8, 2], 'r--', 'LineWidth', 2);
    xlabel('PSF Sigma X'); ylabel('PSF Sigma Y');
    title('PSF Symmetry (points on diagonal = symmetric)');
    axis equal; grid on;
end
```

### Diagnosis

**1. Determine actual PSF shape** - Measure from data:

```matlab
% Fit sigma as free parameter to see true PSF
SMF_diag = SMF;
SMF_diag.Fitting.FitType = 'XYNBSXSY';  % Fit X and Y sigma separately
SMF_diag.Fitting.PSFSigma = 1.5;  % Initial guess

LD = smi_core.LocalizeData(subset, SMF_diag);
SMD_diag = LD.genLocalizations();

fprintf('Measured PSF characteristics:\n');
fprintf('  Sigma X: %.3f +/- %.3f pixels\n', mean(SMD_diag.PSFSigmaX), std(SMD_diag.PSFSigmaX));
fprintf('  Sigma Y: %.3f +/- %.3f pixels\n', mean(SMD_diag.PSFSigmaY), std(SMD_diag.PSFSigmaY));
fprintf('  Asymmetry: %.1f%%\n', ...
    100 * abs(mean(SMD_diag.PSFSigmaX) - mean(SMD_diag.PSFSigmaY)) / ...
    mean([mean(SMD_diag.PSFSigmaX), mean(SMD_diag.PSFSigmaY)]));

% Compare to assumed
fprintf('  Assumed PSFSigma: %.3f pixels\n', SMF.Fitting.PSFSigma);
```

**2. Check for spatial dependence** - Indicates aberrations:

```matlab
% Divide field into quadrants
mid_x = SMD_diag.XSize / 2;
mid_y = SMD_diag.YSize / 2;

quadrants = {
    'Top-Left', SMD_diag.X < mid_x & SMD_diag.Y < mid_y;
    'Top-Right', SMD_diag.X >= mid_x & SMD_diag.Y < mid_y;
    'Bottom-Left', SMD_diag.X < mid_x & SMD_diag.Y >= mid_y;
    'Bottom-Right', SMD_diag.X >= mid_x & SMD_diag.Y >= mid_y
};

fprintf('\nPSF by quadrant:\n');
for i = 1:size(quadrants, 1)
    mask = quadrants{i,2};
    fprintf('  %s: SigmaX=%.3f, SigmaY=%.3f\n', ...
        quadrants{i,1}, mean(SMD_diag.PSFSigmaX(mask)), mean(SMD_diag.PSFSigmaY(mask)));
end
```

**3. Check wavelength and optical parameters** - Theoretical PSF:

```matlab
% Theoretical PSF sigma (pixels)
wavelength_nm = 670;  % Emission wavelength (e.g., Alexa647)
NA = 1.4;  % Objective numerical aperture
pixel_size_nm = SMF.Data.PixelSize * 1000;  % nm

% Diffraction-limited sigma (Gaussian approximation)
sigma_theory_nm = 0.21 * wavelength_nm / NA;  % FWHM / 2.355
sigma_theory_pixels = sigma_theory_nm / pixel_size_nm;

fprintf('Theoretical PSF sigma: %.3f pixels\n', sigma_theory_pixels);
fprintf('Measured PSF sigma: %.3f pixels\n', mean([mean(SMD_diag.PSFSigmaX), mean(SMD_diag.PSFSigmaY)]));

ratio = mean([mean(SMD_diag.PSFSigmaX), mean(SMD_diag.PSFSigmaY)]) / sigma_theory_pixels;
fprintf('Measured / Theoretical: %.2f\n', ratio);

if ratio < 0.7 || ratio > 1.3
    fprintf('WARNING: Measured PSF deviates significantly from theory\n');
    fprintf('  Check focus, optical alignment, or refractive index mismatch\n');
end
```

### Solutions

**Solution 1: Use correct PSFSigma value**

Most common issue:

```matlab
% Use measured value from diagnosis
optimal_sigma_x = median(SMD_diag.PSFSigmaX);
optimal_sigma_y = median(SMD_diag.PSFSigmaY);
optimal_sigma = (optimal_sigma_x + optimal_sigma_y) / 2;

fprintf('Updating PSFSigma: %.3f -> %.3f pixels\n', ...
    SMF.Fitting.PSFSigma, optimal_sigma);

SMF.Fitting.PSFSigma = optimal_sigma;
SMF.Fitting.FitType = 'XYNB';  % Fix sigma at correct value

% Update dependent parameters
SMF.BoxFinding.BoxSize = ceil(4.5 * optimal_sigma);
fprintf('Updated BoxSize: %d pixels\n', SMF.BoxFinding.BoxSize);

% Set reasonable filter bounds
SMF.Thresholding.MinPSFSigma = optimal_sigma * 0.7;
SMF.Thresholding.MaxPSFSigma = optimal_sigma * 1.3;

% Reanalyze
LD = smi_core.LocalizeData(subset, SMF);
SMD_corrected = LD.genLocalizations();

% Check improvement
fprintf('\nResults comparison:\n');
fprintf('  Precision: %.1f nm -> %.1f nm\n', ...
    median(SMD.X_SE) * SMF.Data.PixelSize * 1000, ...
    median(SMD_corrected.X_SE) * SMF.Data.PixelSize * 1000);
fprintf('  Median p-value: %.3f -> %.3f\n', ...
    median(SMD.PValue), median(SMD_corrected.PValue));
```

**Solution 2: Use asymmetric PSF model**

For astigmatism or elliptical PSFs:

```matlab
% If X and Y sigmas differ significantly
asymmetry = abs(optimal_sigma_x - optimal_sigma_y) / optimal_sigma;

if asymmetry > 0.15  % >15%% asymmetry
    fprintf('Using asymmetric PSF model (%.1f%% asymmetry)\n', 100 * asymmetry);

    SMF.Fitting.FitType = 'XYNBSXSY';  % Fit X and Y sigma independently
    SMF.Fitting.PSFSigma = optimal_sigma;  % Initial guess

    % Reanalyze
    LD = smi_core.LocalizeData(subset, SMF);
    SMD_corrected = LD.genLocalizations();

    fprintf('Asymmetric fit results:\n');
    fprintf('  SigmaX: %.3f +/- %.3f pixels\n', ...
        mean(SMD_corrected.PSFSigmaX), std(SMD_corrected.PSFSigmaX));
    fprintf('  SigmaY: %.3f +/- %.3f pixels\n', ...
        mean(SMD_corrected.PSFSigmaY), std(SMD_corrected.PSFSigmaY));
    fprintf('  Median p-value: %.3f\n', median(SMD_corrected.PValue));
else
    fprintf('Symmetric PSF model appropriate (%.1f%% asymmetry)\n', 100 * asymmetry);
end
```

**Solution 3: For 3D imaging with astigmatism**

If using astigmatic PSF for 3D:

```matlab
% smite supports astigmatic PSF for 3D localization
% Requires calibration curve: Z position vs (sigmaX, sigmaY)

% Set up for 3D fitting (if calibrated)
if false  % Enable if you have calibration
    SMF.Fitting.FitType = 'XYNBSXSYZ';  % 3D with astigmatism
    % Need to load calibration data
    % SMF.Fitting.PSFCalibration = ... load calibration file
end
```

**Solution 4: Check optical setup**

If measured PSF differs significantly from theory:

```matlab
fprintf('\nOptical troubleshooting checklist:\n');
fprintf('1. Focus: Is sample in focus?\n');
fprintf('2. Immersion: Correct immersion medium (oil/water/glycerol)?\n');
fprintf('3. Coverslip thickness: #1.5 (0.17mm) required for high-NA objectives\n');
fprintf('4. Refractive index: Sample medium matches objective design?\n');
fprintf('5. Chromatic aberration: Using correct wavelength for focus?\n');
fprintf('6. Optical alignment: Microscope aligned and collimated?\n');
fprintf('7. Dust/dirt: Clean all optical surfaces\n');

% Measured-to-theory ratio gives clues:
if ratio > 1.3
    fprintf('\nPSF too large suggests:\n');
    fprintf('  - Out of focus\n');
    fprintf('  - Wrong immersion medium\n');
    fprintf('  - Spherical aberration\n');
elseif ratio < 0.7
    fprintf('\nPSF too small suggests:\n');
    fprintf('  - Wrong pixel size calibration\n');
    fprintf('  - Magnification error\n');
    fprintf('  - Check camera settings\n');
end
```

### Expected Outcomes

After correcting PSF model:
- Fitted PSF parameters stable across field of view
- P-values improve significantly
- Precision estimates accurate
- Photon estimates unbiased
- Consistent results across different sample regions

## Problem 6: Drift Artifacts

### Symptoms

Stage drift during acquisition causes blurring:

```matlab
% Check if drift correction was performed
if isfield(SMD, 'DriftX') && ~isempty(SMD.DriftX)
    total_drift_x = SMD.DriftX(end,1) - SMD.DriftX(1,1);
    total_drift_y = SMD.DriftY(end,1) - SMD.DriftY(1,1);
    total_drift_pixels = sqrt(total_drift_x^2 + total_drift_y^2);
    total_drift_nm = total_drift_pixels * SMF.Data.PixelSize * 1000;

    fprintf('Drift correction results:\n');
    fprintf('  Total drift: %.1f nm (%.2f pixels)\n', ...
        total_drift_nm, total_drift_pixels);
    fprintf('  X drift: %.1f nm\n', abs(total_drift_x) * SMF.Data.PixelSize * 1000);
    fprintf('  Y drift: %.1f nm\n', abs(total_drift_y) * SMF.Data.PixelSize * 1000);

    if total_drift_nm > 100
        fprintf('WARNING: Significant drift detected (>100 nm)\n');
    end
else
    fprintf('WARNING: No drift correction performed\n');
    fprintf('  Set SMF.DriftCorrection.On = true\n');
end
```

**Key indicators:**
- Blurry or doubled structures in super-resolution image
- Rainbow streaks in drift-colored visualization
- Elongated or smeared features
- Loss of fine structural detail
- Temporal color separation in drift images

**Visual checks:**
```matlab
% Drift trajectory
if isfield(SMD, 'DriftX')
    figure;
    plot(SMD.DriftX(:,1) * SMF.Data.PixelSize * 1000, ...
         SMD.DriftY(:,1) * SMF.Data.PixelSize * 1000, '-o', 'LineWidth', 1.5);
    xlabel('Drift X (nm)'); ylabel('Drift Y (nm)');
    title('Stage Drift Trajectory');
    axis equal; grid on;

    % Drift over time
    figure;
    subplot(2,1,1);
    plot(SMD.DriftX(:,1) * SMF.Data.PixelSize * 1000, 'LineWidth', 1.5);
    ylabel('X Drift (nm)'); xlabel('Frame');
    title('Drift vs Time');
    grid on;

    subplot(2,1,2);
    plot(SMD.DriftY(:,1) * SMF.Data.PixelSize * 1000, 'LineWidth', 1.5);
    ylabel('Y Drift (nm)'); xlabel('Frame');
    grid on;
end

% Drift-colored image
[DriftIm, DriftImRGB] = smi_vis.GenerateImages.driftImage(SMD, 20);
figure;
imshow(DriftImRGB);
title('Drift Image (Blue=early, Red=late)');
% Single-color structures = good
% Rainbow streaks = uncorrected drift
```

### Diagnosis

**1. Quantify drift magnitude and pattern:**

```matlab
if isfield(SMD, 'DriftX')
    % Frame-to-frame drift
    frame_drift_x = diff(SMD.DriftX(:,1));
    frame_drift_y = diff(SMD.DriftY(:,1));
    frame_drift_mag = sqrt(frame_drift_x.^2 + frame_drift_y.^2);

    fprintf('Frame-to-frame drift:\n');
    fprintf('  Mean: %.2f nm/frame\n', mean(frame_drift_mag) * SMF.Data.PixelSize * 1000);
    fprintf('  Max: %.2f nm/frame\n', max(frame_drift_mag) * SMF.Data.PixelSize * 1000);

    % Drift rate
    frame_rate = SMF.Data.FrameRate;
    drift_rate = mean(frame_drift_mag) * SMF.Data.PixelSize * 1000 * frame_rate;
    fprintf('  Drift rate: %.1f nm/second\n', drift_rate);

    % Check for systematic drift (linear) vs random wandering
    x_trend = polyfit(1:length(SMD.DriftX(:,1)), SMD.DriftX(:,1), 1);
    y_trend = polyfit(1:length(SMD.DriftY(:,1)), SMD.DriftY(:,1), 1);

    fprintf('\nDrift pattern:\n');
    fprintf('  X trend: %.3f pixels/frame\n', x_trend(1));
    fprintf('  Y trend: %.3f pixels/frame\n', y_trend(1));

    if abs(x_trend(1)) > 0.01 || abs(y_trend(1)) > 0.01
        fprintf('  -> Systematic drift detected\n');
    else
        fprintf('  -> Random drift (thermal)\n');
    end
else
    fprintf('Cannot diagnose drift - no drift correction was performed\n');
end
```

**2. Check drift correction method:**

```matlab
fprintf('\nDrift correction settings:\n');
fprintf('  Enabled: %d\n', SMF.DriftCorrection.On);
if SMF.DriftCorrection.On
    fprintf('  Method: %s\n', SMF.DriftCorrection.Method);
    fprintf('  FrameRate: %.1f Hz\n', SMF.Data.FrameRate);

    % Method-specific parameters
    switch SMF.DriftCorrection.Method
        case 'DC-KNN'
            fprintf('  K neighbors: %d\n', SMF.DriftCorrection.KNN);
        case 'RCC'
            fprintf('  Frame binning: %d\n', SMF.DriftCorrection.FrameBinning);
        case 'GCC'
            fprintf('  Frame binning: %d\n', SMF.DriftCorrection.FrameBinning);
    end
end
```

**3. Assess correction quality:**

```matlab
% Compare corrected vs uncorrected positions
if isfield(SMD, 'X_SE')
    % Effective precision includes drift
    frame_groups = discretize(SMD.FrameNum, 10);  % Divide into 10 time bins

    fprintf('\nPrecision by time bin:\n');
    for i = 1:10
        mask = (frame_groups == i);
        if sum(mask) > 10
            prec = median(SMD.X_SE(mask)) * SMF.Data.PixelSize * 1000;
            fprintf('  Bin %d: %.1f nm\n', i, prec);
        end
    end
    fprintf('If precision degrades over time -> drift not fully corrected\n');
end
```

### Solutions

**Solution 1: Enable drift correction**

If not enabled:

```matlab
% Enable drift correction
SMF.DriftCorrection.On = true;

% Choose method (DC-KNN is default and robust)
SMF.DriftCorrection.Method = 'DC-KNN';

% DC-KNN parameters
SMF.DriftCorrection.KNN = 10;  % Number of nearest neighbors
SMF.DriftCorrection.PDegree = 1;  % Polynomial degree (1=linear)
SMF.DriftCorrection.FrameBinning = 1;  % Typically 1

fprintf('Enabled drift correction: %s\n', SMF.DriftCorrection.Method);

% Reanalyze
LD = smi_core.LocalizeData(rawData, SMF);
SMD_corrected = LD.genLocalizations();

% Check improvement
if isfield(SMD_corrected, 'DriftX')
    fprintf('Drift after correction: %.1f nm\n', ...
        sqrt(SMD_corrected.DriftX(end,1)^2 + SMD_corrected.DriftY(end,1)^2) * ...
        SMF.Data.PixelSize * 1000);
end
```

**Solution 2: Try different drift correction methods**

If current method inadequate:

```matlab
% Available methods in smite:
methods = {'DC-KNN', 'RCC', 'GCC', 'RDC'};

fprintf('Testing drift correction methods:\n');
for i = 1:length(methods)
    SMF_test = SMF;
    SMF_test.DriftCorrection.On = true;
    SMF_test.DriftCorrection.Method = methods{i};

    % Method-specific tuning
    switch methods{i}
        case 'DC-KNN'
            SMF_test.DriftCorrection.KNN = 10;
        case {'RCC', 'GCC'}
            SMF_test.DriftCorrection.FrameBinning = 50;  % Frames per bin
        case 'RDC'
            % Redundant cross-correlation
            SMF_test.DriftCorrection.FrameBinning = 50;
    end

    LD = smi_core.LocalizeData(subset, SMF_test);
    SMD_test = LD.genLocalizations();

    if isfield(SMD_test, 'DriftX')
        total = sqrt(SMD_test.DriftX(end,1)^2 + SMD_test.DriftY(end,1)^2);
        fprintf('  %s: %.1f nm residual\n', methods{i}, total * SMF.Data.PixelSize * 1000);
    end
end

% Choose best method (lowest residual drift)
fprintf('\nRecommendations:\n');
fprintf('  DC-KNN: Best for sparse, unevenly distributed localizations\n');
fprintf('  RCC: Fast, good for dense uniform data\n');
fprintf('  GCC: Similar to RCC, try if RCC fails\n');
fprintf('  RDC: Robust redundant method, slower\n');
```

**Solution 3: Tune drift correction parameters**

For DC-KNN (most commonly used):

```matlab
% Test different KNN values
knn_values = [5, 10, 20, 50];

fprintf('Tuning DC-KNN parameter:\n');
for i = 1:length(knn_values)
    SMF_test = SMF;
    SMF_test.DriftCorrection.Method = 'DC-KNN';
    SMF_test.DriftCorrection.KNN = knn_values(i);

    LD = smi_core.LocalizeData(subset, SMF_test);
    SMD_test = LD.genLocalizations();

    if isfield(SMD_test, 'DriftX')
        total = sqrt(SMD_test.DriftX(end,1)^2 + SMD_test.DriftY(end,1)^2);
        fprintf('  KNN=%d: %.1f nm residual\n', knn_values(i), ...
            total * SMF.Data.PixelSize * 1000);
    end
end

% Rule of thumb:
%   Sparse data: Lower KNN (5-10)
%   Dense data: Higher KNN (20-50)
%   Start with 10 and adjust

% Optimal value minimizes residual without overfitting
```

For RCC/GCC (correlation-based):

```matlab
% Test different frame binning
binning_values = [20, 50, 100, 200];

fprintf('Tuning RCC frame binning:\n');
for i = 1:length(binning_values)
    SMF_test = SMF;
    SMF_test.DriftCorrection.Method = 'RCC';
    SMF_test.DriftCorrection.FrameBinning = binning_values(i);

    LD = smi_core.LocalizeData(subset, SMF_test);
    SMD_test = LD.genLocalizations();

    if isfield(SMD_test, 'DriftX')
        total = sqrt(SMD_test.DriftX(end,1)^2 + SMD_test.DriftY(end,1)^2);
        fprintf('  Binning=%d: %.1f nm residual\n', binning_values(i), ...
            total * SMF.Data.PixelSize * 1000);
    end
end

% Trade-off:
%   Smaller bins: Better time resolution but noisier
%   Larger bins: Smoother but may miss fast drift
%   Typical: 50-100 frames
```

**Solution 4: Use fiducial markers**

Most accurate drift correction (requires markers in sample):

```matlab
% If your sample has fiducial markers (e.g., gold nanoparticles)
% Track them separately and use for drift correction

% This is an advanced technique requiring:
% 1. Fiducials visible in all frames
% 2. Separate localization of fiducials
% 3. Drift trajectory from fiducial tracks

fprintf('Fiducial-based drift correction:\n');
fprintf('  Requires: Gold nanoparticles or fluorescent beads in sample\n');
fprintf('  Advantage: Most accurate (<10 nm)\n');
fprintf('  Process:\n');
fprintf('    1. Localize fiducials in each frame\n');
fprintf('    2. Track fiducial positions\n');
fprintf('    3. Use tracks as drift trajectory\n');
fprintf('    4. Apply to molecular localizations\n');

% Implementation would require custom code or smite's fiducial tracking
```

**Solution 5: Improve microscope stability**

Prevent drift at source:

```matlab
fprintf('\nMicroscope stability improvements:\n');
fprintf('1. Temperature control:\n');
fprintf('   - Equilibrate microscope (1-2 hours before acquisition)\n');
fprintf('   - Stable room temperature\n');
fprintf('   - Enclosure around microscope\n');
fprintf('2. Mechanical stability:\n');
fprintf('   - Anti-vibration table\n');
fprintf('   - Minimize air currents\n');
fprintf('   - Secure all components\n');
fprintf('3. Focus stability:\n');
fprintf('   - Hardware autofocus system\n');
fprintf('   - Objective heater to prevent thermal drift\n');
fprintf('4. Sample preparation:\n');
fprintf('   - Secure coverslip mounting\n');
fprintf('   - Minimize evaporation (seal chamber)\n');
fprintf('   - Thermal equilibration of sample\n');
```

### Expected Outcomes

After correcting drift:
- Total drift <50 nm (excellent) or 50-100 nm (good)
- Uniform color in drift-colored images (no rainbow streaks)
- Sharp structural features in super-resolution images
- No temporal smearing or doubling of structures
- Consistent precision across all time points

## Summary and Quick Reference

### Decision Tree

```matlab
% Quick diagnostic
load('Results.mat', 'SMD', 'SMF');

% 1. Check precision
precision_nm = median(SMD.X_SE) * SMF.Data.PixelSize * 1000;
if precision_nm > 30
    fprintf('PROBLEM: Poor precision -> See Problem 1\n');
end

% 2. Check detection
avg_per_frame = length(SMD.X) / SMD.NFrames;
if avg_per_frame < expected_per_frame * 0.8  % Know expected density
    fprintf('PROBLEM: Low detection -> See Problem 2\n');
end

% 3. Check false positives
low_photon_fraction = sum(SMD.Photons < 200) / length(SMD.Photons);
if low_photon_fraction > 0.3
    fprintf('PROBLEM: Many false positives -> See Problem 3\n');
end

% 4. Check fit quality
bad_pval_fraction = sum(SMD.PValue < 0.01) / length(SMD.PValue);
if bad_pval_fraction > 0.2
    fprintf('PROBLEM: Poor fit quality -> See Problem 4\n');
end

% 5. Check PSF
if isfield(SMD, 'PSFSigma')
    sigma_cv = std(SMD.PSFSigma) / mean(SMD.PSFSigma);
    if sigma_cv > 0.2
        fprintf('PROBLEM: Wrong PSF model -> See Problem 5\n');
    end
end

% 6. Check drift
if isfield(SMD, 'DriftX')
    total_drift_nm = sqrt(SMD.DriftX(end,1)^2 + SMD.DriftY(end,1)^2) * ...
                     SMF.Data.PixelSize * 1000;
    if total_drift_nm > 100
        fprintf('PROBLEM: Drift artifacts -> See Problem 6\n');
    end
end
```

### Parameter Quick Fixes

| Problem | Most Common Solution | Parameter to Adjust |
|---------|---------------------|---------------------|
| Poor precision | Increase signal | Increase laser power, exposure time |
| Poor precision | Reduce background | Sample prep, background subtraction |
| Poor precision | Fix PSFSigma | Fit with 'XYNBS', use median value |
| Low detection | Lower threshold | Reduce SMF.BoxFinding.MinPhotons |
| Low detection | Relax filtering | Increase SMF.Thresholding.MaxXY_SE |
| False positives | Raise threshold | Increase SMF.BoxFinding.MinPhotons |
| False positives | Enable filtering | SMF.Thresholding.On = true |
| Bad p-values | Correct PSFSigma | Fit and update SMF.Fitting.PSFSigma |
| Bad p-values | Reduce density | Lower emitter concentration |
| Wrong PSF | Measure empirically | FitType = 'XYNBS', find median |
| Wrong PSF | Use asymmetric | FitType = 'XYNBSXSY' |
| Drift | Enable correction | SMF.DriftCorrection.On = true |
| Drift | Try different method | Test DC-KNN, RCC, GCC, RDC |

### Target Quality Metrics

**Excellent results:**
- Precision: <15 nm
- Photons: >800
- Background: <20 photons/pixel
- SBR: >10
- P-value <0.01: <10%
- Detection rate: >95%
- False positive rate: <3%
- Drift: <50 nm

**Good results:**
- Precision: 15-25 nm
- Photons: 400-800
- Background: 20-35 photons/pixel
- SBR: 5-10
- P-value <0.01: 10-20%
- Detection rate: 85-95%
- False positive rate: 3-8%
- Drift: 50-100 nm

**Needs improvement:**
- Precision: >30 nm
- Photons: <300
- Background: >40 photons/pixel
- SBR: <5
- P-value <0.01: >25%
- Detection rate: <80%
- False positive rate: >10%
- Drift: >150 nm

## See Also

- [Understanding Results](../getting-started/understanding-results.md) - Comprehensive guide to result interpretation
- [Tune Parameters](../how-to/tune-parameters.md) - Systematic parameter optimization
- [Threshold Results](../how-to/threshold-results.md) - Quality filtering techniques
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete analysis pipeline
- [SMF Structure](../core-concepts/smf-structure.md) - All parameter documentation
- [SMD Structure](../core-concepts/smd-structure.md) - Results data format
