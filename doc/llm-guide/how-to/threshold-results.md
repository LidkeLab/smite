---
title: "How to Apply Quality Filters to Localization Results"
category: "how-to"
level: "beginner"
tags: ["thresholding", "quality-control", "filtering", "precision", "photons"]
prerequisites: ["localize-molecules.md", "../core-concepts/smd-structure.md"]
related: ["../workflows/smlm-analysis.md", "../core-concepts/smf-structure.md"]
summary: "Guide to applying quality filters to localization results based on precision, photons, fit quality, and PSF parameters"
estimated_time: "12 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# How to Apply Quality Filters to Localization Results

## Purpose

Not all localizations are created equal. Thresholding filters out poor quality fits caused by noise, overlapping molecules, or fitting failures. This guide explains how to configure and apply quality filters based on localization precision, photon counts, PSF parameters, and fit quality metrics.

## Prerequisites

- Understanding of [localization process](localize-molecules.md)
- Familiarity with [SMD structure](../core-concepts/smd-structure.md)
- Basic knowledge of MLE fitting quality metrics

## Overview

Thresholding in smite operates in two stages:

1. **Flag Creation**: Evaluate each localization against quality criteria and encode failures as bits in a `ThreshFlag` field
2. **Application**: Remove flagged localizations from the dataset

This approach allows you to inspect which filters removed which localizations before permanently filtering the data. The `smi_core.Threshold` class handles all thresholding operations.

## Threshold Criteria

smite supports filtering on multiple quality metrics:

| Criterion | SMF Parameter | Typical Range | Purpose |
|-----------|---------------|---------------|---------|
| **Precision (X,Y)** | `MaxXY_SE` | 0.1-0.2 pixels | Remove uncertain fits |
| **Precision (Z)** | `MaxZ_SE` | 0.05-0.5 microns | Remove uncertain Z (3D) |
| **Photons** | `MinPhotons` | 100-200 | Remove dim spots |
| **Background** | `MaxBg` | 20-50 photons/pixel | Remove high background |
| **PSF Sigma** | `MinPSFSigma`, `MaxPSFSigma` | 0.5-2.0 pixels | Remove misfits |
| **P-value** | `MinPValue` | 0.01-0.05 | Remove poor fits |
| **Log-likelihood** | `AutoThreshLogL` | automatic | Alternative to p-value |

### Understanding Each Filter

**Localization Precision (X_SE, Y_SE, Z_SE):**

Precision is the standard error from the Cramer-Rao Lower Bound (CRLB), representing the theoretical uncertainty of the fit. Low photon counts or high background lead to poor precision.

```matlab
% Typical precision threshold: 0.15-0.2 pixels (~15-20 nm)
SMF.Thresholding.MaxXY_SE = 0.2;  % pixels
SMF.Thresholding.MaxZ_SE = 0.05;  % microns (for 3D)
```

**Photon Count:**

Minimum photons ensures sufficient signal strength. Too few photons result in unreliable fits dominated by noise.

```matlab
% Typical range: 100-300 photons
SMF.Thresholding.MinPhotons = 150;
```

**Background:**

Maximum background rejects localizations in regions with excessive autofluorescence or scattered light.

```matlab
% Typical range: 20-50 photons/pixel
SMF.Thresholding.MaxBg = 40;  % photons/pixel
```

**PSF Sigma:**

For fits where PSF sigma is a free parameter (FitType='XYNBS'), this filters out aberrant fits with unrealistic PSF widths.

```matlab
% Typical range: 0.8-1.8 pixels for 1.3 pixel nominal sigma
SMF.Thresholding.MinPSFSigma = 0.8;
SMF.Thresholding.MaxPSFSigma = 1.8;
```

**P-value:**

P-value from chi-squared test of fit residuals. Low p-values indicate poor fit quality, suggesting overlapping molecules or model mismatch.

```matlab
% Typical threshold: 0.01
SMF.Thresholding.MinPValue = 0.01;
```

**Log-likelihood (alternative to p-value):**

Automatic thresholding based on log-likelihood distribution, more robust than fixed p-value cutoffs.

```matlab
% Enable automatic log-likelihood thresholding
SMF.Thresholding.AutoThreshLogL = true;
SMF.Thresholding.AutoThreshPrctile = 1e-4;  % Percentile cutoff
```

## Configuring Thresholds in SMF

### Basic Configuration

Enable thresholding and set standard filters:

```matlab
SMF = smi_core.SingleMoleculeFitting();

% Enable thresholding
SMF.Thresholding.On = true;

% Precision filter
SMF.Thresholding.MaxXY_SE = 0.15;  % 15-20 nm for 100 nm pixel size

% Photon filter
SMF.Thresholding.MinPhotons = 150;

% Fit quality
SMF.Thresholding.MinPValue = 0.01;

% PSF sigma range (if fitting sigma)
SMF.Thresholding.MinPSFSigma = 0.9;
SMF.Thresholding.MaxPSFSigma = 1.7;

% Background limit
SMF.Thresholding.MaxBg = 40;  % photons/pixel
```

### Conservative vs. Permissive Settings

**Conservative (high quality, fewer localizations):**

```matlab
SMF.Thresholding.MaxXY_SE = 0.12;      % Stricter precision
SMF.Thresholding.MinPhotons = 200;     % More photons required
SMF.Thresholding.MinPValue = 0.05;     % Higher p-value threshold
SMF.Thresholding.MaxBg = 30;           % Lower background tolerance
```

**Permissive (more localizations, variable quality):**

```matlab
SMF.Thresholding.MaxXY_SE = 0.25;      % Relaxed precision
SMF.Thresholding.MinPhotons = 100;     % Accept dimmer spots
SMF.Thresholding.MinPValue = 0.01;     # Standard p-value
SMF.Thresholding.MaxBg = 50;           # Higher background tolerance
```

### Disabling Specific Filters

To disable a filter, set it to a non-restrictive value:

```matlab
% Disable background filter
SMF.Thresholding.MaxBg = Inf;

% Disable PSF sigma filter (if you're not fitting sigma anyway)
SMF.Thresholding.MinPSFSigma = 0;
SMF.Thresholding.MaxPSFSigma = Inf;

% Disable photon filter
SMF.Thresholding.MinPhotons = 0;
```

## Applying Thresholds During Localization

Thresholding happens automatically during the localization pipeline when `Thresholding.On = true`:

```matlab
% Configure thresholds
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.15;
SMF.Thresholding.MinPhotons = 150;
SMF.Thresholding.MinPValue = 0.01;

% Localize (thresholding happens internally)
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();

% SMD.ThreshFlag is populated but localizations not yet removed
% Check how many pass
passed = sum(SMD.ThreshFlag == 0);
total = length(SMD.ThreshFlag);
fprintf('Passed thresholds: %d / %d (%.1f%%)\n', passed, total, 100*passed/total);
```

The workflow automatically applies thresholds, but you can also apply them post-hoc.

## Post-Hoc Threshold Application

### Apply Existing ThreshFlag

If `ThreshFlag` already exists, apply it to remove rejected localizations:

```matlab
% Load results
load('Results.mat', 'SMD');

% Apply thresholds
SMD_filtered = smi_core.Threshold.applyThresh(SMD, 1);

fprintf('Before: %d localizations\n', length(SMD.X));
fprintf('After: %d localizations\n', length(SMD_filtered.X));
```

### Create New ThreshFlag

Recalculate thresholds with different criteria:

```matlab
% Load unfiltered results
load('Results.mat', 'SMD', 'SMF');

% Create new threshold criteria
SMF.Thresholding.MaxXY_SE = 0.12;  # Stricter than original
SMF.Thresholding.MinPhotons = 200;

% Create MinMax structure from SMF
MinMax = smi_core.Threshold.setMinMax(SMF);

% Set new ThreshFlag
THR = smi_core.Threshold();
[SMD, TFlag] = THR.setThreshFlag(SMD, MinMax);

% Apply thresholds
SMD_filtered = smi_core.Threshold.applyThresh(SMD, 1);
```

### Selective Filtering

Apply only specific filters by creating custom MinMax:

```matlab
% Filter only on precision and photons (ignore p-value)
MinMax = [];
MinMax.X_SE = [0, 0.15];
MinMax.Y_SE = [0, 0.15];
MinMax.Photons = [150, inf];

% Apply
THR = smi_core.Threshold();
[SMD, TFlag] = THR.setThreshFlag(SMD, MinMax);
SMD_filtered = smi_core.Threshold.applyThresh(SMD, 1);
```

## Understanding ThreshFlag

`ThreshFlag` is a 32-bit unsigned integer where each bit indicates failure of a specific criterion:

| Bit | Field | Meaning |
|-----|-------|---------|
| 1 | X | Position X out of range |
| 2 | Y | Position Y out of range |
| 3 | Z | Position Z out of range |
| 4 | Photons | Photon count outside range |
| 5 | Bg | Background outside range |
| 6 | PSFSigma | PSF sigma outside range |
| 7 | X_SE | X precision too poor |
| 8 | Y_SE | Y precision too poor |
| 9 | Z_SE | Z precision too poor |
| 10 | PValue | P-value too low |
| 11 | LogLikelihood | Log-likelihood outside range |

**Reading ThreshFlag:**

```matlab
% ThreshFlag = 0: Passed all filters
% ThreshFlag > 0: Failed at least one filter

% Example: ThreshFlag = 536
% Binary: 00000000000000000000001000011000
% Failed bits: 4 (Photons), 5 (Bg), 10 (PValue)
dec2bin(536, 32)  % Show binary representation
```

### Translating ThreshFlag

Convert numeric flags to human-readable descriptions:

```matlab
% Translate single flag
flag_value = 536;
[readable, hot_bits] = smi_core.Threshold.translateThreshFlag(flag_value);
fprintf('Flag %d failed: %s\n', flag_value, readable{1});
% Output: "Flag 536 failed: Photons Bg PValue"

% Translate all flags
[readable_all, hot_bits_all] = smi_core.Threshold.translateThreshFlag(SMD.ThreshFlag);

% Show first 10
for i = 1:min(10, length(readable_all))
    if SMD.ThreshFlag(i) > 0
        fprintf('Loc %d failed: %s\n', i, readable_all{i});
    end
end
```

### Analyzing Rejection Reasons

Determine which filters rejected the most localizations:

```matlab
% Count failures by criterion
THR = smi_core.Threshold();
Fields = THR.Fields;

for i = 1:length(Fields)
    bit_failures = sum(bitget(SMD.ThreshFlag, i));
    if bit_failures > 0
        fprintf('%s: %d failures (%.1f%%)\n', Fields{i}, ...
            bit_failures, 100*bit_failures/length(SMD.ThreshFlag));
    end
end
```

## Visualizing Threshold Effects

### Spatial Distribution of Rejected Fits

Visualize where rejected localizations occur:

```matlab
% Create threshold visualization
THR = smi_core.Threshold();
THR.Verbose = 1;
THR.rejectedLocalizations(SMD, 'RNM', '');
% R = individual reasons, N = number of reasons, M = major reasons

% This produces three plots:
% 1. Separate plots for each rejection criterion
% 2. Plot colored by number of rejection reasons
% 3. Plot colored by major rejection reason
```

### Before/After Comparison

Compare localization maps before and after thresholding:

```matlab
% Before thresholding
figure('Position', [100, 100, 1200, 500]);

subplot(1,2,1);
plot(SMD.X, SMD.Y, 'k.', 'MarkerSize', 1);
axis equal; axis tight;
title(sprintf('Before Thresholding: %d localizations', length(SMD.X)));
xlabel('X (pixels)'); ylabel('Y (pixels)');

% After thresholding
SMD_filt = smi_core.Threshold.applyThresh(SMD, 0);

subplot(1,2,2);
plot(SMD_filt.X, SMD_filt.Y, 'k.', 'MarkerSize', 1);
axis equal; axis tight;
title(sprintf('After Thresholding: %d localizations', length(SMD_filt.X)));
xlabel('X (pixels)'); ylabel('Y (pixels)');
```

### Precision Distribution

Compare precision before and after filtering:

```matlab
figure;
histogram(SMD.X_SE * SMD.PixelSize * 1000, 50, 'DisplayName', 'Before');
hold on;
histogram(SMD_filt.X_SE * SMD_filt.PixelSize * 1000, 50, 'DisplayName', 'After');
xlabel('Localization Precision (nm)');
ylabel('Count');
title('Precision Distribution');
legend;
hold off;

fprintf('Median precision before: %.1f nm\n', ...
    median(SMD.X_SE) * SMD.PixelSize * 1000);
fprintf('Median precision after: %.1f nm\n', ...
    median(SMD_filt.X_SE) * SMD_filt.PixelSize * 1000);
```

## Guidelines for Choosing Thresholds

### Start with Data Characteristics

Examine your data distributions first:

```matlab
% Check precision distribution
figure('Position', [100, 100, 1200, 800]);

subplot(2,3,1);
histogram(SMD.X_SE, 50);
xlabel('X Precision (pixels)'); ylabel('Count');
title('X Precision Distribution');
xline(0.15, 'r--', 'LineWidth', 2, 'Label', 'Typical threshold');

subplot(2,3,2);
histogram(SMD.Photons, 50);
xlabel('Photons'); ylabel('Count');
title('Photon Distribution');
xline(150, 'r--', 'LineWidth', 2, 'Label', 'Typical threshold');

subplot(2,3,3);
histogram(SMD.Bg, 50);
xlabel('Background (photons/pixel)'); ylabel('Count');
title('Background Distribution');
xline(40, 'r--', 'LineWidth', 2, 'Label', 'Typical threshold');

subplot(2,3,4);
histogram(SMD.PValue, 50);
xlabel('P-value'); ylabel('Count');
title('P-value Distribution');
xline(0.01, 'r--', 'LineWidth', 2, 'Label', 'Typical threshold');

subplot(2,3,5);
if isfield(SMD, 'PSFSigma') && ~isempty(SMD.PSFSigma)
    histogram(SMD.PSFSigma, 50);
    xlabel('PSF Sigma (pixels)'); ylabel('Count');
    title('PSF Sigma Distribution');
end

subplot(2,3,6);
% SNR calculation
SNR = SMD.Photons ./ sqrt(SMD.Photons + SMD.Bg * pi * 1.3^2);
histogram(SNR, 50);
xlabel('SNR'); ylabel('Count');
title('Signal-to-Noise Ratio');
```

### Empirical Approach

1. **Start permissive**: Begin with loose thresholds to see most data
2. **Identify outliers**: Look at extreme values in distributions
3. **Set thresholds at natural breaks**: Place cutoffs where distributions show clear separation
4. **Iterate**: Tighten gradually and check localization quality

```matlab
% Example iterative approach
thresholds_to_test = [0.1, 0.15, 0.2, 0.25];

for i = 1:length(thresholds_to_test)
    SMF.Thresholding.MaxXY_SE = thresholds_to_test(i);
    MinMax = smi_core.Threshold.setMinMax(SMF);
    THR = smi_core.Threshold();
    [SMD_temp, ~] = THR.setThreshFlag(SMD, MinMax);
    SMD_temp = smi_core.Threshold.applyThresh(SMD_temp, 0);

    fprintf('Threshold %.2f: %d localizations (%.1f%% kept)\n', ...
        thresholds_to_test(i), length(SMD_temp.X), ...
        100*length(SMD_temp.X)/length(SMD.X));
end
```

### Application-Specific Guidelines

**High-precision STORM/PAINT:**
- Precision: MaxXY_SE = 0.12-0.15 pixels (10-15 nm)
- Photons: MinPhotons = 200-300
- P-value: MinPValue = 0.05

**Standard SMLM:**
- Precision: MaxXY_SE = 0.15-0.20 pixels (15-20 nm)
- Photons: MinPhotons = 150-200
- P-value: MinPValue = 0.01

**Challenging conditions (low SNR):**
- Precision: MaxXY_SE = 0.25 pixels (25 nm)
- Photons: MinPhotons = 100
- P-value: MinPValue = 0.01
- Consider AutoThreshLogL instead of fixed p-value

**3D astigmatism:**
- XY Precision: MaxXY_SE = 0.15-0.20 pixels
- Z Precision: MaxZ_SE = 0.03-0.05 microns
- PSF sigma filtering crucial to remove out-of-range Z

## Complete Example

```matlab
% ========== Load Data ==========
load('Results.mat', 'SMD', 'SMF');
fprintf('Loaded %d localizations\n', length(SMD.X));

% ========== Examine Distributions ==========
fprintf('\n=== Quality Metrics ===\n');
fprintf('Precision: %.3f ± %.3f pixels\n', mean(SMD.X_SE), std(SMD.X_SE));
fprintf('Photons: %.0f ± %.0f\n', mean(SMD.Photons), std(SMD.Photons));
fprintf('Background: %.1f ± %.1f photons/pixel\n', mean(SMD.Bg), std(SMD.Bg));
fprintf('P-value: %.3f (median)\n', median(SMD.PValue));

% ========== Configure Thresholds ==========
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.15;    % 15 nm at 100 nm/pixel
SMF.Thresholding.MinPhotons = 150;
SMF.Thresholding.MinPValue = 0.01;
SMF.Thresholding.MinPSFSigma = 0.9;
SMF.Thresholding.MaxPSFSigma = 1.7;
SMF.Thresholding.MaxBg = 40;

% ========== Apply Thresholds ==========
MinMax = smi_core.Threshold.setMinMax(SMF);
THR = smi_core.Threshold();
THR.Verbose = 1;

[SMD, TFlag] = THR.setThreshFlag(SMD, MinMax);
SMD_filtered = smi_core.Threshold.applyThresh(SMD, 1);

% ========== Analyze Results ==========
fprintf('\n=== Thresholding Results ===\n');
fprintf('Original: %d localizations\n', length(SMD.X));
fprintf('Filtered: %d localizations (%.1f%% kept)\n', ...
    length(SMD_filtered.X), 100*length(SMD_filtered.X)/length(SMD.X));

% Individual failure rates
fprintf('\n=== Failure Breakdown ===\n');
fprintf('X_SE failures: %d\n', sum(bitget(SMD.ThreshFlag, 7)));
fprintf('Y_SE failures: %d\n', sum(bitget(SMD.ThreshFlag, 8)));
fprintf('Photon failures: %d\n', sum(bitget(SMD.ThreshFlag, 4)));
fprintf('Bg failures: %d\n', sum(bitget(SMD.ThreshFlag, 5)));
fprintf('PValue failures: %d\n', sum(bitget(SMD.ThreshFlag, 10)));
fprintf('PSFSigma failures: %d\n', sum(bitget(SMD.ThreshFlag, 6)));

% ========== Visualize ==========
% Spatial distribution of rejections
THR.rejectedLocalizations(SMD, 'M', '');  % Major reasons plot

% Quality improvement
figure('Position', [100, 100, 1200, 400]);

subplot(1,3,1);
plot(SMD.X, SMD.Y, '.', 'MarkerSize', 1, 'Color', [0.7 0.7 0.7]);
hold on;
plot(SMD_filtered.X, SMD_filtered.Y, 'k.', 'MarkerSize', 1);
axis equal; axis tight;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Before (gray) vs After (black)');
legend('Rejected', 'Kept');

subplot(1,3,2);
histogram(SMD.X_SE * SMD.PixelSize * 1000, 50, 'FaceAlpha', 0.5);
hold on;
histogram(SMD_filtered.X_SE * SMD_filtered.PixelSize * 1000, 50, 'FaceAlpha', 0.5);
xlabel('Precision (nm)'); ylabel('Count');
title('Precision Distribution');
legend('Before', 'After');

subplot(1,3,3);
histogram(SMD.Photons, 50, 'FaceAlpha', 0.5);
hold on;
histogram(SMD_filtered.Photons, 50, 'FaceAlpha', 0.5);
xlabel('Photons'); ylabel('Count');
title('Photon Distribution');
legend('Before', 'After');

% ========== Save ==========
save('Results_filtered.mat', 'SMD_filtered', 'SMF');
fprintf('\nFiltered results saved to Results_filtered.mat\n');
```

## See Also

- [Localize Molecules](localize-molecules.md) - Where thresholding happens
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete analysis pipeline
- [SMD Structure](../core-concepts/smd-structure.md) - Understanding localization data
- [SMF Structure](../core-concepts/smf-structure.md) - All threshold parameters
