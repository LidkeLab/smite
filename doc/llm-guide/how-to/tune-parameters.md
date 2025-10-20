---
title: "How to Tune SMF Parameters for Optimal Results"
category: "how-to"
level: "intermediate"
tags: ["parameters", "optimization", "tuning", "quality", "smf", "troubleshooting"]
prerequisites: ["create-smf.md", "localize-molecules.md", "threshold-results.md"]
related: ["simulate-data.md", "../workflows/smlm-analysis.md", "../core-concepts/smf-structure.md"]
summary: "Systematic approach to optimizing SMF parameters for best localization quality using test data, evaluation metrics, and iterative refinement"
estimated_time: "25 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# How to Tune SMF Parameters for Optimal Results

## Purpose

Getting excellent localization results requires carefully tuned parameters. Default SMF values provide a starting point, but optimal performance demands systematic parameter optimization tailored to your specific imaging conditions. This guide presents a methodical approach to parameter tuning, from initial testing through iterative refinement, with emphasis on understanding parameter interactions and evaluating result quality.

## Prerequisites

- Ability to [create and configure SMF](create-smf.md)
- Understanding of [localization process](localize-molecules.md)
- Familiarity with [quality filtering](threshold-results.md)
- Basic knowledge of [simulated data](simulate-data.md) (helpful but not required)

## Overview

Parameter tuning follows a systematic workflow:

1. **Prepare Test Data**: Use simulation or representative real data subset
2. **Establish Baseline**: Run with default parameters
3. **Evaluate Quality**: Quantify localization performance
4. **Identify Issues**: Diagnose specific problems
5. **Adjust Parameters**: Make targeted improvements
6. **Iterate**: Repeat evaluation and adjustment
7. **Validate**: Test on independent data

**Key insight**: Parameters are interdependent. Changing one often requires retuning others. The goal is finding the optimal parameter combination for your specific imaging conditions.

## Stage 1: Prepare Test Data

### Using Simulated Data

Simulated data provides ground truth for quantitative validation:

```matlab
% Create realistic test data matching your imaging conditions
obj = smi_sim.SimSMLM();

% Match your imaging parameters
obj.SZ = 128;                  % Image size
obj.NFrames = 200;             % Enough for statistics
obj.EmissionRate = 1000;       % Typical photon count
obj.Bg = 10;                   % Estimated background
obj.PSFSigma = 1.3;            % Known PSF width

% Match your photophysics (DNA-PAINT example)
obj.K_OnToOff = 1;             % Turn-off rate
obj.K_OffToOn = 0.01;          % Turn-on rate
obj.K_OnToBleach = 0.005;      % Bleaching rate
obj.LabelingEfficiency = 1;

% Generate structured pattern for visual assessment
NWings = 16;
obj.simStar(NWings);
obj.applyLabelEffic();
obj.genBlinks('Equib');
[~, testData] = obj.genImageStack();

% Store ground truth
SMD_true = obj.SMD_Model;
ConnectID = obj.SMD_Model.ConnectID;

fprintf('Test data: %d frames, %d localizations expected\n', ...
    obj.NFrames, length(SMD_true.X));
```

**Why simulation?**
- Known ground truth enables quantitative accuracy measurement
- Control over all parameters
- Fast iteration without experimental overhead
- Reproducible baseline for comparison

### Using Real Data Subset

Alternatively, use a small representative subset of real data:

```matlab
% Load full dataset
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'experiment.h5'};
SMF.Data.PixelSize = 0.108;
SMF.Data.FrameRate = 100;
SMF.Data.CameraGain = 2.5;
SMF.Data.CameraOffset = 100;

LD_loader = smi_core.LoadData();
[~, fullData, SMF] = LD_loader.loadRawData(SMF, 1);

% Extract representative subset (e.g., 100 frames)
testData = fullData(:, :, 1:100);

fprintf('Test subset: %d frames from real data\n', size(testData, 3));

% Examine typical frame
figure;
imagesc(testData(:,:,50));
axis image; colormap gray; colorbar;
title('Representative Frame');
```

**Advantages of real data:**
- Tests actual experimental conditions
- Reveals real noise characteristics
- Validates pipeline on target data

**Disadvantages:**
- No ground truth for quantitative validation
- Harder to isolate parameter effects
- Slower iteration

## Stage 2: Establish Baseline

Run localization with sensible defaults to establish performance baseline:

```matlab
% Configure SMF with reasonable starting parameters
SMF = smi_core.SingleMoleculeFitting();

% Camera parameters (match your setup)
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 1;        % Data already in photons (for simulated)
SMF.Data.CameraOffset = 0;
SMF.Data.PixelSize = 0.108;
SMF.Data.FrameRate = 100;

% Detection parameters - conservative defaults
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 200;
SMF.BoxFinding.BoxOverlap = 2;

% Fitting parameters - standard 2D
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';
SMF.Fitting.Iterations = 20;

% Thresholding - moderate filters
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.2;
SMF.Thresholding.MinPhotons = 100;
SMF.Thresholding.MinPValue = 0.01;

% Localize with baseline parameters
LD = smi_core.LocalizeData(testData, SMF);
LD.Verbose = 1;
SMD_baseline = LD.genLocalizations();

% Save baseline for comparison
save('baseline_results.mat', 'SMD_baseline', 'SMF');

fprintf('\n=== Baseline Results ===\n');
fprintf('Localizations: %d\n', length(SMD_baseline.X));
fprintf('Median precision: %.1f nm\n', ...
    median(SMD_baseline.X_SE) * SMF.Data.PixelSize * 1000);
fprintf('Median photons: %.0f\n', median(SMD_baseline.Photons));
```

**Document baseline performance** - you'll compare all subsequent attempts against this.

## Stage 3: Evaluate Quality

### Key Quality Metrics

Quantify localization quality with these metrics:

**1. Localization Precision**

```matlab
% X/Y precision (standard error)
precision_x_nm = SMD_baseline.X_SE * SMF.Data.PixelSize * 1000;  % nm
precision_y_nm = SMD_baseline.Y_SE * SMF.Data.PixelSize * 1000;  % nm

fprintf('Precision (X): %.1f +/- %.1f nm (mean +/- std)\n', ...
    mean(precision_x_nm), std(precision_x_nm));
fprintf('Precision (Y): %.1f +/- %.1f nm (mean +/- std)\n', ...
    mean(precision_y_nm), std(precision_y_nm));

% Target: < 20 nm for high-quality SMLM
% Target: < 50 nm for moderate quality

% Visualize distribution
figure;
histogram(precision_x_nm, 30);
xlabel('Localization Precision (nm)');
ylabel('Count');
title('Precision Distribution');
xline(20, 'r--', 'LineWidth', 2, 'Label', '20 nm target');
```

**2. Detection Rate (with ground truth)**

```matlab
% Match localizations to ground truth
if exist('SMD_true', 'var')
    max_dist = 2;  % pixels
    matched = 0;
    errors = [];

    for i = 1:length(SMD_true.X)
        frame = SMD_true.FrameNum(i);
        same_frame = (SMD_baseline.FrameNum == frame);

        if sum(same_frame) > 0
            dx = SMD_baseline.X(same_frame) - SMD_true.X(i);
            dy = SMD_baseline.Y(same_frame) - SMD_true.Y(i);
            dist = sqrt(dx.^2 + dy.^2);
            [min_dist, ~] = min(dist);

            if min_dist < max_dist
                matched = matched + 1;
                errors = [errors; min_dist];
            end
        end
    end

    detection_rate = 100 * matched / length(SMD_true.X);
    fprintf('Detection rate: %.1f%% (%d/%d)\n', ...
        detection_rate, matched, length(SMD_true.X));
    fprintf('Mean localization error: %.3f pixels (%.1f nm)\n', ...
        mean(errors), mean(errors) * SMF.Data.PixelSize * 1000);

    % Target detection rate: > 90%
    % Target localization error: < 0.2 pixels (< 20 nm)
end
```

**3. False Positive Rate**

```matlab
% Count excess localizations (with ground truth)
if exist('SMD_true', 'var')
    n_detected = length(SMD_baseline.X);
    n_true = length(SMD_true.X);
    n_false_positives = n_detected - matched;
    false_positive_rate = 100 * n_false_positives / n_detected;

    fprintf('False positives: %d (%.1f%% of detections)\n', ...
        n_false_positives, false_positive_rate);

    % Target: < 5% false positives
end
```

**4. Signal-to-Noise Ratio**

```matlab
% Estimate SNR from photons and background
SNR = SMD_baseline.Photons ./ sqrt(SMD_baseline.Photons + ...
    SMD_baseline.Bg .* pi .* SMF.Fitting.PSFSigma^2);

fprintf('SNR: %.1f +/- %.1f (mean +/- std)\n', mean(SNR), std(SNR));
fprintf('SNR range: %.1f to %.1f\n', min(SNR), max(SNR));

% Target SNR: > 5 for reliable localization
% Target SNR: > 10 for high precision

figure;
histogram(SNR, 30);
xlabel('Signal-to-Noise Ratio');
ylabel('Count');
title('SNR Distribution');
xline(5, 'r--', 'Label', 'Min recommended');
```

**5. Fit Quality (P-value)**

```matlab
fprintf('P-value: %.3f +/- %.3f (mean +/- std)\n', ...
    mean(SMD_baseline.PValue), std(SMD_baseline.PValue));
fprintf('P-value < 0.01: %.1f%%\n', ...
    100 * sum(SMD_baseline.PValue < 0.01) / length(SMD_baseline.PValue));

% P-values should be roughly uniform for good fits
% Excess low p-values suggests model mismatch

figure;
histogram(SMD_baseline.PValue, 30);
xlabel('Fit P-value');
ylabel('Count');
title('Fit Quality Distribution');
```

**6. Photon Recovery**

```matlab
fprintf('Photons: %.0f +/- %.0f (mean +/- std)\n', ...
    mean(SMD_baseline.Photons), std(SMD_baseline.Photons));

% Compare to known photons if using simulation
if exist('obj', 'var')
    expected_photons = obj.EmissionRate;
    bias = mean(SMD_baseline.Photons) - expected_photons;
    fprintf('Expected photons: %.0f\n', expected_photons);
    fprintf('Photon bias: %.1f (%.1f%%)\n', ...
        bias, 100 * bias / expected_photons);
end

figure;
histogram(SMD_baseline.Photons, 50);
xlabel('Photons');
ylabel('Count');
title('Photon Distribution');
```

### Create Quality Report

```matlab
% Generate comprehensive quality report
fprintf('\n========== QUALITY REPORT ==========\n');
fprintf('Parameters:\n');
fprintf('  BoxSize: %d\n', SMF.BoxFinding.BoxSize);
fprintf('  MinPhotons (detection): %.0f\n', SMF.BoxFinding.MinPhotons);
fprintf('  PSFSigma: %.2f\n', SMF.Fitting.PSFSigma);
fprintf('  MaxXY_SE (threshold): %.2f\n', SMF.Thresholding.MaxXY_SE);
fprintf('\nResults:\n');
fprintf('  Localizations: %d\n', length(SMD_baseline.X));
fprintf('  Median precision: %.1f nm\n', ...
    median(SMD_baseline.X_SE) * SMF.Data.PixelSize * 1000);
fprintf('  Mean SNR: %.1f\n', mean(SNR));
if exist('detection_rate', 'var')
    fprintf('  Detection rate: %.1f%%\n', detection_rate);
    fprintf('  False positive rate: %.1f%%\n', false_positive_rate);
end
fprintf('====================================\n\n');
```

## Stage 4: Identify Issues

Common problems and their symptoms:

### Issue 1: Poor Precision (> 25 nm)

**Symptoms:**
- Median precision > 0.2 pixels
- Blurry reconstructions
- Wide precision distribution

**Likely causes:**
- Insufficient photons
- High background
- Wrong PSFSigma
- Overlapping emitters

### Issue 2: Low Detection Rate (< 85%)

**Symptoms:**
- Missing expected localizations
- Sparse reconstructions
- Ground truth comparison shows missed emitters

**Likely causes:**
- MinPhotons too high
- BoxSize too small
- Wrong CameraGain/Offset
- Thresholding too strict

### Issue 3: High False Positive Rate (> 10%)

**Symptoms:**
- More localizations than expected
- Background noise detected as molecules
- Random scattered localizations

**Likely causes:**
- MinPhotons too low
- Thresholding too permissive
- Background not properly subtracted

### Issue 4: Biased Photon Estimates

**Symptoms:**
- Systematic over/under-estimation
- Photon distribution shifted from expected

**Likely causes:**
- Wrong CameraGain
- Wrong PSFSigma
- BoxSize too small (loses photons)

### Issue 5: Poor Fit Quality (Low P-values)

**Symptoms:**
- Many localizations with PValue < 0.01
- High LogLikelihood rejection rate

**Likely causes:**
- Wrong PSFSigma
- Overlapping molecules
- Non-Gaussian PSF
- Residual background structure

## Stage 5: Parameter Tuning

### Key Parameters and Tuning Strategy

The most critical parameters for optimization:

**1. PSFSigma** - Single most important parameter

```matlab
% Determine correct PSFSigma empirically
SMF_test = SMF;  % Copy baseline
SMF_test.Fitting.FitType = 'XYNBS';  % Fit sigma
SMF_test.Fitting.PSFSigma = 1.5;     % Initial guess

LD = smi_core.LocalizeData(testData, SMF_test);
SMD_test = LD.genLocalizations();

% Examine fitted sigmas
fitted_sigma = median(SMD_test.PSFSigma);
fprintf('Fitted PSF sigma: %.3f pixels\n', fitted_sigma);

% Visualize distribution
figure;
histogram(SMD_test.PSFSigma, 30);
xlabel('Fitted PSF Sigma (pixels)');
ylabel('Count');
title('PSF Sigma Distribution');

% Use median fitted sigma for subsequent analyses
SMF.Fitting.PSFSigma = fitted_sigma;
SMF.Fitting.FitType = 'XYNB';  % Fix sigma at this value

% Set reasonable thresholds for PSFSigma filter
SMF.Thresholding.MinPSFSigma = fitted_sigma * 0.7;
SMF.Thresholding.MaxPSFSigma = fitted_sigma * 1.3;
```

**Impact**: Wrong PSFSigma causes biased positions, photons, and background estimates. This cascades into poor precision and fit quality.

**2. BoxSize** - Balances signal capture vs. background

```matlab
% Test different box sizes
box_sizes = [5, 7, 9, 11];
results = struct();

for i = 1:length(box_sizes)
    SMF_test = SMF;
    SMF_test.BoxFinding.BoxSize = box_sizes(i);

    LD = smi_core.LocalizeData(testData, SMF_test);
    SMD_test = LD.genLocalizations();

    results(i).BoxSize = box_sizes(i);
    results(i).NLocs = length(SMD_test.X);
    results(i).MedianPrecision = median(SMD_test.X_SE);
    results(i).MedianPhotons = median(SMD_test.Photons);
    results(i).MedianBg = median(SMD_test.Bg);

    fprintf('BoxSize=%d: %d locs, %.1f nm precision, %.0f photons\n', ...
        box_sizes(i), results(i).NLocs, ...
        results(i).MedianPrecision * SMF.Data.PixelSize * 1000, ...
        results(i).MedianPhotons);
end

% Choose BoxSize that minimizes precision while maintaining photon recovery
% Rule of thumb: BoxSize ≈ 4-5 × PSFSigma
optimal_box_size = ceil(4.5 * SMF.Fitting.PSFSigma);
fprintf('Recommended BoxSize: %d\n', optimal_box_size);
SMF.BoxFinding.BoxSize = optimal_box_size;
```

**Impact**: Too small loses photons; too large includes excess background. Both degrade precision.

**3. MinPhotons (Detection)** - Balances sensitivity vs. false positives

```matlab
% Sweep MinPhotons threshold
min_photon_values = [100, 150, 200, 250, 300];
results = struct();

for i = 1:length(min_photon_values)
    SMF_test = SMF;
    SMF_test.BoxFinding.MinPhotons = min_photon_values(i);

    LD = smi_core.LocalizeData(testData, SMF_test);
    SMD_test = LD.genLocalizations();

    results(i).MinPhotons = min_photon_values(i);
    results(i).NDetected = length(SMD_test.X);

    % Calculate detection rate if ground truth available
    if exist('SMD_true', 'var')
        matched = 0;
        for j = 1:length(SMD_true.X)
            frame = SMD_true.FrameNum(j);
            same_frame = (SMD_test.FrameNum == frame);
            if sum(same_frame) > 0
                dist = sqrt((SMD_test.X(same_frame) - SMD_true.X(j)).^2 + ...
                           (SMD_test.Y(same_frame) - SMD_true.Y(j)).^2);
                if min(dist) < 2
                    matched = matched + 1;
                end
            end
        end
        results(i).DetectionRate = 100 * matched / length(SMD_true.X);
        fprintf('MinPhotons=%d: %d detected, %.1f%% detection rate\n', ...
            min_photon_values(i), results(i).NDetected, results(i).DetectionRate);
    else
        fprintf('MinPhotons=%d: %d detected\n', ...
            min_photon_values(i), results(i).NDetected);
    end
end

% Choose MinPhotons that maximizes detection rate while minimizing noise
% Typically: 50-70% of median photon count
recommended_min_photons = 0.6 * median(SMD_baseline.Photons);
fprintf('Recommended MinPhotons: %.0f\n', recommended_min_photons);
SMF.BoxFinding.MinPhotons = recommended_min_photons;
```

**Impact**: Too low detects noise; too high misses dim emitters.

**4. Thresholding Parameters** - Quality control

```matlab
% Tune precision threshold
% Goal: Remove worst ~10-20% of localizations by precision

precision_percentiles = prctile(SMD_baseline.X_SE, [75, 80, 85, 90, 95]);
fprintf('Precision percentiles (75th-95th): %.3f, %.3f, %.3f, %.3f, %.3f\n', ...
    precision_percentiles);

% Choose threshold at ~85th percentile
SMF.Thresholding.MaxXY_SE = precision_percentiles(3);
fprintf('Set MaxXY_SE = %.3f (85th percentile)\n', SMF.Thresholding.MaxXY_SE);

% Tune photon threshold (post-fit)
% Should be lower than detection threshold
photon_percentiles = prctile(SMD_baseline.Photons, [5, 10, 15, 20]);
fprintf('Photon percentiles (5th-20th): %.0f, %.0f, %.0f, %.0f\n', ...
    photon_percentiles);

% Remove bottom ~10%
SMF.Thresholding.MinPhotons = photon_percentiles(2);
fprintf('Set MinPhotons = %.0f (10th percentile)\n', SMF.Thresholding.MinPhotons);

% P-value threshold
% Standard: 0.01 (rejects ~5-10% typically)
SMF.Thresholding.MinPValue = 0.01;
```

**Impact**: Filters poor-quality localizations while retaining good ones.

### Parameter Interdependencies

Parameters don't act independently:

**PSFSigma affects BoxSize:**
```matlab
% Always re-evaluate BoxSize after changing PSFSigma
SMF.BoxFinding.BoxSize = ceil(4.5 * SMF.Fitting.PSFSigma);
```

**BoxSize affects photon estimates:**
```matlab
% Larger boxes recover more photons but include more background
% May need to adjust MinPhotons when changing BoxSize
```

**Thresholding depends on all upstream parameters:**
```matlab
% After changing detection/fitting parameters, re-tune thresholds
% based on new localization distributions
```

## Stage 6: Iterate and Refine

Systematic iteration procedure:

```matlab
% Iteration loop
iteration = 1;
converged = false;

while ~converged && iteration <= 5
    fprintf('\n=== Iteration %d ===\n', iteration);

    % Run localization with current parameters
    LD = smi_core.LocalizeData(testData, SMF);
    LD.Verbose = 0;
    SMD = LD.genLocalizations();

    % Evaluate quality
    precision_nm = median(SMD.X_SE) * SMF.Data.PixelSize * 1000;
    n_locs = length(SMD.X);

    fprintf('Results: %d locs, %.1f nm precision\n', n_locs, precision_nm);

    % Check convergence criteria
    if iteration > 1
        precision_change = abs(precision_nm - prev_precision) / prev_precision;
        locs_change = abs(n_locs - prev_n_locs) / prev_n_locs;

        fprintf('Changes: %.1f%% precision, %.1f%% localizations\n', ...
            100*precision_change, 100*locs_change);

        if precision_change < 0.05 && locs_change < 0.05
            fprintf('Converged (< 5%% change)\n');
            converged = true;
        end
    end

    % Store for comparison
    prev_precision = precision_nm;
    prev_n_locs = n_locs;

    % Adjust parameters based on evaluation
    if ~converged
        % Example adjustment logic
        if precision_nm > 25  % Poor precision
            % Try stricter thresholding
            SMF.Thresholding.MaxXY_SE = SMF.Thresholding.MaxXY_SE * 0.9;
            fprintf('Adjusting: Stricter precision threshold\n');
        end

        if n_locs < 0.8 * length(SMD_true.X)  % Low detection
            % Relax detection threshold
            SMF.BoxFinding.MinPhotons = SMF.BoxFinding.MinPhotons * 0.9;
            fprintf('Adjusting: Lower detection threshold\n');
        end
    end

    iteration = iteration + 1;
end

fprintf('\nFinal optimized parameters saved\n');
save('optimized_SMF.mat', 'SMF');
```

## Stage 7: Validate

Test optimized parameters on independent data:

```matlab
% Load validation dataset (different from tuning data)
[~, validationData, ~] = LD_loader.loadRawData(SMF, 2);

% Apply optimized parameters
LD = smi_core.LocalizeData(validationData, SMF);
SMD_validation = LD.genLocalizations();

% Compare metrics to tuning results
fprintf('\n=== Validation Results ===\n');
fprintf('Tuning precision: %.1f nm\n', ...
    median(SMD.X_SE) * SMF.Data.PixelSize * 1000);
fprintf('Validation precision: %.1f nm\n', ...
    median(SMD_validation.X_SE) * SMF.Data.PixelSize * 1000);

fprintf('\nTuning localizations: %d\n', length(SMD.X));
fprintf('Validation localizations: %d\n', length(SMD_validation.X));

% Parameters are well-tuned if validation metrics are similar to tuning metrics
```

## Complete Tuning Example

End-to-end parameter optimization:

```matlab
%% Setup: Create test data
obj = smi_sim.SimSMLM();
obj.SZ = 128;
obj.NFrames = 200;
obj.EmissionRate = 1000;
obj.Bg = 10;
obj.PSFSigma = 1.3;
obj.K_OnToOff = 1;
obj.K_OffToOn = 0.01;
obj.K_OnToBleach = 0.005;
obj.simStar(16);
obj.applyLabelEffic();
obj.genBlinks('Equib');
[~, testData] = obj.genImageStack();
SMD_true = obj.SMD_Model;

%% Stage 1: Baseline
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraGain = 1;
SMF.Data.CameraOffset = 0;
SMF.Data.PixelSize = 0.108;
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 200;
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

LD = smi_core.LocalizeData(testData, SMF);
SMD_baseline = LD.genLocalizations();

fprintf('Baseline: %d locs, %.1f nm precision\n', ...
    length(SMD_baseline.X), ...
    median(SMD_baseline.X_SE) * SMF.Data.PixelSize * 1000);

%% Stage 2: Optimize PSFSigma
SMF.Fitting.FitType = 'XYNBS';
LD = smi_core.LocalizeData(testData, SMF);
SMD_test = LD.genLocalizations();

optimal_sigma = median(SMD_test.PSFSigma);
fprintf('Optimal PSFSigma: %.3f\n', optimal_sigma);

SMF.Fitting.PSFSigma = optimal_sigma;
SMF.Fitting.FitType = 'XYNB';

%% Stage 3: Optimize BoxSize
SMF.BoxFinding.BoxSize = ceil(4.5 * optimal_sigma);
fprintf('Optimal BoxSize: %d\n', SMF.BoxFinding.BoxSize);

%% Stage 4: Optimize MinPhotons
LD = smi_core.LocalizeData(testData, SMF);
SMD_test = LD.genLocalizations();

optimal_min_photons = 0.6 * median(SMD_test.Photons);
SMF.BoxFinding.MinPhotons = optimal_min_photons;
fprintf('Optimal MinPhotons: %.0f\n', optimal_min_photons);

%% Stage 5: Optimize Thresholding
LD = smi_core.LocalizeData(testData, SMF);
SMD_test = LD.genLocalizations();

SMF.Thresholding.MaxXY_SE = prctile(SMD_test.X_SE, 85);
SMF.Thresholding.MinPhotons = prctile(SMD_test.Photons, 10);
fprintf('Precision threshold: %.3f pixels\n', SMF.Thresholding.MaxXY_SE);

%% Final optimized run
LD = smi_core.LocalizeData(testData, SMF);
SMD_final = LD.genLocalizations();

fprintf('\n=== Final Optimized Results ===\n');
fprintf('Localizations: %d\n', length(SMD_final.X));
fprintf('Precision: %.1f nm\n', ...
    median(SMD_final.X_SE) * SMF.Data.PixelSize * 1000);
fprintf('Photons: %.0f\n', median(SMD_final.Photons));

% Save optimized parameters
save('optimized_parameters.mat', 'SMF');
```

## Tips and Best Practices

### Start with Good Defaults

Don't start from scratch:

```matlab
% Use experiment-specific templates from create-smf.md
% DNA-PAINT, dSTORM, SPT all have reasonable starting points
```

### One Parameter at a Time

```matlab
% Change parameters systematically
% Document each change and its effect
% Don't change multiple parameters simultaneously
```

### Use Visualization

```matlab
% Visual inspection complements metrics
LD.Verbose = 3;  % Show overlays during localization

% Reconstruct image to assess quality
smi_vis.GenerateImages.super_resolution_movie(SMD_final, SMF, ...
    'PixelSize', 10);  % 10 nm pixels
```

### Keep Records

```matlab
% Log all parameter sets tested
diary('tuning_log.txt');
% ... run experiments ...
diary off;

% Save parameter history
parameter_history{iteration} = SMF;
save('parameter_history.mat', 'parameter_history');
```

### Know When to Stop

Parameter tuning has diminishing returns. Stop when:
- Precision < 20 nm (excellent)
- Detection rate > 90%
- False positive rate < 5%
- Metrics stable across datasets

## Troubleshooting

### Optimization Not Converging

Try:
- Smaller parameter adjustments
- Check for data quality issues
- Verify camera calibration
- Test on simpler simulated data first

### Parameters Don't Transfer to Real Data

Causes:
- Camera calibration incorrect
- Simulation doesn't match reality
- Background structure in real data
- PSF shape differs from Gaussian

Solutions:
- Verify CameraGain and CameraOffset empirically
- Add realistic background to simulation
- Consider spatially-varying PSF if necessary

### Results Vary Between Datasets

Indicates:
- Parameters too sensitive
- Data quality varies
- Need robust thresholding (AutoThreshLogL)

## See Also

- [Create SMF](create-smf.md) - Starting parameter templates
- [Localize Molecules](localize-molecules.md) - Understanding the pipeline
- [Threshold Results](threshold-results.md) - Quality filtering details
- [Simulate Data](simulate-data.md) - Creating test data
- [SMF Structure](../core-concepts/smf-structure.md) - Complete parameter reference
- [SMLM Workflow](../workflows/smlm-analysis.md) - Full analysis pipeline
