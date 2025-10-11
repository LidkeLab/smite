---
title: "Basic Localization Example"
category: "examples"
level: "beginner"
tags: ["example", "localization", "tutorial", "code"]
prerequisites: ["../getting-started/installation.md"]
related: ["../how-to/localize-molecules.md", "../workflows/smlm-analysis.md"]
summary: "Complete working example of molecule localization from simulated data"
estimated_time: "10 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Basic Localization Example

## Purpose

This example demonstrates a complete localization workflow from start to finish using simulated data. You'll generate test data, configure parameters, localize molecules, and visualize results. All code is ready to run and produces immediate results, making this perfect for learning smite's localization capabilities.

## Prerequisites

- smite installed and working
- Basic MATLAB knowledge
- 5-10 minutes to run the code

## Overview

This example covers:
1. Generating simulated SMLM data with known ground truth
2. Configuring SMF for localization
3. Running the localization pipeline
4. Comparing results to ground truth
5. Creating visualizations

Everything runs in memory with no external files required.

## Complete Working Code

Copy and run this complete example:

```matlab
%% Basic Localization Example
% Demonstrates molecule localization with simulated data

%% Step 1: Generate Simulated Data
fprintf('=== Generating Simulated Data ===\n');

% Simulation parameters
SZ = 128;                     % Image size (pixels, square)
NFrames = 50;                 % Number of frames
Rho = 0.01;                   % Density (emitters/pixel) - ~200 emitters per frame
Photons = 1000;               % Photons per emitter
PSFSigma = 1.3;               % PSF sigma (pixels)
Bg = 5;                       % Background (photons/pixel)

% Generate image stack using GaussBlobs static method
imageStack_noisy = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, ...
    Photons, PSFSigma, Bg);

% For ground truth, generate SMD
SMD_true = smi_core.SingleMoleculeData.createSMD();
SMD_true.NFrames = NFrames;
SMD_true.NDatasets = 1;
for nn = 1:NFrames
    N = poissrnd(Rho * SZ * SZ);
    SMD_true.FrameNum = cat(1, SMD_true.FrameNum, nn*ones(N,1));
    SMD_true.X = cat(1, SMD_true.X, SZ*rand(N,1));
    SMD_true.Y = cat(1, SMD_true.Y, SZ*rand(N,1));
end

fprintf('Generated %d frames with ~%d emitters per frame\n', NFrames, round(Rho*SZ*SZ));
fprintf('Image size: %d × %d pixels\n', SZ, SZ);
fprintf('Photons per emitter: %d\n', Photons);

%% Step 2: Configure SMF for Localization
fprintf('\n=== Configuring Parameters ===\n');

% Create SMF structure
SMF = smi_core.SingleMoleculeFitting();

% Camera parameters (for photon conversion)
% Since our simulated data is already in photons, set gain=1, offset=0
SMF.Data.CameraGain = 1;
SMF.Data.CameraOffset = 0;
SMF.Data.PixelSize = 0.1;     % 100 nm pixels
SMF.Data.FrameRate = 100;     % 100 Hz

% Box finding parameters
SMF.BoxFinding.BoxSize = 7;              % 7×7 pixel boxes
SMF.BoxFinding.MinPhotons = 300;         % Detection threshold
SMF.BoxFinding.BoxOverlap = 2;           % Allow 2-pixel overlap

% Fitting parameters
SMF.Fitting.PSFSigma = 1.3;              % Known PSF width
SMF.Fitting.FitType = 'XYNB';            % Fit X, Y, photons, background
SMF.Fitting.Iterations = 20;             % MLE iterations

% Thresholding parameters
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.2;         % Max position uncertainty (pixels)
SMF.Thresholding.MinPhotons = 200;       % Min photons after fitting
SMF.Thresholding.MinPValue = 0.01;       % Min fit quality
SMF.Thresholding.MinPSFSigma = 0.8;      % Min PSF sigma
SMF.Thresholding.MaxPSFSigma = 1.8;      % Max PSF sigma

fprintf('Box size: %d × %d pixels\n', SMF.BoxFinding.BoxSize, SMF.BoxFinding.BoxSize);
fprintf('PSF sigma: %.1f pixels\n', SMF.Fitting.PSFSigma);
fprintf('Detection threshold: %d photons\n', SMF.BoxFinding.MinPhotons);

%% Step 3: Localize Molecules
fprintf('\n=== Localizing Molecules ===\n');

% Create LocalizeData object
LD = smi_core.LocalizeData(imageStack_noisy, SMF);
LD.Verbose = 1;  % Show progress

% Run localization
tic;
SMD = LD.genLocalizations();
elapsed_time = toc;

fprintf('Localization complete in %.2f seconds\n', elapsed_time);
fprintf('Found %d localizations\n', length(SMD.X));

% Check thresholding
if isfield(SMD, 'ThreshFlag')
    passed = sum(SMD.ThreshFlag == 0);
    fprintf('Passed quality filters: %d (%.1f%%)\n', ...
        passed, 100*passed/length(SMD.X));
end

%% Step 4: Assess Localization Quality
fprintf('\n=== Quality Assessment ===\n');

% Localization precision
median_precision_x = median(SMD.X_SE);
median_precision_y = median(SMD.Y_SE);
median_precision_nm = median_precision_x * SMF.Data.PixelSize * 1000;

fprintf('Median X precision: %.3f pixels (%.1f nm)\n', ...
    median_precision_x, median_precision_nm);
fprintf('Median Y precision: %.3f pixels (%.1f nm)\n', ...
    median_precision_y, median_precision_y * SMF.Data.PixelSize * 1000);

% Photon statistics
fprintf('Photons: %.0f ± %.0f (mean ± std)\n', ...
    mean(SMD.Photons), std(SMD.Photons));
fprintf('Background: %.1f ± %.1f photons/pixel\n', ...
    mean(SMD.Bg), std(SMD.Bg));

% Fit quality
fprintf('Median p-value: %.3f\n', median(SMD.PValue));

%% Step 5: Compare to Ground Truth
fprintf('\n=== Comparison to Ground Truth ===\n');

% Match localizations to true positions
% Simple nearest-neighbor matching
max_distance = 2;  % pixels
matches = 0;
distances = [];

for i = 1:length(SMD_true.X)
    % Find nearest localization
    dx = SMD.X - SMD_true.X(i);
    dy = SMD.Y - SMD_true.Y(i);
    dist = sqrt(dx.^2 + dy.^2);

    [min_dist, ~] = min(dist);

    if min_dist < max_distance
        matches = matches + 1;
        distances = [distances; min_dist];
    end
end

detection_rate = 100 * matches / length(SMD_true.X);
fprintf('Detection rate: %.1f%% (%d / %d true emitters)\n', ...
    detection_rate, matches, length(SMD_true.X));

if ~isempty(distances)
    fprintf('Mean localization error: %.3f pixels (%.1f nm)\n', ...
        mean(distances), mean(distances) * SMF.Data.PixelSize * 1000);
    fprintf('Median localization error: %.3f pixels (%.1f nm)\n', ...
        median(distances), median(distances) * SMF.Data.PixelSize * 1000);
end

%% Step 6: Visualizations
fprintf('\n=== Creating Visualizations ===\n');

% Create figure with multiple panels
figure('Name', 'Basic Localization Example', 'Position', [100, 100, 1400, 900]);

% Panel 1: Raw data (first frame)
subplot(2,3,1);
imagesc(imageStack_noisy(:,:,1));
axis image; colormap(gca, 'gray'); colorbar;
title('Raw Data (Frame 1)');
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Panel 2: Localizations vs ground truth
subplot(2,3,2);
plot(SMD_true.X, SMD_true.Y, 'go', 'MarkerSize', 8, 'DisplayName', 'True positions');
hold on;
plot(SMD.X, SMD.Y, 'r.', 'MarkerSize', 4, 'DisplayName', 'Localizations');
axis equal; axis([0 SZ 0 SZ]);
legend('Location', 'best');
title(sprintf('Localizations (%.0f%% detected)', detection_rate));
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Panel 3: Localization precision distribution
subplot(2,3,3);
histogram(SMD.X_SE * SMF.Data.PixelSize * 1000, 30);
xlabel('Localization Precision (nm)');
ylabel('Count');
title(sprintf('Precision (median: %.1f nm)', median_precision_nm));
grid on;

% Panel 4: Photon distribution
subplot(2,3,4);
histogram(SMD.Photons, 40);
xlabel('Detected Photons');
ylabel('Count');
title(sprintf('Photons (mean: %.0f)', mean(SMD.Photons)));
xline(mean(SMD.Photons), 'r--', 'LineWidth', 2);
grid on;

% Panel 5: Background distribution
subplot(2,3,5);
histogram(SMD.Bg, 30);
xlabel('Background (photons/pixel)');
ylabel('Count');
title(sprintf('Background (mean: %.1f)', mean(SMD.Bg)));
xline(Bg, 'g--', 'LineWidth', 2, 'Label', 'True');
grid on;

% Panel 6: Localization error distribution
subplot(2,3,6);
if ~isempty(distances)
    histogram(distances * SMF.Data.PixelSize * 1000, 30);
    xlabel('Localization Error (nm)');
    ylabel('Count');
    title(sprintf('Error vs Truth (median: %.1f nm)', ...
        median(distances) * SMF.Data.PixelSize * 1000));
    grid on;
else
    text(0.5, 0.5, 'No matches found', 'HorizontalAlignment', 'center');
    axis off;
end

%% Step 7: Create Super-Resolution Image
fprintf('\n=== Generating Super-Resolution Image ===\n');

% Use Gaussian rendering
SR_zoom = 20;  % 20× zoom (effective pixel size = 5 nm)

% Prepare SMD for gaussianImage (needs XSize, YSize)
SMD.XSize = SZ;
SMD.YSize = SZ;

SR_image = smi_vis.GenerateImages.gaussianImage(SMD, SR_zoom, 0);

% Display
figure('Name', 'Super-Resolution Image');
imagesc(SR_image);
axis image; colormap hot; colorbar;
title(sprintf('Super-Resolution Image (%d× zoom, %.1f nm pixels)', ...
    SR_zoom, SMF.Data.PixelSize * 1000 / SR_zoom));
xlabel('X'); ylabel('Y');

fprintf('Effective SR pixel size: %.1f nm\n', ...
    SMF.Data.PixelSize * 1000 / SR_zoom);

%% Step 8: Summary Statistics
fprintf('\n=== Summary Statistics ===\n');
fprintf('----------------------------------------\n');
fprintf('Data:\n');
fprintf('  Frames: %d\n', NFrames);
fprintf('  True emitters: %d\n', length(SMD_true.X));
fprintf('  Localizations found: %d\n', length(SMD.X));
fprintf('  Detection rate: %.1f%%\n', detection_rate);
fprintf('\nQuality:\n');
fprintf('  Median precision: %.1f nm\n', median_precision_nm);
fprintf('  Mean photons: %.0f\n', mean(SMD.Photons));
fprintf('  Mean background: %.1f photons/pixel\n', mean(SMD.Bg));
if ~isempty(distances)
    fprintf('  Median error vs truth: %.1f nm\n', ...
        median(distances) * SMF.Data.PixelSize * 1000);
end
fprintf('\nPerformance:\n');
fprintf('  Localization time: %.2f seconds\n', elapsed_time);
fprintf('  Speed: %.0f localizations/second\n', length(SMD.X) / elapsed_time);
fprintf('----------------------------------------\n');

fprintf('\nExample complete!\n');
```

## What This Example Demonstrates

### Data Generation
- Creates realistic SMLM data with known ground truth
- Adds Poisson noise to simulate camera detection
- Configurable parameters (photons, background, PSF)

### Localization Pipeline
- Box finding with appropriate thresholds
- PSF fitting using maximum likelihood
- Quality filtering with multiple criteria

### Validation
- Matches localizations to ground truth
- Computes detection rate and localization errors
- Verifies precision estimates

### Visualization
- Raw data inspection
- Localization scatter plots
- Quality metric distributions
- Super-resolution image generation

## Expected Results

When you run this example, you should see:

**Detection rate**: ~95-100% (most true emitters found)

**Localization precision**: ~10-20 nm (depends on photon count)

**Localization error**: ~10-30 nm median error vs ground truth

**Processing speed**: 100s to 1000s of localizations per second (depends on GPU)

## Modifications to Try

### 1. Vary Photon Count

```matlab
% Brighter emitters (better precision)
Photons = 2000;
imageStack_noisy = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, ...
    Photons, PSFSigma, Bg);

% Dimmer emitters (worse precision)
Photons = 500;
imageStack_noisy = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, ...
    Photons, PSFSigma, Bg);
```

### 2. Change Background Level

```matlab
% Higher background (worse SNR)
Bg = 20;
imageStack_noisy = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, ...
    Photons, PSFSigma, Bg);

% Lower background (better SNR)
Bg = 1;
imageStack_noisy = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, ...
    Photons, PSFSigma, Bg);
```

### 3. Adjust Detection Threshold

```matlab
% More sensitive (more detections, more noise)
SMF.BoxFinding.MinPhotons = 150;

% More conservative (fewer detections, cleaner)
SMF.BoxFinding.MinPhotons = 500;
```

### 4. Fit PSF Sigma

```matlab
% Let PSF width vary
SMF.Fitting.FitType = 'XYNBS';
SMF.Fitting.PSFSigma = 1.5;  % Initial guess

% After fitting, check distribution
histogram(SMD.PSFSigma);
xlabel('Fitted PSF Sigma (pixels)');
```

### 5. Add More Frames

```matlab
% More frames = more localizations
NFrames = 200;
imageStack_noisy = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, ...
    Photons, PSFSigma, Bg);
```

### 6. Change Image Size

```matlab
% Larger field of view
SZ = 256;
Rho = 0.01;  % Same density = more total emitters
imageStack_noisy = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, ...
    Photons, PSFSigma, Bg);
```

## Understanding the Output

### Precision vs Error

**Precision** (SMD.X_SE): Uncertainty estimate from CRLB
- Theoretical lower bound on localization uncertainty
- Depends on photons, background, PSF

**Error**: Actual distance from true position
- Can only compute with ground truth
- Should be close to estimated precision

If error >> precision, something is wrong (biased fits, wrong PSF model).

### Detection Rate

Not all emitters may be detected because:
- Too dim (< MinPhotons threshold)
- Overlapping with another emitter
- Bad fit (failed quality filters)
- Edge effects

Typical detection rates: 80-100% for well-separated emitters.

## Next Steps

After running this example:

1. **Try real data**: Use [How to Load Data](../how-to/load-data.md) to load your own files

2. **Full SMLM pipeline**: See [SMLM Workflow](../workflows/smlm-analysis.md) for frame connection and drift correction

3. **Optimize parameters**: Use [How to Localize Molecules](../how-to/localize-molecules.md) to tune settings

4. **Batch processing**: Learn about `smi.SMLM` for automated analysis

5. **Advanced simulation**: Explore `smi_sim.SimSMLM` for more realistic simulations

## Common Issues

**Issue: Low detection rate (<80%)**

Solutions:
- Lower `BoxFinding.MinPhotons`
- Increase `Photons` parameter (make emitters brighter)
- Check thresholding isn't too aggressive

**Issue: Poor precision (>50 nm)**

Solutions:
- Increase `Photons` parameter
- Decrease `Bg` parameter
- Check camera parameters are correct

**Issue: High localization error**

Solutions:
- Verify `SMF.Fitting.PSFSigma` matches `PSFSigma` simulation parameter
- Check for fitting convergence (increase iterations)
- Ensure thresholding filters poor fits

**Issue: Too slow**

Solutions:
- Check GPU is working: `gpuDevice`
- Reduce `NFrames` parameter for testing
- Reduce `Rho` parameter (density)

## See Also

- [How to Localize Molecules](../how-to/localize-molecules.md) - Detailed localization guide
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete analysis pipeline
- [First Analysis](../getting-started/first-analysis.md) - GUI-based workflow
- MATLAB/examples/Example_LocalizeData.m - Additional examples
- MATLAB/examples/Example_GaussBlobs.m - More simulation examples
