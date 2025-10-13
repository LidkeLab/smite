---
title: "How to Calibrate Your Camera"
category: "how-to"
level: "intermediate"
tags: ["camera", "calibration", "gain", "offset", "scmos", "emccd", "photon-conversion"]
prerequisites: ["../core-concepts/smf-structure.md", "create-smf.md"]
related: ["load-data.md", "localize-molecules.md", "../workflows/smlm-analysis.md"]
summary: "Complete guide to camera calibration for accurate photon conversion, including gain/offset measurements and sCMOS variance map generation"
estimated_time: "30 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# How to Calibrate Your Camera

## Purpose

Camera calibration is essential for converting raw camera output (ADU - Analog to Digital Units) to physical photon counts. Accurate calibration enables precise localization, quantitative photon measurements, and proper uncertainty estimates. This guide covers why calibration matters, the differences between EMCCD and sCMOS calibration, how to measure gain and offset, how to generate sCMOS variance maps, and the proper calibration file format for smite.

## Prerequisites

- Understanding of [SMF structure](../core-concepts/smf-structure.md)
- Basic knowledge of camera operation
- Access to camera control software (e.g., MIC - MATLAB Instrument Control)
- For sCMOS: Controllable light source for multi-level measurements

## Overview

Camera calibration determines three critical parameters per pixel:

- **Offset**: Baseline signal with no light (ADU)
- **Gain**: Conversion factor from photoelectrons to ADU (ADU/e-)
- **Variance**: Pixel-wise read noise (ADU^2)

For **EMCCD cameras**, these parameters are uniform across pixels (single scalar values).

For **sCMOS cameras**, each pixel has unique values (2D arrays matching sensor dimensions).

The calibration equation:
```
RawData_ADU = Gain * Photons + Offset + Noise
Photons = (RawData_ADU - Offset) / Gain
```

## Why Camera Calibration Matters

### Impact on Analysis Quality

**Without proper calibration:**
- Localization precision estimates are incorrect
- Photon counts are meaningless (arbitrary units)
- Thresholding on photons doesn't work properly
- Emitter brightness comparisons are invalid
- Statistical tests (likelihood ratios, chi-square) fail
- Results cannot be compared across experiments

**With proper calibration:**
- Accurate photon counts enable quantitative measurements
- Localization precision (Cramer-Rao Lower Bound) is correct
- Brightness-based filtering removes dim/bad localizations
- Results are reproducible and comparable
- Multi-color experiments align properly

### Physical Meaning

When you calibrate correctly, smite reports:
- **Photons**: Physical photons detected (not ADU)
- **Precision (sigma)**: Localization uncertainty in nanometers
- **Intensity**: Actual emitter brightness for comparing molecules

Without calibration, all numbers are relative and experiment-specific.

## EMCCD vs sCMOS Calibration

### EMCCD Camera Characteristics

**Uniform gain and offset:**
- Single gain value applies to entire sensor
- Single offset value for whole sensor
- Read noise is constant across pixels
- Electron multiplication amplifies signal uniformly

**Advantages:**
- Simple calibration (3 scalar values)
- Very low read noise when using EM gain
- Excellent for low-light imaging

**Calibration frequency:**
- Gain: Rarely changes (check if EM gain setting changes)
- Offset: Can drift slightly over time (monthly checks)
- Variance: Stable once measured

### sCMOS Camera Characteristics

**Pixel-wise variation:**
- Every pixel has unique gain
- Every pixel has unique offset
- Read noise varies significantly across sensor (2-10× variation)
- No electron multiplication

**Advantages:**
- Fast readout (100+ fps)
- Large field of view
- Lower cost per pixel

**Calibration requirements:**
- Gain map: 2D array (Y × X pixels)
- Offset map: 2D array (Y × X pixels)
- Variance map: 2D array (Y × X pixels)

**Calibration frequency:**
- Must recalibrate if camera ROI changes
- Recommended every 3-6 months
- Critical after sensor maintenance/replacement

### Decision Guide

**Use existing calibration if:**
- Camera settings unchanged (same ROI, same readout mode)
- Calibration is recent (< 6 months)
- Analysis results look reasonable

**Recalibrate if:**
- First time using camera with smite
- Changed camera ROI or readout mode
- Localization precision looks wrong
- Photon counts seem off by 2-10×
- Camera serviced or settings reset

## EMCCD Camera Calibration

### Measuring Camera Offset

Offset is the mean signal when no light hits the sensor.

**Procedure:**

```matlab
% 1. Set up dark conditions
% - Close shutter or block all light
% - Turn off all illumination sources
% - Wait for thermal equilibrium (5-10 minutes)

% 2. Acquire dark frames
% Use your camera control software or MIC
% Example with generic MATLAB acquisition:
num_frames = 1000;  % More frames = better statistics
dark_frames = % ... acquire 1000 dark frames (Y × X × 1000)

% 3. Calculate offset
offset = mean(dark_frames(:));
fprintf('Camera Offset: %.2f ADU\n', offset);

% Typical EMCCD offset: 100-300 ADU (with preamp offset)
% Should be stable across frames
```

**Verification:**

```matlab
% Check offset stability over time
mean_per_frame = squeeze(mean(mean(dark_frames, 1), 2));
figure;
plot(mean_per_frame);
xlabel('Frame'); ylabel('Mean ADU');
title('Offset Stability Check');

% Standard deviation should be << offset
fprintf('Offset stability: %.3f ADU (%.1f%%)\n', ...
    std(mean_per_frame), 100*std(mean_per_frame)/offset);

% Good stability: < 1% variation
```

**Common issues:**
- **Offset too low (< 50 ADU)**: Risk of negative values after subtraction. Increase camera offset setting.
- **Offset unstable**: Camera warming up or environmental noise. Wait longer, check cooling.
- **Offset varies spatially**: Normal for sCMOS, unexpected for EMCCD. Check camera health.

### Measuring Camera Gain

Gain converts photoelectrons to ADU. For EMCCD with electron multiplication, gain includes both EM gain and preamp gain.

**Method 1: From Specifications**

Check camera manual or calibration certificate:

```matlab
% Example specifications:
% - Preamp gain setting: 2× (ADU/e- without EM)
% - EM gain setting: 300× (electron multiplication)
% - Result: Total gain depends on preamp

% For iXon, typical preamp gains:
% Gain 1: ~3.0 ADU/e-
% Gain 2: ~6.5 ADU/e-
% Gain 3: ~12.0 ADU/e-

% These are WITHOUT EM gain applied
gain_without_EM = 3.0;  % From specs (ADU/e-)
```

**Note:** When EM gain is enabled, you still use the preamp ADU/e- gain value in smite because the EM amplification preserves the Poisson statistics of the detected photoelectrons. The offset and gain correct for the detection chain, while the Poisson statistics remain related to the original photon count.

**Method 2: Photon Transfer Curve**

Measure gain from variance vs mean relationship:

```matlab
% 1. Acquire frames at multiple light levels
% Start dim, increase gradually to ~50% saturation
light_levels = 10;  % Number of different intensities
frames_per_level = 500;

mean_signal = zeros(light_levels, 1);
var_signal = zeros(light_levels, 1);

for i = 1:light_levels
    % Set light level (use your illumination control)
    % lamp.setPower(i * max_power / light_levels);

    % Acquire frames
    frames = % ... acquire frames at this level

    % Calculate mean and variance (after offset subtraction)
    frames_corrected = frames - offset;
    mean_signal(i) = mean(frames_corrected(:));
    var_signal(i) = var(frames_corrected(:));
end

% 2. Fit variance vs mean
% For ideal camera: Variance = Gain * Mean + ReadNoise^2
% Slope = Gain

coeffs = polyfit(mean_signal, var_signal, 1);
gain = coeffs(1);
readnoise_sq = coeffs(2);

fprintf('Camera Gain: %.3f ADU/e-\n', gain);
fprintf('Read Noise: %.2f ADU^2\n', readnoise_sq);

% Plot to verify linear relationship
figure;
plot(mean_signal, var_signal, 'o');
hold on;
plot(mean_signal, gain*mean_signal + readnoise_sq, '-');
xlabel('Mean Signal (ADU)');
ylabel('Variance (ADU^2)');
title('Photon Transfer Curve');
legend('Data', 'Linear Fit');
```

**Typical EMCCD gains:**
- Without EM gain: 1-15 ADU/e- (depends on preamp setting)
- With EM gain: Use preamp value, NOT total amplification
- iXon Ultra: 2-5 ADU/e- (common preamp settings)
- Evolve: 3-8 ADU/e- (typical range)

### Measuring Camera Read Noise

Read noise is the variance in dark frames.

```matlab
% Use the dark frames from offset measurement
read_noise_variance = var(dark_frames(:));
read_noise_std = sqrt(read_noise_variance);

fprintf('Read Noise Variance: %.2f ADU^2\n', read_noise_variance);
fprintf('Read Noise Std Dev: %.2f ADU\n', read_noise_std);

% Store variance (not std dev) for smite
camera_variance = read_noise_variance;

% Typical EMCCD read noise (with EM gain):
% Variance: 5-25 ADU^2
% Std dev: 2-5 ADU
```

### Setting EMCCD Parameters in SMF

```matlab
% Create SMF
SMF = smi_core.SingleMoleculeFitting();

% Camera type
SMF.Data.CameraType = 'EMCCD';

% Calibration parameters (from measurements above)
SMF.Data.CameraGain = gain;           % ADU/e- (e.g., 2.5)
SMF.Data.CameraOffset = offset;       % ADU (e.g., 100)
SMF.Data.CameraReadNoise = camera_variance;  % ADU^2 (e.g., 8)

% Other required parameters
SMF.Data.PixelSize = 0.108;  % micrometers
SMF.Data.FrameRate = 100;    % Hz

fprintf('EMCCD calibration configured in SMF\n');
```

**Quick validation:**

```matlab
% Load a typical frame
test_frame = % ... load one frame of real data
photons_converted = (test_frame - SMF.Data.CameraOffset) / SMF.Data.CameraGain;

% Check converted values are reasonable
fprintf('Photon range: %.0f to %.0f photons/pixel\n', ...
    min(photons_converted(:)), max(photons_converted(:)));

% Background should be ~0 photons
% Bright emitters: 500-5000 photons typical for SMLM
% If numbers seem wrong, recheck gain and offset
```

## sCMOS Camera Calibration

sCMOS calibration requires pixel-wise measurements at multiple light levels.

### Overview of sCMOS Calibration

**What you need:**
1. Controllable light source (LED, lamp, uniform illumination)
2. 1000+ frames at each of 10-20 light levels
3. Dark frames at several exposure times
4. Image analysis to extract per-pixel gain/offset/variance

**Time required:**
- Data collection: 30-60 minutes
- Processing: 5-10 minutes
- Do this once per ROI configuration, save for reuse

**Equipment:**
- sCMOS camera with stable temperature control
- Uniform illumination source (LED or lamp)
- MIC (MATLAB Instrument Control) or equivalent camera control

### Data Collection Procedure

This procedure follows the method described in doc/FileFormats/CalibrationFile.md.

```matlab
%% Setup Equipment
% Connect to sCMOS camera
CameraSCMOS = MIC_HamamatsuCamera();  % Or your camera class
CameraSCMOS.ReturnType = 'matlab';
CameraSCMOS.gui();  % Optional: check camera view

% Connect to controllable light source
Lamp660 = MIC_ThorlabsLED('Dev1', 'ao0');  % Or your lamp
Lamp660.gui();

%% Configure Camera Settings
% Use slow, high-quality readout mode
CameraSCMOS.DefectCorrection = 1;  % Minimal or no correction
CameraSCMOS.ScanMode = 1;          % Slow scan (best quality)
CameraSCMOS.ExpTime_Sequence = 0.01;  % 10 ms exposure
CameraSCMOS.SequenceLength = 1000;    % 1000 frames per level

% Define ROI if needed (or use full sensor)
% ROI format: [XStart, XEnd, YStart, YEnd]
CameraSCMOS.ROI = [897, 1152, 897, 1152];  % Example: 256×256 region
% Or full sensor:
% CameraSCMOS.ROI = [1, 2048, 1, 2048];  % Example for 2048×2048 camera

%% Define Light Levels
% Choose range from dark to ~60% of saturation
% For 12-bit camera (0-4095 ADU), target max ~2500 ADU
% For 16-bit camera (0-65535 ADU), target max ~40000 ADU

% Adjust lamp power range based on your setup
LampPowerRange = linspace(0, 3.6, 20);  % 20 levels from off to 3.6V

%% Collect Calibration Data
Params = [];
MeanLevel = [];   % Mean intensity at each level
VarLevel = [];    % Variance at each level

CameraSCMOS.AcquisitionType = 'sequence';
CameraSCMOS.setup_acquisition();

fprintf('Collecting calibration data...\n');
for ii = 1:numel(LampPowerRange)
    fprintf('Light level %d of %d (power %.2f)\n', ...
        ii, numel(LampPowerRange), LampPowerRange(ii));

    % Set light level
    Lamp660.setPower(LampPowerRange(ii));
    pause(1);  % Allow stabilization

    % Acquire frames
    CameraSCMOS.start_sequence();

    % Calculate mean and variance across frames
    MeanLevel = cat(3, MeanLevel, mean(CameraSCMOS.Data, 3));
    VarLevel = cat(3, VarLevel, var(single(CameraSCMOS.Data), [], 3));
end

% Turn off light
Lamp660.setPower(0);

% Store collected data
Params.MeanLevel = single(MeanLevel);  % Y × X × NumLevels
Params.VarLevel = single(VarLevel);    % Y × X × NumLevels
Params.LampPowerRange = single(LampPowerRange);

fprintf('Data collection complete\n');
```

### Processing Calibration Data

Perform pixel-wise least squares fit to extract gain and offset.

```matlab
%% Fit Gain and Offset Per Pixel
% For each pixel: Variance = Gain * Mean + Offset_variance
% We fit this relationship to get gain and baseline variance

[nY, nX, nLevels] = size(MeanLevel);
Beta = NaN(nY, nX, 2);  % [:,:,1] = offset, [:,:,2] = gain

fprintf('Fitting gain and offset for %d × %d pixels...\n', nY, nX);
for ii = 1:nY
    if mod(ii, 50) == 0
        fprintf('  Row %d of %d\n', ii, nY);
    end
    for jj = 1:nX
        % Extract variance vs mean for this pixel
        x_data = squeeze(MeanLevel(ii, jj, :));
        y_data = squeeze(VarLevel(ii, jj, :));
        weights = 1 ./ y_data;  % Weight by inverse variance

        % Perform weighted least squares fit
        % y = beta(1) + beta(2)*x  =>  beta(1)=offset, beta(2)=gain
        Beta(ii, jj, 1:2) = smi_stat.leastSquaresFit(x_data, y_data, weights);
    end
end

fprintf('Fitting complete\n');

%% Extract Calibration Parameters
% Gain is the slope (Beta(:,:,2))
Params.Gain = single(Beta(:, :, 2));

% Offset is the mean at the lowest light level (dark level)
Params.CCDOffset = single(MeanLevel(:, :, 1));

% Variance is the intercept (variance at zero signal)
Params.CCDVar = single(VarLevel(:, :, 1));

%% Quality Checks
% Visualize gain map
figure;
imagesc(Params.Gain);
axis image; colorbar;
title('Gain Map (ADU/e-)');
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Check gain distribution
figure;
histogram(Params.Gain(:), 50);
xlabel('Gain (ADU/e-)');
ylabel('Pixel Count');
title('Gain Distribution');

fprintf('\nCalibration Statistics:\n');
fprintf('Gain:   %.3f ± %.3f ADU/e- (range: %.3f - %.3f)\n', ...
    mean(Params.Gain(:)), std(Params.Gain(:)), ...
    min(Params.Gain(:)), max(Params.Gain(:)));

fprintf('Offset: %.1f ± %.1f ADU (range: %.1f - %.1f)\n', ...
    mean(Params.CCDOffset(:)), std(Params.CCDOffset(:)), ...
    min(Params.CCDOffset(:)), max(Params.CCDOffset(:)));

fprintf('Variance: %.2f ± %.2f ADU^2 (range: %.2f - %.2f)\n', ...
    mean(Params.CCDVar(:)), std(Params.CCDVar(:)), ...
    min(Params.CCDVar(:)), max(Params.CCDVar(:)));

% Typical sCMOS values:
% Gain: 0.3-3.0 ADU/e- (mean ~1.0)
% Offset: 100-300 ADU
% Variance: 1-20 ADU^2 (varies significantly pixel-to-pixel)
```

**Expected results:**
- Gain map shows slight pixel-to-pixel variation (5-20% coefficient of variation)
- Offset map relatively uniform (unless cosmetic defects present)
- Variance map shows 2-10× variation across sensor (this is normal for sCMOS)

**Red flags:**
- Gain values < 0.1 or > 10 ADU/e-: Check light levels or fitting
- Negative variances: Data collection issue or insufficient light levels
- Extreme outliers in maps: Dead/hot pixels (expected, not a problem)

### Saving sCMOS Calibration File

```matlab
%% Save Calibration to File
% Include ROI information for future reference
Params.CameraObj.ROI = CameraSCMOS.ROI;

% Save with timestamp for version tracking
SaveDir = 'Y:\sCMOS Calibrations\Sequential SR';  % Your calibration directory
if ~exist(SaveDir, 'dir')
    mkdir(SaveDir);
end

FileName = fullfile(SaveDir, ...
    ['GainCalibration-', smi_helpers.genTimeString()]);

% Save as MATLAB v7.3 for large arrays
save(FileName, 'Params', '-v7.3');

fprintf('\nCalibration saved to:\n%s.mat\n', FileName);

%% Verify Saved File
% Load and check structure
CalibrationFile = load([FileName, '.mat']);
fprintf('\nCalibration file structure:\n');
disp(CalibrationFile.Params);

% Should contain:
% Params.Gain       - Y × X single array (ADU/e-)
% Params.CCDOffset  - Y × X single array (ADU)
% Params.CCDVar     - Y × X single array (ADU^2)
% Params.CameraObj.ROI - [XStart, XEnd, YStart, YEnd]
% Params.LampPowerRange - Light levels used
% Params.MeanLevel  - Raw mean data (optional, for reprocessing)
% Params.VarLevel   - Raw variance data (optional, for reprocessing)
```

### Using sCMOS Calibration in SMF

```matlab
% Create SMF
SMF = smi_core.SingleMoleculeFitting();

% Specify sCMOS camera type
SMF.Data.CameraType = 'SCMOS';

% Point to calibration file
SMF.Data.CalibrationFilePath = [FileName, '.mat'];

% Verify calibration loads correctly
if ~exist(SMF.Data.CalibrationFilePath, 'file')
    error('Calibration file not found');
end

% smite will automatically load Gain, CCDOffset, and CCDVar from file
% when you run analysis

% Other required parameters
SMF.Data.PixelSize = 0.1;    % micrometers
SMF.Data.FrameRate = 200;    % Hz

fprintf('sCMOS calibration configured\n');
```

**Important notes:**

1. **ROI consistency**: If you change camera ROI, you must recalibrate or extract the correct subregion from your full-sensor calibration.

2. **CalibrationROI**: smite uses `CalibrationROI` to extract the correct portion of the calibration maps to match your data ROI. This is stored in the calibration file and handled automatically.

3. **Memory**: Full-sensor calibration files can be large (several MB for 2048×2048 sensor). This is normal.

## Calibration File Format

### Required Variables

smite expects calibration files to contain specific variable names:

**For use with DataToPhotons class:**

```matlab
% Variable names that smite looks for:
Params.Gain        % or CameraGain or CCDGain
Params.CCDOffset   % or CameraOffset or Offset
Params.CCDVar      % or CameraReadNoise or CameraNoise
```

The `smi_core.LoadData.loadDataCalibration()` method handles these variations, but the standard naming convention uses:
- `Gain` (or `CCDGain`)
- `CCDOffset`
- `CCDVar`

### Complete Calibration File Structure

Example from the documentation:

```matlab
>> CalibrationFile.Params

ans =

  struct with fields:

            MeanLevel: [256×256×20 single]  % Mean at each light level
             VarLevel: [256×256×20 single]  % Variance at each level
       DarkImagesMean: [256×256×5 single]   % Dark frames (optional)
        DarkImagesVar: [256×256×5 single]   % Dark variance (optional)
    DarkImagesExpTime: [0.01 0.13 0.26 0.38 0.50]  % Exposures (optional)
DarkImagesSequenceLength: 1000              % Frames per dark (optional)
               CCDVar: [256×256 single]     % Read noise variance (required)
              CCDGain: [256×256 single]     % Gain map (required, or "Gain")
            CCDOffset: [256×256 single]     % Offset map (required)
            CameraObj: [1×1 struct]         % Camera info (optional)
       LampPowerRange: [1×20 single]        % Light levels used (optional)
                 Gain: [256×256 single]     % Alternative name for CCDGain
```

### Minimal Calibration File

Minimum required for smite:

```matlab
% For sCMOS
Params.Gain = % ... Y × X single array (ADU/e-)
Params.CCDOffset = % ... Y × X single array (ADU)
Params.CCDVar = % ... Y × X single array (ADU^2)

save('minimal_scmos_calib.mat', 'Params');
```

```matlab
% For EMCCD (if using calibration file instead of SMF.Data fields)
Params.Gain = 2.5;        % scalar
Params.CCDOffset = 100;   % scalar
Params.CCDVar = 8;        % scalar

save('minimal_emccd_calib.mat', 'Params');
```

## Validating Calibration Quality

### Visual Inspection

```matlab
% Load a typical data frame
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'test_data.h5'};
LD = smi_core.LoadData();
[~, raw_data, ~] = LD.loadRawData(SMF, 1);

% Convert to photons manually
if strcmp(SMF.Data.CameraType, 'EMCCD')
    gain = SMF.Data.CameraGain;
    offset = SMF.Data.CameraOffset;
else  % SCMOS
    calib = load(SMF.Data.CalibrationFilePath);
    gain = calib.Params.Gain;
    offset = calib.Params.CCDOffset;
end

% Convert one frame
frame_adu = raw_data(:, :, 1);
frame_photons = (frame_adu - offset) ./ gain;

% Visualize
figure;
subplot(1, 2, 1);
imagesc(frame_adu);
axis image; colorbar;
title('Raw Data (ADU)');

subplot(1, 2, 2);
imagesc(frame_photons);
axis image; colorbar;
title('Converted (Photons)');

% Check values
fprintf('ADU range: %.0f to %.0f\n', min(frame_adu(:)), max(frame_adu(:)));
fprintf('Photon range: %.0f to %.0f\n', min(frame_photons(:)), max(frame_photons(:)));
```

**Expected photon counts:**
- Background: -10 to +10 photons (noise around zero)
- Dim emitters: 50-200 photons
- Typical emitters: 200-1000 photons
- Bright emitters: 1000-5000 photons
- If photon counts seem 10× too high or low, recheck gain

### Quantitative Validation

Compare measured variance to expected Poisson statistics:

```matlab
% For a uniform region (no emitters), variance should equal mean
% Select background region
bg_region = frame_photons(50:100, 50:100);

mean_bg = mean(bg_region(:));
var_bg = var(bg_region(:));

fprintf('Background mean: %.2f photons\n', mean_bg);
fprintf('Background variance: %.2f photons^2\n', var_bg);
fprintf('Variance/Mean ratio: %.2f (should be ~1 for Poisson)\n', var_bg/mean_bg);

% Ratio should be close to 1 for properly calibrated data
% If ratio is far from 1, calibration may be incorrect
```

### Localization Test

Best validation: localize some emitters and check results.

```matlab
% Run localization on calibrated data
LD = smi_core.LocalizeData(raw_data, SMF);
[SMD] = LD.genLocalizations();

% Check photon distribution
figure;
histogram(SMD.Photons, 50);
xlabel('Photons per Localization');
ylabel('Count');
title('Photon Distribution');

fprintf('Photon statistics:\n');
fprintf('  Median: %.0f photons\n', median(SMD.Photons));
fprintf('  Mean: %.0f photons\n', mean(SMD.Photons));
fprintf('  Range: %.0f to %.0f photons\n', min(SMD.Photons), max(SMD.Photons));

% Check precision vs photons relationship
figure;
scatter(SMD.Photons, SMD.X_SE, 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Photons'); ylabel('X Precision (pixels)');
title('Localization Precision vs Brightness');
set(gca, 'XScale', 'log', 'YScale', 'log');

% Precision should improve with photon count (1/sqrt(N) relationship)
% If no relationship, calibration is likely wrong
```

## Troubleshooting

### Issue: Negative Photon Counts

**Symptom:** After conversion, many pixels have negative photon values.

**Causes:**
1. Offset too low (underestimated)
2. Hot pixels or cosmic rays in dark frames
3. Camera baseline changed since calibration

**Solutions:**

```matlab
% Check offset distribution
fprintf('Offset: %.2f ADU\n', mean(offset(:)));
fprintf('Raw data min: %.0f ADU\n', min(raw_data(:)));

% If raw data minimum < offset, offset is too high OR
% your data includes unusually dark pixels

% Option 1: Remeasure offset
% Follow dark frame procedure again

% Option 2: Adjust offset downward slightly
% Only if you know offset is correct but data is dark
SMF.Data.CameraOffset = SMF.Data.CameraOffset - 10;  % Reduce by 10 ADU
```

### Issue: Photon Counts Too High/Low

**Symptom:** Localizations show 10× or 100× expected photon counts.

**Cause:** Gain incorrect (most common calibration error).

**Solutions:**

```matlab
% Check gain reasonableness
fprintf('Current gain: %.3f ADU/e-\n', SMF.Data.CameraGain);

% Typical ranges:
% EMCCD: 1-10 ADU/e- (without considering EM gain)
% sCMOS: 0.3-3.0 ADU/e-

% If photon counts are 5× too high:
SMF.Data.CameraGain = SMF.Data.CameraGain * 5;

% If photon counts are 10× too low:
SMF.Data.CameraGain = SMF.Data.CameraGain / 10;

% Then reprocess data
```

### Issue: sCMOS Calibration Shows Artifacts

**Symptom:** Gain or variance maps show stripes, blocks, or sudden discontinuities.

**Causes:**
1. Non-uniform illumination during calibration
2. Insufficient frames per light level (< 500)
3. Light source instability
4. Vibration or drift during acquisition

**Solutions:**

1. **Check illumination uniformity:**
```matlab
% Examine mean intensity across sensor at mid-level
mid_level = round(size(Params.MeanLevel, 3) / 2);
figure;
imagesc(Params.MeanLevel(:, :, mid_level));
axis image; colorbar;
title('Illumination Uniformity at Mid-Level');

% Should be smooth, gradual variations only
% Stripes or blocks indicate illumination problems
```

2. **Recalibrate with improved setup:**
- Use diffuser for more uniform illumination
- Increase frames per level to 1000+
- Stabilize light source (longer warm-up time)
- Reduce vibration

3. **For minor artifacts in low-signal regions:**
```matlab
% Can proceed if artifacts are in unused regions
% Or if variations are < 20% and smoothly varying
```

### Issue: ROI Mismatch

**Symptom:** Error loading calibration: "CalibrationROI must contain RawDataROI"

**Cause:** Data ROI doesn't match calibration ROI.

**Solution:**

```matlab
% Check ROIs
fprintf('Data ROI: [%d, %d, %d, %d]\n', SMF.Data.DataROI);

calib = load(SMF.Data.CalibrationFilePath);
if isfield(calib.Params, 'CameraObj')
    fprintf('Calibration ROI: [%d, %d, %d, %d]\n', ...
        calib.Params.CameraObj.ROI);
end

% Option 1: Recalibrate with matching ROI
% Option 2: Use full-sensor calibration and let smite extract correct region
% Option 3: Adjust data ROI to match calibration
SMF.Data.DataROI = [];  % Use full sensor
```

## Best Practices

1. **Calibrate regularly:**
   - EMCCD: Every 6 months or after gain changes
   - sCMOS: Every 3-6 months or after ROI changes

2. **Document calibration conditions:**
   - Date and time
   - Camera settings (ROI, readout mode, temperature)
   - Who performed calibration
   - Any unusual observations

3. **Save raw calibration data:**
   - Keep `Params.MeanLevel` and `Params.VarLevel` in file
   - Allows reprocessing if needed
   - Only ~50-200 MB for full sensor

4. **Create calibration library:**
```matlab
% Organize calibrations by camera and ROI
% Y:\Calibrations\
%   Camera1\
%     FullSensor_2024-10-11.mat
%     ROI_256x256_center_2024-10-11.mat
%   Camera2\
%     FullSensor_2024-09-15.mat
```

5. **Validate with known sample:**
   - After calibration, image a well-characterized sample
   - Check photon counts match expected values
   - Compare precision estimates to previous calibrations

6. **Version control:**
```matlab
% Include version info in calibration file
Params.CalibrationVersion = '1.0';
Params.CalibrationDate = datestr(now);
Params.CalibratedBy = 'Username';
Params.Notes = 'First calibration after camera service';
```

## See Also

- [How to Create SMF](create-smf.md) - Setting camera parameters in SMF
- [How to Load Data](load-data.md) - Using calibrated data
- [SMF Structure](../core-concepts/smf-structure.md) - SMF.Data camera fields
- doc/FileFormats/CalibrationFile.md - Detailed calibration file format
- `smi_core.DataToPhotons` - Photon conversion implementation
- `smi_stat.leastSquaresFit` - Fitting algorithm for sCMOS calibration
