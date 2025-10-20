---
title: "How to Create and Configure SMF"
category: "how-to"
level: "beginner"
tags: ["smf", "configuration", "setup", "parameters", "camera"]
prerequisites: ["../core-concepts/smf-structure.md"]
related: ["load-data.md", "localize-molecules.md", "../workflows/smlm-analysis.md"]
summary: "Practical guide for creating and configuring SMF structures for different experiment types and camera setups"
estimated_time: "15 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# How to Create and Configure SMF

## Purpose

The SMF (Single Molecule Fitting) structure is your analysis blueprint in smite. This guide walks you through creating SMF from scratch, configuring it for different camera types, setting common parameters, and saving/loading configurations for reproducible analyses.

## Prerequisites

- Basic understanding of [SMF structure](../core-concepts/smf-structure.md)
- Knowledge of your experimental setup (camera type, pixel size, etc.)
- Raw data files ready for analysis

## Overview

Creating an SMF involves three main steps:

1. **Initialize**: Create SMF object with default values
2. **Configure**: Set parameters for your specific experiment
3. **Save**: Store configuration for reuse and reproducibility

SMF contains all parameters needed to go from raw data to final results. Once configured correctly, you can reuse the same SMF for similar experiments, ensuring consistency and reproducibility.

## Creating a Basic SMF

### From Scratch

```matlab
% Create new SMF with all defaults
SMF = smi_core.SingleMoleculeFitting();

% SMF is now a structure with sub-structures
% SMF.Data          - File and camera information
% SMF.BoxFinding    - Molecule detection parameters
% SMF.Fitting       - PSF fitting parameters
% SMF.Thresholding  - Quality filtering parameters
% SMF.FrameConnection - Linking across frames
% SMF.DriftCorrection - Stage drift compensation
% SMF.Tracking      - Particle tracking parameters
```

### Using the GUI

For interactive configuration:

```matlab
SMF = smi_core.SingleMoleculeFitting();
SMF.gui();  % Opens graphical interface

% The GUI provides:
% - Parameter descriptions and tooltips
% - Valid ranges and types
% - Real-time validation
% - Save/load functionality
% - File/directory browsers
```

The GUI is especially helpful when learning which parameters affect analysis.

### Quick Check

Verify SMF was created:

```matlab
% Check default values
fprintf('Box size: %d pixels\n', SMF.BoxFinding.BoxSize);
fprintf('PSF sigma: %.1f pixels\n', SMF.Fitting.PSFSigma);
fprintf('Camera type: %s\n', SMF.Data.CameraType);

% List all sub-structures
fields = fieldnames(SMF);
fprintf('\nSMF contains %d sub-structures:\n', length(fields));
for i = 1:length(fields)
    fprintf('  %s\n', fields{i});
end
```

## Configuring for EMCCD Camera

EMCCD cameras have uniform gain and offset across all pixels.

### Basic EMCCD Configuration

```matlab
% Create SMF
SMF = smi_core.SingleMoleculeFitting();

% Camera parameters
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;    % ADU per photon (electron multiplication gain)
SMF.Data.CameraOffset = 100;  % Baseline ADU (dark current offset)
SMF.Data.CameraNoise = 5;     % Read noise variance (ADU^2)

% Acquisition parameters
SMF.Data.PixelSize = 0.108;   % Back-projected pixel size (micrometers)
SMF.Data.FrameRate = 100;     % Acquisition rate (Hz)

% File location
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'experiment1.h5'};
SMF.Data.ResultsDir = '/path/to/data/Results';  % Where results are saved
```

### Determining EMCCD Parameters

**Camera Gain:**

From camera specifications or calibration:

```matlab
% Typical EMCCD gains: 1-10 ADU/photon
% Check camera manual or calibration certificate
SMF.Data.CameraGain = 2.5;  % Example value

% If unknown, start with 1 and adjust based on photon counts
```

**Camera Offset:**

From dark frames (no light):

```matlab
% Record dark frames (shutter closed, lights off)
dark_frames = % ... load your dark frames (Y × X × N_frames)

% Calculate offset (mean of dark frames)
SMF.Data.CameraOffset = mean(dark_frames(:));
fprintf('Camera offset: %.1f ADU\n', SMF.Data.CameraOffset);

% Typical range: 50-200 ADU for EMCCD
```

**Camera Noise:**

From dark frame variance:

```matlab
% Calculate variance of dark frames
SMF.Data.CameraNoise = var(dark_frames(:));
fprintf('Camera noise: %.2f ADU^2\n', SMF.Data.CameraNoise);

% Note: This is variance (σ²), not standard deviation
% Typical range: 1-25 ADU^2
```

### Complete EMCCD Example

```matlab
% ========== EMCCD Configuration ==========
SMF = smi_core.SingleMoleculeFitting();

% Camera (Andor iXon or similar)
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;
SMF.Data.CameraOffset = 100;
SMF.Data.CameraNoise = 8;

% Microscope setup
SMF.Data.PixelSize = 0.108;  % 100× objective, 1.6× tube lens, 6.5μm pixels
SMF.Data.FrameRate = 100;    % 10ms exposure

% Data location
SMF.Data.FileDir = '/data/2024-10-10/Cell1';
SMF.Data.FileName = {'PAINT_647.h5'};
SMF.Data.ResultsDir = '/data/2024-10-10/Cell1/Results';
SMF.Data.AnalysisID = 'v1';  % Optional identifier

fprintf('EMCCD configuration complete\n');
```

## Configuring for sCMOS Camera

sCMOS cameras have pixel-wise varying gain, offset, and noise that must be calibrated.

### Basic sCMOS Configuration

```matlab
% Create SMF
SMF = smi_core.SingleMoleculeFitting();

% Camera type
SMF.Data.CameraType = 'SCMOS';

% Point to calibration file
SMF.Data.CalibrationFilePath = '/path/to/calibration/scmos_calib.mat';

% Calibration file MUST contain these variables:
%   CameraGain    - Y × X array (ADU per photon)
%   CameraOffset  - Y × X array (baseline ADU)
%   CameraNoise   - Y × X array (read noise variance, ADU^2)

% Acquisition parameters
SMF.Data.PixelSize = 0.1;     % micrometers
SMF.Data.FrameRate = 200;     % Hz (sCMOS can be faster)

% File location
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'data.h5'};
```

### Creating sCMOS Calibration File

Calibrate your sCMOS camera once, then reuse the file:

```matlab
% ========== Calibration Procedure ==========

% 1. Collect calibration data (follow camera manufacturer's procedure)
%    - Dark frames at different exposure times
%    - Flat-field frames at different intensities
%    - Photon transfer curve measurements

% 2. Process calibration to get pixel-wise parameters
%    (This is camera-specific; consult manufacturer software)

% 3. Save to .mat file
CameraGain = % ... Y × X array from calibration
CameraOffset = % ... Y × X array from calibration
CameraNoise = % ... Y × X array from calibration

% Verify dimensions match camera
[height, width] = size(CameraGain);
fprintf('Calibration for %d × %d camera\n', height, width);

% Check reasonable ranges
fprintf('Gain range: %.2f to %.2f ADU/photon\n', ...
    min(CameraGain(:)), max(CameraGain(:)));
fprintf('Offset range: %.1f to %.1f ADU\n', ...
    min(CameraOffset(:)), max(CameraOffset(:)));
fprintf('Noise range: %.2f to %.2f ADU^2\n', ...
    min(CameraNoise(:)), max(CameraNoise(:)));

% Save calibration file
save('scmos_calibration.mat', 'CameraGain', 'CameraOffset', 'CameraNoise');
fprintf('Calibration saved to scmos_calibration.mat\n');
```

### Complete sCMOS Example

```matlab
% ========== sCMOS Configuration ==========
SMF = smi_core.SingleMoleculeFitting();

% Camera (Hamamatsu Orca Flash 4.0 or similar)
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = '/calibration/orca_flash_2024.mat';

% Verify calibration file exists
if ~exist(SMF.Data.CalibrationFilePath, 'file')
    error('Calibration file not found: %s', SMF.Data.CalibrationFilePath);
end

% Load and verify calibration
calib = load(SMF.Data.CalibrationFilePath);
fprintf('Calibration loaded: %d × %d pixels\n', size(calib.CameraGain));

% Microscope setup
SMF.Data.PixelSize = 0.1;    % 100nm pixels
SMF.Data.FrameRate = 200;    % Fast sCMOS acquisition

% Data location
SMF.Data.FileDir = '/data/sCMOS_experiment';
SMF.Data.FileName = {'fast_acquisition.h5'};
SMF.Data.ResultsDir = '/data/sCMOS_experiment/Results';

fprintf('sCMOS configuration complete\n');
```

## Setting Common Analysis Parameters

Configure detection, fitting, and filtering for your experiment type.

### For DNA-PAINT

High photon count, sparse emitters:

```matlab
% Box finding - high threshold for bright emitters
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 300;  % Strict threshold (bright emitters)
SMF.BoxFinding.BoxOverlap = 2;

% Fitting - standard 2D
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';
SMF.Fitting.Iterations = 20;

% Thresholding - aggressive filtering for high precision
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.15;    % Tight precision requirement
SMF.Thresholding.MinPhotons = 200;   % Post-fit photon filter
SMF.Thresholding.MinPValue = 0.01;
SMF.Thresholding.MinPSFSigma = 0.9;
SMF.Thresholding.MaxPSFSigma = 1.7;

% Frame connection - for blinking emitters
SMF.FrameConnection.On = true;
SMF.FrameConnection.MaxSeparation = 1.0;  % pixels
SMF.FrameConnection.MaxFrameGap = 5;      % frames

% Drift correction - important for long acquisitions
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.Method = 'DC-KNN';
```

### For dSTORM

Dense emitters with blinking:

```matlab
% Box finding - lower threshold (dimmer emitters)
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 150;  % More permissive
SMF.BoxFinding.BoxOverlap = 2;

% Fitting
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';
SMF.Fitting.Iterations = 20;

% Thresholding - balanced
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.2;     % Moderate precision
SMF.Thresholding.MinPhotons = 100;
SMF.Thresholding.MinPValue = 0.01;
SMF.Thresholding.MinPSFSigma = 0.8;
SMF.Thresholding.MaxPSFSigma = 1.8;

% Frame connection - critical for dSTORM
SMF.FrameConnection.On = true;
SMF.FrameConnection.MaxSeparation = 1.0;
SMF.FrameConnection.MaxFrameGap = 10;  # Longer gaps (more blinking)
SMF.FrameConnection.MinNFrameConns = 1;  # Don't filter by connections

% Drift correction
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.Method = 'DC-KNN';
```

### For Single Particle Tracking

Moving particles:

```matlab
% Box finding - moderate threshold
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 100;  # Dimmer particles OK
SMF.BoxFinding.BoxOverlap = 2;

% Fitting
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';
SMF.Fitting.Iterations = 20;

% Thresholding - permissive (tracking handles filtering)
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.3;     # More permissive
SMF.Thresholding.MinPhotons = 50;
SMF.Thresholding.MinPValue = 0.005;

% Tracking - core parameters
SMF.Tracking.D = 0.1;                # Diffusion constant (pixels^2/frame)
SMF.Tracking.MaxDistFF = 5;          # Max frame-to-frame distance
SMF.Tracking.MaxDistGC = 10;         # Max gap closing distance
SMF.Tracking.MaxFrameGap = 5;        # Max frames to bridge
SMF.Tracking.MinTrackLength = 10;    # Minimum trajectory length
SMF.Tracking.TrajwiseD = true;       # Estimate D per trajectory

% Disable SMLM-specific features
SMF.FrameConnection.On = false;      # Use Tracking instead
SMF.DriftCorrection.On = false;      # Or enable if needed
```

### For 3D Astigmatism

Z-position from astigmatic PSF:

```matlab
% Box finding
SMF.BoxFinding.BoxSize = 9;  # Slightly larger for elliptical PSF
SMF.BoxFinding.MinPhotons = 250;
SMF.BoxFinding.BoxOverlap = 2;

% Fitting - 3D astigmatism
SMF.Fitting.FitType = 'XYZNB';
SMF.Fitting.PSFSigma = 1.3;  # Initial guess
SMF.Fitting.Iterations = 20;

% Astigmatism calibration (from calibration procedure)
SMF.Fitting.ZFitStruct.Ax = 266.5;
SMF.Fitting.ZFitStruct.Bx = -1533;
SMF.Fitting.ZFitStruct.Ay = 266.5;
SMF.Fitting.ZFitStruct.By = 1533;
SMF.Fitting.ZFitStruct.Gamma = 0.5;
SMF.Fitting.ZFitStruct.D = 0.5;

% Thresholding - include Z precision
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.2;
SMF.Thresholding.MaxZ_SE = 0.05;  # Z precision (micrometers)
SMF.Thresholding.MinPhotons = 150;
SMF.Thresholding.MinPValue = 0.01;

% Frame connection
SMF.FrameConnection.On = true;
SMF.FrameConnection.MaxSeparation = 1.0;
SMF.FrameConnection.MaxFrameGap = 5;

% Drift correction - 3D aware
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.Method = 'DC-KNN';
SMF.DriftCorrection.PixelSizeZUnit = 0.1;  # XY pixel size for 3D
```

## Calculating Pixel Size

Pixel size depends on your optical setup.

### Formula

```
pixel_size_μm = camera_pixel_μm / total_magnification

total_magnification = objective_mag × additional_mag
```

### Common Setups

**Setup 1: Basic**
- Camera: 6.5 μm pixels
- Objective: 100×
- No additional magnification

```matlab
pixel_size = 6.5 / 100;  % = 0.065 μm = 65 nm
SMF.Data.PixelSize = 0.065;
```

**Setup 2: With tube lens**
- Camera: 6.5 μm pixels
- Objective: 100×
- Tube lens: 1.6× magnification

```matlab
pixel_size = 6.5 / (100 * 1.6);  % = 0.0406 μm = 40.6 nm
SMF.Data.PixelSize = 0.0406;
```

**Setup 3: EMCCD typical**
- Camera: 16 μm pixels (common EMCCD)
- Objective: 100×
- Additional: 1.5×

```matlab
pixel_size = 16 / (100 * 1.5);  % = 0.107 μm = 107 nm
SMF.Data.PixelSize = 0.107;
```

### Empirical Measurement

Verify pixel size using known standard:

```matlab
% Image a calibrated ruler or nanoparticle array
% Measure distance in pixels
distance_pixels = 500;  % Example: 500 pixels

% Known physical distance (micrometers)
distance_um = 50;  % Example: 50 μm

% Calculate pixel size
pixel_size = distance_um / distance_pixels;
fprintf('Measured pixel size: %.4f μm = %.1f nm\n', ...
    pixel_size, pixel_size * 1000);

SMF.Data.PixelSize = pixel_size;
```

## Saving and Loading SMF

Preserve configurations for reproducibility.

### Saving SMF

```matlab
% Save SMF to .mat file
save('SMF_PAINT_config.mat', 'SMF');
fprintf('SMF saved to SMF_PAINT_config.mat\n');

% Save with descriptive name
config_file = sprintf('SMF_%s_%s.mat', ...
    SMF.Data.CameraType, datestr(now, 'yyyy-mm-dd'));
save(config_file, 'SMF');
fprintf('SMF saved to %s\n', config_file);

% SMF is automatically saved with results
% When you run analysis, results files contain both SMD and SMF
```

### Loading SMF

```matlab
% Load previously saved SMF
load('SMF_PAINT_config.mat', 'SMF');
fprintf('SMF loaded\n');

% Verify it loaded correctly
fprintf('Camera type: %s\n', SMF.Data.CameraType);
fprintf('Box size: %d\n', SMF.BoxFinding.BoxSize);

% Modify for new experiment
SMF.Data.FileDir = '/new/experiment/path';
SMF.Data.FileName = {'new_data.h5'};
```

### Loading from Results

Analysis results always include the SMF used:

```matlab
% Load results from previous analysis
load('Results.mat', 'SMF', 'SMD');

% SMF contains exact parameters used
% Use as template for new analysis
SMF_new = SMF;  % Copy
SMF_new.Data.FileDir = '/new/data/path';
SMF_new.Data.FileName = {'new_file.h5'};

% Or examine old SMF to understand analysis
fprintf('Previous analysis used:\n');
fprintf('  PSF sigma: %.2f\n', SMF.Fitting.PSFSigma);
fprintf('  Min photons: %.0f\n', SMF.BoxFinding.MinPhotons);
fprintf('  Precision threshold: %.2f pixels\n', SMF.Thresholding.MaxXY_SE);
```

### Creating Templates

Build a library of configurations:

```matlab
% Template for DNA-PAINT
SMF_PAINT = smi_core.SingleMoleculeFitting();
% ... configure for PAINT ...
save('template_PAINT.mat', 'SMF_PAINT');

% Template for dSTORM
SMF_STORM = smi_core.SingleMoleculeFitting();
% ... configure for STORM ...
save('template_STORM.mat', 'SMF_STORM');

% Template for tracking
SMF_tracking = smi_core.SingleMoleculeFitting();
% ... configure for SPT ...
save('template_tracking.mat', 'SMF_tracking');

% Later: Load appropriate template
load('template_PAINT.mat', 'SMF_PAINT');
SMF = SMF_PAINT;
SMF.Data.FileDir = '/my/experiment';
% ... ready to analyze ...
```

## Validating SMF Configuration

Check configuration before analysis.

### Basic Validation

```matlab
% Check required fields are set
assert(~isempty(SMF.Data.FileDir), 'FileDir not set');
assert(~isempty(SMF.Data.FileName), 'FileName not set');
assert(SMF.Data.PixelSize > 0, 'Invalid PixelSize');
assert(SMF.Data.FrameRate > 0, 'Invalid FrameRate');

fprintf('Basic validation passed\n');
```

### Camera Validation

```matlab
% Validate camera parameters
if strcmp(SMF.Data.CameraType, 'EMCCD')
    % Check scalar parameters
    assert(isscalar(SMF.Data.CameraGain), 'EMCCD gain must be scalar');
    assert(isscalar(SMF.Data.CameraOffset), 'EMCCD offset must be scalar');
    assert(SMF.Data.CameraGain > 0, 'Gain must be positive');
    fprintf('EMCCD parameters validated\n');

elseif strcmp(SMF.Data.CameraType, 'SCMOS')
    % Check calibration file exists
    assert(~isempty(SMF.Data.CalibrationFilePath), ...
        'sCMOS requires CalibrationFilePath');
    assert(exist(SMF.Data.CalibrationFilePath, 'file'), ...
        'Calibration file not found');

    % Verify calibration file contents
    calib = load(SMF.Data.CalibrationFilePath);
    assert(isfield(calib, 'CameraGain'), 'Calibration missing CameraGain');
    assert(isfield(calib, 'CameraOffset'), 'Calibration missing CameraOffset');
    assert(isfield(calib, 'CameraNoise'), 'Calibration missing CameraNoise');
    fprintf('sCMOS calibration validated\n');
end
```

### Parameter Reasonableness

```matlab
% Check parameters are in reasonable ranges
warning_count = 0;

% PSF sigma (typical: 0.8-2.0 pixels)
if SMF.Fitting.PSFSigma < 0.5 || SMF.Fitting.PSFSigma > 3
    warning('PSF sigma %.2f seems unusual (typical: 1-2 pixels)', ...
        SMF.Fitting.PSFSigma);
    warning_count = warning_count + 1;
end

% Box size (typical: 5-11 pixels)
if SMF.BoxFinding.BoxSize < 5 || SMF.BoxFinding.BoxSize > 15
    warning('Box size %d seems unusual (typical: 7-9 pixels)', ...
        SMF.BoxFinding.BoxSize);
    warning_count = warning_count + 1;
end

% Pixel size (typical: 0.05-0.15 μm)
if SMF.Data.PixelSize < 0.04 || SMF.Data.PixelSize > 0.2
    warning('Pixel size %.3f μm seems unusual (typical: 0.05-0.15 μm)', ...
        SMF.Data.PixelSize);
    warning_count = warning_count + 1;
end

if warning_count == 0
    fprintf('All parameters in reasonable ranges\n');
else
    fprintf('%d warnings - review configuration\n', warning_count);
end
```

## Troubleshooting

### Issue: GUI won't open

```matlab
% Check MATLAB version (requires R2021a+)
version_str = version;
fprintf('MATLAB version: %s\n', version_str);

% Try command-line configuration instead
SMF = smi_core.SingleMoleculeFitting();
% Set parameters manually
```

### Issue: Calibration file won't load

```matlab
% Check file exists
calib_path = SMF.Data.CalibrationFilePath;
if ~exist(calib_path, 'file')
    error('Calibration file not found: %s', calib_path);
end

% Check file contents
calib = load(calib_path);
required_fields = {'CameraGain', 'CameraOffset', 'CameraNoise'};
for i = 1:length(required_fields)
    if ~isfield(calib, required_fields{i})
        error('Calibration missing required field: %s', required_fields{i});
    end
end

% Check dimensions match
fprintf('Calibration dimensions: %d × %d\n', size(calib.CameraGain));
```

### Issue: Can't remember parameter values

```matlab
% Check current SMF settings
fprintf('\n=== Current SMF Configuration ===\n');
fprintf('Camera: %s\n', SMF.Data.CameraType);
fprintf('Pixel size: %.4f μm\n', SMF.Data.PixelSize);
fprintf('Frame rate: %.1f Hz\n', SMF.Data.FrameRate);
fprintf('Box size: %d pixels\n', SMF.BoxFinding.BoxSize);
fprintf('PSF sigma: %.2f pixels\n', SMF.Fitting.PSFSigma);
fprintf('Fit type: %s\n', SMF.Fitting.FitType);

% Or use GUI to browse all parameters
SMF.gui();
```

## See Also

- [SMF Structure](../core-concepts/smf-structure.md) - Complete parameter reference
- [How to Load Data](load-data.md) - Using configured SMF to load data
- [How to Localize Molecules](localize-molecules.md) - Using SMF for analysis
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete analysis pipeline
- doc/DataStructures/SMF.md - Technical SMF documentation
