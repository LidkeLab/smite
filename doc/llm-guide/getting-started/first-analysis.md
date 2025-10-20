---
title: "Your First SMLM Analysis"
category: "getting-started"
level: "beginner"
tags: ["smlm", "quickstart", "tutorial", "gui"]
prerequisites: ["installation.md"]
related: ["workflows/smlm-analysis.md", "examples/basic-localization.md"]
summary: "Step-by-step guide to running your first single molecule localization microscopy analysis in smite"
estimated_time: "15 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Your First SMLM Analysis

## Purpose

This guide walks you through your first complete SMLM (Single Molecule Localization Microscopy) analysis using smite's graphical interface. You'll learn how to load data, configure analysis parameters, run a test fit, and generate super-resolution images. This is the fastest way to get meaningful results from smite.

## Prerequisites

- smite installed and working (see [installation.md](installation.md))
- Basic understanding of SMLM concepts
- Access to sample data (we'll show you how to generate test data if needed)

## Overview

smite provides a high-level `smi.SMLM` class that handles the complete SMLM analysis pipeline: loading raw camera data, finding molecules, fitting their positions, connecting localizations across frames, correcting drift, and generating visualizations. The GUI makes this accessible without writing code.

The analysis follows this flow:
1. Load raw image stack (.h5 or .mat file)
2. Configure fitting parameters (box size, PSF model, thresholds)
3. Test on a single dataset to verify parameters
4. Run full analysis on all datasets
5. View results and super-resolution images

## Generate Test Data (Optional)

If you don't have real SMLM data yet, generate a test dataset:

```matlab
% Create temporary directory for test data
SaveDir = smi_helpers.mkSMITETmpDir('first_analysis');

% Generate simulated SMLM movie
SimSMLM = smi_sim.SimSMLM();
SimSMLM.NDatasets = 1;
SimSMLM.NFrames = 100;
SimSMLM.FrameSize = [128, 128];
SimSMLM.Rho = 2e-3;  % Emitter density (emitters/pixel^2)
SimSMLM.Photons = 1000;  % Photons per emitter
[SMD, ~, sequence] = SimSMLM.simulateSMLM();

% Save as .mat file
TestFile = fullfile(SaveDir, 'SMLM_testData.mat');
save(TestFile, 'sequence');

fprintf('Test data saved to: %s\n', TestFile);
```

This creates a 100-frame movie with blinking fluorophores at realistic densities.

## Launch the SMLM GUI

Start MATLAB and create an SMLM object without arguments to open the GUI:

```matlab
SMLMobj = smi.SMLM()
```

The GUI window appears with several sections:
- **File Selection**: Choose your data file
- **SMF Parameters**: Configure analysis settings
- **Actions**: Run test fits or full analysis

## Load Your Data

1. **Set File Directory**: Click "Browse" next to "FileDir" and navigate to your data folder (or the `SaveDir` from above)

2. **Select File**: Click "Browse" next to "FileName" and select your .h5 or .mat file

3. **Verify Settings**:
   - FileType: Should auto-detect as 'h5' or 'mat'
   - DataVariable: For .mat files, verify this is 'sequence' (default)
   - PixelSize: Set to your camera pixel size in micrometers (e.g., 0.1)
   - FrameRate: Set to acquisition rate in Hz (e.g., 100)

Alternatively, set these programmatically:

```matlab
SMLMobj.SMF.Data.FileDir = SaveDir;
SMLMobj.SMF.Data.FileName = {'SMLM_testData.mat'};
SMLMobj.SMF.Data.PixelSize = 0.1;  % micrometers
SMLMobj.SMF.Data.FrameRate = 100;  % Hz
```

## Configure Analysis Parameters

The SMF (Single Molecule Fitting) structure controls all analysis. Key parameters to check:

### BoxFinding Parameters

```matlab
SMLMobj.SMF.BoxFinding.BoxSize = 7;      % Fitting box size (pixels)
SMLMobj.SMF.BoxFinding.MinPhotons = 200; % Minimum photons to detect
```

**BoxSize** should be ~4-5× your PSF sigma. For a typical PSF sigma of 1.3 pixels, BoxSize=7 works well.

### Fitting Parameters

```matlab
SMLMobj.SMF.Fitting.PSFSigma = 1.3;  % Expected PSF width (pixels)
SMLMobj.SMF.Fitting.FitType = 'XYNB'; % Fit X, Y, photons (N), background (B)
```

**PSFSigma** depends on your optical setup. For a 100× objective with 100 nm pixel size, typical values are 1.0-1.5 pixels.

### Thresholding Parameters

```matlab
SMLMobj.SMF.Thresholding.On = true;
SMLMobj.SMF.Thresholding.MaxXY_SE = 0.2;    % Max position uncertainty (pixels)
SMLMobj.SMF.Thresholding.MinPhotons = 100;  % Min photons after fitting
SMLMobj.SMF.Thresholding.MinPValue = 0.01;  % Min p-value for fit quality
```

These filters remove poor quality localizations.

### Frame Connection

```matlab
SMLMobj.SMF.FrameConnection.On = true;
SMLMobj.SMF.FrameConnection.MaxSeparation = 1;  % Max distance (pixels)
SMLMobj.SMF.FrameConnection.MaxFrameGap = 5;    % Max frames to bridge
```

Frame connection links the same emitter across multiple frames, improving precision.

### Drift Correction

```matlab
SMLMobj.SMF.DriftCorrection.On = true;
SMLMobj.SMF.DriftCorrection.Method = 'DC-KNN';  % K-nearest neighbors method
```

Corrects stage drift during acquisition.

## Run a Test Fit

Before analyzing all data, test your parameters on a single dataset:

```matlab
% Test dataset 1 (or specify which dataset number)
SMLMobj.SMF.Data.DatasetList = int32(1);

% Run test with verbose output
SMLMobj.VerboseTest = 3;  % Shows detailed plots
SMLMobj.testFit();
```

The test fit will:
1. Load the dataset
2. Find and fit molecules
3. Display diagnostic plots:
   - Raw data with detected boxes
   - Localizations overlaid on image
   - Histograms of photons, background, precision
   - Super-resolution image

**What to look for:**
- Are boxes finding the bright spots?
- Do photon counts look reasonable (100s to 1000s)?
- Are position uncertainties small (< 0.2 pixels)?
- Does the SR image show clear structure?

If results look poor, adjust parameters and run `testFit()` again.

## Run Full Analysis

Once parameters are satisfactory:

```matlab
% Analyze all datasets (leave DatasetList empty or set to all dataset numbers)
SMLMobj.SMF.Data.DatasetList = int32([]);  % Empty = all datasets

% Run complete analysis
SMLMobj.fullAnalysis();
```

This will:
1. Loop over all datasets
2. Localize molecules in each frame
3. Connect localizations across frames
4. Correct for drift
5. Save results to `FileDir/Results/`
6. Generate all plots and SR images

Progress updates appear in the MATLAB command window.

## View Results

Results are saved in your data directory under `Results/`:

```
FileDir/
├── Results/
│   ├── Results.mat           % SMD and SMF structures
│   ├── Dataset_1/
│   │   ├── GaussIm.png       % Super-resolution image
│   │   ├── HistIm.png        % Histogram image
│   │   ├── Photons_hist.png  % Photon distribution
│   │   ├── X_SE_hist.png     % Precision histogram
│   │   └── ...               % Other diagnostic plots
```

Load and inspect results:

```matlab
% Load results
ResultsFile = fullfile(SMLMobj.SMF.Data.FileDir, 'Results', 'Results.mat');
load(ResultsFile, 'SMD', 'SMF');

% Check number of localizations
fprintf('Total localizations: %d\n', numel(SMD.X));

% Check median precision
median_precision = median(SMD.X_SE);
fprintf('Median X precision: %.3f pixels (%.1f nm)\n', ...
    median_precision, median_precision * SMF.Data.PixelSize * 1000);

% View a super-resolution image
imshow(imread(fullfile(SMLMobj.SMF.Data.FileDir, 'Results', 'Dataset_1', 'GaussIm.png')));
```

The **SMD** structure contains all localizations:
- `SMD.X`, `SMD.Y`: Positions (pixels)
- `SMD.Photons`: Detected photons
- `SMD.X_SE`, `SMD.Y_SE`: Position uncertainties
- `SMD.FrameNum`: Frame number of each localization
- `SMD.ConnectID`: Identifier linking same emitter across frames

## Common Issues and Solutions

**Problem: No localizations found**
- Check `BoxFinding.MinPhotons` - may be too high
- Verify your data actually has bright spots (`imagesc(sequence(:,:,1))`)
- Ensure `PSFSigma` is reasonable for your data

**Problem: Too many localizations (including noise)**
- Increase `BoxFinding.MinPhotons`
- Enable thresholding: `SMF.Thresholding.On = true`
- Increase `Thresholding.MinPhotons`

**Problem: Poor localization precision**
- Emitters may be too dim - check photon counts
- Background may be high - check `Bg` values in SMD
- PSF model may be wrong - verify `PSFSigma` matches your data

**Problem: Super-resolution image looks blurry**
- Not enough localizations - acquire more frames
- Drift correction failed - check drift plots
- Localization precision too poor - improve SNR

## Next Steps

Now that you've completed your first analysis:

1. **Understand the data structures**: Read about [SMF](../core-concepts/smf-structure.md) and [SMD](../core-concepts/smd-structure.md)
2. **Learn the full workflow**: See [SMLM Analysis](../workflows/smlm-analysis.md) for detailed explanation
3. **Customize visualization**: Use `smi_vis.GenerateImages` to create custom SR images
4. **Batch processing**: See `smi.Publish` for analyzing many datasets automatically
5. **Advanced analysis**: Explore clustering, drift correction methods, or tracking

## See Also

- [SMLM Analysis Workflow](../workflows/smlm-analysis.md) - Complete SMLM pipeline details
- [SMF Structure](../core-concepts/smf-structure.md) - All analysis parameters
- [SMD Structure](../core-concepts/smd-structure.md) - Results data format
- [How to Load Data](../how-to/load-data.md) - Data loading details
- [Basic Localization Example](../examples/basic-localization.md) - Script-based approach
