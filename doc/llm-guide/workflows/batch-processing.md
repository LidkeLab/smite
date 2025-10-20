---
title: "Batch Processing with smi.Publish"
category: "workflows"
level: "intermediate"
tags: ["batch", "publish", "workflow", "multi-cell", "automation"]
prerequisites: ["../core-concepts/smf-structure.md", "../core-concepts/smd-structure.md", "smlm-analysis.md"]
related: ["../getting-started/first-analysis.md", "spt-tracking.md"]
summary: "Comprehensive guide to batch processing multiple datasets using the smi.Publish class with standard directory structures"
estimated_time: "20 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Batch Processing with smi.Publish

## Purpose

This document explains how to efficiently process multiple super-resolution datasets in batch mode using smite's `smi.Publish` class. You'll learn how to organize data following the standard naming convention, configure batch processing parameters, analyze multiple cells and labels automatically, and aggregate results across entire experiments. This is essential for high-throughput analysis and ensuring consistency across datasets.

## Prerequisites

- Understanding of [SMF structure](../core-concepts/smf-structure.md)
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Completion of [SMLM analysis workflow](smlm-analysis.md)
- Familiarity with basic SMLM analysis
- Understanding of directory structures and file organization

## Overview

The `smi.Publish` class provides automated batch processing capabilities designed specifically for experiments with:

- **Multiple cells** imaged under the same conditions
- **Multiple labels** per cell (e.g., sequential super-resolution)
- **Consistent directory structure** following the naming convention
- **Common analysis parameters** applied across all datasets

Key benefits:
1. **Automation**: Analyze dozens of datasets with a single command
2. **Consistency**: Apply identical parameters to all datasets
3. **Organization**: Results saved to predictable, organized locations
4. **Error handling**: Continues processing even if individual datasets fail
5. **Multi-channel support**: Generates overlay images for multi-label experiments
6. **Aggregation**: Collects statistics across all cells and labels

## Standard Directory Structure

### Required Organization

The `smi.Publish` class expects data organized in a three-level hierarchy:

```
CoverslipDir/
├── Cell_01/
│   ├── Label_1/
│   │   ├── Data_001.h5
│   │   ├── Data_002.h5
│   │   └── ...
│   ├── Label_2/
│   │   ├── Data_001.h5
│   │   └── ...
│   └── ...
├── Cell_02/
│   ├── Label_1/
│   │   └── Data_001.h5
│   └── Label_2/
│       └── Data_001.h5
├── Cell_03/
│   └── ...
└── Results/           (created by Publish)
    ├── Cell_01/
    │   ├── Label_1/
    │   └── Label_2/
    ├── Cell_02/
    └── ...
```

### Naming Conventions

**CoverslipDir**: Top-level directory containing the entire experiment
- Example: `/data/2024-01-10_DNA-PAINT_Actin/`

**Cell directories**: Must match pattern `Cell*`
- `Cell_01`, `Cell_02`, `Cell_03`, ...
- `Cell01`, `Cell02`, `Cell03`, ...
- Numbers extracted automatically using regular expressions

**Label directories**: Must match pattern `Label*`
- `Label_1`, `Label_2`, `Label_3`, ...
- `Label1`, `Label2`, `Label3`, ...
- Used for different fluorescent labels or sequential imaging

**Data files**: Must match pattern `Data*.h5`
- `Data_001.h5`, `Data_002.h5`, ...
- `Data01.h5`, `Data02.h5`, ...
- Multiple datasets per label are supported

### Special Files

**Bleaching datasets**: Files containing `bleach` in name
- Example: `Data_001_bleach.h5`
- Skipped by default (use `AnalyzeBleaching = true` to include)

**Focus images**: Stored within .h5 files
- Used for brightfield drift correction
- Optional but recommended for multi-channel experiments

## Basic Batch Processing

### Example 1: Simple Batch Analysis

Process all cells and labels in a directory:

```matlab
% Prepare SMF with analysis parameters
SMF = smi_core.SingleMoleculeFitting();

% Data configuration
SMF.Data.PixelSize = 0.108;      % 108 nm/pixel
SMF.Data.FrameRate = 100;        % Hz
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;       % ADU/photon
SMF.Data.CameraOffset = 100;     % ADU

% Box finding
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 250;

% Fitting
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Thresholding
SMF.Thresholding.MaxXY_SE = 0.15;
SMF.Thresholding.MinPhotons = 200;

% Frame connection
SMF.FrameConnection.On = true;
SMF.FrameConnection.MaxSeparation = 1.0;

% Drift correction
SMF.DriftCorrection.On = true;

% Create and configure Publish object
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/2024-01-10_Experiment/';
Pub.Verbose = 1;
Pub.GenerateSR = true;
Pub.SRImageZoom = 20;

% Run batch analysis
Pub.performFullAnalysis();

% Results saved to:
% /data/2024-01-10_Experiment/Results/
```

### Example 2: sCMOS Camera with Calibration

Use pixel-wise calibration for sCMOS cameras:

```matlab
% Prepare SMF
SMF = smi_core.SingleMoleculeFitting();

% sCMOS-specific configuration
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = '/calibrations/sCMOS_2024-01-05.mat';
SMF.Data.PixelSize = 0.0954;  % 95.4 nm/pixel

% Analysis parameters
SMF.Data.AnalysisID = 'dSTORM_run1';
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNBS';  % Fit PSF sigma
SMF.FrameConnection.Method = 'LAP-FC';
SMF.FrameConnection.On = true;
SMF.DriftCorrection.On = true;

% Thresholding for dSTORM
SMF.Thresholding.MaxPSFSigma = 1.5;
SMF.Thresholding.MinPSFSigma = 0.9;
SMF.Thresholding.MinPhotons = 200;
SMF.Thresholding.MaxXY_SE = 0.15;

% Configure Publish
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/dSTORM_experiment/';
Pub.Verbose = 1;
Pub.GenerateSR = true;
Pub.GenerateImagingStats = true;
Pub.SRImageZoom = 20;

% Run analysis
Pub.performFullAnalysis();
```

## Advanced Configuration

### Selective Processing

#### Process Specific Cells

```matlab
% Only analyze Cell_02 and Cell_05
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/experiment/';
Pub.CellList = [2, 5];  % Cell numbers to analyze
Pub.performFullAnalysis();
```

#### Process Specific Labels

```matlab
% Only analyze Label_1 (e.g., first channel only)
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/experiment/';
Pub.LabelID = 1;  % Label numbers to analyze
Pub.performFullAnalysis();

% Analyze Label_1 and Label_3 only
Pub.LabelID = [1, 3];
Pub.performFullAnalysis();
```

#### Combine Cell and Label Selection

```matlab
% Analyze only Label_1 from Cell_01, Cell_03, and Cell_04
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/experiment/';
Pub.CellList = [1, 3, 4];
Pub.LabelID = 1;
Pub.performFullAnalysis();
```

### Different Parameters per Label

Use `SMFperLabel` to apply different parameters to each label:

```matlab
% Prepare SMF for Label_1 (DNA-PAINT with bright signal)
SMF1 = smi_core.SingleMoleculeFitting();
SMF1.Data.PixelSize = 0.108;
SMF1.Data.CameraGain = 2.5;
SMF1.Data.CameraOffset = 100;
SMF1.BoxFinding.MinPhotons = 500;  % Higher threshold
SMF1.Thresholding.MinPhotons = 400;
SMF1.Fitting.PSFSigma = 1.3;
SMF1.FrameConnection.On = true;
SMF1.DriftCorrection.On = true;

% Prepare SMF for Label_2 (dSTORM with dimmer signal)
SMF2 = copy(SMF1);  % Start from SMF1
SMF2.BoxFinding.MinPhotons = 200;  % Lower threshold for dim label
SMF2.Thresholding.MinPhotons = 150;
SMF2.Fitting.PSFSigma = 1.4;  % Slightly different PSF

% NOTE: These fields must match across all SMFs:
% - Data.FileDir, Data.FileName (set by Publish)
% - Data.ResultsDir (set by Publish)
% - Data.AnalysisID
% - Data.PixelSize
% - DriftCorrection.On

% Create Publish object
Pub = smi.Publish(SMF1);  % Use SMF1 as default
Pub.SMFperLabel = {SMF1, SMF2};  % Specify per-label SMFs
Pub.CoverslipDir = '/data/two-color-experiment/';
Pub.Verbose = 1;
Pub.GenerateSR = true;
Pub.GenerateOverlayStats = true;  % Generate two-color overlays
Pub.performFullAnalysis();
```

### Control Output Generation

Fine-tune which outputs are generated:

```matlab
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/experiment/';

% Control what gets generated
Pub.GenerateSR = true;              % Generate SR images (default: true)
Pub.GenerateImagingStats = true;    % Generate imaging statistics (default: true)
Pub.GenerateOverlayStats = false;   % Skip overlay statistics (default: false)
Pub.AnalyzeBleaching = false;       % Skip bleaching datasets (default: false)

% Control SR image zoom factors
Pub.SRImageZoom = 20;               % Standard SR images (default: 20)
Pub.SRCircleImageZoom = 50;         % Circle images showing precision (default: 50)

% Verbosity level
Pub.Verbose = 0;  % Minimal output
Pub.Verbose = 1;  % Standard progress messages (default)
Pub.Verbose = 2;  % Detailed output

Pub.performFullAnalysis();
```

### Custom Results Directory

Specify where results should be saved:

```matlab
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/experiment/';

% Default: Results saved to CoverslipDir/Results/
% Custom location:
Pub.SaveBaseDir = '/data/experiment/Analysis_2024-01-10/';

Pub.performFullAnalysis();
```

## Multi-Channel Analysis

### Two-Color Overlay Generation

For multi-label experiments (e.g., sequential super-resolution):

```matlab
% Configure SMF for two-color experiment
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.PixelSize = 0.108;
SMF.Data.CameraGain = 2.5;
SMF.Data.CameraOffset = 100;
SMF.Data.AnalysisID = 'TwoColor';

% Standard analysis parameters
SMF.BoxFinding.MinPhotons = 300;
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';
SMF.Thresholding.MaxXY_SE = 0.15;
SMF.FrameConnection.On = true;
SMF.DriftCorrection.On = true;

% Configure Publish for overlay generation
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/two-color/';
Pub.GenerateSR = true;
Pub.GenerateOverlayStats = true;  % Enable overlay generation

% Run analysis
Pub.performFullAnalysis();

% Overlays saved to:
% CoverslipDir/Results/Cell_XX_Overlay_*.png
```

### Brightfield Registration

Use brightfield images for registration between channels:

```matlab
% Configure SMF with drift correction
SMF = smi_core.SingleMoleculeFitting();
% ... (configure SMF as above) ...
SMF.DriftCorrection.On = true;  % Enable post-processing DC

% Configure Publish with brightfield DC
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/experiment/';

% Enable brightfield drift correction
Pub.UseBrightfieldDC = true;  % Apply brightfield DC before post-processing DC
Pub.GenerateImagingStats = true;  % Generate brightfield registration stats

% Optional: Define maximum acceptable shift for masking
% Regions with shifts > MaxBrightfieldShift are masked out
Pub.MaxBrightfieldShift = 0.2 / SMF.Data.PixelSize;  % 200 nm in pixels

Pub.performFullAnalysis();

% Registration results saved to:
% CoverslipDir/Results/Cell_XX/Label_Y/Data_ZZZ/AlignReg_*
```

### Registration-Based Coordinate Shifting

Shift coordinates to best-registered dataset:

```matlab
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/experiment/';

% Method 1: Shift based on brightfield registration
Pub.ShiftToReg = true;  % Use brightfield focus images
Pub.GenerateSR = true;

% Method 2: Shift using drift correction registration
Pub.ShiftToRegViaDC = true;  % Use DC-based registration
% Note: Label_1 is used as reference, other labels shifted to match

Pub.performFullAnalysis();
```

## Results Organization

### Output Directory Structure

After running `performFullAnalysis()`, results are organized as:

```
CoverslipDir/Results/
├── Cell_01/
│   ├── Label_1/
│   │   ├── Data_001_AnalysisID/
│   │   │   ├── Data_001_AnalysisID_Results.mat    % SMD + SMF
│   │   │   ├── Data_001_AnalysisID_GaussIm.png   % Gaussian SR image
│   │   │   ├── Data_001_AnalysisID_HistIm.png    % Histogram SR image
│   │   │   ├── Data_001_AnalysisID_CircleIm.png  % Precision circles
│   │   │   ├── Data_001_AnalysisID_Drift.png     % Drift plot
│   │   │   └── Data_001_AnalysisID_*.png         % Various diagnostics
│   │   ├── Data_002_AnalysisID/
│   │   │   └── ...
│   │   └── ...
│   ├── Label_2/
│   │   └── ...
│   └── ...
├── Cell_02/
│   └── ...
├── ResultsStructs/                               % Centralized results
│   ├── Cell_01_Label_1_Data_001_Results.mat
│   ├── Cell_01_Label_2_Data_001_Results.mat
│   ├── Cell_02_Label_1_Data_001_Results.mat
│   └── ...
├── Cell_01_Overlay_*.png                        % Multi-channel overlays
├── Cell_02_Overlay_*.png
├── ResultsStruct.mat                            % Aggregated statistics
└── Log_*.mat                                    % Error log
```

### Key Output Files

**Individual results**: `Cell_XX/Label_Y/Data_ZZZ_AnalysisID/`
- Contains SMD, SMF, and visualizations for each dataset
- Organized by cell, label, and data file

**Centralized results**: `ResultsStructs/`
- Copies of Results.mat files with standardized names
- Easy access to all results without navigating directory tree

**Overlay images**: `Cell_XX_Overlay_*.png`
- Multi-color overlay images for cells with multiple labels
- Automatically generated when `GenerateOverlayStats = true`

**Aggregate statistics**: `ResultsStruct.mat`
- Contains `ResultsStruct` array with statistics from all cells/labels
- Useful for cross-dataset comparisons

**Error log**: `Log_*.mat`
- Records any errors encountered during processing
- Contains `ErrorLog` cell array with details
- Includes start and end times

### Accessing Results

```matlab
% Load aggregate results
load('CoverslipDir/Results/ResultsStruct.mat', 'ResultsStruct');

% ResultsStruct is organized as:
% ResultsStruct(CellNumber, LabelNumber)

% Access specific cell/label
Cell1_Label1 = ResultsStruct(1, 1);

% View available fields
fieldnames(Cell1_Label1)

% Load individual SMD
load('CoverslipDir/Results/ResultsStructs/Cell_01_Label_1_Data_001_Results.mat');
% Now have: SMD, SMF

% Examine localization count
fprintf('Cell 1, Label 1: %d localizations\n', length(SMD.X));
```

## Parallel Processing

### Using MATLAB Parallel Pool

Speed up batch processing with parallel workers:

```matlab
% Start parallel pool with 4 workers
parpool('local', 4);

% Prepare SMF
SMF = smi_core.SingleMoleculeFitting();
% ... (configure SMF) ...

% Run batch analysis
% Note: Publish itself runs serially through cells,
% but individual SMLM analyses can use parallel processing
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/experiment/';
Pub.performFullAnalysis();

% Clean up
delete(gcp('nocreate'));
```

### Process Multiple Coverslips

Process multiple independent coverslip directories:

```matlab
% Define coverslip directories
CoverslipDirs = {
    '/data/2024-01-10_Sample1/'
    '/data/2024-01-11_Sample2/'
    '/data/2024-01-12_Sample3/'
    '/data/2024-01-13_Sample4/'
};

% Prepare common SMF
SMF = smi_core.SingleMoleculeFitting();
% ... (configure SMF) ...

% Process each coverslip
for ii = 1:length(CoverslipDirs)
    fprintf('Processing coverslip %d of %d: %s\n', ...
        ii, length(CoverslipDirs), CoverslipDirs{ii});

    Pub = smi.Publish(SMF);
    Pub.CoverslipDir = CoverslipDirs{ii};
    Pub.Verbose = 1;
    Pub.GenerateSR = true;

    try
        Pub.performFullAnalysis();
        fprintf('Successfully completed %s\n', CoverslipDirs{ii});
    catch ME
        warning('Failed to process %s: %s', CoverslipDirs{ii}, ME.message);
        continue  % Continue with next coverslip
    end
end

fprintf('All coverslips processed!\n');
```

## Workflow Control Flow

### Internal Processing Hierarchy

Understanding the processing flow helps with debugging:

```
performFullAnalysis()
  │
  ├─► Find Cell* directories
  ├─► Filter by CellList (if specified)
  │
  └─► For each Cell:
       │
       processCell(CellName)
         │
         ├─► Find Label* directories
         ├─► Filter by LabelID (if specified)
         │
         └─► For each Label:
              │
              processLabel(CellName, LabelName)
                │
                ├─► Find Data*.h5 files
                ├─► Skip bleaching files (unless AnalyzeBleaching=true)
                │
                └─► For each Data file:
                     │
                     ├─► Load focus images (if present)
                     ├─► Generate brightfield registration stats
                     │   (if GenerateImagingStats=true)
                     │
                     ├─► Run SMLM analysis:
                     │    ├─► Load raw data
                     │    ├─► Localize molecules
                     │    ├─► Frame connection
                     │    ├─► Drift correction
                     │    └─► Generate plots/images
                     │
                     ├─► Apply brightfield DC (if UseBrightfieldDC=true)
                     ├─► Shift coordinates (if ShiftToReg=true)
                     └─► Save results
       │
       └─► Generate multi-label overlay (if NLabels > 1)
  │
  ├─► Generate overlay statistics (if GenerateOverlayStats=true)
  └─► Save aggregate ResultsStruct
```

### Error Handling

The Publish class continues processing even when errors occur:

```matlab
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/experiment/';
Pub.performFullAnalysis();

% Check error log after completion
load([Pub.SaveBaseDir, '/Log_*.mat'], 'ErrorLog');

% ErrorLog is cell array with columns:
% {CellName, LabelName, MException}

if ~isempty(ErrorLog)
    fprintf('Encountered %d errors:\n', size(ErrorLog, 1));
    for ii = 1:size(ErrorLog, 1)
        fprintf('  %s/%s: %s\n', ...
            ErrorLog{ii,1}, ErrorLog{ii,2}, ErrorLog{ii,3}.message);
    end
else
    fprintf('All datasets processed successfully!\n');
end
```

## Complete Examples

### Example 1: DNA-PAINT Batch Processing

```matlab
% DNA-PAINT with high photon counts
SMF = smi_core.SingleMoleculeFitting();

% Camera and data
SMF.Data.PixelSize = 0.108;
SMF.Data.FrameRate = 10;  % 10 Hz for DNA-PAINT
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;
SMF.Data.CameraOffset = 100;
SMF.Data.AnalysisID = 'PAINT_Analysis';

% Box finding (high photon threshold)
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 500;

% Fitting
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Thresholding (strict for DNA-PAINT)
SMF.Thresholding.MaxXY_SE = 0.10;  % 10.8 nm
SMF.Thresholding.MinPhotons = 400;
SMF.Thresholding.AutoThreshLogL = true;

% Frame connection (DNA-PAINT specific)
SMF.FrameConnection.On = true;
SMF.FrameConnection.MaxSeparation = 1.0;
SMF.FrameConnection.MaxFrameGap = 3;
SMF.FrameConnection.MinNFrameConns = 2;  % Must appear at least twice

% Drift correction
SMF.DriftCorrection.On = true;
SMF.DriftCorrection.Method = 'DC-KNN';

% Configure Publish
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/DNA-PAINT_Microtubules/';
Pub.Verbose = 1;
Pub.GenerateSR = true;
Pub.GenerateImagingStats = true;
Pub.SRImageZoom = 20;

% Run analysis
Pub.performFullAnalysis();
```

### Example 2: Multi-Color dSTORM

```matlab
% dSTORM with two colors (sequential imaging)
SMF1 = smi_core.SingleMoleculeFitting();

% Common settings for both channels
SMF1.Data.CameraType = 'SCMOS';
SMF1.Data.CalibrationFilePath = '/cal/sCMOS_cal.mat';
SMF1.Data.PixelSize = 0.0954;
SMF1.Data.AnalysisID = 'TwoColor_dSTORM';

% Channel 1 (e.g., Alexa647 - bright)
SMF1.BoxFinding.MinPhotons = 250;
SMF1.Fitting.PSFSigma = 1.3;
SMF1.Fitting.FitType = 'XYNBS';
SMF1.Thresholding.MinPhotons = 200;
SMF1.Thresholding.MaxXY_SE = 0.15;
SMF1.Thresholding.MinPSFSigma = 0.9;
SMF1.Thresholding.MaxPSFSigma = 1.5;
SMF1.FrameConnection.On = true;
SMF1.DriftCorrection.On = true;

% Channel 2 (e.g., Alexa555 - dimmer)
SMF2 = copy(SMF1);
SMF2.BoxFinding.MinPhotons = 200;  % Lower threshold
SMF2.Thresholding.MinPhotons = 150;

% Configure Publish for two channels
Pub = smi.Publish(SMF1);
Pub.SMFperLabel = {SMF1, SMF2};
Pub.CoverslipDir = '/data/dSTORM_Two_Color/';
Pub.Verbose = 1;
Pub.GenerateSR = true;
Pub.GenerateImagingStats = true;
Pub.GenerateOverlayStats = true;  % Generate two-color overlays
Pub.UseBrightfieldDC = true;      % Use brightfield registration
Pub.ShiftToRegViaDC = true;       % Align channels via drift correction
Pub.SRImageZoom = 20;

% Run analysis
Pub.performFullAnalysis();

% Check for overlay images
overlays = dir([Pub.SaveBaseDir, '/Cell_*_Overlay_*.png']);
fprintf('Generated %d overlay images\n', length(overlays));
```

### Example 3: Selective Reprocessing

Reprocess specific cells with updated parameters:

```matlab
% Load existing SMF and modify parameters
load('/data/experiment/Results/Cell_01/Label_1/Data_001_Old/Results.mat', 'SMF');
SMF = smi_core.SingleMoleculeFitting.reloadSMF(SMF);

% Update analysis parameters
SMF.Data.AnalysisID = 'Reprocessed_v2';
SMF.Thresholding.MaxXY_SE = 0.12;  % Stricter threshold
SMF.Thresholding.MinPhotons = 250;  % Higher photon requirement

% Reprocess only Cell_01 and Cell_03
Pub = smi.Publish(SMF);
Pub.CoverslipDir = '/data/experiment/';
Pub.CellList = [1, 3];  % Only these cells
Pub.Verbose = 1;
Pub.GenerateSR = true;
Pub.SRImageZoom = 20;

% Run reprocessing
Pub.performFullAnalysis();

% New results saved with new AnalysisID:
% Results/Cell_01/Label_1/Data_001_Reprocessed_v2/
```

## Troubleshooting

### Issue: No Cell directories found

**Symptoms:**
```
Publish.performFullAnalysis(): 0 cell directories found
```

**Diagnose:**
```matlab
% Check directory contents
dir(Pub.CoverslipDir)

% Check for Cell* pattern
cellDirs = dir(fullfile(Pub.CoverslipDir, 'Cell*'));
cellDirs.name
```

**Solutions:**
- Verify `CoverslipDir` path is correct
- Ensure cell directories match `Cell*` pattern (case-sensitive on Linux)
- Check directory exists: `isfolder(Pub.CoverslipDir)`

### Issue: No Data files found in Label directory

**Symptoms:**
```
Publish.processLabel(): 0 files found
```

**Diagnose:**
```matlab
labelDir = fullfile(Pub.CoverslipDir, 'Cell_01', 'Label_1');
dataFiles = dir(fullfile(labelDir, 'Data*'));
dataFiles.name
```

**Solutions:**
- Verify data files match `Data*` pattern
- Check file extensions (.h5 files expected)
- Ensure files are not hidden or locked

### Issue: Processing fails for specific datasets

**Symptoms:**
Warning messages during processing, some cells/labels skipped

**Diagnose:**
```matlab
% Load error log after completion
logFiles = dir(fullfile(Pub.SaveBaseDir, 'Log_*.mat'));
load(fullfile(logFiles(1).folder, logFiles(1).name), 'ErrorLog');

% Display errors
for ii = 1:size(ErrorLog, 1)
    fprintf('Error %d: %s/%s\n', ii, ErrorLog{ii,1}, ErrorLog{ii,2});
    fprintf('  Message: %s\n', ErrorLog{ii,3}.message);
    fprintf('  Identifier: %s\n', ErrorLog{ii,3}.identifier);
end
```

**Solutions:**
- Check individual error messages for specific issues
- Verify data files are not corrupted
- Ensure camera calibration files exist
- Check GPU availability for fitting operations

### Issue: Overlay images not generated

**Symptoms:**
No overlay images in Results directory

**Diagnose:**
```matlab
% Check configuration
fprintf('GenerateSR: %d\n', Pub.GenerateSR);
fprintf('GenerateOverlayStats: %d\n', Pub.GenerateOverlayStats);

% Check if multiple labels exist
cellDir = fullfile(Pub.CoverslipDir, 'Cell_01');
labelDirs = dir(fullfile(cellDir, 'Label*'));
fprintf('Number of labels: %d\n', length(labelDirs));
```

**Solutions:**
- Set `Pub.GenerateOverlayStats = true`
- Ensure cells have multiple labels (overlays only for multi-label data)
- Check that all labels processed successfully
- Verify sufficient memory for image generation

### Issue: Brightfield drift correction fails

**Symptoms:**
Warnings about missing FocusImages or AlignReg structures

**Diagnose:**
```matlab
% Check if focus images present in h5 file
filePath = fullfile(Pub.CoverslipDir, 'Cell_01', 'Label_1', 'Data_001.h5');
try
    [~, ~, ~, ~, groups] = smi_core.LoadData.seqH5Data(filePath);
    fprintf('Groups found:\n');
    disp(groups);
catch ME
    fprintf('Error reading file: %s\n', ME.message);
end
```

**Solutions:**
- Set `Pub.UseBrightfieldDC = false` if focus images not available
- Ensure data files were collected with brightfield imaging enabled
- Use standard drift correction only (`SMF.DriftCorrection.On = true`)

## Performance Optimization

### Memory Management

For large experiments (many cells/labels):

```matlab
% Process in chunks
allCells = 1:20;  % 20 cells total
chunkSize = 5;

for ii = 1:chunkSize:length(allCells)
    cellChunk = allCells(ii:min(ii+chunkSize-1, end));

    Pub = smi.Publish(SMF);
    Pub.CoverslipDir = '/data/large_experiment/';
    Pub.CellList = cellChunk;
    Pub.Verbose = 1;
    Pub.GenerateSR = true;

    Pub.performFullAnalysis();

    % Clear variables between chunks
    clear Pub
end
```

### Disk Space Considerations

Monitor disk usage during batch processing:

```matlab
% Estimate output size before running
% Typical outputs per dataset:
% - Results.mat: 10-100 MB (depends on localization count)
% - Images: ~5 MB per image × ~10 images = 50 MB
% - Total: ~100-200 MB per dataset

% Calculate expected space
nCells = 10;
nLabelsPerCell = 2;
nDatasetsPerLabel = 1;
totalDatasets = nCells * nLabelsPerCell * nDatasetsPerLabel;
estimatedSpace_GB = totalDatasets * 0.2;  % 200 MB per dataset

fprintf('Estimated disk space needed: %.1f GB\n', estimatedSpace_GB);

% Check available space before running
[status, result] = system(['df -h ', Pub.SaveBaseDir]);
fprintf('Available space:\n%s\n', result);
```

## See Also

- [SMLM Analysis Workflow](smlm-analysis.md) - Understanding individual analysis
- [SPT Tracking](spt-tracking.md) - Single particle tracking batch processing
- [SMF Structure](../core-concepts/smf-structure.md) - All parameters
- [SMD Structure](../core-concepts/smd-structure.md) - Results format
- MATLAB/+smi/@Publish/README.md - Class documentation
- MATLAB/examples/Example_Publish_generic.m - Generic example script
- MATLAB/examples/Example_Publish_SMF.m - SMF-based example script
