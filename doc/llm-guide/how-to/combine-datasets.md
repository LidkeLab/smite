---
title: "How to Combine Multiple Datasets"
category: "how-to"
level: "intermediate"
tags: ["smd", "tr", "concatenation", "datasets", "drift", "alignment"]
prerequisites: ["../core-concepts/smd-structure.md", "../core-concepts/tr-structure.md"]
related: ["save-and-load-smd.md", "../workflows/smlm-analysis.md"]
summary: "Learn how to properly combine multiple SMD/TR datasets while preserving metadata and handling drift"
estimated_time: "15 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# How to Combine Multiple Datasets

## Purpose

When analyzing multiple datasets together (e.g., combining replicates, pooling channels, or assembling time-series data), you need to concatenate SMD or TR structures while preserving their integrity. This guide explains how to properly combine datasets using smite's concatenation tools, handle the DatasetNum field, manage drift correction across datasets, and align temporal and spatial coordinates.

## Prerequisites

- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Understanding of [TR structure](../core-concepts/tr-structure.md)
- Basic knowledge of drift correction concepts

## Overview

smite provides two main concatenation methods:

1. **`smi_core.SingleMoleculeData.catSMD`**: Concatenates two SMD structures
2. **`smi_core.TrackingResults.catTR`**: Concatenates two TR arrays

Both methods handle field padding, dataset numbering, and special fields (ConnectID, IndSMD, drift arrays) automatically.

## Key Concepts

### DatasetNum Field

The `DatasetNum` field identifies which dataset each localization originated from. When combining datasets, you must decide whether they represent:

- **New data** (`NewDataFlag=true`): Dataset numbers are incremented, treating them as distinct experiments
- **Same data** (`NewDataFlag=false`): Dataset numbers are preserved, treating them as parts of the same experiment

### Drift Arrays

Drift arrays (`DriftX`, `DriftY`, `DriftZ`) have dimensions `[NFrames × NDatasets]`. When concatenating:

- Drift arrays are concatenated horizontally (column-wise)
- Each dataset maintains its own drift trajectory
- Post-concatenation, drift can be referenced by dataset number

### ConnectID and IndSMD

These fields maintain indices that track frame-connected localizations:

- **ConnectID**: Unique identifier for each emitter across frames
- **IndSMD**: Indices in pre-frame-connection SMD corresponding to combined localizations

When concatenating, these indices are automatically offset to maintain uniqueness.

## Combining SMD Structures

### Basic Concatenation

```matlab
% Load two datasets
SMD1 = load('dataset1_Results.mat', 'SMD').SMD;
SMD2 = load('dataset2_Results.mat', 'SMD').SMD;

% Concatenate treating them as new distinct datasets
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2);

% Check results
fprintf('SMD1: %d localizations, %d datasets\n', length(SMD1.X), SMD1.NDatasets);
fprintf('SMD2: %d localizations, %d datasets\n', length(SMD2.X), SMD2.NDatasets);
fprintf('Combined: %d localizations, %d datasets\n', ...
    length(SMD_combined.X), SMD_combined.NDatasets);
```

Default behavior:
- `NewDataFlag = true`: SMD2 dataset numbers are incremented by SMD1.NDatasets
- `ShowWarnings = true`: Warns about field mismatches

### NewDataFlag: New vs. Same Data

**Treating as new data (default):**

```matlab
% Combining separate experiments
SMD1.NDatasets = 3;  % Datasets 1, 2, 3
SMD2.NDatasets = 2;  % Originally datasets 1, 2

SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2, true);

% Result: SMD2 datasets are renumbered to 4, 5
% SMD_combined.NDatasets = 5
% SMD_combined.DatasetNum contains [1,2,3,4,5]
```

This is appropriate when:
- Combining data from different days
- Pooling biological replicates
- Merging experiments with independent drift

**Treating as same data:**

```matlab
% Combining subsets of the same experiment
% (e.g., after processing in chunks)
SMD1.DatasetNum = [1; 1; 2; 2];  % 4 localizations from datasets 1, 2
SMD2.DatasetNum = [1; 2];        % 2 localizations from datasets 1, 2

SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2, false);

% Result: Dataset numbers preserved
% SMD_combined.DatasetNum = [1; 1; 2; 2; 1; 2]
% SMD_combined.NDatasets = 2
```

This is appropriate when:
- Recombining filtered subsets
- Merging ROIs from same dataset
- Concatenating pre/post processed data

### Handling Field Mismatches

If SMD structures have different fields, catSMD handles this automatically:

```matlab
% SMD1 has field 'CustomField', SMD2 doesn't
SMD1.CustomField = ones(length(SMD1.X), 1);

% Concatenate with warnings
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2);
% Warning: 'The following SMD1 fields aren't present in SMD2: CustomField'

% Result: SMD_combined.CustomField = SMD1.CustomField (only from SMD1)

% Suppress warnings if expected
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2, true, false);
```

### Combining Multiple Datasets

To combine more than two datasets:

```matlab
% Method 1: Sequential concatenation
SMD_all = smi_core.SingleMoleculeData.createSMD();
file_list = {'data1.mat', 'data2.mat', 'data3.mat', 'data4.mat'};

for ii = 1:length(file_list)
    data = load(file_list{ii}, 'SMD');
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, data.SMD, true);
    fprintf('Added dataset %d: now %d localizations\n', ii, length(SMD_all.X));
end

fprintf('\nFinal combined dataset:\n');
fprintf('  Total localizations: %d\n', length(SMD_all.X));
fprintf('  Total datasets: %d\n', SMD_all.NDatasets);
fprintf('  Dataset range: %d to %d\n', min(SMD_all.DatasetNum), max(SMD_all.DatasetNum));
```

```matlab
% Method 2: Pairwise combination
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2);
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD_combined, SMD3);
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD_combined, SMD4);
```

Both methods produce identical results. Method 1 is more flexible for variable numbers of files.

## Combining TR Arrays

TR (Tracking Results) structures are arrays where each element represents one trajectory. Concatenation is simpler than SMD but requires field compatibility.

### Basic TR Concatenation

```matlab
% Load tracking results from two experiments
TR1 = load('tracking1_Results.mat', 'TR').TR;
TR2 = load('tracking2_Results.mat', 'TR').TR;

% Concatenate
TR_combined = smi_core.TrackingResults.catTR(TR1, TR2);

fprintf('TR1: %d trajectories\n', length(TR1));
fprintf('TR2: %d trajectories\n', length(TR2));
fprintf('Combined: %d trajectories\n', length(TR_combined));
```

### Field Padding

If TR structures have different fields, `catTR` automatically pads missing fields:

```matlab
% TR1 has field 'DiffusionCoeff', TR2 doesn't
TR1(1).DiffusionCoeff = 0.5;

TR_combined = smi_core.TrackingResults.catTR(TR1, TR2);

% Result: TR2 elements get DiffusionCoeff = [] (empty)
% All elements now share the same fields
```

This uses `padTR` internally:

```matlab
% Manual padding (rarely needed)
TR2_padded = smi_core.TrackingResults.padTR(TR2, TR1);
% TR2_padded now has all fields from TR1
```

### Orientation Handling

`catTR` ensures both inputs are column arrays:

```matlab
TR1 = TR_row.';  % Force column orientation
TR2 = TR_col;    % Already column

TR_combined = smi_core.TrackingResults.catTR(TR1, TR2, true);
% Result: Always column array

% Disable automatic orientation check
TR_combined = smi_core.TrackingResults.catTR(TR1, TR2, false);
% Use only if you're certain of orientations (faster)
```

### Combining Multiple TR Arrays

```matlab
% Initialize with empty TR
TR_all = smi_core.TrackingResults.createTR();
tracking_files = dir('Tracking_*.mat');

for ii = 1:length(tracking_files)
    data = load(tracking_files(ii).name, 'TR');
    TR_all = smi_core.TrackingResults.catTR(TR_all, data.TR, false);
    fprintf('Loaded %s: %d trajectories\n', tracking_files(ii).name, length(data.TR));
end

fprintf('\nTotal trajectories: %d\n', length(TR_all));
```

## Drift Correction Across Datasets

### Understanding Drift Arrays

Drift arrays store the drift correction for each frame and dataset:

```matlab
% Example SMD with 2 datasets, 100 frames each
SMD.NDatasets = 2;
SMD.NFrames = 100;
SMD.DriftX = zeros(100, 2);  % [NFrames × NDatasets]
SMD.DriftY = zeros(100, 2);

% Each column represents one dataset's drift trajectory
% SMD.DriftX(:, 1) = drift for dataset 1
% SMD.DriftX(:, 2) = drift for dataset 2
```

### Concatenating Drift Arrays

When using `catSMD`, drift arrays are concatenated horizontally:

```matlab
% SMD1: 2 datasets, 100 frames
% SMD2: 3 datasets, 100 frames

SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2, true);

% Result: SMD_combined.DriftX is [100 × 5]
% Columns 1-2: from SMD1
% Columns 3-5: from SMD2 (renumbered to datasets 3-5)

fprintf('Combined drift array size: %d frames × %d datasets\n', ...
    size(SMD_combined.DriftX, 1), size(SMD_combined.DriftX, 2));
```

### Per-Dataset Drift Correction

When combining datasets, each maintains independent drift:

```matlab
% Apply intra-dataset drift correction before concatenation
SMD_all = smi_core.SingleMoleculeData.createSMD();

for ii = 1:NDatasets
    % Load dataset
    SMD_i = load(sprintf('dataset%d.mat', ii), 'SMD').SMD;

    % Perform intra-dataset drift correction
    DC = smi_core.DriftCorrection(SMF, SMD_i);
    [SMD_corrected, ~] = DC.driftCorrectKNNIntra(SMD_i, ii, ii);

    % Concatenate
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_corrected, false);
end

% Each dataset has its own drift trajectory in SMD_all
```

### Inter-Dataset Drift Correction

For global alignment across datasets:

```matlab
% First, perform intra-dataset drift correction
SMD_intra = smi_core.SingleMoleculeData.createSMD();
DC = smi_core.DriftCorrection(SMF, SMD_all);

for ii = 1:SMD_all.NDatasets
    [SMD_intra_i, ~] = DC.driftCorrectKNNIntra(SMD_all, ii, ii);
    SMD_intra = smi_core.SingleMoleculeData.catSMD(SMD_intra, SMD_intra_i, false);
end

% Then, perform inter-dataset drift correction
[SMD_final, Statistics] = DC.driftCorrectKNNInter(SMD_intra);

% SMD_final is now aligned across all datasets
fprintf('Inter-dataset drift correction complete\n');
fprintf('Final RMSE: %.3f pixels\n', Statistics.RMSE);
```

### Extracting Specific Datasets

After concatenation, you can extract subsets:

```matlab
% Extract datasets 2 and 5 from combined SMD
datasets_to_extract = [2, 5];
[SMD_subset, keep_bool] = smi_core.SingleMoleculeData.extractDatasets(...
    SMD_combined, datasets_to_extract);

% Optionally compress dataset numbers (2,5 → 1,2)
[SMD_subset, keep_bool] = smi_core.SingleMoleculeData.extractDatasets(...
    SMD_combined, datasets_to_extract, true);

fprintf('Extracted %d localizations from %d datasets\n', ...
    length(SMD_subset.X), SMD_subset.NDatasets);
```

## Temporal and Spatial Alignment

### Frame Number Consistency

Frame numbers are preserved independently for each dataset:

```matlab
% SMD1: Dataset 1, frames 1-100
% SMD2: Dataset 1, frames 1-80

SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2, false);

% Both still have frames 1-100 (SMD1) and 1-80 (SMD2)
% FrameNum is per-dataset, not global

% To get global frame indices:
global_frame = SMD_combined.NFrames * (SMD_combined.DatasetNum - 1) + SMD_combined.FrameNum;
```

### Spatial Coordinate Systems

Concatenation assumes coordinates are in the same reference frame:

```matlab
% If datasets have different coordinate systems, align first
% Example: Register SMD2 to SMD1 coordinate system

% Apply transformation
transform_matrix = load('registration_transform.mat');
SMD2.X = transform_matrix.tform.transformPointsForward([SMD2.X, SMD2.Y]);
SMD2.Y = SMD2_transformed(:, 2);
SMD2.X = SMD2_transformed(:, 1);

% Then concatenate
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2);
```

For multi-channel registration:

```matlab
% Load channel registration transform
CR = load(SMF.Data.RegistrationFilePath);

% Transform SMD2 (channel 2) to SMD1 (channel 1) frame
SMD2_transformed = smi_core.ChannelRegistration.transformSMD(SMD2, CR.tform);
SMD2_transformed.IsTransformed = true;

% Combine channels
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2_transformed);
```

### Preserving Metadata

Scalar fields (PixelSize, FrameRate, etc.) must be consistent:

```matlab
% Check consistency before concatenation
if SMD1.PixelSize ~= SMD2.PixelSize
    warning('Pixel sizes differ! SMD1: %.3f, SMD2: %.3f', ...
        SMD1.PixelSize, SMD2.PixelSize);

    % Option 1: Rescale coordinates
    scale_factor = SMD2.PixelSize / SMD1.PixelSize;
    SMD2.X = SMD2.X * scale_factor;
    SMD2.Y = SMD2.Y * scale_factor;
    SMD2.PixelSize = SMD1.PixelSize;

    % Option 2: Proceed anyway (catSMD will warn)
end

SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2);
```

If scalar fields differ, catSMD stores them as cell arrays:

```matlab
% If pixel sizes differ:
% SMD_combined.PixelSize = {0.108; 0.100}  (cell array, not scalar)

% Convert back to uniform if acceptable
if iscell(SMD_combined.PixelSize)
    SMD_combined.PixelSize = SMD_combined.PixelSize{1};
    warning('Using first PixelSize value for combined SMD');
end
```

## Complete Examples

### Example 1: Combining Biological Replicates

```matlab
% Combine three replicate experiments
replicate_files = {'Replicate1_Results.mat', ...
                   'Replicate2_Results.mat', ...
                   'Replicate3_Results.mat'};

SMD_all = smi_core.SingleMoleculeData.createSMD();

for ii = 1:length(replicate_files)
    data = load(replicate_files{ii}, 'SMD');
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, data.SMD, true);
    fprintf('Replicate %d: %d localizations, %d datasets\n', ...
        ii, length(data.SMD.X), data.SMD.NDatasets);
end

fprintf('\nCombined: %d localizations from %d datasets\n', ...
    length(SMD_all.X), SMD_all.NDatasets);

% Save combined results
save('All_Replicates_Results.mat', 'SMD_all', '-v7.3');
```

### Example 2: Batch Tracking with Concatenation

```matlab
% Track multiple files and concatenate results
obj = smi.SPT(SMF);
obj.FilePattern = 'Cell*_Data.mat';
obj.FindFiles = true;

% Batch track
[TR, SMD, SMDPreThresh, FileList, ~] = obj.batchTrack();

% TR and SMD are cell arrays, one per file
% obj.SMDBatch contains concatenated SMD across all files

fprintf('Tracked %d files\n', length(FileList));
fprintf('Combined SMD has %d localizations from %d datasets\n', ...
    length(obj.SMDBatch.X), obj.SMDBatch.NDatasets);

% Combine TR arrays
TR_all = smi_core.TrackingResults.createTR();
for ii = 1:length(TR)
    TR_all = smi_core.TrackingResults.catTR(TR_all, TR{ii});
end

fprintf('Total trajectories: %d\n', length(TR_all));
```

### Example 3: Drift Correction Workflow

```matlab
% Load multiple datasets
SMD_all = smi_core.SingleMoleculeData.createSMD();
file_list = dir('Dataset*.mat');

for ii = 1:length(file_list)
    data = load(file_list(ii).name, 'SMD');
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, data.SMD, true);
end

fprintf('Loaded %d datasets with %d localizations\n', ...
    SMD_all.NDatasets, length(SMD_all.X));

% Intra-dataset drift correction
SMF = smi_core.SingleMoleculeFitting();
DC = smi_core.DriftCorrection(SMF, SMD_all);
SMD_intra = smi_core.SingleMoleculeData.createSMD();

for ii = 1:SMD_all.NDatasets
    [SMD_intra_i, stats] = DC.driftCorrectKNNIntra(SMD_all, ii, ii);
    SMD_intra = smi_core.SingleMoleculeData.catSMD(SMD_intra, SMD_intra_i, false);
    fprintf('Dataset %d: RMSE = %.3f pixels\n', ii, stats.RMSE);
end

% Inter-dataset drift correction
[SMD_final, stats_inter] = DC.driftCorrectKNNInter(SMD_intra);

fprintf('\nInter-dataset drift correction:\n');
fprintf('  RMSE: %.3f pixels\n', stats_inter.RMSE);
fprintf('  Max drift: %.2f pixels\n', max(abs(SMD_final.DriftX(:))));

% Save drift-corrected results
save('Combined_DriftCorrected.mat', 'SMD_final', 'SMD_all', 'SMD_intra', '-v7.3');

% Visualize drift
DC.plotDriftCorrection(SMD_final, 'both');
```

## Common Pitfalls

### 1. Incorrect NewDataFlag

```matlab
% WRONG: Treating replicates as same data
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2, false);
% Dataset numbers overlap, drift arrays misaligned

% CORRECT: Treat as new data
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2, true);
```

### 2. Forgetting Empty Initialization

```matlab
% WRONG: First concatenation with undefined SMD
SMD_all = smi_core.SingleMoleculeData.catSMD(SMD1, SMD2);  % SMD1 must exist

% CORRECT: Initialize empty
SMD_all = smi_core.SingleMoleculeData.createSMD();
for ii = 1:N
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_list{ii});
end
```

### 3. Ignoring Coordinate Systems

```matlab
% WRONG: Combining without registration
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD_ch1, SMD_ch2);
% Channels misaligned

% CORRECT: Register first
SMD_ch2_reg = smi_core.ChannelRegistration.transformSMD(SMD_ch2, tform);
SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD_ch1, SMD_ch2_reg);
```

## See Also

- [SMD Structure](../core-concepts/smd-structure.md) - Detailed field descriptions
- [TR Structure](../core-concepts/tr-structure.md) - Tracking results format
- [Save and Load SMD](save-and-load-smd.md) - File operations
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete pipeline
