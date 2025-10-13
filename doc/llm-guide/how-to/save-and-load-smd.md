---
title: "How to Save and Load SMD/TR Results"
category: "how-to"
level: "beginner"
tags: ["saving", "loading", "persistence", "smd", "tr", "smf", "export", "csv", "hdf5"]
prerequisites: ["../core-concepts/smd-structure.md", "localize-molecules.md"]
related: ["visualize-results.md", "../workflows/smlm-analysis.md", "../workflows/spt-tracking.md"]
summary: "Complete guide to persisting and loading SMD/TR/SMF results including partial loading and export to other formats"
estimated_time: "12 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# How to Save and Load SMD/TR Results

## Purpose

After completing localization or tracking analysis, you need to save results for later use, sharing, or publication. This guide explains how to save and load SMD (Single Molecule Data), TR (Tracking Results), and SMF (parameter) structures using MATLAB's native .mat format, how to organize results directories, efficiently load partial data, and export to other formats like CSV and HDF5.

## Prerequisites

- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Completed localization or tracking analysis
- Basic MATLAB file I/O knowledge

## Overview

smite uses MATLAB's .mat format for saving results because it:
- Preserves all data types and structures
- Supports compression (-v7.3 format for large files)
- Allows partial loading of specific variables
- Is MATLAB's native format with fast read/write

Results are organized in a structured directory hierarchy:
```
DataDirectory/
  Results/
    DataFileName_Results.mat          (SMLM results)
    DataFileName_AnalysisID_Results.mat
    DataFileName/                     (plots and images)
    TrackingData_Results.mat          (SPT results)
```

## Automatic Saving from Workflows

### SMLM Analysis Results

The `smi.SMLM` workflow automatically saves results:

```matlab
% Configure SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'experiment1.h5'};
SMF.Data.ResultsDir = fullfile(SMF.Data.FileDir, 'Results');
SMF.Data.AnalysisID = 'test1';  % Optional identifier

% Run analysis (saves automatically)
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();

% Results saved to:
% /path/to/data/Results/experiment1_test1_Results.mat
```

**What gets saved:**
- `SMD`: Final filtered localization results
- `SMF`: All analysis parameters (packaged structure)

**File naming convention:**
- Without AnalysisID: `FileName_Results.mat`
- With AnalysisID: `FileName_AnalysisID_Results.mat`

### SPT Tracking Results

The `smi.SPT` workflow saves tracking results:

```matlab
% Configure and run tracking
SPTobj = smi.SPT();
SPTobj.SMF.Data.FileDir = '/path/to/data';
SPTobj.SMF.Data.FileName = {'tracking_data.mat'};
SPTobj.SMF.Data.ResultsDir = fullfile(SPTobj.SMF.Data.FileDir, 'Results');
SPTobj.performFullAnalysis();

% Results saved to:
% /path/to/data/Results/tracking_data_Results.mat
```

**What gets saved:**
- `SMD`: Individual localizations from tracking
- `SMDPreThresh`: Pre-threshold localizations
- `TR`: Tracking Results (trajectories organized by track)
- `SMF`: All analysis parameters

**Optional pre-registration results:**
If channel registration was performed, pre-registration results are also saved:
- `tracking_data_Results_PreCR.mat`: Contains `SMD`, `SMDPreThresh`, `TR`, `SMF` before registration

## Manual Saving

### Save SMD from LocalizeData

After localization without a full workflow:

```matlab
% Generate localizations
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();

% Save manually
results_dir = fullfile(SMF.Data.FileDir, 'Results');
if ~isfolder(results_dir)
    mkdir(results_dir);
end

% Package SMF for saving
SMF_packed = SMF.packageSMF();

% Save
save(fullfile(results_dir, 'my_results.mat'), 'SMD', 'SMF_packed', '-v7.3');
fprintf('Saved to: %s\n', fullfile(results_dir, 'my_results.mat'));
```

### Save Multiple Structures

Save complete analysis results:

```matlab
% After full analysis
SMD = % ... your final SMD
SMD_prethresh = % ... pre-threshold results
SMF_packed = SMF.packageSMF();

% Save all structures
save('complete_results.mat', 'SMD', 'SMD_prethresh', 'SMF_packed', '-v7.3');
```

### Save with Compression

For large datasets (>2GB or millions of localizations):

```matlab
% Use MATLAB's HDF5-based format with compression
save('large_results.mat', 'SMD', 'SMF_packed', '-v7.3', '-nocompression');
% or with compression (default):
save('large_results.mat', 'SMD', 'SMF_packed', '-v7.3');
```

**Format options:**
- `-v7.3`: Required for files >2GB, supports partial loading
- `-v7`: Legacy format, faster but 2GB limit
- `-nocompression`: Faster save/load but larger file size

## Loading Results

### Basic Loading

Load complete results file:

```matlab
% Load all variables
results = load('experiment1_Results.mat');
SMD = results.SMD;
SMF = results.SMF;

% Or load specific variables
load('experiment1_Results.mat', 'SMD');

% Check what's in file
vars = whos('-file', 'experiment1_Results.mat');
for i = 1:length(vars)
    fprintf('%s: %s [%s]\n', vars(i).name, mat2str(vars(i).size), vars(i).class);
end
```

### Load from Results Directory

Locate and load results from standard directory:

```matlab
% Construct path
data_dir = '/path/to/data';
file_name = 'experiment1';
analysis_id = 'test1';

if isempty(analysis_id)
    results_file = sprintf('%s_Results.mat', file_name);
else
    results_file = sprintf('%s_%s_Results.mat', file_name, analysis_id);
end

results_path = fullfile(data_dir, 'Results', results_file);

% Load
if isfile(results_path)
    load(results_path, 'SMD', 'SMF');
    fprintf('Loaded: %d localizations\n', length(SMD.X));
else
    error('Results file not found: %s', results_path);
end
```

### Load SPT Results

Load tracking results with all components:

```matlab
% Load SPT results
load('tracking_data_Results.mat', 'SMD', 'TR', 'SMF');

fprintf('Loaded %d localizations\n', length(SMD.X));
fprintf('Loaded %d trajectories\n', length(TR));

% Access trajectory information
for i = 1:min(5, length(TR))
    fprintf('Track %d: %d points, %d frames\n', ...
        i, length(TR(i).X), length(TR(i).FrameNum));
end
```

## Partial Loading (Memory Efficient)

For very large results files, load only what you need:

### Load Specific Fields

```matlab
% Load only positions (not photons, backgrounds, etc.)
m = matfile('large_results.mat');  % Create handle without loading

% Access specific fields
X = m.SMD.X;  % Only loads X positions
Y = m.SMD.Y;

fprintf('Loaded %d positions without loading full SMD\n', length(X));

% Check file contents without loading
info = whos(m);
for i = 1:length(info)
    fprintf('%s: %.1f MB\n', info(i).name, info(i).bytes/1e6);
end
```

### Load Subset of Localizations

```matlab
% Load results handle
m = matfile('large_results.mat');

% Load first 10000 localizations only
subset_size = 10000;
SMD_subset = struct();
SMD_subset.X = m.SMD.X(1:subset_size);
SMD_subset.Y = m.SMD.Y(1:subset_size);
SMD_subset.Photons = m.SMD.Photons(1:subset_size);
SMD_subset.FrameNum = m.SMD.FrameNum(1:subset_size);

% Copy scalar fields
SMD_subset.PixelSize = m.SMD.PixelSize;
SMD_subset.NFrames = m.SMD.NFrames;

fprintf('Loaded subset: %d localizations\n', length(SMD_subset.X));
```

### Load by Frame Range

```matlab
% Load handle
m = matfile('large_results.mat');

% Determine frames to load
frame_start = 100;
frame_end = 500;

% Find indices in this frame range
frame_nums = m.SMD.FrameNum(:);  % Load frame numbers
indices = find(frame_nums >= frame_start & frame_nums <= frame_end);

fprintf('Found %d localizations in frames %d-%d\n', ...
    length(indices), frame_start, frame_end);

% Load only those localizations
SMD_windowed = struct();
SMD_windowed.X = m.SMD.X(indices);
SMD_windowed.Y = m.SMD.Y(indices);
SMD_windowed.Photons = m.SMD.Photons(indices);
SMD_windowed.FrameNum = m.SMD.FrameNum(indices);
```

### Load by Spatial ROI

```matlab
% Load positions first
m = matfile('large_results.mat');
X = m.SMD.X(:);
Y = m.SMD.Y(:);

% Define ROI [X_min, X_max, Y_min, Y_max] in pixels
roi = [50, 150, 50, 150];

% Find localizations in ROI
in_roi = (X >= roi(1)) & (X <= roi(2)) & (Y >= roi(3)) & (Y <= roi(4));
indices = find(in_roi);

fprintf('Found %d localizations in ROI\n', length(indices));

% Load only those localizations
SMD_roi = struct();
SMD_roi.X = X(indices);
SMD_roi.Y = Y(indices);
SMD_roi.Photons = m.SMD.Photons(indices);
SMD_roi.FrameNum = m.SMD.FrameNum(indices);
```

## Results Directory Organization

### Standard Directory Structure

```
MyExperiment/
  Cell1_Data.h5                    # Raw data
  Cell2_Data.h5
  Results/                         # Created automatically
    Cell1_Data_Results.mat         # Localization results
    Cell1_Data/                    # Plots and images
      SR_Image.png
      Histogram_Photons.fig
      DriftCorrection.fig
    Cell1_Data_test1_Results.mat   # With AnalysisID
    Cell1_Data_test1/
      SR_Image.png
      ...
    Cell2_Data_Results.mat
    Cell2_Data/
      ...
```

### Configuring Results Directory

```matlab
% Default: Results subdirectory of data location
SMF.Data.FileDir = '/path/to/data';
SMF.Data.ResultsDir = fullfile(SMF.Data.FileDir, 'Results');

% Custom location
SMF.Data.ResultsDir = '/path/to/custom/results';

% Ensure directory exists
if ~isfolder(SMF.Data.ResultsDir)
    mkdir(SMF.Data.ResultsDir);
end
```

### Working with Analysis IDs

Use Analysis IDs to save multiple analyses of the same data:

```matlab
% First analysis
SMF.Data.AnalysisID = 'conservative';
SMF.Thresholding.MaxXY_SE = 0.12;
% Saves to: DataFile_conservative_Results.mat

% Second analysis
SMF.Data.AnalysisID = 'permissive';
SMF.Thresholding.MaxXY_SE = 0.25;
% Saves to: DataFile_permissive_Results.mat

% Load specific analysis
load('DataFile_conservative_Results.mat', 'SMD');
SMD_conservative = SMD;
load('DataFile_permissive_Results.mat', 'SMD');
SMD_permissive = SMD;

fprintf('Conservative: %d localizations\n', length(SMD_conservative.X));
fprintf('Permissive: %d localizations\n', length(SMD_permissive.X));
```

## Exporting to Other Formats

### Export to CSV

Export localization table for external analysis:

```matlab
% Load results
load('Results.mat', 'SMD');

% Create table with key fields
T = table(SMD.FrameNum, SMD.X, SMD.Y, SMD.Photons, ...
    SMD.X_SE, SMD.Y_SE, SMD.Bg, ...
    'VariableNames', {'Frame', 'X_pixels', 'Y_pixels', 'Photons', ...
    'X_SE_pixels', 'Y_SE_pixels', 'Background'});

% Convert to physical units
T.X_nm = SMD.X * SMD.PixelSize * 1000;
T.Y_nm = SMD.Y * SMD.PixelSize * 1000;
T.X_SE_nm = SMD.X_SE * SMD.PixelSize * 1000;
T.Y_SE_nm = SMD.Y_SE * SMD.PixelSize * 1000;

% Save to CSV
writetable(T, 'localizations.csv');
fprintf('Exported %d localizations to CSV\n', height(T));
```

### Export Trajectories to CSV

Export SPT trajectories:

```matlab
% Load tracking results
load('tracking_Results.mat', 'TR');

% Convert TR to table (one row per localization)
all_data = [];
for i = 1:length(TR)
    track_id = i * ones(length(TR(i).X), 1);
    track_data = [track_id, TR(i).FrameNum, TR(i).X, TR(i).Y, ...
                  TR(i).Photons, TR(i).X_SE, TR(i).Y_SE];
    all_data = [all_data; track_data];
end

% Create table
T = array2table(all_data, 'VariableNames', ...
    {'TrackID', 'Frame', 'X', 'Y', 'Photons', 'X_SE', 'Y_SE'});

writetable(T, 'trajectories.csv');
fprintf('Exported %d trajectories (%d points)\n', length(TR), height(T));
```

### Export Specific Columns

Export minimal data for external tools:

```matlab
% Minimal format for ThunderSTORM, etc.
T_minimal = table(SMD.X, SMD.Y, SMD.Photons, SMD.FrameNum, ...
    'VariableNames', {'x', 'y', 'intensity', 'frame'});

writetable(T_minimal, 'localizations_minimal.csv');
```

### Export to HDF5

Save results in HDF5 format for Python/other tools:

```matlab
% Create HDF5 file
h5_file = 'results.h5';
if isfile(h5_file)
    delete(h5_file);
end

% Write datasets
h5create(h5_file, '/X', size(SMD.X));
h5write(h5_file, '/X', SMD.X);

h5create(h5_file, '/Y', size(SMD.Y));
h5write(h5_file, '/Y', SMD.Y);

h5create(h5_file, '/Photons', size(SMD.Photons));
h5write(h5_file, '/Photons', SMD.Photons);

h5create(h5_file, '/FrameNum', size(SMD.FrameNum));
h5write(h5_file, '/FrameNum', double(SMD.FrameNum));

% Write attributes (metadata)
h5writeatt(h5_file, '/', 'PixelSize', SMD.PixelSize);
h5writeatt(h5_file, '/', 'NFrames', SMD.NFrames);
h5writeatt(h5_file, '/', 'NLocalizations', length(SMD.X));

fprintf('Exported to HDF5: %s\n', h5_file);

% Verify
h5disp(h5_file);
```

### Export to JSON (for web applications)

```matlab
% Create structure with key data
export_struct = struct();
export_struct.positions = [SMD.X, SMD.Y];
export_struct.photons = SMD.Photons;
export_struct.frames = SMD.FrameNum;
export_struct.metadata.pixel_size = SMD.PixelSize;
export_struct.metadata.n_frames = SMD.NFrames;

% Convert to JSON
json_text = jsonencode(export_struct);

% Save
fid = fopen('results.json', 'w');
fprintf(fid, '%s', json_text);
fclose(fid);

fprintf('Exported to JSON\n');
```

## Complete Workflow Examples

### Example 1: Save and Reload SMLM Results

```matlab
% ===== Initial Analysis =====
% Load and analyze data
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'cell1.h5'};
SMF.Data.ResultsDir = fullfile(SMF.Data.FileDir, 'Results');

SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();

fprintf('Analysis complete. Results saved automatically.\n');

% ===== Later: Reload Results =====
% Clear workspace
clear all;

% Load saved results
data_dir = '/path/to/data';
load(fullfile(data_dir, 'Results', 'cell1_Results.mat'), 'SMD', 'SMF');

fprintf('Reloaded: %d localizations\n', length(SMD.X));

% Continue analysis
% Apply additional filtering, visualization, etc.
```

### Example 2: Batch Analysis with Organized Saving

```matlab
% Process multiple files
file_dir = '/path/to/data';
files = {'cell1.h5', 'cell2.h5', 'cell3.h5'};

for i = 1:length(files)
    fprintf('Processing %s...\n', files{i});

    % Configure
    SMF = smi_core.SingleMoleculeFitting();
    SMF.Data.FileDir = file_dir;
    SMF.Data.FileName = files(i);
    SMF.Data.ResultsDir = fullfile(file_dir, 'Results');
    SMF.Data.AnalysisID = sprintf('batch_%d', i);

    % Analyze
    SMLMobj = smi.SMLM(SMF);
    SMLMobj.Verbose = 1;
    SMLMobj.fullAnalysis();

    fprintf('Saved: %s\n', fullfile(SMF.Data.ResultsDir, ...
        sprintf('%s_batch_%d_Results.mat', files{i}(1:end-3), i)));
end

fprintf('\nBatch processing complete.\n');
```

### Example 3: Load, Filter, and Re-save

```matlab
% Load original results
load('experiment_Results.mat', 'SMD', 'SMF');
fprintf('Original: %d localizations\n', length(SMD.X));

% Apply stricter post-hoc filtering
SMF.Thresholding.MaxXY_SE = 0.12;  % Stricter
MinMax = smi_core.Threshold.setMinMax(SMF);
THR = smi_core.Threshold();
[SMD, ~] = THR.setThreshFlag(SMD, MinMax);
SMD_filtered = smi_core.Threshold.applyThresh(SMD, 1);

fprintf('Filtered: %d localizations\n', length(SMD_filtered.X));

% Save filtered version
SMF_packed = SMF.packageSMF();
save('experiment_filtered_Results.mat', 'SMD_filtered', 'SMF_packed', '-v7.3');

% Also export to CSV for external analysis
T = table(SMD_filtered.X, SMD_filtered.Y, SMD_filtered.Photons, ...
    'VariableNames', {'X', 'Y', 'Photons'});
writetable(T, 'experiment_filtered.csv');

fprintf('Saved filtered results in .mat and .csv formats\n');
```

### Example 4: Memory-Efficient Loading of Large Dataset

```matlab
% Work with large file without loading everything
results_file = 'large_experiment_Results.mat';

% Create file handle
m = matfile(results_file);

% Get size information
fprintf('File contains %d localizations\n', length(m.SMD.X));

% Process in chunks
chunk_size = 100000;
n_chunks = ceil(length(m.SMD.X) / chunk_size);

all_densities = [];
for chunk = 1:n_chunks
    % Load this chunk
    idx_start = (chunk-1)*chunk_size + 1;
    idx_end = min(chunk*chunk_size, length(m.SMD.X));

    X_chunk = m.SMD.X(idx_start:idx_end);
    Y_chunk = m.SMD.Y(idx_start:idx_end);

    % Process chunk (example: compute local density)
    density_chunk = compute_density(X_chunk, Y_chunk);
    all_densities = [all_densities; density_chunk];

    fprintf('Processed chunk %d/%d\n', chunk, n_chunks);
end

% Save processed results
save('density_results.mat', 'all_densities');
fprintf('Processing complete\n');

function density = compute_density(X, Y)
    % Placeholder for density calculation
    density = ones(size(X));
end
```

## Troubleshooting

### Issue: File too large (>2GB error)

```matlab
% Use -v7.3 format
save('results.mat', 'SMD', 'SMF', '-v7.3');

% Or reduce data size by filtering
SMD_subset = smi_core.SingleMoleculeData.isolateSubSMD(SMD, 1:100000);
save('results_subset.mat', 'SMD_subset', 'SMF');
```

### Issue: Loading takes too long

```matlab
% Use partial loading
m = matfile('large_results.mat');

% Load only what you need
X = m.SMD.X;
Y = m.SMD.Y;

% Instead of:
% load('large_results.mat');  % Loads everything
```

### Issue: Variable not found in file

```matlab
% Check file contents
vars = whos('-file', 'results.mat');
disp(vars);

% Load with error handling
if ismember('SMD', {vars.name})
    load('results.mat', 'SMD');
else
    error('SMD not found in file');
end
```

### Issue: Out of memory when loading

```matlab
% Use matfile for memory-mapped access
m = matfile('results.mat');

% Access directly without loading to workspace
n_locs = length(m.SMD.X);
mean_photons = mean(m.SMD.Photons);

fprintf('File has %d localizations\n', n_locs);
fprintf('Mean photons: %.1f\n', mean_photons);
```

### Issue: Corrupted or incomplete save

```matlab
% Verify save succeeded
save('results.mat', 'SMD', 'SMF', '-v7.3');

% Verify by loading
try
    test = load('results.mat', 'SMD');
    fprintf('Save verified: %d localizations\n', length(test.SMD.X));
catch ME
    error('Save failed: %s', ME.message);
end
```

## Best Practices

1. **Always use -v7.3 format** for results files (supports large data and partial loading)
2. **Use descriptive AnalysisID** to distinguish different analyses
3. **Save SMF with results** so you know analysis parameters later
4. **Organize with Results directory** to keep data separate from results
5. **Use matfile for large datasets** instead of loading everything
6. **Export to CSV for sharing** with collaborators using other tools
7. **Verify saves succeeded** before deleting intermediate data
8. **Document custom exports** if using non-standard formats

## See Also

- [SMD Structure](../core-concepts/smd-structure.md) - Understanding result fields
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete analysis with saving
- [SPT Workflow](../workflows/spt-tracking.md) - Tracking analysis with saving
- [Visualize Results](visualize-results.md) - Creating figures from loaded data
- [Threshold Results](threshold-results.md) - Post-hoc filtering of loaded data
