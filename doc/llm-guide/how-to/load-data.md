---
title: "How to Load Data"
category: "how-to"
level: "beginner"
tags: ["data-loading", "files", "h5", "mat", "roi"]
prerequisites: ["../core-concepts/smf-structure.md"]
related: ["localize-molecules.md", "../workflows/smlm-analysis.md"]
summary: "Practical guide to loading raw camera data from .h5 and .mat files into smite"
estimated_time: "10 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# How to Load Data

## Purpose

Loading data correctly is the first critical step in any smite analysis. This guide shows you how to load raw camera images from .h5 and .mat files, configure camera parameters, apply ROIs, select specific datasets, and verify that data loaded correctly.

## Prerequisites

- Basic understanding of [SMF structure](../core-concepts/smf-structure.md)
- Raw data file (.h5 or .mat format)

## Overview

smite supports two primary data formats:
- **.h5 (HDF5)**: Preferred format, supports multiple datasets, metadata
- **.mat (MATLAB)**: Simple format, requires variable name specification

The `smi_core.LoadData` class handles all loading operations, controlled by `SMF.Data` parameters. Data is returned as a 3D array: `(Y_pixels × X_pixels × Frames)`.

## Basic Data Loading

### From .h5 File

```matlab
% Configure SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'experiment1.h5'};
SMF.Data.FileType = 'h5';  % Auto-detected from extension

% Load first dataset
LD = smi_core.LoadData();
[~, sequence, SMF] = LD.loadRawData(SMF, 1);

% sequence is Y × X × Frames array
fprintf('Loaded: %d × %d × %d\n', size(sequence));
```

**.h5 files** can contain multiple datasets. Specify which to load:

```matlab
% Load dataset 3
[~, sequence, SMF] = LD.loadRawData(SMF, 3);

% Load all datasets
SMF.Data.DatasetList = int32([]);  % Empty = all
```

### From .mat File

```matlab
% Configure SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'experiment1.mat'};
SMF.Data.FileType = 'mat';
SMF.Data.DataVariable = 'sequence';  % Variable name in .mat file

% Load
LD = smi_core.LoadData();
[~, sequence, SMF] = LD.loadRawData(SMF, 1);
```

**Common variable names:**
- `sequence` (default)
- `data`
- `imageStack`

If your .mat file uses a different name, set `DataVariable`:

```matlab
SMF.Data.DataVariable = 'myData';
```

### Direct Loading (Alternative)

For simple cases, load directly:

```matlab
% For .mat files
data_struct = load('experiment1.mat');
sequence = data_struct.sequence;  % Or your variable name

% For .h5 files
sequence = h5read('experiment1.h5', '/Data/dataset_1');
```

But using `LoadData` is recommended because it:
- Handles file format differences automatically
- Updates SMF with metadata
- Applies ROI and camera corrections
- Validates data integrity

## Configuring Camera Parameters

Camera parameters convert raw ADU (analog-to-digital units) to photons.

### EMCCD Camera

Typically uniform gain and offset:

```matlab
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;    % ADU per photon
SMF.Data.CameraOffset = 100;  % Baseline ADU
SMF.Data.CameraNoise = 5;     % Read noise (electrons RMS)

% Photon conversion happens automatically during analysis:
% photons = (sequence - CameraOffset) / CameraGain
```

**Finding your camera parameters:**

1. **Gain**: From camera specs or calibration (typical: 1-10 ADU/photon for EMCCD)
2. **Offset**: Mean value of dark frame (lights off, shutter closed)
3. **Noise**: Standard deviation of dark frame

```matlab
% If you have dark frames
dark_frames = % ... load your dark frames ...
SMF.Data.CameraOffset = mean(dark_frames(:));
SMF.Data.CameraNoise = std(dark_frames(:));
```

### sCMOS Camera

Pixel-wise varying gain, offset, and noise:

```matlab
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = '/path/to/calibration.mat';

% Calibration file must contain:
%   CameraGain    (Y × X array, ADU/photon)
%   CameraOffset  (Y × X array, ADU)
%   CameraNoise   (Y × X array, electrons RMS)
```

**Creating calibration file:**

```matlab
% From calibration measurements
CameraGain = % ... Y × X array ...
CameraOffset = % ... Y × X array ...
CameraNoise = % ... Y × X array ...

save('camera_calibration.mat', 'CameraGain', 'CameraOffset', 'CameraNoise');
```

## Working with Multiple Datasets

.h5 files can contain multiple datasets (e.g., different imaging conditions, time points).

### List Available Datasets

```matlab
% Check what datasets are in the file
file_path = fullfile(SMF.Data.FileDir, SMF.Data.FileName{1});
info = h5info(file_path);

% Display datasets
for i = 1:length(info.Groups(1).Groups)
    dataset_name = info.Groups(1).Groups(i).Name;
    fprintf('Dataset %d: %s\n', i, dataset_name);
end
```

### Select Specific Datasets

```matlab
% Analyze only datasets 1, 3, and 5
SMF.Data.DatasetList = int32([1, 3, 5]);

% Or use inclusion/exclusion
SMF.Data.DatasetMods = {1:10; [2,4,7]};  % Include 1-10, exclude 2,4,7
% Results in: 1,3,5,6,8,9,10
```

**How DatasetMods works:**

```matlab
% {inclusion_list; exclusion_list}
SMF.Data.DatasetMods = {[]; []};        % All datasets
SMF.Data.DatasetMods = {1:5; []};       % Datasets 1-5
SMF.Data.DatasetMods = {[]; [3,4]};     % All except 3 and 4
SMF.Data.DatasetMods = {1:10; [5,6]};   % Datasets 1-10 except 5,6
```

`DatasetList` is computed from `DatasetMods` automatically by `LoadData`.

### Loop Over Datasets

```matlab
% Load and process each dataset
for dataset_num = 1:5
    [~, sequence, ~] = LD.loadRawData(SMF, dataset_num);

    % Process this dataset
    fprintf('Dataset %d: %d frames\n', dataset_num, size(sequence, 3));

    % Your analysis here...
end
```

## Applying Region of Interest (ROI)

Analyze only a portion of the image:

```matlab
% Define ROI: [YStart, XStart, YEnd, XEnd] (1-indexed)
SMF.Data.DataROI = [50, 50, 200, 200];

% Load data - ROI applied automatically
[~, sequence, ~] = LD.loadRawData(SMF, 1);

% sequence is now (200-50+1) × (200-50+1) × Frames
fprintf('ROI size: %d × %d\n', size(sequence, 1), size(sequence, 2));
```

**Interactive ROI selection:**

```matlab
% Load full image first
SMF.Data.DataROI = [];
[~, sequence_full, ~] = LD.loadRawData(SMF, 1);

% Display first frame
imagesc(sequence_full(:,:,1));
axis image; colorbar;
title('Click to select ROI');

% Use MATLAB's ROI tools
h = drawrectangle;  % Or imrect, etc.
pos = round(h.Position);  % [X, Y, Width, Height]

% Convert to smite format [YStart, XStart, YEnd, XEnd]
SMF.Data.DataROI = [pos(2), pos(1), pos(2)+pos(4)-1, pos(1)+pos(3)-1];

fprintf('ROI: [%d, %d, %d, %d]\n', SMF.Data.DataROI);

% Reload with ROI
[~, sequence_roi, ~] = LD.loadRawData(SMF, 1);
```

**Why use ROI:**
- Faster processing (smaller data)
- Focus on area of interest
- Reduce memory usage
- Exclude bad regions (dust, dead pixels)

## Configuring Acquisition Parameters

Essential metadata for analysis:

```matlab
% Pixel size (micrometers)
SMF.Data.PixelSize = 0.108;  % 108 nm for typical 100× objective

% Frame rate (Hz)
SMF.Data.FrameRate = 100;  % 100 frames per second

% These are used for:
% - Converting positions to physical units
% - Drift correction
% - Diffusion analysis
% - Proper MSD calculations
```

**Finding pixel size:**

```
pixel_size_μm = magnification × camera_pixel_μm / (objective_mag × tube_lens_mag)
```

For example:
- Camera: 6.5 μm pixels
- Objective: 100×
- Additional magnification: 1.6×

```
pixel_size = 6.5 / (100 × 1.6) = 0.0406 μm = 40.6 nm
```

## Verifying Loaded Data

Always verify data loaded correctly:

### Check Dimensions

```matlab
[~, sequence, SMF] = LD.loadRawData(SMF, 1);

fprintf('Data dimensions:\n');
fprintf('  Y size: %d pixels\n', size(sequence, 1));
fprintf('  X size: %d pixels\n', size(sequence, 2));
fprintf('  Frames: %d\n', size(sequence, 3));

% Should match expected values
assert(size(sequence, 3) > 0, 'No frames loaded!');
```

### Visualize Data

```matlab
% View first frame
figure;
imagesc(sequence(:,:,1));
axis image; colorbar;
title('First Frame');
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Check intensity range
fprintf('Intensity range: %d to %d ADU\n', min(sequence(:)), max(sequence(:)));

% View middle frame
mid_frame = round(size(sequence, 3) / 2);
figure;
imagesc(sequence(:,:,mid_frame));
axis image; colorbar;
title(sprintf('Frame %d', mid_frame));
```

### Check for Issues

```matlab
% Check for all-zero frames
zero_frames = squeeze(sum(sum(sequence, 1), 2)) == 0;
if any(zero_frames)
    warning('Found %d all-zero frames', sum(zero_frames));
    fprintf('Zero frames: %s\n', mat2str(find(zero_frames)));
end

% Check for hot pixels
max_per_frame = squeeze(max(max(sequence, [], 1), [], 2));
if any(max_per_frame > 10000)  % Adjust threshold
    warning('Possible hot pixels detected');
end

% Check for saturation
saturated_pixels = sum(sequence(:) >= 4095);  % 12-bit camera
if saturated_pixels > 0
    warning('%d saturated pixels (%.2f%%)', ...
        saturated_pixels, 100*saturated_pixels/numel(sequence));
end
```

## Common Data Loading Patterns

### Pattern 1: Simple Single File

```matlab
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = pwd;  % Current directory
SMF.Data.FileName = {'data.mat'};
SMF.Data.DataVariable = 'sequence';
SMF.Data.PixelSize = 0.1;
SMF.Data.FrameRate = 100;

LD = smi_core.LoadData();
[~, sequence, SMF] = LD.loadRawData(SMF, 1);
```

### Pattern 2: Multiple Files, Same Directory

```matlab
SMF.Data.FileDir = '/data/2024-01-10';
SMF.Data.FileName = {'Cell1.h5', 'Cell2.h5', 'Cell3.h5'};

% Process each file
for file_idx = 1:length(SMF.Data.FileName)
    SMF.Data.FileName = SMF.Data.FileName(file_idx);
    [~, sequence, ~] = LD.loadRawData(SMF, 1);

    % Analyze...
end
```

### Pattern 3: Batch with ROI

```matlab
% Define ROI once
SMF.Data.DataROI = [50, 50, 200, 200];

% Load multiple datasets
for ds = 1:5
    [~, sequence, ~] = LD.loadRawData(SMF, ds);
    fprintf('Dataset %d: %d × %d × %d\n', ds, size(sequence));
end
```

### Pattern 4: Memory-Efficient (Load Frame by Frame)

For very large files:

```matlab
% Get dimensions first
[~, sequence_small, ~] = LD.loadRawData(SMF, 1);
[nY, nX, nFrames] = size(sequence_small);
clear sequence_small;

% Process frame by frame
for frame = 1:nFrames
    % Load single frame
    % (Requires custom implementation or subsampling)

    % For now, smite loads full dataset
    % Consider using DataROI to reduce size
end
```

## Troubleshooting

### Issue: File not found

```matlab
% Check file exists
file_path = fullfile(SMF.Data.FileDir, SMF.Data.FileName{1});
if ~exist(file_path, 'file')
    error('File not found: %s', file_path);
end

% Check directory
if ~isfolder(SMF.Data.FileDir)
    error('Directory not found: %s', SMF.Data.FileDir);
end
```

### Issue: Wrong variable name in .mat file

```matlab
% List variables in .mat file
vars = whos('-file', 'data.mat');
fprintf('Variables in file:\n');
for i = 1:length(vars)
    fprintf('  %s: %s\n', vars(i).name, vars(i).class);
end

% Set correct variable
SMF.Data.DataVariable = 'correct_name';
```

### Issue: Data dimensions wrong

```matlab
% After loading
[nY, nX, nFrames] = size(sequence);

% Check if dimensions make sense
if nFrames < nY  % Might be transposed
    warning('Data might be incorrectly oriented');
    fprintf('Size: %d × %d × %d\n', nY, nX, nFrames);
end

% Verify against expected
expected_Y = 256;
expected_X = 256;
if nY ~= expected_Y || nX ~= expected_X
    warning('Unexpected image size');
end
```

### Issue: Data is in photons, not ADU

If data already converted to photons:

```matlab
% Set gain = 1, offset = 0
SMF.Data.CameraGain = 1;
SMF.Data.CameraOffset = 0;

% Or skip conversion in analysis
```

### Issue: ROI out of bounds

```matlab
% Validate ROI
[~, sequence_full, ~] = LD.loadRawData(SMF, 1);
[nY, nX, ~] = size(sequence_full);

roi = SMF.Data.DataROI;
if any(roi < 1) || roi(3) > nY || roi(4) > nX
    error('ROI [%d,%d,%d,%d] exceeds image bounds (%d×%d)', ...
        roi, nY, nX);
end
```

## See Also

- [SMF Structure](../core-concepts/smf-structure.md) - Data parameters reference
- [How to Localize Molecules](localize-molecules.md) - Next step after loading
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete analysis pipeline
- doc/FileFormats/HDF5.md - .h5 file format details
- doc/FileFormats/CalibrationFile.md - Camera calibration format
