---
title: "How to Export Results to Other Formats and Tools"
category: "how-to"
level: "beginner"
tags: ["export", "csv", "hdf5", "json", "interoperability", "thunderstorm", "visp", "python"]
prerequisites: ["save-and-load-smd.md", "../core-concepts/smd-structure.md"]
related: ["visualize-results.md", "localize-molecules.md", "../workflows/smlm-analysis.md"]
summary: "Complete guide to exporting SMD/TR data to CSV, HDF5, JSON, and other tool-specific formats for analysis in external applications"
estimated_time: "15 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# How to Export Results to Other Formats and Tools

## Purpose

While smite uses MATLAB's .mat format for native storage, you often need to export results to other formats for collaboration, publication, or analysis in external tools. This guide covers exporting SMD (Single Molecule Data) and TR (Tracking Results) to widely-used formats including CSV for spreadsheets, HDF5 for Python/other languages, JSON for web applications, and tool-specific formats for programs like ViSP, ThunderSTORM, and rapidSTORM.

## Prerequisites

- SMD or TR structure with localization/tracking data
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Knowledge of [saving and loading](save-and-load-smd.md) basics
- Target application or format requirements

## Overview

Export formats serve different purposes:

- **CSV**: Human-readable, spreadsheet-compatible, universally supported
- **HDF5**: Efficient binary format for large datasets, Python-friendly
- **JSON**: Web applications, JavaScript, readable structured data
- **Tool-specific**: Pre-configured formats for SMLM analysis software

This guide shows how to export coordinates, metadata, uncertainties, and other fields while maintaining data integrity and proper units.

## Export to CSV

CSV (Comma-Separated Values) is the most universal format for data exchange.

### Basic Localization Export

Export essential localization data for spreadsheet analysis:

```matlab
% Load results
load('Results.mat', 'SMD');

% Create table with core fields
T = table(SMD.FrameNum, SMD.X, SMD.Y, SMD.Photons, ...
    SMD.X_SE, SMD.Y_SE, SMD.Bg, ...
    'VariableNames', {'Frame', 'X_pixels', 'Y_pixels', 'Photons', ...
    'X_SE_pixels', 'Y_SE_pixels', 'Background'});

% Write to CSV
writetable(T, 'localizations.csv');
fprintf('Exported %d localizations to CSV\n', height(T));
```

**Output file structure:**
```
Frame,X_pixels,Y_pixels,Photons,X_SE_pixels,Y_SE_pixels,Background
1,45.32,67.89,1234,0.08,0.09,123
1,123.45,89.12,2345,0.07,0.08,145
2,46.01,68.23,1189,0.09,0.10,118
...
```

### Export with Physical Units

Include both pixel and nanometer coordinates:

```matlab
% Load results
load('Results.mat', 'SMD');

% Convert to physical units (nm)
pixel_size_nm = SMD.PixelSize * 1000;  % Convert microns to nm

% Create comprehensive table
T = table();
T.Frame = SMD.FrameNum;
T.X_pixels = SMD.X;
T.Y_pixels = SMD.Y;
T.X_nm = SMD.X * pixel_size_nm;
T.Y_nm = SMD.Y * pixel_size_nm;
T.Photons = SMD.Photons;
T.Background = SMD.Bg;
T.X_Precision_pixels = SMD.X_SE;
T.Y_Precision_pixels = SMD.Y_SE;
T.X_Precision_nm = SMD.X_SE * pixel_size_nm;
T.Y_Precision_nm = SMD.Y_SE * pixel_size_nm;

% Add metadata as header comment
fid = fopen('localizations_with_units.csv', 'w');
fprintf(fid, '# PixelSize: %.4f microns\n', SMD.PixelSize);
fprintf(fid, '# NFrames: %d\n', SMD.NFrames);
fprintf(fid, '# NLocalizations: %d\n', length(SMD.X));
fprintf(fid, '# Camera: %s\n', SMD.CameraType);
fclose(fid);

% Append data
writetable(T, 'localizations_with_units.csv', 'WriteMode', 'Append');
fprintf('Exported with metadata and dual units\n');
```

### Export Tracking Trajectories to CSV

Export SPT results with trajectory IDs:

```matlab
% Load tracking results
load('tracking_Results.mat', 'TR');

% Convert TR array to table (one row per localization)
data_cell = cell(length(TR), 1);
for i = 1:length(TR)
    n_points = length(TR(i).X);
    track_data = table();
    track_data.TrackID = repmat(i, n_points, 1);
    track_data.Frame = TR(i).FrameNum;
    track_data.X = TR(i).X;
    track_data.Y = TR(i).Y;
    track_data.Photons = TR(i).Photons;
    track_data.X_SE = TR(i).X_SE;
    track_data.Y_SE = TR(i).Y_SE;
    track_data.ConnectID = TR(i).ConnectID;
    data_cell{i} = track_data;
end

% Concatenate all tracks
T = vertcat(data_cell{:});

% Write to CSV
writetable(T, 'trajectories.csv');
fprintf('Exported %d trajectories (%d points total)\n', length(TR), height(T));
```

### Export Trajectory Statistics

Export per-trajectory summary statistics:

```matlab
% Load tracking results
load('tracking_Results.mat', 'TR');

% Compute statistics for each trajectory
n_tracks = length(TR);
T = table();
T.TrackID = (1:n_tracks)';
T.NPoints = arrayfun(@(x) length(x.X), TR)';
T.Duration_frames = arrayfun(@(x) max(x.FrameNum) - min(x.FrameNum) + 1, TR)';
T.MeanPhotons = arrayfun(@(x) mean(x.Photons), TR)';
T.TotalDisplacement = arrayfun(@(x) sqrt((x.X(end)-x.X(1))^2 + (x.Y(end)-x.Y(1))^2), TR)';

% Compute MSD if available
if isfield(TR, 'MSD')
    T.MSD_slope = arrayfun(@(x) x.MSD(1), TR)';
end

writetable(T, 'trajectory_statistics.csv');
fprintf('Exported statistics for %d trajectories\n', n_tracks);
```

### Selective Field Export

Export only specific columns for external tools:

```matlab
% Minimal format for general SMLM tools
T_minimal = table(SMD.X, SMD.Y, SMD.Photons, SMD.FrameNum, ...
    'VariableNames', {'x', 'y', 'intensity', 'frame'});
writetable(T_minimal, 'localizations_minimal.csv');

% Precision-weighted format
T_weighted = table(SMD.X, SMD.Y, 1./(SMD.X_SE.^2), SMD.FrameNum, ...
    'VariableNames', {'x', 'y', 'weight', 'frame'});
writetable(T_weighted, 'localizations_weighted.csv');

% Channel-specific export (if dual-channel)
if isfield(SMD, 'Channel') && max(SMD.Channel) > 1
    for ch = 1:max(SMD.Channel)
        idx = SMD.Channel == ch;
        T_ch = table(SMD.X(idx), SMD.Y(idx), SMD.Photons(idx), ...
            SMD.FrameNum(idx), 'VariableNames', {'x', 'y', 'intensity', 'frame'});
        writetable(T_ch, sprintf('localizations_channel%d.csv', ch));
    end
end
```

## Export to HDF5

HDF5 is an efficient binary format ideal for large datasets and Python integration.

### Basic HDF5 Export

Export core localization data:

```matlab
% Load results
load('Results.mat', 'SMD');

% Create HDF5 file
h5_file = 'results.h5';
if isfile(h5_file)
    delete(h5_file);  % Remove existing file
end

% Write localization arrays
h5create(h5_file, '/X', size(SMD.X));
h5write(h5_file, '/X', SMD.X);

h5create(h5_file, '/Y', size(SMD.Y));
h5write(h5_file, '/Y', SMD.Y);

h5create(h5_file, '/Photons', size(SMD.Photons));
h5write(h5_file, '/Photons', SMD.Photons);

h5create(h5_file, '/FrameNum', size(SMD.FrameNum));
h5write(h5_file, '/FrameNum', double(SMD.FrameNum));

h5create(h5_file, '/X_SE', size(SMD.X_SE));
h5write(h5_file, '/X_SE', SMD.X_SE);

h5create(h5_file, '/Y_SE', size(SMD.Y_SE));
h5write(h5_file, '/Y_SE', SMD.Y_SE);

% Write metadata as attributes
h5writeatt(h5_file, '/', 'PixelSize', SMD.PixelSize);
h5writeatt(h5_file, '/', 'PixelSize_units', 'microns');
h5writeatt(h5_file, '/', 'NFrames', SMD.NFrames);
h5writeatt(h5_file, '/', 'NLocalizations', length(SMD.X));
h5writeatt(h5_file, '/', 'DatasetName', SMD.DatasetName);

fprintf('Exported to HDF5: %s\n', h5_file);

% Verify structure
h5disp(h5_file);
```

### Hierarchical HDF5 Structure

Organize data in groups for complex exports:

```matlab
% Load results
load('Results.mat', 'SMD');

h5_file = 'results_organized.h5';
if isfile(h5_file)
    delete(h5_file);
end

% Create groups
% Localizations group
loc_group = '/localizations';
h5create(h5_file, [loc_group '/positions_x'], size(SMD.X));
h5write(h5_file, [loc_group '/positions_x'], SMD.X);

h5create(h5_file, [loc_group '/positions_y'], size(SMD.Y));
h5write(h5_file, [loc_group '/positions_y'], SMD.Y);

h5create(h5_file, [loc_group '/frame_numbers'], size(SMD.FrameNum));
h5write(h5_file, [loc_group '/frame_numbers'], double(SMD.FrameNum));

% Photometry group
phot_group = '/photometry';
h5create(h5_file, [phot_group '/photons'], size(SMD.Photons));
h5write(h5_file, [phot_group '/photons'], SMD.Photons);

h5create(h5_file, [phot_group '/background'], size(SMD.Bg));
h5write(h5_file, [phot_group '/background'], SMD.Bg);

% Precision group
prec_group = '/precision';
h5create(h5_file, [prec_group '/x_stderr'], size(SMD.X_SE));
h5write(h5_file, [prec_group '/x_stderr'], SMD.X_SE);

h5create(h5_file, [prec_group '/y_stderr'], size(SMD.Y_SE));
h5write(h5_file, [prec_group '/y_stderr'], SMD.Y_SE);

% Metadata group attributes
h5writeatt(h5_file, '/', 'version', '1.0');
h5writeatt(h5_file, '/', 'software', 'smite');
h5writeatt(h5_file, loc_group, 'units', 'pixels');
h5writeatt(h5_file, loc_group, 'pixel_size_um', SMD.PixelSize);

fprintf('Exported hierarchical HDF5 structure\n');
```

### Export Trajectories to HDF5

Export tracking results preserving trajectory structure:

```matlab
% Load tracking results
load('tracking_Results.mat', 'TR');

h5_file = 'trajectories.h5';
if isfile(h5_file)
    delete(h5_file);
end

% Create trajectory group for each track
for i = 1:length(TR)
    track_group = sprintf('/track_%04d', i);

    % Write trajectory data
    h5create(h5_file, [track_group '/x'], size(TR(i).X));
    h5write(h5_file, [track_group '/x'], TR(i).X);

    h5create(h5_file, [track_group '/y'], size(TR(i).Y));
    h5write(h5_file, [track_group '/y'], TR(i).Y);

    h5create(h5_file, [track_group '/frame'], size(TR(i).FrameNum));
    h5write(h5_file, [track_group '/frame'], double(TR(i).FrameNum));

    h5create(h5_file, [track_group '/photons'], size(TR(i).Photons));
    h5write(h5_file, [track_group '/photons'], TR(i).Photons);

    % Track metadata
    h5writeatt(h5_file, track_group, 'track_length', length(TR(i).X));
    h5writeatt(h5_file, track_group, 'start_frame', min(TR(i).FrameNum));
    h5writeatt(h5_file, track_group, 'end_frame', max(TR(i).FrameNum));
end

% Global metadata
h5writeatt(h5_file, '/', 'n_trajectories', length(TR));

fprintf('Exported %d trajectories to HDF5\n', length(TR));
```

### Python-Compatible HDF5 Export

Format optimized for Python/pandas reading:

```matlab
% Load results
load('Results.mat', 'SMD');

h5_file = 'results_python.h5';
if isfile(h5_file)
    delete(h5_file);
end

% Create dataset group (pandas-compatible)
dset = '/data';

% Store as columns (easier for pandas DataFrame)
h5create(h5_file, [dset '/frame'], size(SMD.FrameNum));
h5write(h5_file, [dset '/frame'], double(SMD.FrameNum));

h5create(h5_file, [dset '/x'], size(SMD.X));
h5write(h5_file, [dset '/x'], SMD.X);

h5create(h5_file, [dset '/y'], size(SMD.Y));
h5write(h5_file, [dset '/y'], SMD.Y);

h5create(h5_file, [dset '/photons'], size(SMD.Photons));
h5write(h5_file, [dset '/photons'], SMD.Photons);

h5create(h5_file, [dset '/sigma_x'], size(SMD.X_SE));
h5write(h5_file, [dset '/sigma_x'], SMD.X_SE);

h5create(h5_file, [dset '/sigma_y'], size(SMD.Y_SE));
h5write(h5_file, [dset '/sigma_y'], SMD.Y_SE);

% Store metadata separately
h5writeatt(h5_file, dset, 'pixel_size_um', SMD.PixelSize);
h5writeatt(h5_file, dset, 'n_frames', SMD.NFrames);

fprintf('Created Python-compatible HDF5 file\n');
```

**Python reading example:**
```python
import h5py
import pandas as pd

# Read HDF5 file
with h5py.File('results_python.h5', 'r') as f:
    data = pd.DataFrame({
        'frame': f['/data/frame'][:],
        'x': f['/data/x'][:],
        'y': f['/data/y'][:],
        'photons': f['/data/photons'][:],
        'sigma_x': f['/data/sigma_x'][:],
        'sigma_y': f['/data/sigma_y'][:]
    })

    # Read metadata
    pixel_size = f['/data'].attrs['pixel_size_um']

print(f"Loaded {len(data)} localizations")
print(f"Pixel size: {pixel_size} µm")
```

## Export to JSON

JSON is ideal for web applications and JavaScript-based tools.

### Basic JSON Export

Export for web visualization:

```matlab
% Load results
load('Results.mat', 'SMD');

% Downsample for web display if needed
max_points = 10000;
if length(SMD.X) > max_points
    indices = randperm(length(SMD.X), max_points);
    fprintf('Downsampling from %d to %d points for web display\n', ...
        length(SMD.X), max_points);
else
    indices = 1:length(SMD.X);
end

% Create structure for export
export_struct = struct();
export_struct.data = struct();
export_struct.data.x = SMD.X(indices);
export_struct.data.y = SMD.Y(indices);
export_struct.data.photons = SMD.Photons(indices);
export_struct.data.frames = SMD.FrameNum(indices);

% Metadata
export_struct.metadata = struct();
export_struct.metadata.pixel_size = SMD.PixelSize;
export_struct.metadata.pixel_size_units = 'microns';
export_struct.metadata.n_frames = SMD.NFrames;
export_struct.metadata.n_localizations = length(indices);
export_struct.metadata.total_localizations = length(SMD.X);
export_struct.metadata.dataset = SMD.DatasetName;

% Convert to JSON with pretty formatting
json_text = jsonencode(export_struct, 'PrettyPrint', true);

% Save
fid = fopen('results.json', 'w');
fprintf(fid, '%s', json_text);
fclose(fid);

fprintf('Exported to JSON: %d localizations\n', length(indices));
```

### Trajectory Export to JSON

Export tracks with time-series structure:

```matlab
% Load tracking results
load('tracking_Results.mat', 'TR');

% Convert trajectories to JSON-friendly format
tracks = cell(length(TR), 1);
for i = 1:length(TR)
    tracks{i} = struct();
    tracks{i}.id = i;
    tracks{i}.length = length(TR(i).X);
    tracks{i}.positions = struct('x', TR(i).X, 'y', TR(i).Y);
    tracks{i}.frames = TR(i).FrameNum;
    tracks{i}.photons = TR(i).Photons;
    tracks{i}.start_frame = min(TR(i).FrameNum);
    tracks{i}.end_frame = max(TR(i).FrameNum);
end

% Create export structure
export_struct = struct();
export_struct.trajectories = tracks;
export_struct.metadata = struct('n_trajectories', length(TR));

% Save
json_text = jsonencode(export_struct);
fid = fopen('trajectories.json', 'w');
fprintf(fid, '%s', json_text);
fclose(fid);

fprintf('Exported %d trajectories to JSON\n', length(TR));
```

## Tool-Specific Exports

### ThunderSTORM Format

Export in ThunderSTORM CSV format:

```matlab
% Load results
load('Results.mat', 'SMD');

% ThunderSTORM expects specific column names
% Units: positions in nm, sigma in nm, intensity in photons
T = table();
T.id = (1:length(SMD.X))';
T.frame = SMD.FrameNum;
T.x_nm = SMD.X * SMD.PixelSize * 1000;
T.y_nm = SMD.Y * SMD.PixelSize * 1000;
T.sigma_nm = (SMD.X_SE + SMD.Y_SE) / 2 * SMD.PixelSize * 1000;
T.intensity = SMD.Photons;
T.offset = SMD.Bg;
T.bkgstd = ones(size(SMD.X)) * 50;  % Estimate if not available
T.uncertainty_nm = T.sigma_nm;

writetable(T, 'thunderstorm_localizations.csv');
fprintf('Exported %d localizations in ThunderSTORM format\n', length(SMD.X));
```

### rapidSTORM Format

Export in rapidSTORM text format:

```matlab
% Load results
load('Results.mat', 'SMD');

% rapidSTORM format: x y frame intensity
% Units: x,y in nm, intensity in photons
T = table();
T.x = SMD.X * SMD.PixelSize * 1000;  % Convert to nm
T.y = SMD.Y * SMD.PixelSize * 1000;
T.frame = SMD.FrameNum;
T.intensity = SMD.Photons;

% Write with tab delimiter (rapidSTORM preference)
writetable(T, 'rapidstorm_localizations.txt', 'Delimiter', '\t');
fprintf('Exported to rapidSTORM format\n');
```

### ViSP (Visualization and Statistical analysis of molecular Position) Format

Export for ViSP analysis:

```matlab
% Load results
load('Results.mat', 'SMD');

% ViSP expects: X(nm), Y(nm), Frame, Photons, Precision(nm)
T = table();
T.X_nm = SMD.X * SMD.PixelSize * 1000;
T.Y_nm = SMD.Y * SMD.PixelSize * 1000;
T.Frame = SMD.FrameNum;
T.Photons = SMD.Photons;
T.Precision_nm = sqrt(SMD.X_SE.^2 + SMD.Y_SE.^2) * SMD.PixelSize * 1000;

% Optional columns
if isfield(SMD, 'Z')
    T.Z_nm = SMD.Z * SMD.PixelSize * 1000;
end

writetable(T, 'visp_localizations.csv');
fprintf('Exported for ViSP analysis\n');
```

### PALMTracer Format

Export trajectories for PALMTracer:

```matlab
% Load tracking results
load('tracking_Results.mat', 'TR');

% PALMTracer format: TrackID, Frame, X(nm), Y(nm), Intensity
all_data = [];
for i = 1:length(TR)
    n_points = length(TR(i).X);
    track_data = [repmat(i, n_points, 1), ...
                  TR(i).FrameNum, ...
                  TR(i).X * SMD.PixelSize * 1000, ...
                  TR(i).Y * SMD.PixelSize * 1000, ...
                  TR(i).Photons];
    all_data = [all_data; track_data];
end

T = array2table(all_data, 'VariableNames', ...
    {'TrackID', 'Frame', 'X_nm', 'Y_nm', 'Intensity'});

writetable(T, 'palmtracer_trajectories.csv');
fprintf('Exported %d tracks for PALMTracer\n', length(TR));
```

## Advanced Export Techniques

### Batch Export Multiple Files

Export multiple result files at once:

```matlab
% Define results directory
results_dir = 'C:\Data\Experiment\Results';
result_files = dir(fullfile(results_dir, '*_Results.mat'));

% Export each file to CSV
for i = 1:length(result_files)
    % Load
    file_path = fullfile(result_files(i).folder, result_files(i).name);
    load(file_path, 'SMD');

    % Create export filename
    [~, name, ~] = fileparts(result_files(i).name);
    csv_name = [name '.csv'];

    % Export
    T = table(SMD.FrameNum, SMD.X, SMD.Y, SMD.Photons, ...
        'VariableNames', {'Frame', 'X', 'Y', 'Photons'});
    writetable(T, fullfile(results_dir, csv_name));

    fprintf('Exported %s: %d localizations\n', csv_name, length(SMD.X));
end

fprintf('Batch export complete: %d files\n', length(result_files));
```

### Export with Custom Filtering

Apply filters during export:

```matlab
% Load results
load('Results.mat', 'SMD');

% Define filters
photon_min = 500;
photon_max = 10000;
precision_max_nm = 20;

% Apply filters
pixel_size_nm = SMD.PixelSize * 1000;
precision_nm = sqrt(SMD.X_SE.^2 + SMD.Y_SE.^2) * pixel_size_nm;

valid = (SMD.Photons >= photon_min) & ...
        (SMD.Photons <= photon_max) & ...
        (precision_nm <= precision_max_nm);

fprintf('Filtering: %d -> %d localizations\n', length(SMD.X), sum(valid));

% Export filtered data
T = table(SMD.X(valid), SMD.Y(valid), SMD.Photons(valid), ...
    SMD.FrameNum(valid), precision_nm(valid), ...
    'VariableNames', {'X', 'Y', 'Photons', 'Frame', 'Precision_nm'});

writetable(T, 'localizations_filtered.csv');
```

### Export ROI Subsets

Export only localizations within regions of interest:

```matlab
% Load results
load('Results.mat', 'SMD');

% Define ROIs [x_min, x_max, y_min, y_max] in pixels
rois = [
    50, 100, 50, 100;   % ROI 1
    150, 200, 150, 200; % ROI 2
    100, 150, 200, 250  % ROI 3
];

% Export each ROI separately
for i = 1:size(rois, 1)
    roi = rois(i, :);
    in_roi = (SMD.X >= roi(1)) & (SMD.X <= roi(2)) & ...
             (SMD.Y >= roi(3)) & (SMD.Y <= roi(4));

    fprintf('ROI %d: %d localizations\n', i, sum(in_roi));

    T = table(SMD.X(in_roi), SMD.Y(in_roi), SMD.Photons(in_roi), ...
        'VariableNames', {'X', 'Y', 'Photons'});
    writetable(T, sprintf('roi_%d_localizations.csv', i));
end
```

### Export with Coordinate Transformations

Apply transformations during export:

```matlab
% Load results
load('Results.mat', 'SMD');

% Define coordinate transformation (rotation and shift)
theta = 15 * pi/180;  % 15 degree rotation
shift_x = 10;  % pixels
shift_y = 5;

% Apply transformation
R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
coords = [SMD.X(:), SMD.Y(:)];
coords_transformed = (R * coords')';

T = table();
T.X_original = SMD.X;
T.Y_original = SMD.Y;
T.X_transformed = coords_transformed(:, 1) + shift_x;
T.Y_transformed = coords_transformed(:, 2) + shift_y;
T.Frame = SMD.FrameNum;
T.Photons = SMD.Photons;

writetable(T, 'localizations_transformed.csv');
fprintf('Exported with coordinate transformation\n');
```

## Complete Workflow Examples

### Example 1: Multi-Format Export Pipeline

Export results to all common formats:

```matlab
% Load results
results_file = 'experiment1_Results.mat';
load(results_file, 'SMD');

fprintf('Exporting %d localizations...\n', length(SMD.X));

% 1. CSV for spreadsheets
T = table(SMD.FrameNum, SMD.X, SMD.Y, SMD.Photons, SMD.X_SE, SMD.Y_SE);
writetable(T, 'export_basic.csv');
fprintf('  [1/3] CSV export complete\n');

% 2. HDF5 for Python
h5_file = 'export_data.h5';
if isfile(h5_file), delete(h5_file); end
h5create(h5_file, '/x', size(SMD.X));
h5write(h5_file, '/x', SMD.X);
h5create(h5_file, '/y', size(SMD.Y));
h5write(h5_file, '/y', SMD.Y);
h5create(h5_file, '/photons', size(SMD.Photons));
h5write(h5_file, '/photons', SMD.Photons);
h5writeatt(h5_file, '/', 'pixel_size', SMD.PixelSize);
fprintf('  [2/3] HDF5 export complete\n');

% 3. JSON for web
export_struct.data = struct('x', SMD.X, 'y', SMD.Y, 'photons', SMD.Photons);
export_struct.metadata = struct('pixel_size', SMD.PixelSize);
fid = fopen('export_data.json', 'w');
fprintf(fid, '%s', jsonencode(export_struct));
fclose(fid);
fprintf('  [3/3] JSON export complete\n');

fprintf('All exports complete!\n');
```

### Example 2: Export for Collaborative Analysis

Prepare data package for collaborators:

```matlab
% Load results
load('Results.mat', 'SMD', 'SMF');

% Create export directory
export_dir = 'Collaborative_Data_Export';
if ~isfolder(export_dir)
    mkdir(export_dir);
end

% 1. Full data in CSV with metadata
T_full = table();
T_full.Frame = SMD.FrameNum;
T_full.X_pixels = SMD.X;
T_full.Y_pixels = SMD.Y;
T_full.X_nm = SMD.X * SMD.PixelSize * 1000;
T_full.Y_nm = SMD.Y * SMD.PixelSize * 1000;
T_full.Photons = SMD.Photons;
T_full.Background = SMD.Bg;
T_full.X_Precision_nm = SMD.X_SE * SMD.PixelSize * 1000;
T_full.Y_Precision_nm = SMD.Y_SE * SMD.PixelSize * 1000;

writetable(T_full, fullfile(export_dir, 'full_localizations.csv'));

% 2. Metadata file
fid = fopen(fullfile(export_dir, 'metadata.txt'), 'w');
fprintf(fid, 'Dataset: %s\n', SMD.DatasetName);
fprintf(fid, 'Pixel Size: %.4f microns\n', SMD.PixelSize);
fprintf(fid, 'Total Frames: %d\n', SMD.NFrames);
fprintf(fid, 'Total Localizations: %d\n', length(SMD.X));
fprintf(fid, 'Camera: %s\n', SMD.CameraType);
fprintf(fid, 'Mean Photons: %.1f\n', mean(SMD.Photons));
fprintf(fid, 'Median Precision: %.2f nm\n', ...
    median(SMD.X_SE) * SMD.PixelSize * 1000);
fclose(fid);

% 3. README
fid = fopen(fullfile(export_dir, 'README.txt'), 'w');
fprintf(fid, 'SMLM Data Export Package\n');
fprintf(fid, '========================\n\n');
fprintf(fid, 'Files included:\n');
fprintf(fid, '- full_localizations.csv: Complete localization table\n');
fprintf(fid, '- metadata.txt: Dataset parameters\n\n');
fprintf(fid, 'Column descriptions:\n');
fprintf(fid, '- Frame: Acquisition frame number\n');
fprintf(fid, '- X_pixels, Y_pixels: Position in camera pixels\n');
fprintf(fid, '- X_nm, Y_nm: Position in nanometers\n');
fprintf(fid, '- Photons: Detected photons per localization\n');
fprintf(fid, '- X_Precision_nm, Y_Precision_nm: Localization precision\n');
fclose(fid);

fprintf('Export package created in: %s\n', export_dir);
```

### Example 3: Export for Publication

Create publication-ready data files with documentation:

```matlab
% Load results
load('Results.mat', 'SMD');

% High-quality CSV with complete information
T = table();
T.Frame = SMD.FrameNum;
T.X_nm = SMD.X * SMD.PixelSize * 1000;
T.Y_nm = SMD.Y * SMD.PixelSize * 1000;
T.Photons = round(SMD.Photons);
T.Background = round(SMD.Bg);
T.Localization_Precision_nm = round(sqrt(SMD.X_SE.^2 + SMD.Y_SE.^2) * ...
    SMD.PixelSize * 1000, 2);

% Add header with full metadata
filename = 'publication_dataset_S1.csv';
fid = fopen(filename, 'w');
fprintf(fid, '# Supplementary Dataset S1\n');
fprintf(fid, '# Single Molecule Localizations\n');
fprintf(fid, '#\n');
fprintf(fid, '# Dataset: %s\n', SMD.DatasetName);
fprintf(fid, '# Acquisition frames: %d\n', SMD.NFrames);
fprintf(fid, '# Pixel size: %.4f µm (%.1f nm)\n', SMD.PixelSize, SMD.PixelSize*1000);
fprintf(fid, '# Camera: %s\n', SMD.CameraType);
fprintf(fid, '# Total localizations: %d\n', length(SMD.X));
fprintf(fid, '# Analysis software: smite (https://github.com/LidkeLab/smite)\n');
fprintf(fid, '#\n');
fprintf(fid, '# Column definitions:\n');
fprintf(fid, '# Frame: Acquisition frame number (1-indexed)\n');
fprintf(fid, '# X_nm, Y_nm: Molecular position in nanometers\n');
fprintf(fid, '# Photons: Detected photon count\n');
fprintf(fid, '# Background: Local background level (photons)\n');
fprintf(fid, '# Localization_Precision_nm: Position uncertainty (nm, 1σ)\n');
fprintf(fid, '#\n');
fclose(fid);

% Append data table
writetable(T, filename, 'WriteMode', 'Append');

fprintf('Publication dataset created: %s\n', filename);
fprintf('  %d localizations with full documentation\n', length(SMD.X));
```

## Troubleshooting

### Issue: CSV file too large for Excel

Excel has 1,048,576 row limit. For larger datasets:

```matlab
% Option 1: Split into multiple files
chunk_size = 1000000;
n_chunks = ceil(length(SMD.X) / chunk_size);

for i = 1:n_chunks
    idx_start = (i-1)*chunk_size + 1;
    idx_end = min(i*chunk_size, length(SMD.X));
    indices = idx_start:idx_end;

    T = table(SMD.X(indices), SMD.Y(indices), SMD.Photons(indices));
    writetable(T, sprintf('localizations_part%d.csv', i));
end

% Option 2: Downsample
sample_rate = 0.1;  % Keep 10%
indices = randperm(length(SMD.X), round(length(SMD.X) * sample_rate));
T = table(SMD.X(indices), SMD.Y(indices), SMD.Photons(indices));
writetable(T, 'localizations_downsampled.csv');
```

### Issue: HDF5 export fails with large arrays

```matlab
% Use chunking for large datasets
chunk_size = [100000, 1];  % Chunk along first dimension

h5create(h5_file, '/X', size(SMD.X), 'ChunkSize', chunk_size, ...
    'Deflate', 5);  % Add compression
h5write(h5_file, '/X', SMD.X);
```

### Issue: JSON file too large

```matlab
% Export summary statistics instead of raw data
summary = struct();
summary.n_localizations = length(SMD.X);
summary.mean_photons = mean(SMD.Photons);
summary.bounds = [min(SMD.X), max(SMD.X), min(SMD.Y), max(SMD.Y)];

% Or compress data by rounding
export_struct.data.x = round(SMD.X, 2);  % 2 decimal places
export_struct.data.y = round(SMD.Y, 2);
```

### Issue: Coordinate units confusion

Always document units clearly:

```matlab
% Include units in column names
T = table();
T.x_pixels = SMD.X;
T.y_pixels = SMD.Y;
T.x_nm = SMD.X * SMD.PixelSize * 1000;
T.y_nm = SMD.Y * SMD.PixelSize * 1000;

% Add metadata comment
fid = fopen('localizations_with_units.csv', 'w');
fprintf(fid, '# IMPORTANT: Pixel size = %.4f microns\n', SMD.PixelSize);
fprintf(fid, '# Columns with _pixels are in camera pixel units\n');
fprintf(fid, '# Columns with _nm are in nanometers\n');
fclose(fid);
writetable(T, 'localizations_with_units.csv', 'WriteMode', 'Append');
```

## Best Practices

1. **Always include metadata**: Pixel size, frame count, and units
2. **Document column meanings**: Use clear variable names and README files
3. **Test imports**: Verify exported files open correctly in target application
4. **Preserve precision**: Don't over-round position data
5. **Include uncertainties**: Export localization precision when available
6. **Use compression for HDF5**: Saves storage for large datasets
7. **Validate exports**: Check row counts and value ranges match original
8. **Version control**: Include software version and export date

## See Also

- [Save and Load SMD](save-and-load-smd.md) - Native MATLAB format persistence
- [Visualize Results](visualize-results.md) - Creating publication figures
- [SMD Structure](../core-concepts/smd-structure.md) - Understanding field meanings
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete analysis pipeline
