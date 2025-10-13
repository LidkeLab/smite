---
title: "Troubleshooting Memory Issues"
category: "troubleshooting"
level: "intermediate"
tags: ["troubleshooting", "memory", "gpu", "performance", "optimization", "errors"]
prerequisites: ["../core-concepts/smf-structure.md", "../how-to/use-gpu.md"]
related: ["../how-to/localize-molecules.md", "../workflows/smlm-analysis.md", "../workflows/batch-processing.md"]
summary: "Diagnose and resolve memory-related errors including out-of-memory errors, GPU exhaustion, and memory leaks"
estimated_time: "15-20 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Troubleshooting Memory Issues

## Purpose

Memory-related errors are among the most common issues when processing large SMLM and SPT datasets in smite. This guide provides comprehensive diagnosis and solutions for memory problems including MATLAB out-of-memory errors, GPU memory exhaustion, memory leaks, and strategies for handling large datasets. Understanding memory management is critical for processing high-resolution, multi-frame datasets efficiently.

## Prerequisites

- Understanding of [SMF structure](../core-concepts/smf-structure.md)
- Familiarity with [GPU acceleration](../how-to/use-gpu.md)
- Basic knowledge of MATLAB memory management
- Understanding of your system's memory capacity

## Overview of Memory in smite

smite uses two types of memory:

**System RAM (CPU memory):**
- Loading raw data from disk
- Creating image stacks
- Storing SMD results
- Frame connection operations
- Drift correction calculations
- Visualization and plotting

**GPU memory:**
- Box finding (FindROI)
- Gaussian MLE fitting (GaussMLE)
- PSF calculations
- Temporary buffers for CUDA operations

Most memory issues occur during localization, where both CPU and GPU memory are heavily utilized. smite includes automatic memory chunking for GPU operations, but large datasets can still cause problems.

## Problem 1: MATLAB Out of Memory Errors

### Symptoms

```matlab
Error using zeros
Out of memory. The likely cause is insufficient memory available in the system.

Error in smi_core.LocalizeData>genLocalizations
```

Or:

```matlab
Error: Requested 12000x12000x5000 (2.7GB) array exceeds maximum array size preference
```

### Diagnosis

Check available system memory:

```matlab
% Check memory status
[user, sys] = memory;

fprintf('System memory available: %.2f GB\n', sys.PhysicalMemory.Available / 1e9);
fprintf('MATLAB memory limit: %.2f GB\n', user.MaxPossibleArrayBytes / 1e9);
fprintf('Memory currently used: %.2f GB\n', (sys.PhysicalMemory.Total - sys.PhysicalMemory.Available) / 1e9);
```

Check data size:

```matlab
% Check size of your data
data_info = whos('sequence');
fprintf('Data size: %.2f GB\n', data_info.bytes / 1e9);

% Or for loaded h5 file
h5_info = h5info('Data_001.h5', '/sequence');
data_size_bytes = prod(h5_info.Dataspace.Size) * 4;  % 4 bytes for single precision
fprintf('H5 data size: %.2f GB\n', data_size_bytes / 1e9);
```

### Solutions

#### Solution 1: Process Data in Frame Batches

Process large datasets in smaller frame chunks:

```matlab
% Load and process in batches
filename = 'Data_001.h5';
h5_info = h5info(filename, '/sequence');
total_frames = h5_info.Dataspace.Size(3);
batch_size = 1000;  % Adjust based on available memory

% Initialize SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileName = filename;
% ... set other parameters ...

% Process in batches
SMD_all = smi_core.SingleMoleculeData.createSMD();

for batch_idx = 1:ceil(total_frames / batch_size)
    start_frame = (batch_idx - 1) * batch_size + 1;
    end_frame = min(batch_idx * batch_size, total_frames);

    fprintf('Processing frames %d to %d of %d...\n', start_frame, end_frame, total_frames);

    % Load only this batch
    sequence_batch = h5read(filename, '/sequence', [1, 1, start_frame], ...
        [h5_info.Dataspace.Size(1), h5_info.Dataspace.Size(2), end_frame - start_frame + 1]);

    % Process batch
    LD = smi_core.LocalizeData(sequence_batch, SMF);
    SMD_batch = LD.genLocalizations();

    % Adjust frame numbers to absolute frame indices
    SMD_batch.FrameNum = SMD_batch.FrameNum + start_frame - 1;

    % Combine with previous results
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_batch);

    % Clear batch data to free memory
    clear sequence_batch SMD_batch LD;

    fprintf('Batch %d complete. Total localizations: %d\n', batch_idx, length(SMD_all.X));
end

% Continue with drift correction, frame connection, etc.
SMD_final = SMD_all;
```

**Key points:**
- Batch size of 500-1000 frames typically works well
- Adjust based on your frame size (256x256 vs 512x512 vs larger)
- Frame numbers must be adjusted after batch processing
- Clear variables after each batch

#### Solution 2: Use ROI-Based Processing

Process only a region of interest instead of entire frames:

```matlab
% Define ROI (pixels)
ROI_XStart = 100;
ROI_YStart = 100;
ROI_Width = 256;
ROI_Height = 256;

% Load only ROI
sequence_roi = h5read('Data_001.h5', '/sequence', ...
    [ROI_YStart, ROI_XStart, 1], ...
    [ROI_Height, ROI_Width, Inf]);  % Inf loads all frames

% Process ROI
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileName = 'Data_001.h5';
% ... set parameters ...

LD = smi_core.LocalizeData(sequence_roi, SMF);
SMD_roi = LD.genLocalizations();

% Adjust coordinates to full frame
SMD_roi.X = SMD_roi.X + ROI_XStart - 1;
SMD_roi.Y = SMD_roi.Y + ROI_YStart - 1;
```

**Multiple ROIs approach:**

```matlab
% Process multiple non-overlapping ROIs
ROI_list = [
    1,   1,   256, 256;
    257, 1,   256, 256;
    1,   257, 256, 256;
    257, 257, 256, 256
];

SMD_all = smi_core.SingleMoleculeData.createSMD();

for roi_idx = 1:size(ROI_list, 1)
    fprintf('Processing ROI %d of %d...\n', roi_idx, size(ROI_list, 1));

    roi = ROI_list(roi_idx, :);
    sequence_roi = h5read('Data_001.h5', '/sequence', ...
        [roi(2), roi(1), 1], [roi(4), roi(3), Inf]);

    LD = smi_core.LocalizeData(sequence_roi, SMF);
    SMD_roi = LD.genLocalizations();

    % Adjust coordinates
    SMD_roi.X = SMD_roi.X + roi(1) - 1;
    SMD_roi.Y = SMD_roi.Y + roi(2) - 1;

    % Combine
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_roi);

    clear sequence_roi SMD_roi LD;
end
```

#### Solution 3: Downsample Data

For initial parameter optimization or when full resolution isn't critical:

```matlab
% Load and downsample
sequence = h5read('Data_001.h5', '/sequence');

% Temporal downsampling (use every Nth frame)
temporal_stride = 2;  % Use every 2nd frame
sequence_downsampled = sequence(:, :, 1:temporal_stride:end);

% Or spatial binning (2x2 binning)
sequence_binned = zeros(size(sequence, 1)/2, size(sequence, 2)/2, size(sequence, 3), 'single');
for i = 1:size(sequence, 3)
    sequence_binned(:, :, i) = (sequence(1:2:end, 1:2:end, i) + ...
                                  sequence(2:2:end, 1:2:end, i) + ...
                                  sequence(1:2:end, 2:2:end, i) + ...
                                  sequence(2:2:end, 2:2:end, i)) / 4;
end

% Adjust parameters for binned data
SMF.Data.PixelSize = SMF.Data.PixelSize * 2;  % Effective pixel size doubled
SMF.Fitting.PSFSigma = SMF.Fitting.PSFSigma / 2;  % PSF sigma in pixels halved

% Process downsampled data
LD = smi_core.LocalizeData(sequence_binned, SMF);
SMD = LD.genLocalizations();

% Scale coordinates back if needed
SMD.X = SMD.X * 2;
SMD.Y = SMD.Y * 2;
```

#### Solution 4: Increase MATLAB Memory Limit (Windows)

On Windows, increase Java heap memory:

```matlab
% Check current limit
maxmem = java.lang.Runtime.getRuntime.maxMemory / 1e9;
fprintf('MATLAB max memory: %.2f GB\n', maxmem);

% Increase in preferences
% Home tab > Preferences > MATLAB > General > Java Heap Memory
% Restart MATLAB after changing
```

#### Solution 5: Free Unused Memory

Explicitly clear variables and force garbage collection:

```matlab
% Clear large temporary variables
clear sequence raw_data temp_stack;

% Force MATLAB to free memory
pack;  % Memory defragmentation (may be slow)

% Or just garbage collection
java.lang.System.gc();
```

## Problem 2: GPU Out of Memory Errors

### Symptoms

```matlab
Error using gpuArray
Out of memory on device.

Error using parallel.gpu.CUDAKernel/feval
CUDA_ERROR_OUT_OF_MEMORY
```

Or during FindROI or GaussMLE:

```matlab
Error in smi_core.FindROI>findROI
An unexpected error occurred during CUDA execution.
```

### Diagnosis

Check GPU memory status:

```matlab
% Check GPU memory
gpu = gpuDevice();
fprintf('GPU: %s\n', gpu.Name);
fprintf('Total memory: %.2f GB\n', gpu.TotalMemory / 1e9);
fprintf('Available memory: %.2f GB\n', gpu.AvailableMemory / 1e9);
fprintf('Memory used: %.2f GB\n', (gpu.TotalMemory - gpu.AvailableMemory) / 1e9);

% Check what's using GPU memory
% Any gpuArray variables in workspace
vars = whos;
gpu_vars = vars(strcmp({vars.class}, 'gpuArray'));
if ~isempty(gpu_vars)
    fprintf('\nGPU variables in workspace:\n');
    for i = 1:length(gpu_vars)
        fprintf('  %s: %.2f MB\n', gpu_vars(i).name, gpu_vars(i).bytes / 1e6);
    end
end
```

Estimate memory requirements:

```matlab
% For FindROI
image_size = size(sequence);
NCopies = 8;  % Out-of-place operations
bytes_needed = 4 * NCopies * prod(image_size);  % 4 bytes per float
fprintf('FindROI needs ~%.2f GB GPU memory\n', bytes_needed / 1e9);

% For GaussMLE
BoxSize = SMF.BoxFinding.BoxSize;
NParams = 4;  % XYNB fit type
n_localizations = 10000;  % Estimate
bytes_per_fit = 4 * (BoxSize^2 + NParams + NParams + 1);
bytes_needed = bytes_per_fit * n_localizations;
fprintf('GaussMLE needs ~%.2f GB for %d fits\n', bytes_needed / 1e9, n_localizations);
```

### Solutions

#### Solution 1: Reset GPU Memory

Clear GPU memory and start fresh:

```matlab
% Reset GPU (clears all GPU memory)
gpu = gpuDevice();
reset(gpu);

% Wait for operations to complete
wait(gpu);

% Verify memory cleared
fprintf('Available memory after reset: %.2f GB\n', gpu.AvailableMemory / 1e9);
```

**Reset between processing batches:**

```matlab
for batch = 1:N_batches
    % Process batch
    LD = smi_core.LocalizeData(sequence_batch, SMF);
    SMD_batch = LD.genLocalizations();

    % Save or accumulate results
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_batch);

    % Clear GPU memory
    clear LD SMD_batch;
    reset(gpuDevice());

    fprintf('Batch %d complete. GPU memory cleared.\n', batch);
end
```

#### Solution 2: Reduce Box Size

Smaller boxes require less GPU memory for fitting:

```matlab
% Default box size
SMF.BoxFinding.BoxSize = 7;  % 7x7 = 49 pixels per box

% Reduce to 5x5 if memory constrained
SMF.BoxFinding.BoxSize = 5;  % 5x5 = 25 pixels per box (almost 50% less memory)

% Only do this if PSF fits within smaller box
% Rule of thumb: BoxSize should be >= 4 * PSFSigma
min_box_size = ceil(4 * max(SMF.Fitting.PSFSigma));
fprintf('Minimum recommended box size: %d\n', min_box_size);
```

**Note:** Smaller boxes may reduce fitting accuracy if PSF extends beyond box edges. Verify fit quality after reducing box size.

#### Solution 3: Process Fewer Frames at Once

smite automatically chunks data for GPU, but you can pre-chunk for more control:

```matlab
% Process in smaller frame batches
frames_per_batch = 500;  % Reduce if still having issues

h5_info = h5info('Data_001.h5', '/sequence');
total_frames = h5_info.Dataspace.Size(3);

SMD_all = smi_core.SingleMoleculeData.createSMD();

for batch_start = 1:frames_per_batch:total_frames
    batch_end = min(batch_start + frames_per_batch - 1, total_frames);

    % Load batch
    sequence_batch = h5read('Data_001.h5', '/sequence', ...
        [1, 1, batch_start], [Inf, Inf, batch_end - batch_start + 1]);

    % Reset GPU before processing
    reset(gpuDevice());

    % Process
    LD = smi_core.LocalizeData(sequence_batch, SMF);
    SMD_batch = LD.genLocalizations();
    SMD_batch.FrameNum = SMD_batch.FrameNum + batch_start - 1;

    % Combine and clear
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_batch);
    clear sequence_batch LD SMD_batch;

    fprintf('Processed frames %d-%d\n', batch_start, batch_end);
end
```

#### Solution 4: Close Other GPU Applications

Free GPU memory by closing other applications:

```matlab
% On command line (Windows):
% nvidia-smi
% Shows all processes using GPU

% On Linux/Mac:
% nvidia-smi
% kill <pid> to terminate processes
```

Close:
- Other MATLAB sessions using GPU
- Deep learning frameworks (TensorFlow, PyTorch)
- Video rendering applications
- GPU-accelerated browsers (Chrome, Firefox with GPU acceleration)
- Games or 3D modeling software

#### Solution 5: Use Smaller Image Regions

Process spatially smaller regions:

```matlab
% Divide large frame into quadrants
frame_size = [512, 512];
quadrant_size = frame_size / 2;

quadrants = {
    [1, quadrant_size(1), 1, quadrant_size(2)]  % top-left
    [quadrant_size(1)+1, frame_size(1), 1, quadrant_size(2)]  % top-right
    [1, quadrant_size(1), quadrant_size(2)+1, frame_size(2)]  % bottom-left
    [quadrant_size(1)+1, frame_size(1), quadrant_size(2)+1, frame_size(2)]  % bottom-right
};

SMD_all = smi_core.SingleMoleculeData.createSMD();

for q = 1:length(quadrants)
    rect = quadrants{q};

    % Load quadrant
    sequence_quad = h5read('Data_001.h5', '/sequence', ...
        [rect(3), rect(1), 1], ...
        [rect(4)-rect(3)+1, rect(2)-rect(1)+1, Inf]);

    % Process
    LD = smi_core.LocalizeData(sequence_quad, SMF);
    SMD_quad = LD.genLocalizations();

    % Adjust coordinates
    SMD_quad.X = SMD_quad.X + rect(1) - 1;
    SMD_quad.Y = SMD_quad.Y + rect(3) - 1;

    % Combine
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_quad);

    clear sequence_quad LD SMD_quad;
    reset(gpuDevice());
end
```

## Problem 3: Memory Leaks

### Symptoms

Memory usage gradually increases over time:

```matlab
% During batch processing or loops
Initial memory: 2.5 GB
After batch 1: 3.2 GB
After batch 2: 3.9 GB
After batch 3: 4.6 GB
...eventually crashes
```

Or:

```matlab
% MATLAB becomes progressively slower
% System becomes unresponsive
% "Out of memory" after processing multiple files
```

### Diagnosis

Monitor memory over time:

```matlab
% Track memory usage during processing
memory_log = [];

for batch = 1:N_batches
    % Get memory before
    [~, sys] = memory;
    mem_before = (sys.PhysicalMemory.Total - sys.PhysicalMemory.Available) / 1e9;

    % Process batch
    LD = smi_core.LocalizeData(sequence_batch, SMF);
    SMD_batch = LD.genLocalizations();

    % Get memory after
    [~, sys] = memory;
    mem_after = (sys.PhysicalMemory.Total - sys.PhysicalMemory.Available) / 1e9;

    memory_log(batch, :) = [mem_before, mem_after, mem_after - mem_before];

    fprintf('Batch %d: Memory %.2f GB -> %.2f GB (delta: %.2f GB)\n', ...
        batch, mem_before, mem_after, mem_after - mem_before);
end

% Plot memory usage
figure;
plot(memory_log(:, 2));
xlabel('Batch number');
ylabel('Memory usage (GB)');
title('Memory usage over time');
```

### Solutions

#### Solution 1: Explicit Variable Clearing

Clear variables in every loop iteration:

```matlab
for batch = 1:N_batches
    % Load and process
    sequence_batch = loadBatch(batch);
    LD = smi_core.LocalizeData(sequence_batch, SMF);
    SMD_batch = LD.genLocalizations();

    % Save or accumulate results
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_batch);

    % CRITICAL: Clear all temporary variables
    clear sequence_batch LD SMD_batch;

    % Also clear GPU memory
    if batch mod 5 == 0  % Every 5 batches
        reset(gpuDevice());
    end
end
```

#### Solution 2: Avoid Growing Arrays in Loops

Don't dynamically grow arrays:

```matlab
% BAD: Array grows each iteration (slow and memory inefficient)
SMD_array = [];
for i = 1:N
    SMD_batch = processBatch(i);
    SMD_array = [SMD_array; SMD_batch];  % BAD
end

% GOOD: Use catSMD which is optimized
SMD_all = smi_core.SingleMoleculeData.createSMD();
for i = 1:N
    SMD_batch = processBatch(i);
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_batch);
    clear SMD_batch;
end

% BETTER: Preallocate if size known
% (not always possible with localization data)
```

#### Solution 3: Process in Fresh MATLAB Sessions

For very large batch jobs, restart MATLAB periodically:

```matlab
% Main processing script: process_all.m
batch_files = dir('Data*.h5');
batches_per_session = 10;

for session = 1:ceil(length(batch_files) / batches_per_session)
    start_idx = (session - 1) * batches_per_session + 1;
    end_idx = min(session * batches_per_session, length(batch_files));

    % Create batch processing script
    fid = fopen('temp_batch.m', 'w');
    fprintf(fid, 'addpath ''~/Documents/MATLAB/smite/MATLAB'';\n');
    fprintf(fid, 'setupSMITE;\n');
    fprintf(fid, 'process_batch(%d, %d);\n', start_idx, end_idx);
    fprintf(fid, 'exit;\n');
    fclose(fid);

    % Run in fresh MATLAB
    system('matlab -nodisplay -r "run(''temp_batch.m'')"');

    fprintf('Completed session %d\n', session);
end
```

#### Solution 4: Use Persistent Variables Carefully

Avoid persistent variables that accumulate data:

```matlab
% BAD: Persistent variable keeps growing
function result = process_with_history(data)
    persistent history;
    if isempty(history)
        history = [];
    end
    history = [history; data];  % Grows forever
    result = analyze(history);
end

% GOOD: Clear persistent when needed
function result = process_with_history(data, clear_flag)
    persistent history;
    if nargin > 1 && clear_flag
        history = [];
        return;
    end
    if isempty(history)
        history = [];
    end
    history = data;  % Replace, don't append
    result = analyze(history);
end
```

#### Solution 5: Monitor and Profile

Use MATLAB's memory profiler:

```matlab
% Profile memory usage
profile -memory on;

% Run your analysis
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();

% View results
profile viewer;
% Look for functions with high memory allocation
```

## Problem 4: Large Dataset Strategies

### Working with Very Large Datasets (>100 GB)

#### Strategy 1: Streaming Processing

Process data without loading entire dataset:

```matlab
% Process frame by frame (slowest but most memory efficient)
h5_info = h5info('huge_data.h5', '/sequence');
frame_size = h5_info.Dataspace.Size(1:2);
total_frames = h5_info.Dataspace.Size(3);

% Process in chunks for efficiency
chunk_size = 100;  % 100 frames at a time

SMD_all = smi_core.SingleMoleculeData.createSMD();

for chunk_start = 1:chunk_size:total_frames
    chunk_end = min(chunk_start + chunk_size - 1, total_frames);

    % Load only this chunk
    sequence_chunk = h5read('huge_data.h5', '/sequence', ...
        [1, 1, chunk_start], [Inf, Inf, chunk_end - chunk_start + 1]);

    % Process
    LD = smi_core.LocalizeData(sequence_chunk, SMF);
    SMD_chunk = LD.genLocalizations();
    SMD_chunk.FrameNum = SMD_chunk.FrameNum + chunk_start - 1;

    % Save intermediate results
    if mod(chunk_start, 1000) == 1
        save(sprintf('intermediate_results_%06d.mat', chunk_start), 'SMD_all');
    end

    % Combine
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_chunk);

    clear sequence_chunk LD SMD_chunk;
    reset(gpuDevice());

    fprintf('Processed %d/%d frames\n', chunk_end, total_frames);
end
```

#### Strategy 2: Parallel Processing Across Multiple Machines

Distribute processing to multiple computers:

```matlab
% On each machine, process a subset of frames
machine_id = 1;  % Change for each machine
n_machines = 4;

h5_info = h5info('huge_data.h5', '/sequence');
total_frames = h5_info.Dataspace.Size(3);

frames_per_machine = ceil(total_frames / n_machines);
start_frame = (machine_id - 1) * frames_per_machine + 1;
end_frame = min(machine_id * frames_per_machine, total_frames);

% Process only this machine's frames
sequence_subset = h5read('huge_data.h5', '/sequence', ...
    [1, 1, start_frame], [Inf, Inf, end_frame - start_frame + 1]);

LD = smi_core.LocalizeData(sequence_subset, SMF);
SMD_subset = LD.genLocalizations();
SMD_subset.FrameNum = SMD_subset.FrameNum + start_frame - 1;

% Save machine-specific results
save(sprintf('SMD_machine_%d.mat', machine_id), 'SMD_subset', 'SMF');

% Later, combine all results
% On main machine:
SMD_combined = smi_core.SingleMoleculeData.createSMD();
for m = 1:n_machines
    load(sprintf('SMD_machine_%d.mat', m));
    SMD_combined = smi_core.SingleMoleculeData.catSMD(SMD_combined, SMD_subset);
end
```

#### Strategy 3: Prefilter Frames

Process only frames with signal:

```matlab
% Identify frames with activity
h5_info = h5info('Data_001.h5', '/sequence');
total_frames = h5_info.Dataspace.Size(3);

active_frames = [];
batch_size = 100;

% Scan through to find active frames
for batch_start = 1:batch_size:total_frames
    batch_end = min(batch_start + batch_size - 1, total_frames);
    sequence_batch = h5read('Data_001.h5', '/sequence', ...
        [1, 1, batch_start], [Inf, Inf, batch_end - batch_start + 1]);

    % Calculate mean intensity per frame
    mean_intensity = squeeze(mean(mean(sequence_batch, 1), 2));

    % Identify frames above threshold
    threshold = median(mean_intensity) + 2 * std(mean_intensity);
    active_in_batch = find(mean_intensity > threshold) + batch_start - 1;
    active_frames = [active_frames; active_in_batch(:)];

    clear sequence_batch;
end

fprintf('Found %d active frames out of %d total (%.1f%%)\n', ...
    length(active_frames), total_frames, 100 * length(active_frames) / total_frames);

% Now process only active frames
sequence_active = zeros(h5_info.Dataspace.Size(1), ...
    h5_info.Dataspace.Size(2), length(active_frames), 'single');

for i = 1:length(active_frames)
    frame = h5read('Data_001.h5', '/sequence', [1, 1, active_frames(i)], [Inf, Inf, 1]);
    sequence_active(:, :, i) = frame;
end

LD = smi_core.LocalizeData(sequence_active, SMF);
SMD = LD.genLocalizations();

% Map back to original frame numbers
SMD.FrameNum = active_frames(SMD.FrameNum);
```

## Prevention Best Practices

### 1. Estimate Memory Before Processing

```matlab
function checkMemoryRequirements(filename, SMF)
    % Estimate memory needed for processing

    h5_info = h5info(filename, '/sequence');
    dims = h5_info.Dataspace.Size;

    % Data size
    data_size = prod(dims) * 4;  % 4 bytes per float

    % Estimate localizations (rough: 1 per 10x10 pixels per frame)
    est_locs = prod(dims) / 100;
    smd_size = est_locs * 8 * 20;  % 20 fields, 8 bytes per double

    % GPU memory (from FindROI)
    gpu_needed = data_size * 8;  % 8 copies for processing

    % Total system memory
    total_needed = data_size + smd_size + data_size * 2;  % 2x for overhead

    fprintf('=== Memory Requirements ===\n');
    fprintf('Data size: %.2f GB\n', data_size / 1e9);
    fprintf('Estimated SMD size: %.2f GB\n', smd_size / 1e9);
    fprintf('GPU memory needed: %.2f GB\n', gpu_needed / 1e9);
    fprintf('Total system memory: %.2f GB\n', total_needed / 1e9);

    % Check available
    [~, sys] = memory;
    gpu = gpuDevice();

    fprintf('\n=== Available ===\n');
    fprintf('System memory: %.2f GB\n', sys.PhysicalMemory.Available / 1e9);
    fprintf('GPU memory: %.2f GB\n', gpu.AvailableMemory / 1e9);

    % Recommend batching
    if total_needed > sys.PhysicalMemory.Available
        recommended_batches = ceil(total_needed / sys.PhysicalMemory.Available);
        fprintf('\n*** RECOMMENDATION: Process in %d batches ***\n', recommended_batches);
    else
        fprintf('\n*** Sufficient memory available ***\n');
    end
end
```

### 2. Monitor Memory During Processing

```matlab
function SMD = processWithMonitoring(sequence, SMF)
    % Process with memory monitoring

    fprintf('Starting processing...\n');
    [~, sys] = memory;
    gpu = gpuDevice();
    fprintf('Initial: System %.2f GB, GPU %.2f GB\n', ...
        sys.PhysicalMemory.Available / 1e9, gpu.AvailableMemory / 1e9);

    LD = smi_core.LocalizeData(sequence, SMF);
    SMD = LD.genLocalizations();

    [~, sys] = memory;
    gpu = gpuDevice();
    fprintf('After processing: System %.2f GB, GPU %.2f GB\n', ...
        sys.PhysicalMemory.Available / 1e9, gpu.AvailableMemory / 1e9);
end
```

### 3. Configure Conservative Memory Settings

```matlab
% Set up SMF with memory-conscious parameters
SMF = smi_core.SingleMoleculeFitting();

% Use smaller box size if appropriate
SMF.BoxFinding.BoxSize = 5;  % Instead of 7

% Disable verbose output that creates large displays
SMF.Data.DataVariable = 'sequence';

% Don't save intermediate results
% (handle saving manually when needed)
```

## See Also

- [GPU Usage Guide](../how-to/use-gpu.md) - GPU memory management details
- [Localize Molecules](../how-to/localize-molecules.md) - Core localization workflow
- [Batch Processing](../workflows/batch-processing.md) - Processing multiple datasets
- [SMLM Analysis](../workflows/smlm-analysis.md) - Complete SMLM workflow

## Quick Reference

**Common memory-saving approaches:**

| Problem | Solution | Memory Savings |
|---------|----------|----------------|
| Large dataset | Frame batching | 50-90% |
| Large frames | ROI processing | Proportional to ROI size |
| Many localizations | Process and save incrementally | 50-80% |
| GPU exhaustion | Reset between batches | Clears 100% GPU |
| Memory leaks | Explicit clearing | Prevents accumulation |

**Rule of thumb memory requirements:**

| Dataset | System RAM | GPU Memory |
|---------|------------|------------|
| 256x256, 10k frames | 2.5 GB | 2 GB |
| 512x512, 10k frames | 10 GB | 8 GB |
| 1024x1024, 10k frames | 40 GB | 32 GB |
| 2048x2048, 10k frames | 160 GB | Requires batching |

**Batch size recommendations:**

| Available RAM | Available GPU | Recommended Batch (256x256) |
|---------------|---------------|----------------------------|
| 8 GB | 2 GB | 500 frames |
| 16 GB | 4 GB | 1000 frames |
| 32 GB | 8 GB | 2000 frames |
| 64 GB | 16 GB | Full dataset (if <5k frames) |
