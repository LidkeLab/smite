---
title: "Troubleshooting GPU Problems"
category: "troubleshooting"
level: "intermediate"
tags: ["gpu", "cuda", "hardware", "errors", "performance", "debugging"]
prerequisites: ["../how-to/use-gpu.md"]
related: ["../getting-started/installation.md", "common-errors.md"]
summary: "Diagnose and resolve GPU-related issues including detection failures, CUDA errors, memory problems, and performance degradation"
estimated_time: "15 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Troubleshooting GPU Problems

## Purpose

GPU acceleration is essential for smite's core localization operations. When GPU issues occur, they can prevent analysis entirely or severely degrade performance. This guide provides systematic diagnosis and solutions for common GPU problems including detection failures, CUDA kernel errors, out-of-memory issues, PTX compilation problems, performance degradation, and multi-GPU complications.

## Prerequisites

- Understanding of [GPU acceleration in smite](../how-to/use-gpu.md)
- NVIDIA GPU with CUDA support
- MATLAB Parallel Computing Toolbox

## Overview

GPU problems in smite typically fall into six categories:

1. **GPU not detected**: MATLAB cannot find or access your GPU
2. **CUDA kernel errors**: PTX files missing or incompatible
3. **Out of memory errors**: Insufficient GPU memory for dataset
4. **PTX compilation failures**: Problems building CUDA kernels
5. **Performance degradation**: GPU available but running slowly
6. **Multiple GPU issues**: Problems selecting or using specific GPUs

Each section provides symptoms, diagnostic steps, and solutions.

## Problem 1: GPU Not Detected

### Symptoms

```matlab
>> gpuDevice()
Error using gpuDevice (line 26)
No supported GPU was detected.
```

Or during localization:

```matlab
>> SMD = LD.genLocalizations();
Error: GPU device not available for CUDA operations.
```

### Diagnosis

Run systematic checks to identify the root cause:

**Step 1: Check if GPU is physically present**

```matlab
% Check if MATLAB can see any GPU
if parallel.gpu.GPUDevice.isAvailable
    fprintf('MATLAB can access GPU\n');
else
    fprintf('MATLAB cannot access GPU\n');
end

% Try to get GPU count
try
    N = gpuDeviceCount;
    fprintf('GPU count: %d\n', N);
catch ME
    fprintf('Error getting GPU count: %s\n', ME.message);
end
```

**Step 2: Check GPU from system level**

```bash
# Windows/Linux command line
nvidia-smi

# Should show GPU name, driver version, CUDA version
# If this fails, problem is at driver level
```

**Step 3: Check Parallel Computing Toolbox**

```matlab
% Verify toolbox is available
if license('test', 'Distrib_Computing_Toolbox')
    fprintf('Parallel Computing Toolbox: Licensed\n');
else
    error('Parallel Computing Toolbox not available - required for GPU operations');
end

% Check toolbox version
ver('parallel')
```

**Step 4: Check GPU type**

```bash
# Linux/Windows
lspci | grep -i nvidia  # Linux
wmic path win32_VideoController get name  # Windows
```

smite requires NVIDIA GPUs. AMD, Intel, or Apple Silicon GPUs are not supported.

### Solutions

**Solution 1A: GPU driver not installed or outdated**

If `nvidia-smi` fails, GPU driver is the problem:

```bash
# Download latest driver from:
# https://www.nvidia.com/Download/index.aspx

# After installation, verify:
nvidia-smi
```

Then restart MATLAB:

```matlab
% Close and reopen MATLAB
% Test again:
gpuDevice()
```

**Solution 1B: Non-NVIDIA GPU**

smite requires NVIDIA GPU for CUDA operations:

```matlab
% Check what GPU you have
% If AMD, Intel, or Apple Silicon:
error('smite requires NVIDIA GPU with CUDA support. Your GPU is not compatible.');
```

**Options:**
- Use a different computer with NVIDIA GPU
- Install NVIDIA GPU if system supports it
- Use remote compute server with NVIDIA GPU

**Solution 1C: MATLAB cannot access GPU (permissions)**

On Linux, user may need to be added to video group:

```bash
# Add user to video group
sudo usermod -a -G video $USER

# Log out and log back in for changes to take effect
```

On Windows, ensure MATLAB has permissions:

```powershell
# Run MATLAB as administrator (right-click, "Run as administrator")
# Test GPU access
```

**Solution 1D: GPU disabled in BIOS**

Some systems allow disabling discrete GPU in BIOS:

1. Restart computer
2. Enter BIOS/UEFI (usually Del, F2, or F10 during boot)
3. Look for GPU settings
4. Ensure discrete GPU is enabled
5. Save and restart

**Solution 1E: Driver conflict after CUDA Toolkit installation**

If GPU worked before but stopped after installing CUDA Toolkit:

```bash
# The CUDA Toolkit can install an incompatible driver
# Uninstall CUDA Toolkit (keep GPU driver)
# Or reinstall driver:
# https://www.nvidia.com/Download/index.aspx
```

Then verify:

```matlab
gpu = gpuDevice();
fprintf('Driver version: %s\n', gpu.DriverVersion);
fprintf('Toolkit version: %s\n', gpu.ToolkitVersion);
```

## Problem 2: CUDA Kernel Errors

### Symptoms

```matlab
>> SMD = LD.genLocalizations();
Error: Unable to find PTX file: smi_cuda_gaussMLEv2.ptx

Or:

Error using parallel.gpu.CUDAKernel
The cu file could not be compiled.

Or:

Error: CUDA kernel launch failure
```

### Diagnosis

**Step 1: Check if PTX files exist**

```matlab
% List required PTX files
required_files = {
    'smi_cuda_FindROI.ptx'
    'smi_cuda_gaussMLEv2.ptx'
    'smi_cuda_gaussBlobROIStack.ptx'
    'smi_cuda_PSFSample3DBlob.ptx'
};

fprintf('Checking for PTX files:\n');
for i = 1:length(required_files)
    file_path = which(required_files{i});
    if ~isempty(file_path)
        fprintf('  %s: FOUND at %s\n', required_files{i}, file_path);
    else
        fprintf('  %s: NOT FOUND\n', required_files{i});
    end
end
```

**Step 2: Check if PTX directory is on path**

```matlab
% PTX files should be in smite/MATLAB/ptx
ptx_dir = fullfile(fileparts(which('setupSMITE')), 'ptx');
fprintf('Expected PTX directory: %s\n', ptx_dir);
fprintf('Directory exists: %d\n', isfolder(ptx_dir));

% Check if on path
path_cell = strsplit(path, pathsep);
on_path = any(strcmp(path_cell, ptx_dir));
fprintf('PTX directory on path: %d\n', on_path);
```

**Step 3: Test kernel loading**

```matlab
% Try to load a kernel
try
    ptx_file = fullfile(ptx_dir, 'smi_cuda_gaussMLEv2.ptx');
    kernel = parallel.gpu.CUDAKernel(ptx_file, 'smi_cuda_gaussMLEv2.cu');
    fprintf('Kernel loaded successfully\n');
catch ME
    fprintf('Kernel loading failed: %s\n', ME.message);
end
```

### Solutions

**Solution 2A: PTX files missing**

If PTX files don't exist:

```matlab
% Option 1: Reinstall smite
% Download fresh copy from GitHub
% Precompiled PTX files should be included

% Option 2: Recompile PTX files (requires CUDA Toolkit)
cd ~/Documents/MATLAB/smite/MATLAB/source/cuda
cuda_Make
```

**Solution 2B: Path not set correctly**

If PTX files exist but not found:

```matlab
% Re-run setup
setupSMITE

% Verify PTX directory added
ptx_dir = fullfile(fileparts(which('setupSMITE')), 'ptx');
addpath(ptx_dir);
savepath  % Save for future sessions

% Test again
which smi_cuda_gaussMLEv2.ptx
```

**Solution 2C: Compute capability mismatch**

If kernel fails to load due to architecture:

```matlab
% Check GPU compute capability
gpu = gpuDevice();
fprintf('Your GPU compute capability: %.1f\n', gpu.ComputeCapability);
fprintf('smite requires: >= 5.0\n');

if gpu.ComputeCapability < 5.0
    error('GPU compute capability too low. Upgrade GPU or use different machine.');
end
```

PTX files are compiled for compute capability 5.0+. If your GPU is older:

```matlab
% Must recompile with lower architecture
% Edit smite/MATLAB/source/cuda/cuda_Make.m
% Change: nvcc_cmd = 'nvcc -arch=sm_50 -ptx %s -o %s\n';
% To: nvcc_cmd = 'nvcc -arch=sm_35 -ptx %s -o %s\n';  % For compute 3.5
% Then recompile (requires CUDA Toolkit)
```

**Solution 2D: Kernel file corruption**

If PTX files present but kernel loading fails:

```matlab
% Check file integrity
ptx_file = which('smi_cuda_gaussMLEv2.ptx');
file_info = dir(ptx_file);
fprintf('File size: %d bytes\n', file_info.bytes);

if file_info.bytes < 1000  % PTX files should be several KB
    warning('PTX file appears corrupted or incomplete');
    fprintf('Recompile PTX files or reinstall smite\n');
end

% Delete and recompile
delete(ptx_file);
cd ~/Documents/MATLAB/smite/MATLAB/source/cuda
cuda_Make
```

## Problem 3: Out of Memory Errors

### Symptoms

```matlab
>> SMD = LD.genLocalizations();
Error: Out of memory on device.

Or:

Error using gpuArray
Out of memory on device.

Or:

CUDA error: out of memory (error code = 2)
```

### Diagnosis

**Step 1: Check available GPU memory**

```matlab
% Check memory before analysis
gpu = gpuDevice();
fprintf('Total GPU memory: %.2f GB\n', gpu.TotalMemory / 1e9);
fprintf('Available memory: %.2f GB\n', gpu.AvailableMemory / 1e9);
fprintf('Memory used: %.2f GB (%.1f%%)\n', ...
    (gpu.TotalMemory - gpu.AvailableMemory) / 1e9, ...
    100 * (1 - gpu.AvailableMemory / gpu.TotalMemory));

% Check if other processes using GPU
if gpu.AvailableMemory < 0.5e9  % Less than 500 MB available
    warning('Very little GPU memory available. Close other GPU applications.');
end
```

**Step 2: Estimate memory requirements**

```matlab
% For your dataset
[nY, nX, nFrames] = size(sequence);
fprintf('Data size: %d × %d × %d frames\n', nY, nX, nFrames);

% Estimate memory needed for FindROI
NCopies = 8;  % Number of temporary arrays
NBytesPerPixel = 4;  % single float
memory_findROI = NCopies * NBytesPerPixel * nY * nX / 1e9;
fprintf('FindROI needs ~%.2f GB\n', memory_findROI);

% Estimate memory for fitting (per 1000 localizations)
BoxSize = SMF.BoxFinding.BoxSize;
NParams = 5;  % For XYNB fit
memory_per_fit = 4 * (BoxSize^2 + NParams + NParams + 1);
fits_per_GB = 1e9 / memory_per_fit;
fprintf('Can fit ~%d localizations per GB\n', round(fits_per_GB));
```

**Step 3: Check for memory leaks**

```matlab
% Monitor memory over time
for i = 1:5
    gpu = gpuDevice();
    fprintf('Iteration %d: %.2f GB available\n', i, gpu.AvailableMemory / 1e9);
    pause(1);
end

% If memory decreases without operations, may have leak
```

### Solutions

**Solution 3A: Clear GPU memory**

```matlab
% Clear all GPU arrays
clear all  % Clears workspace including GPU arrays

% Reset GPU device
gpu = gpuDevice();
reset(gpu);

% Verify memory freed
fprintf('Available memory after reset: %.2f GB\n', gpu.AvailableMemory / 1e9);
```

**Solution 3B: Process smaller batches**

For large datasets, process in chunks:

```matlab
% Batch processing approach
NFrames = size(sequence, 3);
BatchSize = 500;  % Adjust based on available memory
NBatches = ceil(NFrames / BatchSize);

SMD_all = smi_core.SingleMoleculeData.createSMD();

for batch = 1:NBatches
    start_frame = (batch - 1) * BatchSize + 1;
    end_frame = min(batch * BatchSize, NFrames);

    fprintf('Processing batch %d/%d (frames %d-%d)\n', ...
        batch, NBatches, start_frame, end_frame);

    % Extract batch
    batch_data = sequence(:, :, start_frame:end_frame);

    % Process batch
    LD_batch = smi_core.LocalizeData(batch_data, SMF);
    SMD_batch = LD_batch.genLocalizations();

    % Adjust frame numbers
    SMD_batch.FrameNum = SMD_batch.FrameNum + start_frame - 1;

    % Combine results
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_batch);

    % Clear GPU memory between batches
    reset(gpuDevice());

    fprintf('  Found %d localizations\n', length(SMD_batch.X));
end

fprintf('Total localizations: %d\n', length(SMD_all.X));
```

**Solution 3C: Reduce spatial dimensions with ROI**

```matlab
% Use ROI to process smaller area
SMF.Data.DataROI = [50, 50, 200, 200];  % [YStart, XStart, YEnd, XEnd]

% This reduces memory by factor of:
full_pixels = nY * nX;
roi_pixels = (200-50+1) * (200-50+1);
reduction_factor = full_pixels / roi_pixels;
fprintf('Memory reduction: %.1fx\n', reduction_factor);
```

**Solution 3D: Reduce box size**

```matlab
% Smaller boxes use less memory
% Default is 7, try 5 if PSF permits
SMF.BoxFinding.BoxSize = 5;  % Instead of 7

% Memory savings:
% 7×7 = 49 pixels per fit
% 5×5 = 25 pixels per fit
% ~2× reduction in fitting memory
```

**Solution 3E: Close other GPU applications**

```bash
# Check what's using GPU
nvidia-smi

# Shows processes using GPU memory
# Close unnecessary applications:
# - Other MATLAB sessions
# - Deep learning frameworks (PyTorch, TensorFlow)
# - Video encoding/gaming applications
# - Desktop compositing (on Linux, may need to stop X server)
```

**Solution 3F: Upgrade GPU**

If dataset routinely exceeds GPU memory:

```matlab
% Calculate needed memory
peak_memory_GB = 4;  % Your typical requirement

% Your current GPU
gpu = gpuDevice();
fprintf('Current GPU: %s (%.1f GB)\n', gpu.Name, gpu.TotalMemory / 1e9);

if gpu.TotalMemory < peak_memory_GB * 1e9
    fprintf('Consider GPU with >= %.0f GB memory\n', peak_memory_GB);
    fprintf('Recommended GPUs:\n');
    fprintf('  RTX 3090 (24 GB)\n');
    fprintf('  RTX A5000 (24 GB)\n');
    fprintf('  RTX A6000 (48 GB)\n');
end
```

## Problem 4: PTX Compilation Failures

### Symptoms

```matlab
>> cd ~/Documents/MATLAB/smite/MATLAB/source/cuda
>> cuda_Make
'nvcc' is not recognized as an internal or external command.

Or:

nvcc fatal : Unsupported gpu architecture 'compute_50'

Or:

error: identifier "atomicAdd" is undefined (on older GPUs)
```

### Diagnosis

**Step 1: Check if nvcc is available**

```bash
# Command line
nvcc --version

# Should show CUDA compiler version
# If not found, CUDA Toolkit not installed or not on PATH
```

**Step 2: Check CUDA Toolkit compatibility**

```matlab
% Check GPU compute capability
gpu = gpuDevice();
fprintf('GPU compute capability: %.1f\n', gpu.ComputeCapability);

% Check MATLAB's supported CUDA version
% See: https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html
```

**Step 3: Check Visual Studio (Windows only)**

```powershell
# Windows command line
cl

# Should show Microsoft C/C++ Compiler
# If not found, Visual Studio C++ compiler not installed or not on PATH
```

### Solutions

**Solution 4A: Install CUDA Toolkit**

```bash
# Download from NVIDIA:
# https://developer.nvidia.com/cuda-downloads

# Linux installation:
sudo sh cuda_*_linux.run

# Windows: Run installer

# Verify installation:
nvcc --version
```

**Solution 4B: Add nvcc to PATH**

If CUDA Toolkit installed but `nvcc` not found:

```matlab
% Edit smite/MATLAB/source/cuda/cuda_Make.m

% Linux/MacOS - add line:
setenv('PATH', ['/usr/local/cuda/bin:' getenv('PATH')]);

% Windows - add line:
setenv('PATH', ['C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\bin;' getenv('PATH')]);

% Adjust version (v12.2) to match your installation
```

**Solution 4C: Install Visual Studio (Windows)**

Windows requires Visual Studio C++ compiler:

```powershell
# Download Visual Studio Community (free):
# https://visualstudio.microsoft.com/downloads/

# During installation, select:
# "Desktop development with C++"

# Add cl.exe to PATH in cuda_Make.m:
setenv('PATH', ['C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.35.32215\bin\Hostx64\x64;' getenv('PATH')]);

# Adjust path to match your VS installation
```

**Solution 4D: Use precompiled PTX files**

Avoid compilation if possible:

```matlab
% smite includes precompiled PTX files for:
% - Linux (64-bit)
% - MacOS (64-bit, Intel)
% - Windows (64-bit)

% If your platform matches, just use included PTX files
% Ensure smite/MATLAB/ptx is on MATLAB path:
addpath('~/Documents/MATLAB/smite/MATLAB/ptx');
savepath;

% Verify files present:
dir('~/Documents/MATLAB/smite/MATLAB/ptx/*.ptx')
```

**Solution 4E: Adjust compute capability**

If compilation fails due to architecture mismatch:

```matlab
% Edit smite/MATLAB/source/cuda/cuda_Make.m
% Current line:
nvcc_cmd = 'nvcc -arch=sm_50 -ptx %s -o %s\n';

% For newer GPUs (compute capability 8.0+):
nvcc_cmd = 'nvcc -arch=sm_80 -ptx %s -o %s\n';

% For older GPUs (compute capability 3.5+):
nvcc_cmd = 'nvcc -arch=sm_35 -ptx %s -o %s\n';

% Then recompile:
cuda_Make
```

Note: Using older architecture (sm_35) on newer GPUs may reduce performance.

## Problem 5: Slow GPU Performance

### Symptoms

```matlab
% GPU operations much slower than expected
% Processing 100 frames takes minutes instead of seconds
% GPU utilization shown as low in nvidia-smi
```

### Diagnosis

**Step 1: Benchmark actual performance**

```matlab
% Create test dataset
ImageSize = 256;
NFrames = 100;
sequence = randn(ImageSize, ImageSize, NFrames, 'single') * 10 + 100;

% Add synthetic molecules
for frame = 1:NFrames
    for mol = 1:50
        x = randi([20, ImageSize-20]);
        y = randi([20, ImageSize-20]);
        sequence(:, :, frame) = sequence(:, :, frame) + ...
            1000 * exp(-((meshgrid(1:ImageSize) - x).^2 + (meshgrid(1:ImageSize)' - y).^2) / (2*1.3^2));
    end
end

% Time GPU processing
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraGain = 1;
SMF.Data.CameraOffset = 100;

tic;
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();
elapsed = toc;

fprintf('Time: %.2f seconds\n', elapsed);
fprintf('Frames per second: %.1f\n', NFrames / elapsed);
fprintf('Expected: >10 fps for typical GPU\n');
```

**Step 2: Check GPU utilization**

```bash
# In separate terminal, monitor GPU
watch -n 1 nvidia-smi

# While running analysis, check:
# - GPU utilization (should be >50% during processing)
# - GPU memory usage
# - GPU temperature
# - Power draw
```

**Step 3: Check data transfer bottleneck**

```matlab
% Time data loading separately
tic;
[~, sequence, SMF] = LD.loadRawData(SMF, 1);
load_time = toc;
fprintf('Data loading time: %.2f seconds\n', load_time);

% If load_time > processing_time, disk I/O is bottleneck
```

**Step 4: Check for thermal throttling**

```bash
# Check GPU temperature
nvidia-smi --query-gpu=temperature.gpu --format=csv

# If temperature >80°C, GPU may be throttling
```

### Solutions

**Solution 5A: Improve cooling**

If GPU overheating (>85°C):

- Clean dust from GPU fans and heatsink
- Improve case airflow (add fans)
- Increase GPU fan speed (using nvidia-settings or GPU-Z)
- Consider aftermarket GPU cooler
- Move workstation to cooler room

**Solution 5B: Close background applications**

```bash
# Check what's using GPU
nvidia-smi

# Close:
# - Web browsers with hardware acceleration
# - Video players
# - Other MATLAB sessions
# - Machine learning frameworks
# - Mining software
# - Desktop effects (Linux: switch to lightweight WM)
```

**Solution 5C: Optimize data loading**

```matlab
% Use SSD instead of HDD for data files
% Preload data into RAM if possible
sequence = loadYourData();  % Load once

% Then process without reloading
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();
```

**Solution 5D: Process multiple datasets consecutively**

GPU initialization has overhead. Batch processing is more efficient:

```matlab
files = {'data1.h5', 'data2.h5', 'data3.h5'};

% GPU initialized once at start
for i = 1:length(files)
    SMF.Data.FileName = files(i);
    [~, sequence, SMF] = LD.loadRawData(SMF, 1);

    LD = smi_core.LocalizeData(sequence, SMF);
    SMD = LD.genLocalizations();

    % Save results
    save(sprintf('results_%d.mat', i), 'SMD');
end
```

**Solution 5E: Check PCIe slot configuration**

GPU may be in wrong PCIe slot:

```bash
# Linux - check PCIe bandwidth
lspci -vv | grep -A 30 VGA

# Look for "LnkSta: Speed 8GT/s, Width x16"
# Should be x16 for best performance
# If x8 or x4, GPU in suboptimal slot
```

Move GPU to PCIe x16 slot closest to CPU (usually top slot).

**Solution 5F: Verify GPU is actually being used**

```matlab
% Explicitly select GPU
gpuDevice(1);

% Monitor during processing
gpu = gpuDevice();
fprintf('Before: %.1f%% utilization\n', gpu.Utilization);

% Run processing
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();

gpu = gpuDevice();
fprintf('After: %.1f%% utilization\n', gpu.Utilization);

% If utilization stays at 0%, GPU not being used
% Check if PTX files loaded correctly
```

## Problem 6: Multiple GPU Issues

### Symptoms

```matlab
% Wrong GPU selected
% "Device 0 out of range" errors
% Inconsistent performance across GPUs
% Parallel pool doesn't distribute across GPUs
```

### Diagnosis

**Step 1: List available GPUs**

```matlab
% Count GPUs
N = gpuDeviceCount;
fprintf('Number of GPUs: %d\n', N);

% List each GPU
for i = 1:N
    gpu = gpuDevice(i);
    fprintf('\nGPU %d:\n', i);
    fprintf('  Name: %s\n', gpu.Name);
    fprintf('  Compute capability: %.1f\n', gpu.ComputeCapability);
    fprintf('  Total memory: %.1f GB\n', gpu.TotalMemory / 1e9);
    fprintf('  Available memory: %.1f GB\n', gpu.AvailableMemory / 1e9);
    fprintf('  Utilization: %.1f%%\n', gpu.Utilization);
end
```

**Step 2: Check which GPU is currently selected**

```matlab
% Current GPU
gpu = gpuDevice();
fprintf('Currently using GPU %d: %s\n', gpu.Index, gpu.Name);
```

**Step 3: Test each GPU individually**

```matlab
for i = 1:gpuDeviceCount
    % Select GPU
    gpu = gpuDevice(i);
    fprintf('Testing GPU %d: %s\n', i, gpu.Name);

    % Run simple test
    try
        A = gpuArray(rand(1000, 'single'));
        B = A * A';
        C = gather(B);
        fprintf('  Status: OK\n');
    catch ME
        fprintf('  Status: FAILED - %s\n', ME.message);
    end
end
```

### Solutions

**Solution 6A: Select specific GPU**

```matlab
% Explicitly select GPU before analysis
gpuDevice(2);  % Use GPU 2

% Verify selection
gpu = gpuDevice();
fprintf('Using: %s\n', gpu.Name);

% Run analysis
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();
```

**Solution 6B: Use fastest GPU**

```matlab
% Find GPU with most memory and highest compute capability
N = gpuDeviceCount;
best_gpu = 1;
best_score = 0;

for i = 1:N
    gpu = gpuDevice(i);
    score = gpu.ComputeCapability * gpu.TotalMemory;
    fprintf('GPU %d score: %.2e\n', i, score);

    if score > best_score
        best_score = score;
        best_gpu = i;
    end
end

fprintf('Selected GPU %d\n', best_gpu);
gpuDevice(best_gpu);
```

**Solution 6C: Distribute work across multiple GPUs**

```matlab
% Setup parallel pool with multiple GPUs
N_GPUs = gpuDeviceCount;
parpool('local', N_GPUs);

% Process multiple files in parallel
files = {'data1.h5', 'data2.h5', 'data3.h5', 'data4.h5'};

parfor i = 1:length(files)
    % Each worker uses different GPU automatically
    gpu = gpuDevice();
    fprintf('Worker %d using GPU %d: %s\n', ...
        i, gpu.Index, gpu.Name);

    % Load and process
    SMF_local = SMF;  % Copy for each worker
    SMF_local.Data.FileName = files(i);
    [~, sequence, ~] = smi_core.LoadData().loadRawData(SMF_local, 1);

    LD = smi_core.LocalizeData(sequence, SMF_local);
    SMD = LD.genLocalizations();

    % Save results
    save(sprintf('results_%d.mat', i), 'SMD');
end

delete(gcp('nocreate'));
```

**Solution 6D: Handle heterogeneous GPUs**

When GPUs have different capabilities:

```matlab
% Check each GPU and assign appropriate workload
for i = 1:gpuDeviceCount
    gpu = gpuDevice(i);

    if gpu.TotalMemory > 8e9  % >8 GB
        fprintf('GPU %d: Use for large datasets\n', i);
        % Process large files on this GPU
    else
        fprintf('GPU %d: Use for small datasets\n', i);
        % Process smaller files on this GPU
    end
end
```

**Solution 6E: Reset all GPUs**

If GPUs in inconsistent state:

```matlab
% Reset each GPU
for i = 1:gpuDeviceCount
    gpu = gpuDevice(i);
    fprintf('Resetting GPU %d: %s\n', i, gpu.Name);
    reset(gpu);
    fprintf('  Available memory: %.2f GB\n', gpu.AvailableMemory / 1e9);
end

% Select GPU 1 for next operation
gpuDevice(1);
```

## General Diagnostic Workflow

When encountering any GPU problem, follow this systematic approach:

### Step 1: Verify basics

```matlab
% Run comprehensive check
function checkGPU()
    fprintf('=== GPU Diagnostic ===\n\n');

    % 1. Check GPU count
    try
        N = gpuDeviceCount;
        fprintf('GPU count: %d\n', N);
    catch ME
        error('Cannot detect GPUs: %s', ME.message);
    end

    % 2. Check each GPU
    for i = 1:N
        fprintf('\n--- GPU %d ---\n', i);
        try
            gpu = gpuDevice(i);
            fprintf('Name: %s\n', gpu.Name);
            fprintf('Compute capability: %.1f\n', gpu.ComputeCapability);
            fprintf('Memory: %.1f GB (%.1f GB available)\n', ...
                gpu.TotalMemory/1e9, gpu.AvailableMemory/1e9);
            fprintf('Driver: %s\n', gpu.DriverVersion);
        catch ME
            fprintf('Error: %s\n', ME.message);
        end
    end

    % 3. Check PTX files
    fprintf('\n--- PTX Files ---\n');
    ptx_files = {'smi_cuda_FindROI.ptx', 'smi_cuda_gaussMLEv2.ptx'};
    for i = 1:length(ptx_files)
        path = which(ptx_files{i});
        if isempty(path)
            fprintf('%s: NOT FOUND\n', ptx_files{i});
        else
            fprintf('%s: OK\n', ptx_files{i});
        end
    end

    fprintf('\n=== Diagnostic Complete ===\n');
end

checkGPU();
```

### Step 2: Test with simple operation

```matlab
% Minimal GPU test
try
    gpu = gpuDevice();
    A = gpuArray(rand(100, 'single'));
    B = A * A';
    C = gather(B);
    fprintf('Basic GPU operations: OK\n');
catch ME
    fprintf('Basic GPU operations FAILED: %s\n', ME.message);
end
```

### Step 3: Test smite-specific operations

```matlab
% Test FindROI kernel
try
    % Small test image
    test_img = randn(128, 128, 10, 'single') * 10 + 100;

    SMF = smi_core.SingleMoleculeFitting();
    SMF.Data.CameraGain = 1;
    SMF.Data.CameraOffset = 100;

    LD = smi_core.LocalizeData(test_img, SMF);
    SMD = LD.genLocalizations();

    fprintf('smite GPU operations: OK\n');
    fprintf('Found %d localizations\n', length(SMD.X));
catch ME
    fprintf('smite GPU operations FAILED: %s\n', ME.message);
    fprintf('Check PTX files and GPU compatibility\n');
end
```

## See Also

- [How to Use GPU Acceleration](../how-to/use-gpu.md) - Comprehensive GPU usage guide
- [Installation Guide](../getting-started/installation.md) - Setup requirements
- [doc/mex+CUDA.md](../../mex+CUDA.md) - Compilation details
- [MATLAB GPU Support](https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html) - Compatibility matrix
- [NVIDIA CUDA Toolkit](https://developer.nvidia.com/cuda-downloads) - Download page
- [NVIDIA Drivers](https://www.nvidia.com/Download/index.aspx) - Driver downloads
