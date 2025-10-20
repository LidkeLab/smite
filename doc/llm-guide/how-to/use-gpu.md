---
title: "How to Use GPU Acceleration"
category: "how-to"
level: "intermediate"
tags: ["gpu", "cuda", "performance", "optimization", "hardware"]
prerequisites: ["../core-concepts/smf-structure.md"]
related: ["localize-molecules.md", "../troubleshooting/common-errors.md"]
summary: "Guide to GPU acceleration in smite for high-performance molecule localization"
estimated_time: "10 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# How to Use GPU Acceleration

## Purpose

smite uses NVIDIA GPU acceleration for computationally intensive operations, particularly molecule localization via GaussMLE fitting and box finding. GPU acceleration provides 10-100x speedup over CPU-only implementations, enabling real-time processing of large SMLM datasets. This guide explains GPU requirements, how to verify GPU availability, optimize GPU performance, troubleshoot GPU issues, and work with multiple GPUs.

## Prerequisites

- Understanding of [SMF structure](../core-concepts/smf-structure.md)
- NVIDIA GPU with CUDA support
- MATLAB Parallel Computing Toolbox

## Overview

GPU acceleration in smite is **required** for core localization operations, not optional. The following operations use GPU:

1. **GaussMLE fitting**: 2D Gaussian PSF fitting (all fit types: XYNB, XYNBS, XYNBSXSY, XYZNB)
2. **Box finding**: Gaussian filtering and local maxima detection
3. **PSF calculations**: Point spread function generation (optional)

Without a compatible GPU, localization will fail. CPU-only operation is not supported.

## GPU Requirements

### Hardware Requirements

**Minimum requirements:**

- NVIDIA GPU (AMD, Intel, Apple Silicon GPUs not supported)
- CUDA compute capability 5.0 or higher
- 2 GB GPU memory (minimum)
- 4+ GB GPU memory (recommended for typical datasets)
- 8+ GB GPU memory (recommended for large datasets)

**Compute capability reference:**

| Compute Capability | GPU Architecture | Example GPUs |
|-------------------|------------------|--------------|
| 5.0 - 5.3 | Maxwell | GTX 750 Ti, GTX 980, Quadro M4000 |
| 6.0 - 6.2 | Pascal | GTX 1080, Tesla P100, Quadro P6000 |
| 7.0 - 7.5 | Volta, Turing | RTX 2080, Tesla V100, Quadro RTX 5000 |
| 8.0 - 8.9 | Ampere | RTX 3090, A100, RTX A6000 |
| 9.0+ | Hopper | H100 |

**Note:** Compute capability must be supported by your MATLAB version. See [MATLAB GPU support by release](https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html).

### Software Requirements

**Required:**

- MATLAB R2021a or later
- MATLAB Parallel Computing Toolbox
- NVIDIA GPU driver (latest recommended)

**For CUDA compilation (only if modifying source):**

- NVIDIA CUDA Toolkit (nvcc compiler)
- Windows only: Visual Studio C++ compiler (cl.exe)

**Note:** Precompiled CUDA files (`.ptx`) are included for Linux, MacOS, and Windows, so CUDA Toolkit installation is not needed unless you modify CUDA source code.

## Checking GPU Availability

### Basic GPU Check

Verify GPU is detected and meets requirements:

```matlab
% Check if GPU is available
try
    gpu = gpuDevice();
    fprintf('GPU detected: %s\n', gpu.Name);
    fprintf('Compute capability: %.1f\n', gpu.ComputeCapability);
    fprintf('Total memory: %.2f GB\n', gpu.TotalMemory / 1e9);
    fprintf('Available memory: %.2f GB\n', gpu.AvailableMemory / 1e9);
    fprintf('Driver version: %s\n', gpu.DriverVersion);
    fprintf('Toolkit version: %s\n', gpu.ToolkitVersion);
catch ME
    error('No compatible GPU found: %s', ME.message);
end
```

### Verify Compute Capability

```matlab
gpu = gpuDevice();
if gpu.ComputeCapability < 5.0
    error('GPU compute capability is %.1f, but smite requires >= 5.0', ...
        gpu.ComputeCapability);
else
    fprintf('GPU compute capability %.1f meets requirements (>= 5.0)\n', ...
        gpu.ComputeCapability);
end
```

### Check MATLAB Compatibility

```matlab
% Verify Parallel Computing Toolbox
if license('test', 'Distrib_Computing_Toolbox')
    fprintf('Parallel Computing Toolbox: Available\n');
else
    error('Parallel Computing Toolbox required but not available');
end

% Check CUDA support
if parallel.gpu.GPUDevice.isAvailable
    fprintf('CUDA support: Available\n');
else
    error('CUDA support not available');
end
```

### Full System Check

```matlab
function checkGPUForSMITE()
    % Comprehensive GPU check for smite

    fprintf('=== GPU System Check for smite ===\n\n');

    % 1. Check GPU availability
    try
        gpu = gpuDevice();
    catch ME
        error('GPU not detected: %s\nsmite requires NVIDIA GPU with CUDA support.', ME.message);
    end

    % 2. Display GPU info
    fprintf('GPU: %s\n', gpu.Name);
    fprintf('Compute Capability: %.1f\n', gpu.ComputeCapability);
    fprintf('Total Memory: %.2f GB\n', gpu.TotalMemory / 1e9);
    fprintf('Available Memory: %.2f GB\n', gpu.AvailableMemory / 1e9);

    % 3. Check compute capability
    if gpu.ComputeCapability < 5.0
        error('Compute capability %.1f is too low. smite requires >= 5.0', ...
            gpu.ComputeCapability);
    else
        fprintf('Status: PASS (compute capability >= 5.0)\n');
    end

    % 4. Check memory
    if gpu.TotalMemory < 2e9
        warning('GPU memory (%.2f GB) is low. May have issues with large datasets.', ...
            gpu.TotalMemory / 1e9);
    else
        fprintf('Status: PASS (sufficient memory)\n');
    end

    % 5. Check CUDA kernel files
    fprintf('\n=== Checking CUDA kernel files ===\n');

    kernels = {
        'smi_cuda_FindROI.ptx'
        'smi_cuda_gaussMLEv2.ptx'
    };

    for i = 1:length(kernels)
        if exist(kernels{i}, 'file')
            fprintf('%s: FOUND\n', kernels{i});
        else
            error('%s: NOT FOUND\nReinstall smite or recompile CUDA files.', kernels{i});
        end
    end

    % 6. Test basic GPU operation
    fprintf('\n=== Testing GPU operations ===\n');
    try
        A = gpuArray(rand(100, 100, 'single'));
        B = A * A';
        C = gather(B);
        fprintf('GPU array operations: PASS\n');
    catch ME
        error('GPU operations failed: %s', ME.message);
    end

    fprintf('\n=== GPU System Check Complete ===\n');
    fprintf('Your system is ready to use smite with GPU acceleration.\n');
end
```

Run this check before starting analysis:

```matlab
checkGPUForSMITE();
```

## GPU Memory Management

### Understanding GPU Memory Usage

smite automatically manages GPU memory by chunking data to fit available memory:

**Box finding (FindROI):**

```matlab
% Memory usage calculation (from FindROI.m)
NCopies = 8;  % For out-of-place operations
NBytesPerPixel = 4;  % single float
NElem = numel(Data);
gpu = gpuDevice;
NLoops = ceil(NBytesPerPixel * NCopies * NElem / gpu.AvailableMemory);
```

**GaussMLE fitting:**

```matlab
% Memory usage calculation (from GaussMLE.m)
BytesPerFloat = 4;
MemoryPerFit = BytesPerFloat * (BoxSize^2 + NParams + NParams + 1);
MaxFits = floor(gpu.AvailableMemory / MemoryPerFit / 2);  % Factor of 2 for safety
```

### Monitoring GPU Memory

Check memory during analysis:

```matlab
% Before starting
gpu = gpuDevice();
fprintf('Available memory before: %.2f GB\n', gpu.AvailableMemory / 1e9);

% Run localization
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();

% After completion
gpu = gpuDevice();
fprintf('Available memory after: %.2f GB\n', gpu.AvailableMemory / 1e9);
```

### Clearing GPU Memory

If encountering memory issues:

```matlab
% Clear GPU memory
gpu = gpuDevice();
reset(gpu);

% Or wait for automatic cleanup
wait(gpu);
```

### Optimizing Memory Usage

**For large datasets:**

```matlab
% Process in batches
NFrames = size(sequence, 3);
BatchSize = 1000;  % Adjust based on available memory
NBatches = ceil(NFrames / BatchSize);

SMD_all = smi_core.SingleMoleculeData.createSMD();

for batch = 1:NBatches
    start_frame = (batch - 1) * BatchSize + 1;
    end_frame = min(batch * BatchSize, NFrames);

    % Process batch
    batch_data = sequence(:, :, start_frame:end_frame);
    LD = smi_core.LocalizeData(batch_data, SMF);
    SMD_batch = LD.genLocalizations();

    % Adjust frame numbers
    SMD_batch.FrameNum = SMD_batch.FrameNum + start_frame - 1;

    % Combine
    SMD_all = smi_core.SingleMoleculeData.catSMD(SMD_all, SMD_batch);

    % Clear GPU memory
    reset(gpuDevice());

    fprintf('Processed batch %d/%d\n', batch, NBatches);
end
```

**Memory-constrained systems:**

```matlab
% Reduce box size (if appropriate for PSF)
SMF.BoxFinding.BoxSize = 5;  % Instead of 7

% Process fewer frames at once
% smite automatically chunks data, but you can pre-chunk
```

## GPU Performance Optimization

### CPU vs GPU Performance Comparison

Typical performance gains with GPU acceleration:

| Operation | CPU Time | GPU Time | Speedup |
|-----------|----------|----------|---------|
| Box finding (1000 frames, 256x256) | ~30 sec | ~0.5 sec | 60x |
| GaussMLE fitting (10,000 boxes) | ~120 sec | ~2 sec | 60x |
| Complete SMLM pipeline | ~10 min | ~15 sec | 40x |

**Benchmark your system:**

```matlab
% Create test data
NFrames = 100;
ImageSize = 256;
sequence = randn(ImageSize, ImageSize, NFrames, 'single') * 10 + 100;

% Add synthetic emitters
for frame = 1:NFrames
    for emitter = 1:50
        x = randi([20, ImageSize-20]);
        y = randi([20, ImageSize-20]);
        amplitude = 1000;
        sigma = 1.3;
        [X, Y] = meshgrid(1:ImageSize, 1:ImageSize);
        sequence(:, :, frame) = sequence(:, :, frame) + ...
            amplitude * exp(-((X - x).^2 + (Y - y).^2) / (2 * sigma^2));
    end
end

% Setup SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraGain = 1;
SMF.Data.CameraOffset = 100;
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 200;
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Benchmark
tic;
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();
gpu_time = toc;

fprintf('GPU processing time: %.2f seconds\n', gpu_time);
fprintf('Frames per second: %.1f\n', NFrames / gpu_time);
fprintf('Localizations: %d\n', length(SMD.X));
```

### Maximizing GPU Utilization

**1. Minimize data transfers between CPU and GPU**

smite handles this automatically, but if writing custom code:

```matlab
% Good: Keep data on GPU
gpu_data = gpuArray(sequence);
% ... multiple GPU operations ...
result = gather(gpu_data);  % Transfer once at end

% Bad: Repeated transfers
for i = 1:N
    gpu_data = gpuArray(data_i);  % Transfer to GPU
    result_i = gather(process(gpu_data));  % Transfer to CPU
end
```

**2. Process multiple datasets consecutively**

GPU initialization has overhead. Process multiple files in one session:

```matlab
files = {'data1.h5', 'data2.h5', 'data3.h5'};

% GPU initialized once
for i = 1:length(files)
    SMF.Data.FileName = files(i);
    [~, sequence, SMF] = LD_loader.loadRawData(SMF, 1);
    LD = smi_core.LocalizeData(sequence, SMF);
    SMD = LD.genLocalizations();
    % Save results...
end
```

**3. Use appropriate data types**

smite uses `single` precision (float32) for GPU operations, which is optimal for GPU performance and memory usage.

## Multiple GPU Support

### Detecting Multiple GPUs

```matlab
% List all GPUs
N = gpuDeviceCount;
fprintf('Number of GPUs: %d\n', N);

for i = 1:N
    gpu = gpuDevice(i);
    fprintf('GPU %d: %s (%.1f GB, compute capability %.1f)\n', ...
        i, gpu.Name, gpu.TotalMemory / 1e9, gpu.ComputeCapability);
end
```

### Selecting a GPU

```matlab
% Select specific GPU
gpu = gpuDevice(2);  % Use GPU 2
fprintf('Selected GPU: %s\n', gpu.Name);

% Run analysis on selected GPU
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();
```

### Parallel Processing with Multiple GPUs

Process different datasets on different GPUs using parallel pool:

```matlab
% Setup parallel pool with multiple GPUs
N_GPUs = gpuDeviceCount;
parpool('local', N_GPUs);

% Process multiple files in parallel
files = {'data1.h5', 'data2.h5', 'data3.h5', 'data4.h5'};

parfor i = 1:length(files)
    % Each worker automatically uses different GPU
    gpu = gpuDevice();
    fprintf('Worker %d using GPU: %s\n', ...
        get(getCurrentTask(), 'ID'), gpu.Name);

    % Load and process
    SMF_local = SMF;  % Copy SMF for each worker
    SMF_local.Data.FileName = files(i);
    [~, sequence, SMF_local] = smi_core.LoadData().loadRawData(SMF_local, 1);

    LD = smi_core.LocalizeData(sequence, SMF_local);
    SMD = LD.genLocalizations();

    % Save results
    save(sprintf('results_%d.mat', i), 'SMD', 'SMF_local');
end

delete(gcp('nocreate'));
```

## Troubleshooting GPU Issues

### GPU Not Detected

**Problem:** `gpuDevice()` fails or returns error

**Solutions:**

```matlab
% Check if GPU is available
if ~parallel.gpu.GPUDevice.isAvailable
    % Possible causes:
    % 1. No NVIDIA GPU installed
    % 2. GPU driver not installed
    % 3. GPU driver too old
    % 4. MATLAB cannot access GPU

    % Check system:
    % - Device Manager (Windows): Display adapters
    % - nvidia-smi (Linux/Windows command line)
    % - System Preferences > Displays (Mac)
end

% Update GPU driver
% Download from: https://www.nvidia.com/Download/index.aspx

% Restart MATLAB after driver update
```

### Compute Capability Too Low

**Problem:** GPU compute capability < 5.0

**Solution:**

```matlab
gpu = gpuDevice();
if gpu.ComputeCapability < 5.0
    fprintf('GPU compute capability: %.1f\n', gpu.ComputeCapability);
    fprintf('smite requires: >= 5.0\n');
    fprintf('Upgrade to newer GPU or use different machine.\n');
end
```

Compatible GPUs: Maxwell architecture (2014) or newer.

### Out of Memory Errors

**Problem:** `Out of memory on device` or `CUDA_ERROR_OUT_OF_MEMORY`

**Solutions:**

```matlab
% 1. Check available memory
gpu = gpuDevice();
fprintf('Available: %.2f GB\n', gpu.AvailableMemory / 1e9);

% 2. Clear GPU memory
reset(gpu);

% 3. Process smaller batches
% (see GPU Memory Management section above)

% 4. Reduce box size if appropriate
SMF.BoxFinding.BoxSize = 5;  % Instead of 7

% 5. Close other GPU applications
% - Check with nvidia-smi (command line)
% - Close other MATLAB sessions using GPU
% - Close deep learning applications, games, etc.
```

### CUDA Kernel Not Found

**Problem:** `Unable to find PTX file` or kernel loading fails

**Solutions:**

```matlab
% Check if PTX files exist
required_files = {
    'smi_cuda_FindROI.ptx'
    'smi_cuda_gaussMLEv2.ptx'
};

for i = 1:length(required_files)
    if exist(required_files{i}, 'file')
        fprintf('%s: Found at %s\n', required_files{i}, which(required_files{i}));
    else
        fprintf('%s: NOT FOUND\n', required_files{i});
    end
end

% If not found:
% 1. Check setupSMITE.m was run
% 2. Check ptx directory is on path
% 3. Recompile CUDA files if needed (see Compilation section)
```

### GPU Performance Degraded

**Problem:** GPU processing slower than expected

**Diagnostics:**

```matlab
% Check if GPU is being used
gpu = gpuDevice();
fprintf('GPU utilization before: %.1f%%\n', gpu.Utilization);

% Monitor during processing
tic;
LD = smi_core.LocalizeData(sequence, SMF);
SMD = LD.genLocalizations();
elapsed = toc;

fprintf('Processing time: %.2f seconds\n', elapsed);
fprintf('GPU utilization after: %.1f%%\n', gpu.Utilization);

% If utilization is low (<50%), possible causes:
% - Data transfer bottleneck (slow disk I/O)
% - Small dataset (GPU initialization overhead dominates)
% - GPU throttling due to temperature
% - Background processes using GPU
```

**Solutions:**

```matlab
% 1. Check GPU temperature (use nvidia-smi or GPU-Z)
% If overheating, improve cooling

% 2. Close background GPU applications

% 3. Check data loading speed
tic;
[~, sequence, SMF] = LD_loader.loadRawData(SMF, 1);
load_time = toc;
fprintf('Data loading time: %.2f seconds\n', load_time);

% 4. Use SSD for data storage (faster than HDD)
```

## Compiling CUDA Files (Advanced)

Compilation is only needed if:

1. Modifying CUDA source code (`.cu` files)
2. Precompiled PTX files don't work on your system
3. Targeting specific GPU architecture for optimization

### Prerequisites

**All platforms:**

- NVIDIA CUDA Toolkit installed
- `nvcc` compiler on system PATH

**Windows only:**

- Visual Studio with C++ compiler
- `cl.exe` compiler on system PATH

### Compilation Steps

```matlab
% Navigate to CUDA source directory
cd ~/Documents/MATLAB/smite/MATLAB/source/cuda

% Run compilation script
cuda_Make

% Output goes to ~/Documents/MATLAB/smite/MATLAB/ptx/
```

### Troubleshooting Compilation

**nvcc not found:**

Edit `cuda_Make.m` to add nvcc to PATH:

```matlab
% MacOS/Linux
setenv('PATH', ['/usr/local/cuda/bin:' getenv('PATH')]);

% Windows
setenv('PATH', ['C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.8\bin;' getenv('PATH')]);
```

**cl.exe not found (Windows):**

Add Visual Studio compiler to PATH:

```matlab
setenv('PATH', ['C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.35.32215\bin\Hostx64\x64;' getenv('PATH')]);
```

Adjust path based on your Visual Studio installation.

## See Also

- [Localize Molecules](localize-molecules.md) - GPU is used during localization
- [SMF Structure](../core-concepts/smf-structure.md) - All fitting parameters
- [Troubleshooting](../troubleshooting/common-errors.md) - Error solutions
- [doc/mex+CUDA.md](../../mex+CUDA.md) - Detailed compilation instructions
- [MATLAB GPU Support](https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html) - Version compatibility
