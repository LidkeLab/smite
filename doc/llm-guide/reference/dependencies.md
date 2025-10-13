---
title: "Software Dependencies Reference"
category: "reference"
level: "beginner"
tags: ["dependencies", "requirements", "matlab", "toolboxes", "cuda", "gpu", "ffmpeg"]
prerequisites: []
related: ["../getting-started/installation.md", "../troubleshooting/installation-issues.md", "../troubleshooting/gpu-problems.md"]
summary: "Complete reference for all software dependencies required and optional for smite"
estimated_time: "10-15 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Software Dependencies Reference

## Purpose

This document provides a comprehensive reference for all software dependencies required to run smite. Understanding these dependencies helps you:
- Determine if your system is compatible
- Identify which features require specific toolboxes
- Check what's already installed
- Troubleshoot missing dependency errors
- Plan installations for new systems

## Dependency Overview

smite has three tiers of dependencies:

1. **Core Requirements** - Essential for basic functionality
2. **Required for Full Functionality** - Needed for standard workflows
3. **Optional** - Enable specific features or performance enhancements

## Core Requirements

These are absolutely essential to run smite at all.

### Operating System

**Supported:**
- Linux (64-bit)
- MacOS (64-bit, Intel or Apple Silicon via Rosetta 2)
- Windows (64-bit)

**Version Requirements:**
- Linux: Any modern distribution (Ubuntu 18.04+, CentOS 7+, etc.)
- MacOS: 10.14 (Mojave) or later
- Windows: Windows 10 or later

**Why:** smite includes platform-specific precompiled binary files (mex and ptx) for these systems.

### MATLAB

**Minimum Version:** R2021a

**Recommended:** R2023a or later

**Why R2021a minimum:**
- Modern GPU support APIs
- Improved class and object handling
- Better parallel computing features
- Function argument validation support

**Check your version:**
```matlab
version
```

**Output format:**
```
9.10.0.1602886 (R2021a)
```

The version number format is: `<major>.<minor>.<patch>.<build> (R<year><release>)`

### MATLAB Toolboxes (Required)

Three toolboxes are required for core smite functionality:

#### 1. Image Processing Toolbox

**Why required:**
- Image filtering and preprocessing
- Morphological operations
- Region of interest (ROI) operations
- Image registration functions

**Used by:**
- All image loading and processing
- ROI analysis tools
- Drift correction algorithms
- Image visualization

#### 2. Parallel Computing Toolbox

**Why required:**
- GPU computing for localization (CUDA operations)
- Parallel processing for batch operations
- Memory management for large datasets
- GPUArray data structures

**Used by:**
- `smi_core.GaussMLE` (GPU-accelerated fitting)
- All CUDA-based localization
- Batch processing workflows
- GPU device management

**Critical note:** This toolbox is REQUIRED even if you don't have a GPU, because smite's core fitting algorithms rely on CUDA. GPU acceleration isn't optional - it's fundamental to how smite performs localization.

#### 3. Statistics and Machine Learning Toolbox

**Why required:**
- Statistical analysis of results
- Distribution fitting
- Hypothesis testing
- Clustering algorithms
- Random number generation with specific distributions

**Used by:**
- Frame connection algorithms
- Drift correction methods
- Clustering analysis
- Statistical result validation
- Simulation tools

## Optional Toolboxes

These toolboxes enable specific features but aren't required for basic SMLM workflows.

### Curve Fitting Toolbox

**Enables:**
- Advanced PSF fitting models
- Fourier Ring Correlation (FRC) resolution estimation
- Diffusion coefficient estimation
- Custom fitting functions

**Specific modules requiring it:**
- `smi_cluster.*` (some clustering algorithms)
- `smi_core.FRC` (resolution estimation)
- `smi_stat.DiffusionEstimator`

**Check if installed:**
```matlab
license('test', 'Curve_Fitting_Toolbox')
```

Returns `1` if available, `0` if not.

### Optimization Toolbox

**Enables:**
- Pair correlation analysis
- Advanced diffusion analysis
- Parameter optimization in fitting
- Constrained optimization problems

**Specific modules requiring it:**
- `smi_cluster.PairCorrelation`
- `smi_stat.DiffusionEstimator`

**Check if installed:**
```matlab
license('test', 'Optimization_Toolbox')
```

### Signal Processing Toolbox

**Enables:**
- Fourier Ring Correlation (FRC)
- Advanced filtering operations
- Frequency domain analysis

**Specific modules requiring it:**
- `smi_core.FRC`

**Check if installed:**
```matlab
license('test', 'Signal_Processing_Toolbox')
```

## GPU Requirements

### NVIDIA GPU

**Required Specifications:**
- NVIDIA GPU (AMD/Intel GPUs not supported)
- CUDA Compute Capability ≥ 5.0
- Minimum 2 GB VRAM (4+ GB recommended)
- Supported by your MATLAB version

**Why NVIDIA only:**
MATLAB's Parallel Computing Toolbox only supports CUDA (NVIDIA's GPU computing platform). Other GPU types (AMD, Intel, Apple Silicon) cannot be used for GPU acceleration.

**Compute Capability explained:**
Compute Capability is NVIDIA's versioning system for GPU architectures. Version 5.0+ includes:
- Maxwell architecture (GTX 900 series, GTX 10 series)
- Pascal architecture (GTX 10 series, Quadro P series)
- Volta architecture (Tesla V100)
- Turing architecture (RTX 20 series, Quadro RTX)
- Ampere architecture (RTX 30 series, A series)
- Ada Lovelace architecture (RTX 40 series)
- Hopper architecture (H100)

**Check your GPU:**
```matlab
gpuDevice
```

**Example output:**
```
ans =

  CUDADevice with properties:

                      Name: 'NVIDIA GeForce RTX 3080'
                     Index: 1
         ComputeCapability: '8.6'
            SupportsDouble: 1
     GraphicsDriverVersion: '535.183'
               DriverModel: 'WDDM'
            ToolkitVersion: 11.8000
        MaxThreadsPerBlock: 1024
          MaxShmemPerBlock: 49152 (49.15 KB)
        MaxThreadBlockSize: [1024 1024 64]
               MaxGridSize: [2.1475e+09 65535 65535]
                 SIMDWidth: 32
               TotalMemory: 10200186880 (1.02e+10 bytes)
           AvailableMemory: 8434274304 (8.43e+09 bytes)
       MultiprocessorCount: 68
              ClockRateKHz: 1710000
               ComputeMode: Default
      GPUOverlapsTransfers: 1
    KernelExecutionTimeout: 1
          CanMapHostMemory: 1
           DeviceSupported: 1
           DeviceAvailable: 1
            DeviceSelected: 1
```

**Key fields to check:**
- `ComputeCapability`: Must be ≥ '5.0'
- `DeviceSupported`: Must be 1
- `DeviceAvailable`: Must be 1

**If no GPU available:**
```
Error using gpuDevice
No supported GPU device found.
```

This means either:
- No GPU is installed
- GPU drivers not installed/updated
- GPU not NVIDIA
- GPU compute capability too old
- Parallel Computing Toolbox not installed

### CUDA Compatibility

MATLAB versions support different CUDA versions and GPU architectures:

**Version compatibility:**
- R2021a: CUDA 11.0-11.2, Compute Capability 3.5-8.6
- R2021b: CUDA 11.0-11.4, Compute Capability 3.5-8.6
- R2022a: CUDA 11.0-11.6, Compute Capability 3.7-8.6
- R2022b: CUDA 11.0-11.8, Compute Capability 3.7-8.6
- R2023a: CUDA 11.0-12.0, Compute Capability 3.7-9.0
- R2023b: CUDA 11.2-12.2, Compute Capability 3.7-9.0
- R2024a: CUDA 11.2-12.4, Compute Capability 5.0-9.0

**Reference:** [MathWorks GPU Support by Release](https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html)

**Why this matters:**
If your GPU is newer than what your MATLAB version supports, CUDA compilation may fail. Solution: Update MATLAB or use precompiled ptx files.

## External Software

### ffmpeg (Linux Only)

**Required for:** Video generation during localization (when `LocalizeData.Verbose >= 3`)

**Platform:** Linux only (not needed for Windows/MacOS)

**Install:**
```bash
# Ubuntu/Debian
sudo apt-get install ffmpeg

# CentOS/RHEL
sudo yum install ffmpeg

# Arch
sudo pacman -S ffmpeg
```

**Verify installation:**
```bash
ffmpeg -version
```

**Why only Linux:**
Video generation during localization is only implemented for Linux. Windows and MacOS users can still generate videos through other smite visualization tools, but the real-time overlay videos during `LocalizeData.genLocalizations()` won't work.

**Alternative without ffmpeg:**
You can still use smite on Linux without ffmpeg - just keep `LocalizeData.Verbose < 3` to skip video generation.

## Compilation Dependencies

Most users don't need these - precompiled binaries are included. Only required if:
- Precompiled files don't work
- You modified source code
- You're developing new features

### For Mex Files

**C/C++ Compiler:**
- **Windows:** MinGW-w64 (via MATLAB Add-On) or Visual Studio
- **Linux:** GCC (usually pre-installed)
- **MacOS:** Xcode Command Line Tools

**Setup:**
```matlab
mex -setup
```

**Should output:**
```
MEX configured to use 'gcc' for C language compilation.
```

**Install compiler if needed:**

**Windows (MinGW-w64):**
```matlab
% In MATLAB command window
% Go to: Add-Ons → Get Add-Ons → Search "MinGW"
% Install "MATLAB Support for MinGW-w64 C/C++ Compiler"
```

**Linux (GCC):**
```bash
sudo apt-get install build-essential  # Ubuntu/Debian
sudo yum groupinstall "Development Tools"  # CentOS/RHEL
```

**MacOS (Xcode):**
```bash
xcode-select --install
```

### For CUDA Files

**NVIDIA CUDA Toolkit:**
- Download from: https://developer.nvidia.com/cuda-toolkit
- Version must be compatible with your MATLAB version (see CUDA Compatibility above)
- `nvcc` compiler must be on system PATH

**Windows also requires:**
- Visual Studio (Community Edition is free)
- `cl.exe` compiler on system PATH

**Verify installation:**
```bash
nvcc --version
```

**Should output:**
```
nvcc: NVIDIA (R) Cuda compiler driver
Copyright (c) 2005-2023 NVIDIA Corporation
Built on Wed_Feb__8_05:53:42_Coordinated_Universal_Time_2023
Cuda compilation tools, release 12.1, V12.1.66
Build cuda_12.1.r12.1/compiler.32415258_0
```

## Checking What's Installed

### Check MATLAB and Toolboxes

**List all installed products:**
```matlab
ver
```

**Example output:**
```
MATLAB Version: 9.14.0.2206163 (R2023a)
MATLAB                                                Version 9.14        (R2023a)
Curve Fitting Toolbox                                 Version 3.9         (R2023a)
Image Processing Toolbox                              Version 11.7        (R2023a)
Optimization Toolbox                                  Version 9.5         (R2023a)
Parallel Computing Toolbox                            Version 7.8         (R2023a)
Signal Processing Toolbox                             Version 9.2         (R2023a)
Statistics and Machine Learning Toolbox               Version 12.5        (R2023a)
```

**Check specific toolbox:**
```matlab
% Returns 1 if installed, 0 if not
license('test', 'Image_Toolbox')
license('test', 'Distrib_Computing_Toolbox')  % Parallel Computing
license('test', 'Statistics_Toolbox')
license('test', 'Curve_Fitting_Toolbox')
license('test', 'Optimization_Toolbox')
license('test', 'Signal_Toolbox')
```

**Check MATLAB version programmatically:**
```matlab
% Get version info
v = ver('MATLAB');
fprintf('MATLAB %s (%s)\n', v.Version, v.Release);

% Check if version >= R2021a
if verLessThan('matlab', '9.10')
    warning('MATLAB R2021a or later required');
end
```

### Check GPU

**Basic check:**
```matlab
gpuDevice
```

**Check if GPU operations work:**
```matlab
try
    A = gpuArray(rand(100));
    B = A * A;
    fprintf('GPU operations working\n');
catch ME
    fprintf('GPU error: %s\n', ME.message);
end
```

**Get detailed GPU info:**
```matlab
gpu = gpuDevice;
fprintf('GPU: %s\n', gpu.Name);
fprintf('Compute Capability: %s\n', gpu.ComputeCapability);
fprintf('Memory: %.2f GB total, %.2f GB available\n', ...
    gpu.TotalMemory/1e9, gpu.AvailableMemory/1e9);
fprintf('Driver: %s\n', gpu.GraphicsDriverVersion);
fprintf('Toolkit: %.1f\n', gpu.ToolkitVersion);
```

### Check External Software

**ffmpeg (Linux):**
```bash
which ffmpeg
ffmpeg -version
```

**CUDA toolkit:**
```bash
which nvcc
nvcc --version
```

**Check CUDA driver:**
```bash
nvidia-smi
```

## Comprehensive Dependency Check Script

Save this as `checkSMITEDeps.m`:

```matlab
function checkSMITEDeps()
%checkSMITEDeps Check all smite dependencies
%
% Runs through all dependencies and reports status

    fprintf('=== SMITE Dependency Check ===\n\n');

    % MATLAB version
    fprintf('--- MATLAB ---\n');
    v = ver('MATLAB');
    fprintf('Version: %s (%s)\n', v.Version, v.Release);
    if verLessThan('matlab', '9.10')
        fprintf('  ⚠ WARNING: R2021a or later recommended\n');
    else
        fprintf('  ✓ OK\n');
    end
    fprintf('\n');

    % Required toolboxes
    fprintf('--- Required Toolboxes ---\n');
    checkToolbox('Image Processing', 'Image_Toolbox', true);
    checkToolbox('Parallel Computing', 'Distrib_Computing_Toolbox', true);
    checkToolbox('Statistics and Machine Learning', 'Statistics_Toolbox', true);
    fprintf('\n');

    % Optional toolboxes
    fprintf('--- Optional Toolboxes ---\n');
    checkToolbox('Curve Fitting', 'Curve_Fitting_Toolbox', false);
    checkToolbox('Optimization', 'Optimization_Toolbox', false);
    checkToolbox('Signal Processing', 'Signal_Toolbox', false);
    fprintf('\n');

    % GPU
    fprintf('--- GPU ---\n');
    try
        gpu = gpuDevice;
        fprintf('Name: %s\n', gpu.Name);
        fprintf('Compute Capability: %s\n', gpu.ComputeCapability);
        cap = str2double(gpu.ComputeCapability);
        if cap >= 5.0
            fprintf('  ✓ Compatible\n');
        else
            fprintf('  ⚠ WARNING: Compute Capability < 5.0\n');
        end
        fprintf('Memory: %.2f GB\n', gpu.TotalMemory/1e9);
        fprintf('Driver: %s\n', gpu.GraphicsDriverVersion);
    catch ME
        fprintf('  ✗ ERROR: %s\n', ME.message);
    end
    fprintf('\n');

    % smite installed
    fprintf('--- SMITE ---\n');
    if exist('setupSMITE', 'file')
        fprintf('Found setupSMITE: %s\n', which('setupSMITE'));
        try
            ver_str = smi_helpers.versionSMITE();
            fprintf('Version: %s\n', ver_str);
            fprintf('  ✓ Installed\n');
        catch
            fprintf('  ⚠ WARNING: setupSMITE not run yet\n');
        end
    else
        fprintf('  ✗ ERROR: setupSMITE not found\n');
        fprintf('  Add smite/MATLAB to path and run setupSMITE\n');
    end
    fprintf('\n');

    fprintf('=== Check Complete ===\n');
end

function checkToolbox(name, license_name, required)
    if license('test', license_name)
        fprintf('%s: ✓ Installed\n', name);
    else
        if required
            fprintf('%s: ✗ MISSING (REQUIRED)\n', name);
        else
            fprintf('%s: ✗ Not installed (optional)\n', name);
        end
    end
end
```

**Usage:**
```matlab
checkSMITEDeps()
```

## Dependency Resolution

### Missing Required Toolbox

**Symptom:**
```
Error: Undefined function or variable 'imfilter'
```

**Diagnosis:**
```matlab
license('test', 'Image_Toolbox')  % Returns 0
```

**Solution:**
Install missing toolbox via MATLAB Add-Ons or license administrator.

### GPU Not Working

**Symptom:**
```
Error: No supported GPU device found
```

**Diagnosis steps:**
1. Check if GPU exists: `gpuDevice`
2. Check driver: Run `nvidia-smi` in terminal
3. Check Parallel Computing Toolbox: `license('test', 'Distrib_Computing_Toolbox')`

**Solutions:**
- Install/update NVIDIA drivers
- Install Parallel Computing Toolbox
- Verify GPU is NVIDIA (not AMD/Intel)
- Check compute capability ≥ 5.0

### Compilation Errors

**Symptom:**
```
Error: mex compiler not configured
```

**Solution:**
```matlab
mex -setup
```

Install compiler if needed (see Compilation Dependencies above).

### Wrong MATLAB Version

**Symptom:**
```
Undefined function or variable 'arguments'
```

**Diagnosis:**
```matlab
version  % Check if < R2021a
```

**Solution:**
Update MATLAB to R2021a or later, or use older smite release compatible with your version.

## Platform-Specific Notes

### Windows

**Typical configuration:**
- Windows 10/11
- MATLAB R2023a+
- NVIDIA GPU with latest drivers
- No additional software needed

**Path separator:** Backslash (`\`)

### Linux

**Typical configuration:**
- Ubuntu 20.04/22.04
- MATLAB R2023a+
- NVIDIA GPU with proprietary drivers
- ffmpeg installed

**Install drivers:**
```bash
# Ubuntu
sudo ubuntu-drivers autoinstall
```

**Path separator:** Forward slash (`/`)

### MacOS

**Intel Macs:**
- Works normally if NVIDIA GPU present (rare)
- Most Intel Macs have AMD GPUs (won't work)

**Apple Silicon (M1/M2/M3/M4):**
- MATLAB must run via Rosetta 2
- No GPU acceleration (no NVIDIA GPU)
- CPU processing only
- All non-GPU features work fine

**Path separator:** Forward slash (`/`)

## Summary Matrix

| Component | Windows | Linux | MacOS Intel | MacOS Silicon |
|-----------|---------|-------|-------------|---------------|
| MATLAB R2021a+ | ✓ | ✓ | ✓ | ✓ (Rosetta 2) |
| Image Processing | ✓ Required | ✓ Required | ✓ Required | ✓ Required |
| Parallel Computing | ✓ Required | ✓ Required | ✓ Required | ✓ Required |
| Statistics & ML | ✓ Required | ✓ Required | ✓ Required | ✓ Required |
| NVIDIA GPU | ✓ Required | ✓ Required | ✓ Required (rare) | ✗ Not available |
| ffmpeg | Not needed | ✓ Optional | Not needed | Not needed |
| C/C++ Compiler | Optional | Optional | Optional | Optional |
| CUDA Toolkit | Optional | Optional | Optional | Not applicable |

## See Also

- [Installation Guide](../getting-started/installation.md) - Step-by-step installation
- [Installation Issues](../troubleshooting/installation-issues.md) - Common problems
- [GPU Problems](../troubleshooting/gpu-problems.md) - GPU-specific troubleshooting
- [Compilation Errors](../troubleshooting/compilation-errors.md) - Mex/CUDA compilation help

## Quick Reference

**Minimum to run smite:**
- MATLAB R2021a
- Image Processing Toolbox
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox
- NVIDIA GPU (Compute Capability ≥ 5.0)

**Check everything:**
```matlab
checkSMITEDeps()  % Run the script above
```

**Most common issue:**
Missing GPU or incompatible GPU. Verify with `gpuDevice`.
