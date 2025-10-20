---
title: "Installation Guide"
category: "getting-started"
level: "beginner"
tags: ["installation", "setup", "dependencies", "gpu", "cuda"]
prerequisites: []
related: ["quickstart.md", "first-analysis.md", "../troubleshooting/installation-issues.md"]
summary: "Complete installation guide for smite including dependencies, GPU setup, and compilation"
estimated_time: "10-15 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Installation Guide

## Purpose

This guide provides complete installation instructions for smite, including:
- Installing MATLAB and required toolboxes
- Setting up smite
- GPU/CUDA configuration (optional)
- Compiling mex and CUDA files (when needed)
- Troubleshooting common issues

## Prerequisites

- **MATLAB R2021a or later**
- **Administrator/sudo access** (for some installations)
- **(Optional) NVIDIA GPU** with CUDA compute capability ≥5.0

## System Requirements

### Required

- **Operating System:** Linux, MacOS, or Windows (64-bit)
- **MATLAB:** R2021a or later
- **MATLAB Toolboxes:**
  - Image Processing Toolbox
  - Parallel Computing Toolbox
  - Statistics and Machine Learning Toolbox

###Optional (for specific features)

- **MATLAB Toolboxes:**
  - Curve Fitting Toolbox (for smi_cluster, smi_core.FRC, smi_stat.DiffusionEstimator)
  - Optimization Toolbox (for smi_cluster.PairCorrelation, smi_stat.DiffusionEstimator)
  - Signal Processing Toolbox (for smi_core.FRC)

- **GPU:** NVIDIA GPU with CUDA compute capability ≥5.0 for GPU acceleration

- **External Software:**
  - ffmpeg (Linux only, for video generation when LocalizeData.Verbose ≥ 3)

## Installation Steps

### 1. Install MATLAB and Toolboxes

Download MATLAB from: https://www.mathworks.com/products/matlab.html

Verify required toolboxes are installed:

```matlab
ver
```

Look for:
- Image Processing Toolbox
- Parallel Computing Toolbox
- Statistics and Machine Learning Toolbox

### 2. Download smite

**Option A: Clone with Git (Recommended)**

```bash
cd ~/Documents/MATLAB
git clone https://github.com/LidkeLab/smite.git
```

This gets the latest development version.

**Option B: Download Release**

Download from: https://github.com/LidkeLab/smite/releases

Extract to: `~/Documents/MATLAB/smite/`

### 3. Configure MATLAB Path

Add these lines to your `startup.m` file:

**Location of startup.m:**
- Linux/Mac: `~/Documents/MATLAB/startup.m`
- Windows: `C:\Users\YourName\Documents\MATLAB\startup.m`

**Add to startup.m:**

```matlab
% Setup smite
addpath '~/Documents/MATLAB/smite/MATLAB'
setupSMITE
```

**Windows example:**
```matlab
addpath 'C:\Users\YourName\Documents\MATLAB\smite\MATLAB'
setupSMITE
```

**If startup.m doesn't exist,** create it with just these lines.

### 4. Restart MATLAB

Close and reopen MATLAB to run the startup script.

**Expected Output:**
```
SMITE version: 1.x.x
```

### 5. Verify Installation

```matlab
% Test that smite is accessible
SMF = smi_core.SingleMoleculeFitting()

% Run a quick test
smi_core.LocalizeData.unitTest()
```

If these commands run without errors, installation succeeded!

## GPU Setup (Optional but Recommended)

GPU acceleration significantly speeds up localization. smite includes precompiled CUDA files for common configurations.

### Check GPU Compatibility

```matlab
% Check if you have a compatible GPU
gpuDevice
```

**Expected Output (if GPU present):**
```
ans =

  CUDADevice with properties:

                      Name: 'NVIDIA GeForce RTX 3080'
           ComputeCapability: '8.6'
```

**Requirement:** ComputeCapability must be ≥ 5.0

**If no GPU or incompatible:**
GPU is required for core fitting operations (GaussMLE). Some features like data loading and visualization will work without GPU, but localization/fitting requires CUDA.

### Verify GPU Works

```matlab
% Test GPU localization
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraType = 'EMCCD';
B = smi_sim.GaussBlobs.genRandomBlobImage();
LD = smi_core.LocalizeData(poissrnd(B), SMF);
[SMD] = LD.genLocalizations();
```

If this runs without GPU errors, GPU is working!

## Compilation (Advanced - Usually Not Needed)

smite includes precompiled files for Linux, MacOS, and Windows. You only need to compile if:
- Precompiled files don't work on your system
- You modified source code
- Your GPU architecture is newer than precompiled versions

### When to Compile

**Mex files:** If you get errors like "Invalid MEX-file"
**CUDA files:** If GPU acceleration doesn't work despite having compatible GPU

### Compiling Mex Files

**Prerequisites:**
- C/C++ compiler configured with MATLAB

**Check compiler:**
```matlab
mex -setup
```

If no compiler, install one:
- **Windows:** Install MinGW-w64 via MATLAB Add-Ons
- **Linux:** `sudo apt-get install build-essential`
- **Mac:** Install Xcode command line tools

**Compile:**
```matlab
cd ~/Documents/MATLAB/smite/MATLAB/source/c
mex_Make
```

**Expected Output:**
```
Building with 'gcc'
MEX completed successfully
```

Compiled files go to: `smite/MATLAB/mex/`

### Compiling CUDA Files

**Prerequisites:**
- NVIDIA GPU Computing Toolkit (CUDA) installed
- `nvcc` compiler on PATH
- (Windows only) Visual Studio with `cl.exe` compiler

**Check CUDA installation:**
```bash
nvcc --version
```

**Compile:**
```matlab
cd ~/Documents/MATLAB/smite/MATLAB/source/cuda
cuda_Make
```

**Expected Output:**
```
ptx compilation successful
```

Compiled files go to: `smite/MATLAB/ptx/`

**Troubleshooting compilation:** See [Compilation Errors](../troubleshooting/compilation-errors.md)

## Testing Installation

### Quick Test
```matlab
% Should complete without errors
smi.SMLM.unitTest()
```

### Full Test Suite
```matlab
% Run all unit tests (takes several minutes)
run_tests
```

Tests save results to: `tempdir/smite/unitTest/`

**Note:** Tests are stochastic, so outputs may vary slightly between runs.

## Common Issues

### "setupSMITE not found"

**Solution:**
```matlab
% Check if path is correct
which setupSMITE

% If not found, add path manually:
addpath '~/Documents/MATLAB/smite/MATLAB'
setupSMITE
```

### "Invalid MEX-file" errors

**Solution:** Recompile mex files (see Compiling Mex Files above)

### GPU not detected

**Check:**
```matlab
gpuDevice
```

**Solutions:**
- Install/update NVIDIA drivers
- Verify ComputeCapability ≥ 5.0
- Check that Parallel Computing Toolbox is installed
- Note: MacOS with non-NVIDIA GPUs won't work (Apple Silicon M1/M2/M3 use different architecture)

### MATLAB version too old

smite requires R2021a or later.

**Check version:**
```matlab
version
```

**Solution:** Update MATLAB or use an older smite release

## Platform-Specific Notes

### Linux

**Install ffmpeg (for video generation):**
```bash
sudo apt-get install ffmpeg
```

### MacOS

**Apple Silicon (M1/M2/M3):**
- MATLAB must run under Rosetta 2
- GPU acceleration not available (no NVIDIA GPU)
- CPU processing works fine

### Windows

**Paths:** Use backslashes in Windows paths:
```matlab
addpath 'C:\Users\YourName\Documents\MATLAB\smite\MATLAB'
```

## Updating smite

### If installed via Git:
```bash
cd ~/Documents/MATLAB/smite
git pull
```

Then restart MATLAB.

### If installed from release:
Download new release and extract, overwriting old files.

## Uninstalling

1. Remove from `startup.m`
2. Delete directory: `~/Documents/MATLAB/smite/`
3. Restart MATLAB

## See Also

- [Quick Start](quickstart.md) - Get running in 5 minutes
- [First Analysis](first-analysis.md) - Run your first complete workflow
- [Installation Issues](../troubleshooting/installation-issues.md) - Detailed troubleshooting
- [GPU Problems](../troubleshooting/gpu-problems.md) - GPU-specific issues
- [Compilation Errors](../troubleshooting/compilation-errors.md) - Mex/CUDA compilation help

## Next Steps

After successful installation:
1. **Quick test:** Try [Quick Start Guide](quickstart.md)
2. **Learn concepts:** Read [Architecture](../core-concepts/architecture.md)
3. **Real analysis:** Follow [First Analysis](first-analysis.md)
4. **Explore features:** Check [SMLM Workflow](../workflows/smlm-analysis.md)
