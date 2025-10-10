---
title: "Quick Start Guide"
category: "getting-started"
level: "beginner"
tags: ["quickstart", "setup", "first-steps", "getting-started"]
prerequisites: []
related: ["installation.md", "first-analysis.md"]
summary: "Get started with smite in 5 minutes - install, verify, and run your first analysis"
estimated_time: "5-10 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Quick Start Guide

## Purpose

This guide gets you up and running with **smite** (Single Molecule Imaging Toolbox Extraordinaire) in just a few minutes. You'll install smite, verify it works, and run a simple localization example.

## Prerequisites

- MATLAB R2021a or later installed
- Basic familiarity with MATLAB
- (Optional) NVIDIA GPU with CUDA capability ≥5.0 for GPU acceleration

## Overview

smite is a MATLAB toolbox for analyzing fluorescence single molecule imaging data, with emphasis on:
- **SMLM** (Single Molecule Localization Microscopy) - super-resolution imaging
- **SPT** (Single Particle Tracking) - tracking particles over time

The toolbox is organized around two key data structures:
- **SMF** (Single Molecule Fitting) - contains all analysis parameters
- **SMD** (Single Molecule Data) - contains all results

## Installation (5 minutes)

### 1. Clone the Repository

```matlab
% In your terminal/command prompt:
cd ~/Documents/MATLAB
git clone https://github.com/LidkeLab/smite.git
```

Or download the latest release from: https://github.com/LidkeLab/smite/releases

### 2. Add to MATLAB Path

Add these lines to your `~/Documents/MATLAB/startup.m` file:

```matlab
addpath '~/Documents/MATLAB/smite/MATLAB'
setupSMITE
```

**Windows users:** Use Windows paths, e.g., `C:\Users\YourName\Documents\MATLAB\smite\MATLAB`

### 3. Restart MATLAB

Close and reopen MATLAB to run the startup script.

**Expected Output:**
```
SMITE version: 1.x.x
```

If you see the version number, installation succeeded!

## Verify Installation (2 minutes)

Test that smite is working:

```matlab
% Create a default SMF (parameter structure)
SMF = smi_core.SingleMoleculeFitting()
```

**Expected Output:**
```
SMF =

  SingleMoleculeFitting with properties:

           Data: [1×1 struct]
     BoxFinding: [1×1 struct]
        Fitting: [1×1 struct]
   Thresholding: [1×1 struct]
FrameConnection: [1×1 struct]
DriftCorrection: [1×1 struct]
       Tracking: [1×1 struct]
```

If you see this structure, smite is properly installed!

## Your First Analysis (3 minutes)

Let's localize molecules in a simulated image:

```matlab
% 1. Generate a test image with random blobs
B = smi_sim.GaussBlobs.genRandomBlobImage();

% 2. Add noise to make it realistic
B = poissrnd(B);

% 3. Create parameter structure with defaults
SMF = smi_core.SingleMoleculeFitting();

% 4. Create localization object
LD = smi_core.LocalizeData(B, SMF);

% 5. Find the localizations!
[SMD] = LD.genLocalizations();

% 6. Display results
fprintf('Found %d localizations\n', numel(SMD.X));
```

**Expected Output:**
```
Found ~50-100 localizations
```

**Explanation:**
- `GaussBlobs` generates simulated fluorescent spots
- `poissrnd()` adds Poisson noise (realistic for photon counting)
- `LocalizeData` does the heavy lifting of finding and fitting blobs
- `SMD` (Single Molecule Data) contains the results: positions (X, Y), photon counts, uncertainties, etc.

## Visualize Results (Bonus)

See your localizations with a colored overlay:

```matlab
% Set verbose mode for visualization
LD.Verbose = 3;

% Localize again to see the overlay
[SMD] = LD.genLocalizations();
```

This creates a color overlay showing:
- **Red**: Raw data
- **Green**: Fitted localizations
- **Blue**: Residuals (what wasn't fit)

## Common Issues

**"Undefined variable or class 'smi_core'"**
- Solution: Restart MATLAB to ensure `setupSMITE` ran
- Or manually run: `setupSMITE` from the MATLAB command line

**"setupSMITE not found"**
- Solution: Check that you added the correct path to `startup.m`
- Verify path with: `which setupSMITE`

**GPU warnings**
- These are normal if you don't have an NVIDIA GPU
- Note: GPU (NVIDIA CUDA) is required for core fitting operations. Some features like visualization and data loading will work without GPU.

## See Also

- [Installation Guide](installation.md) - Detailed setup including GPU, dependencies, and troubleshooting
- [First Analysis](first-analysis.md) - Complete walkthrough of analyzing real data
- [Architecture](../core-concepts/architecture.md) - Understanding SMF and SMD structures

## Next Steps

Now that smite is running:

1. **Learn the concepts**: Read [Architecture](../core-concepts/architecture.md) to understand SMF/SMD paradigm
2. **Try real data**: Follow [First Analysis](first-analysis.md) for a complete workflow
3. **Explore workflows**: Check [SMLM Analysis](../workflows/smlm-analysis.md) for super-resolution imaging

Welcome to smite! 🔬
