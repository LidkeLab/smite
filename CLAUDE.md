# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Repository Overview

**smite** (Single Molecule Imaging Toolbox Extraordinaire) is a MATLAB-based toolbox for fluorescence single molecule imaging analysis, with emphasis on single molecule localization microscopy (SMLM) and single particle tracking (SPT).

## Core Architecture

### Data Structure Philosophy

smite is designed around two fundamental data structures that define all analysis:

- **SMF (Single Molecule Fitting)**: A structure of structures containing all parameters required to go from raw data to results. This completely defines the data analysis pipeline.
- **SMD (Single Molecule Data)**: Contains all localization results including positions, photons, uncertainties, frame numbers, etc.
- **TR (Tracking Results)**: A variation of SMD specifically for tracking data.

These structures are inputs/outputs throughout the codebase. While implemented as classes, they function primarily as structures with fixed fields.

### Namespace Organization

smite uses MATLAB's namespace system with `+` directories:

- **+smi**: Highest-level entry point with main workflows (SMLM, SPT, BaGoL, Publish)
- **+smi_core**: Core functionality (LocalizeData, DriftCorrection, FrameConnection, LoadData, etc.)
- **+smi_sim**: Simulation tools (GaussBlobs, SimSMLM, SimSPT)
- **+smi_cluster**: Clustering algorithms (DBSCAN, Voronoi, H-SET, PairCorrelation)
- **+smi_stat**: Statistical methods (HMM, ChangeDetection, DiffusionEstimator)
- **+smi_vis**: Visualization tools (GenerateImages, GenerateMovies)
- **+smi_psf**: Point spread function tools
- **+smi_helpers**: Helper utilities (Filters, ROITools)

### Main Workflows

1. **SMLM (smi.SMLM)**: Processes 2D super-resolution data from .h5/.mat files
2. **SPT (smi.SPT)**: Single particle tracking analysis
3. **BaGoL (smi.BaGoL)**: Bayesian Grouping of Localizations for determining emitter positions
4. **Publish (smi.Publish)**: Batch processing with standard naming convention (CoverslipDir/Cell*/Label*/Data*.h5)

### Image Coordinate System

- Images follow MATLAB's column-major format
- Coordinate (1,1) = center of top-left pixel
- Coordinate (2,1) = one pixel down from top, leftmost column

## Development Commands

### Setup and Installation

After cloning, add to `startup.m`:
```matlab
addpath '~/Documents/MATLAB/smite/MATLAB'
setupSMITE
```

This runs `setupSMITE.m` which:
- Adds required paths (MATLAB, mex, ptx directories)
- Sets up external software dependencies
- Reports smite version

### Running Tests

Run complete test suite:
```matlab
run_tests
```

This executes unit tests for all major components. Results saved to `tempdir/smite/unitTest/name_of_test/`. Expected results available in `MATLAB/ExpectedResults/`.

Key tests (run these first as indicators of proper installation):
- `smi.SMLM.unitTest` - Tests major SMLM functionality
- `smi.SPT.unitTestFFGC` - Tests frame-to-frame and gap closing processes

Run individual class unit tests:
```matlab
smi_core.LocalizeData.unitTest()
smi_cluster.Clustering.unitTest()
smi_psf.PointSpreadFunction.unitTest()
```

### Compilation (When Needed)

**Mex files** (if source files modified or mex files don't work):
```matlab
cd ~/Documents/MATLAB/smite/MATLAB/source/c
mex_Make
```
Requires C/C++ compiler setup (`mex -setup`). Outputs to `MATLAB/mex/` with platform-specific extensions (.mexa64, .mexmaci64, .mexw64).

**CUDA files** (if source files modified or ptx files don't work):
```matlab
cd ~/Documents/MATLAB/smite/MATLAB/source/cuda
cuda_Make
```
Requires `nvcc` compiler on PATH (and `cl.exe` on Windows). Outputs to `MATLAB/ptx/`. GPU compute capability must be ≥5.0.

## Key Design Patterns

### Working with SMF

```matlab
% Create SMF object with defaults
SMF = smi_core.SingleMoleculeFitting()

% Access properties
B = SMF.BoxFinding.BoxOverlap

% Set properties
SMF.BoxFinding.BoxOverlap = 0

% Interactive GUI
SMF.gui()
```

### Typical Analysis Pipeline

```matlab
% 1. Create/configure SMF
SMF = smi_core.SingleMoleculeFitting()

% 2. Load and localize data
LD = smi_core.LocalizeData(imageStack, SMF)
[SMD] = LD.genLocalizations()

% 3. Optional: Frame connection
if SMF.FrameConnection.On
    % Frame connection happens in pipeline
end

% 4. Optional: Drift correction
if SMF.DriftCorrection.On
    % Drift correction happens in pipeline
end
```

### GUI-Based Workflow

For SMLM analysis with GUI:
```matlab
SMLMobj = smi.SMLM()  % Opens GUI when no arguments
% Use GUI to navigate to data, set parameters, run analysis
```

## File Formats

- **Raw data**: .h5 (preferred) or .mat files with default variable name 'sequence'
- **Camera calibration**: Required for full functionality (see doc/FileFormats/CalibrationFile.md)
- **Channel registration**: Optional .mat files for multi-channel data
- **Results**: Saved to `Data.ResultsDir` (default: FileDir/Results)

## Testing Philosophy

Tests are stochastic by nature, so outputs may vary between runs. Tests generate their own data when needed. Large result files are excluded from ExpectedResults to keep repository size manageable.

## Dependencies

- MATLAB R2021a or later
- Required toolboxes: Image Processing, Parallel Computing, Statistics and Machine Learning
- Optional toolboxes: Curve Fitting, Optimization, Signal Processing (for specific features)
- NVIDIA GPU with CUDA compute capability ≥5.0 (supported by MATLAB version)
- Linux: ffmpeg for video generation (when LocalizeData.Verbose ≥ 3)

## Important Notes

- Precompiled mex and CUDA files are included for 64-bit Linux, MacOS, and Windows
- Recompilation only needed if source files change or precompiled files fail
- MATLAB doesn't support non-NVIDIA GPUs for CUDA operations
- The main branch is used for both development and releases
- Issues and support requests go through GitHub Issues
- New contributions should be on separate branches with Pull Requests to main
