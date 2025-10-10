---
title: "smite Architecture Overview"
category: "core-concepts"
level: "intermediate"
tags: ["architecture", "design", "namespaces", "structure"]
prerequisites: ["getting-started/installation.md"]
related: ["smf-structure.md", "smd-structure.md"]
summary: "Understanding smite's organizational structure, namespaces, and design philosophy"
estimated_time: "10 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# smite Architecture Overview

## Purpose

Understanding smite's architecture helps you navigate the codebase, choose the right tools for your analysis, and extend smite with custom functionality. This document explains the organizational philosophy, namespace structure, and key design patterns that make smite both powerful and accessible.

## Prerequisites

- Familiarity with MATLAB basics
- Understanding of object-oriented programming concepts
- Completion of [installation](../getting-started/installation.md)

## Overview

smite (Single Molecule Imaging Toolbox Extraordinaire) is built on a philosophy of **complete parametric definition**: every analysis is fully specified by parameter structures (SMF) and produces standardized result structures (SMD/TR). This design ensures:

- **Reproducibility**: Same SMF always produces same results
- **Modularity**: Components can be used independently or combined
- **Extensibility**: New methods integrate seamlessly
- **Transparency**: All parameters are explicit and accessible

The codebase uses MATLAB's namespace system (`+directory` notation) to organize functionality into logical groups, making it easy to discover and use the right tools.

## Core Design Philosophy

### Data Structure-Centric Design

smite revolves around three fundamental structures:

1. **SMF (Single Molecule Fitting)**: A structure of structures containing ALL parameters needed to define an analysis pipeline from raw data to results

2. **SMD (Single Molecule Data)**: A structure containing ALL localization results (positions, uncertainties, photon counts, etc.)

3. **TR (Tracking Results)**: An array of SMD structures, one per trajectory, organizing tracking data

**Key insight**: These are implemented as classes to provide helper methods and GUIs, but they function primarily as structures with fixed fields. You work with them like structures:

```matlab
% Create SMF with all defaults
SMF = smi_core.SingleMoleculeFitting();

% Access fields like a structure
boxSize = SMF.BoxFinding.BoxSize;

% Modify fields
SMF.Fitting.PSFSigma = 1.3;
SMF.Thresholding.MaxXY_SE = 0.2;
```

This approach means:
- **No hidden state**: Everything is in the structure
- **Easy to save/load**: Just save the structure
- **Simple to modify**: Direct field access
- **Transparent**: You can see all parameters

### Namespace Organization

smite uses MATLAB's namespace system where `+directory` creates a namespace. This keeps related functionality together and prevents name conflicts.

## Namespace Hierarchy

### High-Level Entry Points: `+smi`

The `+smi` namespace contains user-facing workflows for complete analyses:

```matlab
% Main workflows
smi.SMLM        % Single Molecule Localization Microscopy
smi.SPT         % Single Particle Tracking
smi.BaGoL       % Bayesian Grouping of Localizations
smi.Publish     % Batch processing with standard file structure
```

**When to use**: When you want a complete, end-to-end analysis with minimal setup.

**Example**:
```matlab
% Complete SMLM analysis with GUI
SMLMobj = smi.SMLM();
SMLMobj.fullAnalysis();
```

### Core Functionality: `+smi_core`

The `+smi_core` namespace contains fundamental building blocks used by high-level workflows:

```matlab
% Data structures
smi_core.SingleMoleculeFitting   % SMF parameter structure
smi_core.SingleMoleculeData      % SMD results structure
smi_core.TrackingResults         % TR trajectory structure

% Core processing
smi_core.LoadData                % Load raw data from files
smi_core.LocalizeData            % Find and fit molecules
smi_core.FrameConnection         % Connect localizations across frames
smi_core.DriftCorrection         % Correct stage drift
smi_core.ChannelRegistration     % Multi-color registration

% Lower-level tools
smi_core.FindROI                 % Find regions of interest (boxes)
smi_core.GaussMLE                % Maximum likelihood Gaussian fitting
smi_core.ThresholdFits           % Filter localizations by quality
```

**When to use**: When you want fine-grained control or to build custom pipelines.

**Example**:
```matlab
% Custom pipeline with specific steps
SMF = smi_core.SingleMoleculeFitting();
LD = smi_core.LocalizeData(imageStack, SMF);
SMD = LD.genLocalizations();
```

### Simulation: `+smi_sim`

Generate synthetic data for testing and validation:

```matlab
smi_sim.GaussBlobs      % Generate random Gaussian blobs
smi_sim.SimSMLM         % Simulate complete SMLM datasets
smi_sim.SimSPT          % Simulate particle tracking data
```

**When to use**: Testing algorithms, validating analysis, generating training data.

**Example**:
```matlab
% Create test data
SimSMLM = smi_sim.SimSMLM();
SimSMLM.NFrames = 100;
[SMD, ~, sequence] = SimSMLM.simulateSMLM();
```

### Clustering: `+smi_cluster`

Spatial analysis and clustering algorithms:

```matlab
smi_cluster.Clustering          % Hierarchical clustering
smi_cluster.DBSCAN              % Density-based clustering
smi_cluster.Voronoi             % Voronoi tessellation analysis
smi_cluster.PairCorrelation     % Pair correlation functions
smi_cluster.StatisticsClustering % H-SET clustering
```

**When to use**: Analyzing spatial organization of localizations.

### Statistics: `+smi_stat`

Statistical analysis methods for SPT and dynamics:

```matlab
smi_stat.HMM                    % Hidden Markov Model analysis
smi_stat.DiffusionEstimator     % Estimate diffusion coefficients
smi_stat.ChangeDetection        % Detect state changes in trajectories
```

**When to use**: Analyzing dynamic behavior, state transitions, diffusion.

### Visualization: `+smi_vis`

Tools for creating images and movies:

```matlab
smi_vis.GenerateImages          % Create super-resolution images
smi_vis.GenerateMovies          % Create movies from tracking data
```

**When to use**: Creating publication-quality visualizations.

### Point Spread Function: `+smi_psf`

PSF modeling and manipulation:

```matlab
smi_psf.PointSpreadFunction     % PSF modeling
smi_psf.Zernike                 % Zernike polynomial tools
```

**When to use**: Advanced PSF modeling, 3D localization.

### Helpers: `+smi_helpers`

Utility functions used throughout smite:

```matlab
smi_helpers.Filters             % Filtering operations
smi_helpers.ROITools            % ROI manipulation
smi_helpers.mkSMITETmpDir       % Create temporary directories
```

**When to use**: Supporting operations in custom scripts.

## Typical Usage Patterns

### Pattern 1: High-Level GUI Workflow

For quick analyses with visual feedback:

```matlab
% SMLM with GUI
SMLMobj = smi.SMLM();  % No args = GUI opens
% Use GUI to set parameters
SMLMobj.fullAnalysis();

% SPT with GUI
SPTobj = smi.SPT();
SPTobj.gui();
```

### Pattern 2: Script-Based Workflow

For reproducible, automated analyses:

```matlab
% Configure SMF
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'data.h5'};
SMF.Fitting.PSFSigma = 1.3;

% Run analysis
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();
```

### Pattern 3: Low-Level Component Usage

For custom pipelines or method development:

```matlab
% Load data manually
LD = smi_core.LoadData();
[~, sequence, SMF] = LD.loadRawData(SMF, 1);

% Localize with custom settings
LDObj = smi_core.LocalizeData(sequence, SMF);
LDObj.Verbose = 3;
SMD = LDObj.genLocalizations();

% Apply custom processing
% ... your custom analysis ...

% Save results manually
save('my_results.mat', 'SMD', 'SMF');
```

### Pattern 4: Batch Processing

For analyzing many datasets with consistent parameters:

```matlab
% Setup Publish object for batch processing
Pub = smi.Publish();
Pub.SMF.Data.FileDir = '/path/to/CoverslipDir';
% Files should be organized as: Cell*/Label*/Data*.h5
Pub.analyzeCells();  % Analyzes all matching files
```

## Data Flow Through smite

A typical SMLM analysis flows through these components:

```
Raw Data (.h5/.mat)
    ↓
LoadData → Image Stack (ADU)
    ↓
DataToPhotons → Image Stack (photons)
    ↓
LocalizeData → SMD (all localizations)
    ↓
ThresholdFits → SMD (filtered)
    ↓
FrameConnection → SMD (connected emitters)
    ↓
DriftCorrection → SMD (drift-corrected)
    ↓
GenerateImages → Super-Resolution Images
```

Each step is controlled by parameters in SMF and produces or modifies SMD.

## Image Coordinate Convention

smite uses MATLAB's image coordinate system:

- **Origin**: (1, 1) is the center of the top-left pixel
- **Axes**: First index is Y (row, vertical), second is X (column, horizontal)
- **Movement**: (2, 1) is one pixel down from (1, 1)

```matlab
% Image indexing
pixel_value = image(y, x);  % Note: Y first, X second

% Localization coordinates in SMD
SMD.X  % X position (horizontal, column)
SMD.Y  % Y position (vertical, row)
```

This is important when plotting or transforming coordinates.

## Extension Points

smite is designed for extension. Common customization points:

### Add Custom Fitting Methods

Extend `smi_core.GaussMLE` or create new fitting classes that accept SMF and return SMD.

### Add New Drift Correction Methods

Add methods to `smi_core.DriftCorrection` following the existing pattern.

### Add Analysis Tools

Create new classes in appropriate namespaces (`+smi_cluster`, `+smi_stat`) that work with SMD structures.

### Custom Visualization

Extend `smi_vis.GenerateImages` or create standalone visualization functions.

## File Organization

Within the repository:

```
smite/
├── MATLAB/
│   ├── +smi/              # High-level workflows
│   ├── +smi_core/         # Core functionality
│   ├── +smi_sim/          # Simulation
│   ├── +smi_cluster/      # Clustering
│   ├── +smi_stat/         # Statistics
│   ├── +smi_vis/          # Visualization
│   ├── +smi_psf/          # PSF tools
│   ├── +smi_helpers/      # Utilities
│   ├── examples/          # Example scripts
│   ├── mex/               # Compiled MEX files
│   ├── ptx/               # CUDA PTX files
│   └── source/            # C/CUDA source code
└── doc/                   # Documentation
```

## Performance Considerations

### GPU Requirement

Fitting operations require NVIDIA CUDA GPU:

```matlab
% Requires CUDA GPU (compute capability ≥5.0)
LD = smi_core.LocalizeData(sequence, SMF);
% GaussMLE fitting requires GPU - CPU-only mode not supported
```

### MEX Acceleration

Performance-critical loops are implemented in MEX (C/C++):

- Frame connection algorithms
- Some fitting routines
- Drift correction calculations

Precompiled binaries are included for Linux, macOS, and Windows.

### Parallel Processing

Some operations can use MATLAB's parallel processing:

```matlab
% Enable parallel pools for large batch jobs
parpool('local', 4);  % 4 workers
Pub = smi.Publish();
Pub.analyzeCells();  % May use parallel workers
```

## See Also

- [SMF Structure](smf-structure.md) - Complete parameter reference
- [SMD Structure](smd-structure.md) - Results data format
- [SMLM Workflow](../workflows/smlm-analysis.md) - End-to-end SMLM pipeline
- [SPT Workflow](../workflows/spt-tracking.md) - Tracking pipeline
- CLAUDE.md (in repository root) - Developer-focused architecture guide
