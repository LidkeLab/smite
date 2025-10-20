---
title: "API Reference: +smi_stat Namespace"
category: "api-reference"
level: "advanced"
tags: ["api", "smi_stat", "statistics", "diffusion", "hmm", "change-detection", "trajectory-analysis"]
prerequisites: ["../core-concepts/architecture.md", "../core-concepts/tr-structure.md", "../workflows/spt-tracking.md"]
related: ["../examples/tracking-diffusion.md", "../workflows/spt-tracking.md", "smi-core.md"]
summary: "Complete API reference for the +smi_stat namespace covering statistical analysis methods for single molecule data, including diffusion estimation, Hidden Markov Models, and change point detection"
estimated_time: "35 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# API Reference: +smi_stat Namespace

## Purpose

The +smi_stat namespace provides advanced statistical analysis tools for single molecule data, with particular emphasis on trajectory analysis. These classes enable quantitative characterization of molecular dynamics, state transitions, and temporal patterns in single particle tracking and single molecule imaging data. Understanding these tools is essential for extracting biophysical parameters from SPT experiments and performing rigorous statistical analyses of molecular behavior.

## Prerequisites

- Understanding of [smite architecture](../core-concepts/architecture.md)
- Familiarity with [TR structure](../core-concepts/tr-structure.md)
- Knowledge of [SPT workflow](../workflows/spt-tracking.md)
- Basic understanding of statistical methods (MSD, maximum likelihood estimation)
- Familiarity with Hidden Markov Models (for HMM class)

## Overview

The +smi_stat namespace provides three primary analysis classes and supporting utilities:

**Core Statistical Classes:**
- `DiffusionEstimator` - Estimate diffusion coefficients from trajectory data using MSD, CDF, or MLE methods
- `HMM` - Hidden Markov Model analysis for detecting state transitions in multi-color tracking data
- `ChangeDetection` - Bayesian change point detection in intensity time series

**Supporting Functions:**
- Image registration utilities (findOffset, findCoordAffine)
- Bootstrap fitting methods (bootstrapFit, bootstrapFitCon)
- Colocalization metrics (mandersSplitCoefs, pearsonCorrCoef)

These tools extend beyond basic localization and tracking to provide quantitative biophysical insights from single molecule data.

---

## DiffusionEstimator

### Description

`DiffusionEstimator` estimates diffusion coefficients from single particle tracking data using multiple analysis methods. It computes Mean Squared Displacement (MSD), cumulative distribution functions of displacements, or maximum likelihood estimates from jump statistics. The class supports both single-trajectory and ensemble analysis, multi-component diffusion models, and automatic standard error estimation.

### Key Concept

Diffusion analysis characterizes the random motion of molecules. The diffusion coefficient D relates to the MSD through `<r^2> = 4Dt` in 2D. Different methods (MSD fitting, CDF fitting, MLE) offer trade-offs between computational speed, statistical efficiency, and robustness to heterogeneity.

### Class Definition

```matlab
classdef DiffusionEstimator < handle
```

### Constructor

```matlab
DE = smi_stat.DiffusionEstimator()
DE = smi_stat.DiffusionEstimator(TR, SMF)
DE = smi_stat.DiffusionEstimator(TR, SMF, Verbose, AutoRun)
```

**Parameters:**
- `TR`: Tracking Results structure array
- `SMF`: Single Molecule Fitting structure (for pixel size, frame rate)
- `Verbose`: Verbosity level (0-3)
- `AutoRun`: Boolean flag to immediately run analysis (default: false)

**Returns:**
- `DE`: DiffusionEstimator object instance

**Example:**
```matlab
% Basic creation
DE = smi_stat.DiffusionEstimator();
DE.TR = myTrajectories;
DE.SMF = mySMF;

% Auto-run on creation
[DE, Results] = smi_stat.DiffusionEstimator(TR, SMF, 1, true);
```

### Properties

#### Algorithm Parameters

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `DiffusionModel` | string | 'Brownian' | Diffusion model type |
| `NComponents` | integer | 1 | Number of diffusing populations |
| `FitMethod` | string | 'WeightedLS' | Fitting method ('WeightedLS' or 'LS') |
| `FitTarget` | string | 'MSD' | Data to fit ('MSD', 'CDFOfJumps', 'LikelihoodOfJumps') |
| `FrameLagRange` | array | [1, 5] | Range of frame lags for analysis |
| `NFitPoints` | integer | 5 | Number of MSD points to fit |
| `NDimensions` | integer | 2 | Spatial dimensions (2 or 3) |

#### Analysis Control

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `EstimateSEs` | logical | true | Compute standard errors (may be slow) |
| `FitIndividualTrajectories` | logical | true | Fit each trajectory separately |
| `UnitFlag` | logical | false | Output units (false=pixels/frames, true=micrometers/seconds) |
| `Verbose` | integer | 0 | Verbosity level (0=silent, 1=progress, 2=detailed, 3=debug) |

#### Input/Output

| Property | Type | Description |
|----------|------|-------------|
| `TR` | struct | Tracking Results input |
| `SMF` | struct | SMF structure for calibration |
| `SaveDir` | string | Directory for saving results |
| `BaseSaveName` | string | Base filename for saved results |

#### Results (SetAccess = protected)

| Property | Type | Description |
|----------|------|-------------|
| `DiffusionStruct` | struct | Main results structure with D estimates |
| `MSDSingleTraj` | struct | Trajectory-wise MSD data |
| `MSDEnsemble` | struct | Ensemble-averaged MSD |

### Methods

#### Primary Analysis Method

##### estimateDiffusionConstant

```matlab
DiffusionStruct = obj.estimateDiffusionConstant()
```

Main analysis method that computes diffusion coefficients from tracking data.

**Workflow:**
1. Computes MSD for all trajectories and ensemble
2. Depending on `FitTarget`:
   - 'MSD': Fits linear model to MSD vs. lag time
   - 'CDFOfJumps': Computes and fits cumulative distribution of displacements
   - 'LikelihoodOfJumps': Maximum likelihood estimation from displacement distribution
3. Estimates standard errors (if `EstimateSEs = true`)
4. Converts units (if `UnitFlag = true`)
5. Returns both trajectory-wise and ensemble results

**Returns:**
- `DiffusionStruct`: Structure array with fields:
  - `Name`: 'trajectory' or 'ensemble'
  - `Units`: {JumpUnit; TimeUnit}
  - `FitParams`: Fit parameters (depends on method)
  - `FitParamsSE`: Standard errors of fit parameters
  - `DiffusionConstant`: Estimated D value(s)
  - `DiffusionConstantSE`: Standard error(s) of D
  - `PopulationRatios`: Fraction in each component (multi-component only)
  - `PopulationRatiosSE`: Standard errors of ratios
  - `NFitPoints`: Number of points used in fit

**Example:**
```matlab
% Load tracking results
load('TrackingResults.mat', 'TR', 'SMF');

% Setup diffusion estimator
DE = smi_stat.DiffusionEstimator();
DE.TR = TR;
DE.SMF = SMF;
DE.FitTarget = 'MSD';
DE.FrameLagRange = [1, 5];
DE.NFitPoints = 5;
DE.UnitFlag = true;  % Get results in physical units
DE.Verbose = 1;

% Estimate diffusion
Results = DE.estimateDiffusionConstant();

% Extract ensemble results
D = Results(2).DiffusionConstant;
D_SE = Results(2).DiffusionConstantSE;
fprintf('Ensemble diffusion: %.3f +/- %.3f um^2/s\n', D, D_SE);

% Extract trajectory-wise results
D_traj = Results(1).DiffusionConstant;
fprintf('Mean trajectory D: %.3f um^2/s\n', mean(D_traj));
fprintf('Median trajectory D: %.3f um^2/s\n', median(D_traj));
```

##### saveResults

```matlab
obj.saveResults()
obj.saveResults(SaveParams)
```

Saves diffusion analysis results to disk.

**Parameters:**
- `SaveParams`: Structure with optional fields:
  - `SaveDir`: Save directory (default: obj.SaveDir)
  - `BaseSaveName`: Base filename (default: obj.BaseSaveName)

**Saves:**
- `[BaseSaveName]_DiffusionResults.mat` containing:
  - `DiffusionStruct`: Main results
  - `MSDSingleTraj`: Trajectory MSDs
  - `MSDEnsemble`: Ensemble MSD
- MSD plots (if generated)

**Example:**
```matlab
DE.SaveDir = '/results/diffusion';
DE.BaseSaveName = 'Cell1_Label1';
DE.saveResults();
```

### Static Methods

##### computeMSD

```matlab
[MSDSingleTraj, MSDEnsemble] = smi_stat.DiffusionEstimator.computeMSD(TR, FrameLagRange, Verbose)
```

Computes Mean Squared Displacement for trajectories.

**Parameters:**
- `TR`: Tracking Results structure
- `FrameLagRange`: [MinLag, MaxLag] in frames
- `Verbose`: Verbosity level

**Returns:**
- `MSDSingleTraj`: Structure array (one per trajectory) with fields:
  - `MSD`: MSD values at each lag
  - `FrameLags`: Lag times in frames
  - `NPoints`: Number of displacement pairs at each lag
- `MSDEnsemble`: Structure with ensemble-averaged MSD:
  - `MSD`: Ensemble MSD
  - `FrameLags`: Lag times
  - `SquaredDisplacement`: All squared displacements
  - `FrameLagsAll`: Corresponding lags for each displacement

**Example:**
```matlab
% Compute MSD for lag times 1-10 frames
[MSDTraj, MSDEns] = smi_stat.DiffusionEstimator.computeMSD(TR, [1, 10], 1);

% Plot ensemble MSD
figure;
errorbar(MSDEns.FrameLags, MSDEns.MSD, MSDEns.MSD_SE);
xlabel('Lag time (frames)');
ylabel('MSD (pixels^2)');
```

##### fitMSD

```matlab
[FitParams, FitParamsSE] = smi_stat.DiffusionEstimator.fitMSD(MSDStruct, FitMethod, NFitPoints, DiffusionModel, Verbose)
```

Fits MSD data to extract diffusion coefficient.

**Parameters:**
- `MSDStruct`: MSD structure from `computeMSD()`
- `FitMethod`: 'WeightedLS' or 'LS'
- `NFitPoints`: Number of initial points to fit
- `DiffusionModel`: 'Brownian' (only option currently)
- `Verbose`: Verbosity level

**Returns:**
- `FitParams`: [Intercept, Slope] where D = Slope / (4 * NDim)
- `FitParamsSE`: Standard errors of fit parameters

**Model:** `MSD(t) = 4Dt + sigma^2` (2D) where sigma^2 is localization uncertainty

##### fitCDFOfJumps

```matlab
[FitParams, FitParamsSE] = smi_stat.DiffusionEstimator.fitCDFOfJumps(MSDStruct, FitMethod, NComponents, DiffusionModel, Verbose)
```

Fits cumulative distribution function of jump distances.

**Parameters:**
- `MSDStruct`: MSD structure with CDFOfJumps computed
- `FitMethod`: 'WeightedLS' or 'LS'
- `NComponents`: Number of diffusing populations (1, 2, 3, ...)
- `DiffusionModel`: 'Brownian'
- `Verbose`: Verbosity level

**Returns:**
- `FitParams`: [D1, D2, ..., DN, p1, p2, ...] where Di are diffusion constants and pi are population fractions
- `FitParamsSE`: Standard errors

**Use case:** Detecting heterogeneous diffusion (multiple populations)

##### mleOfJumps

```matlab
[MLEParams, MLEParamsSE] = smi_stat.DiffusionEstimator.mleOfJumps(MSDStruct, NComponents, DiffusionModel, Verbose)
```

Maximum likelihood estimation from jump distribution.

**Parameters:**
- `MSDStruct`: MSD structure
- `NComponents`: Number of components
- `DiffusionModel`: 'Brownian'
- `Verbose`: Verbosity level

**Returns:**
- `MLEParams`: MLE parameter estimates
- `MLEParamsSE`: Standard errors from observed Fisher information

**Advantage:** Most statistically efficient method, optimal for multi-component models

##### plotEnsembleMSD

```matlab
PlotAxes = smi_stat.DiffusionEstimator.plotEnsembleMSD(PlotAxes, MSDEnsemble, DiffusionStruct, DiffusionModel, UnitFlag)
```

Plots ensemble MSD with fit overlay.

**Parameters:**
- `PlotAxes`: Axes handle ([] for new figure)
- `MSDEnsemble`: Ensemble MSD structure
- `DiffusionStruct`: Results from fitting
- `DiffusionModel`: 'Brownian'
- `UnitFlag`: Use physical units

**Returns:**
- `PlotAxes`: Axes handle

##### plotEnsembleCDFOfJumps

```matlab
PlotAxes = smi_stat.DiffusionEstimator.plotEnsembleCDFOfJumps(PlotAxes, MSDEnsemble, DiffusionStruct, UnitFlag)
```

Plots CDF of jumps with fit overlay.

### Usage Patterns

#### Pattern 1: Basic MSD Analysis

```matlab
% Load tracking data
load('Tracking_TR.mat', 'TR');
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.PixelSize = 0.108;  % micrometers
SMF.Data.FrameRate = 50;     % Hz

% Setup diffusion analysis
DE = smi_stat.DiffusionEstimator();
DE.TR = TR;
DE.SMF = SMF;
DE.FitTarget = 'MSD';
DE.FrameLagRange = [1, 5];
DE.NFitPoints = 4;
DE.UnitFlag = true;  % Physical units

% Run analysis
Results = DE.estimateDiffusionConstant();

% Report ensemble results
D = Results(2).DiffusionConstant;
D_SE = Results(2).DiffusionConstantSE;
fprintf('Diffusion coefficient: %.3f +/- %.3f um^2/s\n', D, D_SE);

% Plot
DE.plotEnsembleMSD([], DE.MSDEnsemble, Results(2), 'Brownian', true);
```

#### Pattern 2: Multi-Component Analysis

```matlab
% Detect heterogeneous diffusion populations
DE = smi_stat.DiffusionEstimator();
DE.TR = TR;
DE.SMF = SMF;
DE.FitTarget = 'CDFOfJumps';  % Better for multi-component
DE.NComponents = 2;            % Test 2 populations
DE.FrameLagRange = [1, 1];     % Single-frame lags
DE.UnitFlag = true;

% Estimate
Results = DE.estimateDiffusionConstant();

% Extract ensemble results
D = Results(2).DiffusionConstant;
D_SE = Results(2).DiffusionConstantSE;
PopRatio = Results(2).PopulationRatios;
PopRatio_SE = Results(2).PopulationRatiosSE;

fprintf('Component 1: D = %.3f +/- %.3f um^2/s (%.1f%%)\n', ...
    D(1), D_SE(1), 100*PopRatio(1));
fprintf('Component 2: D = %.3f +/- %.3f um^2/s (%.1f%%)\n', ...
    D(2), D_SE(2), 100*(1-PopRatio(1)));

% Plot CDF with fit
DE.plotEnsembleCDFOfJumps([], DE.MSDEnsemble, Results(2), true);
```

#### Pattern 3: Trajectory-Wise Analysis

```matlab
% Analyze diffusion on per-trajectory basis
DE = smi_stat.DiffusionEstimator();
DE.TR = TR;
DE.SMF = SMF;
DE.FitTarget = 'MSD';
DE.FitIndividualTrajectories = true;
DE.UnitFlag = true;

Results = DE.estimateDiffusionConstant();

% Get trajectory-wise D values
D_traj = Results(1).DiffusionConstant;
D_traj_SE = Results(1).DiffusionConstantSE;

% Histogram of diffusion coefficients
figure;
histogram(D_traj, 20);
xlabel('Diffusion coefficient (um^2/s)');
ylabel('Count');
title(sprintf('Mean D = %.3f um^2/s', mean(D_traj)));

% Find fast and slow trajectories
D_threshold = median(D_traj);
FastTraj = find(D_traj > D_threshold);
SlowTraj = find(D_traj <= D_threshold);
fprintf('%d fast, %d slow trajectories\n', length(FastTraj), length(SlowTraj));
```

#### Pattern 4: Maximum Likelihood Estimation

```matlab
% Most statistically efficient method
DE = smi_stat.DiffusionEstimator();
DE.TR = TR;
DE.SMF = SMF;
DE.FitTarget = 'LikelihoodOfJumps';
DE.NComponents = 1;
DE.EstimateSEs = true;  % Get standard errors from Fisher information
DE.UnitFlag = true;

Results = DE.estimateDiffusionConstant();

% MLE provides tightest confidence intervals
D_MLE = Results(2).DiffusionConstant;
D_MLE_SE = Results(2).DiffusionConstantSE;
fprintf('MLE: D = %.4f +/- %.4f um^2/s\n', D_MLE, D_MLE_SE);
fprintf('95%% CI: [%.4f, %.4f] um^2/s\n', ...
    D_MLE - 1.96*D_MLE_SE, D_MLE + 1.96*D_MLE_SE);
```

### See Also

- [SPT Tracking Workflow](../workflows/spt-tracking.md)
- [Tracking Example](../examples/tracking-diffusion.md)
- `smi.SPT` - Tracking analysis
- `smi_core.TrackingResults` - TR structure

---

## HMM

### Description

`HMM` performs Hidden Markov Model analysis on multi-color single particle tracking data to detect state transitions and quantify interaction kinetics. Designed for analyzing dimer formation between fluorescently labeled proteins, the class identifies transitions between "free" and "dimer" states from inter-channel separation statistics and estimates association/dissociation rate constants.

### Key Concept

Hidden Markov Models assume observed data (inter-channel separations) arise from a hidden state sequence (free vs. bound). By modeling emission probabilities (separation distributions per state) and transition rates between states, HMM analysis infers the most likely state sequence and quantifies kinetic parameters like dimerization rates.

### Class Definition

```matlab
classdef HMM < handle
```

### Constructor

```matlab
H = smi_stat.HMM()
H = smi_stat.HMM(TRArray, SMF)
```

**Parameters:**
- `TRArray`: N × 2 array of TR structures (dimer candidate pairs)
- `SMF`: Single Molecule Fitting structure(s)

**Returns:**
- `H`: HMM object instance

**Example:**
```matlab
% Basic creation
H = smi_stat.HMM();

% With data
H = smi_stat.HMM(DimerCandidates, SMF);
```

### Properties

#### Model Parameters

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `ModelSpecifier` | string | 'DF' | Pre-built model ('DF'=dimer/free, 'DDF'=dimer/domain/free) |
| `DimerSeparation` | numeric | 0.5 | Expected dimer separation (pixels) |
| `DomainSeparation` | numeric | 2.0 | Typical domain size (pixels) |
| `MaxSeparation` | numeric | 5.0 | Max separation for dimer candidates (pixels) |
| `RegErrorInflation` | numeric | 0.0 | Additional registration error (pixels) |
| `PDFHandles` | cell | [] | Custom emission PDF function handles |
| `RateParametersGuess` | array | [0.01; 0.01] | Initial guess for rate parameters (1/frames) |

#### Data Organization

| Property | Type | Description |
|----------|------|-------------|
| `TRArray` | struct | N × 2 array of trajectory pairs |
| `SMF` | struct | SMF structure(s) for pixel size, frame rate |
| `ChannelNames` | cell | {'Channel 1'; 'Channel 2'} |
| `StateNames` | cell | {'Dimer'; 'Free'} |

#### Output Control

| Property | Type | Default | Description |
|----------|------|---------|-------------|
| `SaveResults` | logical | true | Save results to disk |
| `GenerateMovies` | logical | false | Generate trajectory movies |
| `GeneratePlots` | logical | [true; true] | [durations histogram; summary plots] |
| `MovieParams` | struct | - | Movie generation parameters |
| `PlotParams` | struct | - | Plotting parameters |
| `SaveDir` | string | '' | Save directory |
| `ConditionLabel` | string | '' | Label for experimental condition |
| `UnitFlag` | logical | false | Output in physical units |
| `Verbose` | integer | 1 | Verbosity level |

#### Results (SetAccess = protected)

| Property | Type | Description |
|----------|------|-------------|
| `RateParameters` | array | Estimated transition rates |
| `RateParametersSE` | array | Standard errors of rates |
| `TRArrayTrunc` | struct | Pre-processed TR array used in analysis |

### Methods

#### Primary Analysis Method

##### performFullAnalysis

```matlab
[RateParameters, RateParametersSE, LogLikelihood] = obj.performFullAnalysis()
```

Main analysis method that performs complete HMM analysis.

**Workflow:**
1. Generates emission PDFs for each state
2. Computes pair-wise separations
3. Prepares emission matrices for all trajectory pairs
4. Estimates rate parameters via maximum likelihood
5. Computes Viterbi path (most likely state sequence) for each pair
6. Calculates dimer durations
7. Saves results (if `SaveResults = true`)

**Returns:**
- `RateParameters`: Estimated transition rates (1/frames)
  - For 'DF' model: [k_on, k_off]
  - For 'DDF' model: [k_on_dimer, k_off_dimer, k_on_domain, k_off_domain]
- `RateParametersSE`: Standard errors of rate parameters
- `LogLikelihood`: Log-likelihood of data given model

**Example:**
```matlab
% Load tracking data
load('TwoColorTracking.mat', 'TRArray', 'SMF');

% Setup HMM analysis
H = smi_stat.HMM(TRArray, SMF);
H.ModelSpecifier = 'DF';
H.DimerSeparation = 0.5;  % pixels
H.MaxSeparation = 5.0;    % pixels
H.UnitFlag = true;        % Physical units
H.Verbose = 2;

% Run analysis
[Rates, RatesSE, LogL] = H.performFullAnalysis();

% Report results
fprintf('Association rate: %.3f +/- %.3f (1/s)\n', Rates(1), RatesSE(1));
fprintf('Dissociation rate: %.3f +/- %.3f (1/s)\n', Rates(2), RatesSE(2));
fprintf('Log-likelihood: %.2f\n', LogL);

% Calculate derived quantities
k_on = Rates(1);
k_off = Rates(2);
K_d = k_off / k_on;  % Equilibrium dissociation constant
tau_dimer = 1 / k_off;  % Mean dimer lifetime
fprintf('K_d: %.3f\n', K_d);
fprintf('Mean dimer lifetime: %.2f s\n', tau_dimer);
```

##### saveResults

```matlab
obj.saveResults()
```

Saves HMM analysis results and generates plots.

**Saves:**
- `HMMResults.mat` containing:
  - `TRArray`: Full trajectory array with state sequences
  - `RateParameters`, `RateParametersSE`
  - Analysis metadata
- Dimer duration histograms (if `GeneratePlots(1) = true`)
- Summary plots (if `GeneratePlots(2) = true`)

### Static Methods

##### findDimerCandidates

```matlab
TRArray = smi_stat.HMM.findDimerCandidates(TR1, TR2, MaxDimerSeparation, MaxSeparation, MinValidPoints, MinPhotons, BorderPadding)
```

Finds trajectory pairs that may represent dimers.

**Parameters:**
- `TR1`, `TR2`: Tracking results from channels 1 and 2
- `MaxDimerSeparation`: Max separation for dimer state (pixels)
- `MaxSeparation`: Max separation for candidate consideration (pixels)
- `MinValidPoints`: Minimum overlapping localizations
- `MinPhotons`: Minimum photon count threshold
- `BorderPadding`: Pixels to exclude near image border

**Returns:**
- `TRArray`: N × 2 array of candidate trajectory pairs

**Example:**
```matlab
% Load two-channel tracking results
load('Channel1_TR.mat', 'TR');
TR1 = TR;
load('Channel2_TR.mat', 'TR');
TR2 = TR;

% Find candidates
TRArray = smi_stat.HMM.findDimerCandidates(...
    TR1, TR2, ...
    0.5, ...  % MaxDimerSeparation
    5.0, ...  % MaxSeparation
    10, ...   % MinValidPoints
    500, ...  % MinPhotons
    5);       % BorderPadding

fprintf('Found %d dimer candidates\n', size(TRArray, 1));
```

##### findDimerCandidatesFromFiles

```matlab
[TRArray, SMFArray, FileList] = smi_stat.HMM.findDimerCandidatesFromFiles(FileDir, FilePatterns, MaxDimerSeparation, MaxSeparation, MinValidPoints, MinPhotons, BorderPadding, Verbose)
```

Batch processing: finds candidates from multiple files.

**Parameters:**
- `FileDir`: Directory containing TR files
- `FilePatterns`: Cell array of filename patterns for each channel
- Other parameters: Same as `findDimerCandidates`
- `Verbose`: Verbosity level

**Returns:**
- `TRArray`: Combined candidate array
- `SMFArray`: Corresponding SMF structures
- `FileList`: List of processed files

##### generateEmissionPDFs

```matlab
PDFHandles = smi_stat.HMM.generateEmissionPDFs(ModelSpecifier)
```

Generates emission probability density functions for specified model.

**Parameters:**
- `ModelSpecifier`: 'DF' or 'DDF'

**Returns:**
- `PDFHandles`: Cell array of function handles for each state's emission PDF

##### estimateRateParameters

```matlab
[RateParameters, RateParametersSE, LogLikelihood] = smi_stat.HMM.estimateRateParameters(EmissionPDFCell, DeltaT, RateParametersGuess, SearchOptions, Verbose)
```

Estimates transition rate parameters via maximum likelihood.

**Parameters:**
- `EmissionPDFCell`: Cell array of emission matrices
- `DeltaT`: Cell array of time intervals between observations
- `RateParametersGuess`: Initial parameter guess
- `SearchOptions`: Optimization options structure
- `Verbose`: Verbosity level

**Returns:**
- `RateParameters`: MLE estimates
- `RateParametersSE`: Standard errors
- `LogLikelihood`: Log-likelihood at optimum

##### computeViterbiPath

```matlab
StateSequence = smi_stat.HMM.computeViterbiPath(StateSpace, InitialProbability, TransitionMatrix, EmissionMatrix)
```

Computes most likely state sequence using Viterbi algorithm.

**Parameters:**
- `StateSpace`: Vector of state identifiers
- `InitialProbability`: Initial state probabilities
- `TransitionMatrix`: State transition probability matrix (or series)
- `EmissionMatrix`: Emission probabilities for each observation

**Returns:**
- `StateSequence`: Most likely state at each time point

##### computeDimerDurations

```matlab
DimerDurations = smi_stat.HMM.computeDimerDurations(StateSequence, FrameNum, DimerState)
```

Extracts dimer event durations from state sequence.

**Parameters:**
- `StateSequence`: State sequence from Viterbi
- `FrameNum`: Frame numbers for each state
- `DimerState`: State ID corresponding to dimer (default: 1)

**Returns:**
- `DimerDurations`: Duration of each dimer event (in frames)

##### computePairSeparation

```matlab
TRArray = smi_stat.HMM.computePairSeparation(TRArray)
```

Computes inter-channel separation for trajectory pairs.

**Parameters:**
- `TRArray`: N × 2 trajectory pair array

**Returns:**
- `TRArray`: Updated with `Separation` field in each trajectory

### Usage Patterns

#### Pattern 1: Basic Dimer Analysis

```matlab
% Find dimer candidates from two-channel tracking
TR1 = load('Channel1_TR.mat');
TR2 = load('Channel2_TR.mat');

% Find candidates
TRArray = smi_stat.HMM.findDimerCandidates(TR1, TR2, 0.5, 5.0, 10, 500, 5);

% Setup HMM
H = smi_stat.HMM(TRArray, SMF);
H.ModelSpecifier = 'DF';
H.DimerSeparation = 0.5;
H.SaveDir = '/results/hmm';
H.ConditionLabel = 'Control';

% Analyze
[Rates, RatesSE] = H.performFullAnalysis();

% Interpret
k_on = Rates(1);   % Association rate
k_off = Rates(2);  % Dissociation rate
tau = 1 / k_off;   % Mean dimer lifetime
fprintf('Dimer lifetime: %.2f frames\n', tau);
```

#### Pattern 2: Batch Processing Multiple Conditions

```matlab
% Compare HMM results across conditions
Conditions = {'Control', 'Drug1', 'Drug2'};
Results = cell(length(Conditions), 1);

for ii = 1:length(Conditions)
    % Load data
    FileDir = fullfile('/data', Conditions{ii});

    % Find candidates
    [TRArray, ~, ~] = smi_stat.HMM.findDimerCandidatesFromFiles(...
        FileDir, {'*_Ch1_TR.mat', '*_Ch2_TR.mat'}, ...
        0.5, 5.0, 10, 500, 5, 1);

    % Analyze
    H = smi_stat.HMM(TRArray, SMF);
    H.ConditionLabel = Conditions{ii};
    H.SaveDir = fullfile('/results', Conditions{ii});
    [Rates, RatesSE] = H.performFullAnalysis();

    % Store
    Results{ii}.Condition = Conditions{ii};
    Results{ii}.Rates = Rates;
    Results{ii}.RatesSE = RatesSE;
    Results{ii}.k_off = Rates(2);
    Results{ii}.k_off_SE = RatesSE(2);
end

% Compare dissociation rates
figure;
k_offs = cellfun(@(x) x.k_off, Results);
k_offs_SE = cellfun(@(x) x.k_off_SE, Results);
errorbar(1:length(Conditions), k_offs, k_offs_SE, 'o');
set(gca, 'XTick', 1:length(Conditions), 'XTickLabel', Conditions);
ylabel('Dissociation rate (1/frame)');
title('Off-rate Comparison');
```

#### Pattern 3: Three-State Model (Dimer/Domain/Free)

```matlab
% Analyze with domain state
H = smi_stat.HMM(TRArray, SMF);
H.ModelSpecifier = 'DDF';  % Dimer/Domain/Free
H.DimerSeparation = 0.5;   % True dimer
H.DomainSeparation = 2.0;  % Domain colocalization
H.StateNames = {'Dimer', 'Domain', 'Free'};

[Rates, RatesSE] = H.performFullAnalysis();

% Interpret three-state results
% Rates = [k_D_to_free, k_free_to_D, k_Dom_to_free, k_free_to_Dom, k_D_to_Dom, k_Dom_to_D]
% (Exact ordering depends on implementation - check documentation)
fprintf('Dimer off-rate: %.3f\n', Rates(2));
fprintf('Domain off-rate: %.3f\n', Rates(4));
```

### See Also

- [SPT Workflow](../workflows/spt-tracking.md)
- `smi.SPT` - Tracking analysis
- `smi_core.ChannelRegistration` - Multi-channel alignment

**Citation:**
Nitta, C. F., Green, E. W., Jhamba, E. D., et al. (2021). EGFR transactivates RON to drive oncogenic crosstalk. eLife, 10:e63678. https://doi.org/10.7554/eLife.63678

---

## ChangeDetection

### Description

`ChangeDetection` identifies change points in integer-valued intensity time series using Bayesian model comparison. The method detects discrete times when the mean intensity of a Poisson process changes, making it ideal for analyzing photobleaching steps, blinking events, or clustering intensity traces in single molecule data.

### Key Concept

Change point detection asks: given a sequence of observations, when did the underlying intensity change? The Bayesian approach computes the probability that a change occurred at each time point by comparing models with vs. without a change point, using Bayes factors to determine significance.

### Class Definition

```matlab
classdef ChangeDetection < handle
```

### Constructor

```matlab
CD = smi_stat.ChangeDetection(Data, LogBayesThreshold)
```

**Parameters:**
- `Data`: Vector of integer intensity values (Poisson-distributed counts)
- `LogBayesThreshold`: Log Bayes factor threshold for accepting change points (positive scalar, e.g., 10-100)

**Returns:**
- `CD`: ChangeDetection object with analysis results

**Note:** Analysis runs automatically in constructor

**Example:**
```matlab
% Analyze photobleaching trace
IntensityTrace = [1000, 985, 1012, 502, 498, 510, 0, 0, 0];  % Step bleaching
CD = smi_stat.ChangeDetection(IntensityTrace, 10);

fprintf('Detected %d change points\n', CD.NchangePoints);
fprintf('Change points at frames: %s\n', num2str(CD.ChangePoints));
fprintf('Intensities: %s\n', num2str(CD.Intensity));
```

### Properties (SetAccess = protected)

#### Input

| Property | Type | Description |
|----------|------|-------------|
| `Data` | array | Input intensity time series |
| `LogBayesThreshold` | numeric | Threshold for accepting change points |
| `Nobservations` | integer | Length of data |

#### Output

| Property | Type | Description |
|----------|------|-------------|
| `NchangePoints` | integer | Number of detected change points |
| `ChangePoints` | array | Frame indices of change points |
| `Intensity` | array | Mean intensity in each segment |
| `IntensityModel` | array | Reconstructed step function intensity trace |
| `LogBayesFactors` | array | Log Bayes factors for accepted change points |
| `NrejectedChangePoints` | integer | Number of rejected candidate change points |
| `RejectedChangePoints` | array | Frame indices of rejected candidates |
| `RejectedLogBayesFactors` | array | Log Bayes factors for rejected points |

### Methods

##### setData

```matlab
obj.setData(Data, LogBayesThreshold)
```

Sets new data and re-runs analysis.

**Parameters:**
- `Data`: New intensity time series
- `LogBayesThreshold`: New threshold (optional)

**Example:**
```matlab
CD = smi_stat.ChangeDetection([], 10);  % Create empty
CD.setData(MyTrace, 20);  % Analyze with threshold 20
```

##### plotIntensityEstimate

```matlab
FigHandle = obj.plotIntensityEstimate()
```

Visualizes detected change points.

**Creates:**
- Subplot 1: Data with estimated step function overlay
- Subplot 2: Log probability of change point at each frame

**Returns:**
- `FigHandle`: Figure handle

**Example:**
```matlab
CD = smi_stat.ChangeDetection(IntensityTrace, 10);
CD.plotIntensityEstimate();
```

### Static Methods

##### simulate

```matlab
Data = smi_stat.ChangeDetection.simulate(NObservations, ChangePoints, Intensity)
```

Simulates intensity trace with specified change points.

**Parameters:**
- `NObservations`: Length of time series
- `ChangePoints`: Frame indices of change points
- `Intensity`: Mean intensities for each segment (length = NChangePoints + 1)

**Returns:**
- `Data`: Simulated Poisson-distributed intensity trace

**Example:**
```matlab
% Simulate three-level bleaching
NFrames = 100;
ChangePoints = [30, 70];
Intensities = [1000, 500, 100];
Data = smi_stat.ChangeDetection.simulate(NFrames, ChangePoints, Intensities);

% Test detection
CD = smi_stat.ChangeDetection(Data, 10);
fprintf('True changes at: %s\n', num2str(ChangePoints));
fprintf('Detected changes at: %s\n', num2str(CD.ChangePoints));
```

##### plotSimulatedEstimate

```matlab
[CD, FigHandle] = smi_stat.ChangeDetection.plotSimulatedEstimate(NObservations, ChangePoints, Intensity, LogBayesThreshold)
```

Simulates and analyzes data with visualization.

**Parameters:**
- Same as `simulate()` plus `LogBayesThreshold`

**Returns:**
- `CD`: ChangeDetection object
- `FigHandle`: Figure handle

##### plotRandSimulatedEstimate

```matlab
[CD, FigHandle] = smi_stat.ChangeDetection.plotRandSimulatedEstimate(NObservations, NChangePoints, meanIntensity, LogBayesThreshold)
```

Random simulation for testing algorithm performance.

**Parameters:**
- `NObservations`: Trace length
- `NChangePoints`: Number of change points to simulate
- `meanIntensity`: Mean intensity level
- `LogBayesThreshold`: Detection threshold

**Returns:**
- `CD`: ChangeDetection object
- `FigHandle`: Figure handle

##### fitClusterIntensity

```matlab
CD = smi_stat.ChangeDetection.fitClusterIntensity(Data, LogBayesThreshold, PlotFlag)
```

Convenience wrapper for analyzing cluster intensity traces.

**Parameters:**
- `Data`: Intensity time series (zeros at end are removed)
- `LogBayesThreshold`: Detection threshold
- `PlotFlag`: Generate plot (default: false)

**Returns:**
- `CD`: ChangeDetection object

**Use case:** Analyzing intensity of localization clusters over time

### Usage Patterns

#### Pattern 1: Photobleaching Step Detection

```matlab
% Analyze single molecule photobleaching
load('PhotonTrace.mat', 'Photons');  % Frame-by-frame photon counts

% Detect bleaching steps
CD = smi_stat.ChangeDetection(round(Photons), 50);

% Report results
fprintf('%d bleaching steps detected\n', CD.NchangePoints);
fprintf('Step locations (frames): %s\n', num2str(CD.ChangePoints));
fprintf('Intensities (photons): %s\n', num2str(round(CD.Intensity)));

% Visualize
CD.plotIntensityEstimate();

% Estimate number of fluorophores
NSteps = CD.NchangePoints;
fprintf('Estimated number of fluorophores: %d\n', NSteps);
```

#### Pattern 2: Optimizing Detection Threshold

```matlab
% Test different thresholds
IntensityTrace = load('Trace.mat');
Thresholds = [10, 20, 50, 100, 200];

figure;
for ii = 1:length(Thresholds)
    CD = smi_stat.ChangeDetection(IntensityTrace, Thresholds(ii));

    subplot(length(Thresholds), 1, ii);
    plot(IntensityTrace, 'k.');
    hold on;
    plot(CD.IntensityModel, 'r-', 'LineWidth', 2);
    ylabel('Intensity');
    title(sprintf('Threshold = %d, %d changes', ...
        Thresholds(ii), CD.NchangePoints));
end
xlabel('Frame');

% Rule of thumb: threshold ~ mean intensity / 10
OptimalThreshold = mean(IntensityTrace) / 10;
fprintf('Suggested threshold: %.1f\n', OptimalThreshold);
```

#### Pattern 3: Clustering Analysis

```matlab
% Analyze intensity of clustered localizations
load('ClusterData.mat', 'ClusterIntensities');  % N clusters × T frames

% Analyze each cluster
Results = struct();
for ii = 1:size(ClusterIntensities, 1)
    Trace = ClusterIntensities(ii, :);

    % Skip if insufficient data
    if sum(Trace > 0) < 5
        continue;
    end

    % Detect changes
    CD = smi_stat.ChangeDetection(round(Trace), 30);

    % Store
    Results(ii).ClusterID = ii;
    Results(ii).NChanges = CD.NchangePoints;
    Results(ii).ChangePoints = CD.ChangePoints;
    Results(ii).Intensities = CD.Intensity;
end

% Summary statistics
NChanges = [Results.NChanges];
fprintf('Mean changes per cluster: %.1f\n', mean(NChanges));
fprintf('Max changes: %d\n', max(NChanges));
```

#### Pattern 4: Testing Algorithm with Simulations

```matlab
% Validate detection performance
NTrials = 100;
TruePositives = 0;
FalsePositives = 0;
FalseNegatives = 0;

for ii = 1:NTrials
    % Simulate single-step bleaching
    TrueChangePoint = randi([30, 70]);
    Data = smi_stat.ChangeDetection.simulate(100, TrueChangePoint, [1000, 100]);

    % Detect
    CD = smi_stat.ChangeDetection(Data, 50);

    % Evaluate
    if CD.NchangePoints == 1
        DetectedChange = CD.ChangePoints;
        if abs(DetectedChange - TrueChangePoint) <= 2  % Within 2 frames
            TruePositives = TruePositives + 1;
        else
            FalsePositives = FalsePositives + 1;
        end
    elseif CD.NchangePoints == 0
        FalseNegatives = FalseNegatives + 1;
    else
        FalsePositives = FalsePositives + 1;
    end
end

% Report performance
fprintf('True positives: %d / %d (%.1f%%)\n', ...
    TruePositives, NTrials, 100*TruePositives/NTrials);
fprintf('False positives: %d\n', FalsePositives);
fprintf('False negatives: %d\n', FalseNegatives);
```

### See Also

- Photobleaching analysis guides
- Clustering analysis workflows

**Citation:**
Ensign, D. L., & Pande, V. S. (2010). Bayesian Detection of Intensity Changes in Single Molecule and Molecular Dynamics Trajectories. J. Phys. Chem. B, 114(1):280-292. https://doi.org/10.1021/jp906786b

---

## Supporting Functions

The +smi_stat namespace includes several supporting functions for image registration, statistical fitting, and colocalization analysis:

### Image Registration

##### findOffset

```matlab
[XOffset, YOffset, C] = smi_stat.findOffset(Image1, Image2, BoxSize, Threshold)
```

Finds translational offset between two images using cross-correlation.

**Parameters:**
- `Image1`, `Image2`: Images to register
- `BoxSize`: Size of correlation region
- `Threshold`: Correlation threshold

**Returns:**
- `XOffset`, `YOffset`: Pixel offsets
- `C`: Cross-correlation matrix

##### findOffsetIter

```matlab
[XOffset, YOffset] = smi_stat.findOffsetIter(Image1, Image2, NIterations)
```

Iterative offset finding with progressive refinement.

##### findCoordAffine

```matlab
Transform = smi_stat.findCoordAffine(Coords1, Coords2)
```

Computes affine transformation between coordinate sets.

**Parameters:**
- `Coords1`, `Coords2`: N × 2 coordinate arrays

**Returns:**
- `Transform`: Affine transformation structure

### Statistical Fitting

##### bootstrapFit

```matlab
[ParamMean, ParamStd] = smi_stat.bootstrapFit(X, Y, FitFunc, NBootstrap)
```

Bootstrap resampling for fit parameter uncertainty estimation.

**Parameters:**
- `X`, `Y`: Data to fit
- `FitFunc`: Function handle for fitting
- `NBootstrap`: Number of bootstrap samples

**Returns:**
- `ParamMean`: Mean parameter values
- `ParamStd`: Standard deviations

##### leastSquaresFit

```matlab
[Params, Residuals] = smi_stat.leastSquaresFit(X, Y, Model)
```

Least squares fitting to specified model.

### Colocalization Analysis

##### pearsonCorrCoef

```matlab
[R, P] = smi_stat.pearsonCorrCoef(Image1, Image2)
```

Computes Pearson correlation coefficient for image colocalization.

**Parameters:**
- `Image1`, `Image2`: Images to compare

**Returns:**
- `R`: Pearson correlation coefficient
- `P`: P-value for significance test

##### mandersSplitCoefs

```matlab
[M1, M2] = smi_stat.mandersSplitCoefs(Image1, Image2, Threshold1, Threshold2)
```

Computes Manders' colocalization coefficients.

**Parameters:**
- `Image1`, `Image2`: Images
- `Threshold1`, `Threshold2`: Intensity thresholds

**Returns:**
- `M1`: Manders' coefficient for Image1
- `M2`: Manders' coefficient for Image2

### Example: Multi-Channel Registration

```matlab
% Register two-channel images
load('Channel1_Image.mat', 'Im1');
load('Channel2_Image.mat', 'Im2');

% Find offset
[XOff, YOff, C] = smi_stat.findOffset(Im1, Im2, 64, 0.5);
fprintf('Offset: (%.2f, %.2f) pixels\n', XOff, YOff);

% Apply offset
Im2_Registered = smi_stat.shiftImage(Im2, XOff, YOff);

% Verify registration
[R_before, ~] = smi_stat.pearsonCorrCoef(Im1, Im2);
[R_after, ~] = smi_stat.pearsonCorrCoef(Im1, Im2_Registered);
fprintf('Correlation before: %.3f\n', R_before);
fprintf('Correlation after: %.3f\n', R_after);

% Compute colocalization
[M1, M2] = smi_stat.mandersSplitCoefs(Im1, Im2_Registered, 100, 100);
fprintf('Manders M1: %.3f, M2: %.3f\n', M1, M2);
```

---

## Performance Considerations

### DiffusionEstimator Speed

Different `FitTarget` options have different computational costs:

```matlab
% Fastest: MSD fitting (milliseconds per trajectory)
DE.FitTarget = 'MSD';

% Medium: CDF fitting (seconds per trajectory)
DE.FitTarget = 'CDFOfJumps';

% Slowest: MLE (seconds to minutes, especially with EstimateSEs = true)
DE.FitTarget = 'LikelihoodOfJumps';
```

For large datasets:
```matlab
% Disable SE estimation for speed
DE.EstimateSEs = false;

% Fit ensemble only
DE.FitIndividualTrajectories = false;
```

### HMM Memory Usage

For large trajectory arrays:
```matlab
% Process in batches
BatchSize = 100;
NTotal = size(TRArray, 1);
for ii = 1:BatchSize:NTotal
    EndIdx = min(ii + BatchSize - 1, NTotal);
    H = smi_stat.HMM(TRArray(ii:EndIdx, :), SMF);
    H.performFullAnalysis();
end
```

### ChangeDetection Scaling

Algorithm scales as O(N^2) for N observations:
```matlab
% For long traces, consider windowing
TraceLength = length(Data);
WindowSize = 200;
Overlap = 50;

for ii = 1:WindowSize-Overlap:TraceLength-WindowSize
    Window = Data(ii:ii+WindowSize-1);
    CD = smi_stat.ChangeDetection(Window, Threshold);
    % Process results...
end
```

---

## Troubleshooting

### DiffusionEstimator Issues

**Problem:** Negative diffusion estimates

**Solution:**
```matlab
% Check for tracking errors
Lengths = arrayfun(@(x) length(x.FrameNum), TR);
histogram(Lengths);  % Are trajectories long enough?

% Increase minimum trajectory length
GoodTraj = Lengths >= 10;
TR_Filtered = TR(GoodTraj);

% Check for localization precision
MeanSE = mean([TR.X_SE]);
fprintf('Mean localization SE: %.2f pixels\n', MeanSE);
% If too large, diffusion analysis may be unreliable
```

**Problem:** Multi-component fit fails

**Solution:**
```matlab
% Reduce number of components
DE.NComponents = 2;  % Start with 2 instead of 3

% Use MLE instead of CDF fitting
DE.FitTarget = 'LikelihoodOfJumps';

% Provide better initial guess
DE.FitParams0 = [0.05, 0.2, 0.5];  % [D1, D2, fraction1]
```

### HMM Issues

**Problem:** Rate parameter estimation fails

**Solution:**
```matlab
% Provide better initial guess
H.RateParametersGuess = [0.001; 0.01];  % Slower rates

% Increase maximum separation for more candidates
H.MaxSeparation = 8.0;

% Check if sufficient candidates
fprintf('%d dimer candidates found\n', size(H.TRArray, 1));
% Need at least ~20 for reliable estimation
```

**Problem:** State sequences look unrealistic

**Solution:**
```matlab
% Adjust dimer separation parameter
H.DimerSeparation = 0.3;  % Tighter threshold

% Include registration error
H.RegErrorInflation = 0.1;  % Add 0.1 pixels uncertainty

% Visualize emission PDFs to check model assumptions
H.plotSepDistribs([], H.TRArray, H.PDFHandles, {}, struct());
```

### ChangeDetection Issues

**Problem:** Too many false positives

**Solution:**
```matlab
% Increase threshold
CD = smi_stat.ChangeDetection(Data, 100);  % Higher threshold

% Check data quality
fprintf('Mean intensity: %.1f\n', mean(Data));
fprintf('Std: %.1f\n', std(Data));
% SNR = mean/sqrt(mean) for Poisson
% Low SNR → unreliable detection
```

**Problem:** Missing obvious change points

**Solution:**
```matlab
% Lower threshold
CD = smi_stat.ChangeDetection(Data, 10);

% Check for intensity near zero
MinIntensity = min(Data(Data > 0));
if MinIntensity < 10
    warning('Low intensity may affect detection sensitivity');
end
```

---

## Summary

The +smi_stat namespace provides advanced statistical tools for quantitative single molecule analysis:

**DiffusionEstimator:** Extracts diffusion coefficients from trajectories using MSD, CDF, or MLE methods. Supports heterogeneous populations and provides rigorous uncertainty estimates.

**HMM:** Detects state transitions in multi-color tracking data using Hidden Markov Models. Quantifies interaction kinetics and identifies transient binding events.

**ChangeDetection:** Identifies intensity change points using Bayesian model comparison. Useful for photobleaching analysis and temporal clustering studies.

**Supporting Functions:** Image registration, colocalization metrics, and statistical fitting utilities.

These tools transform trajectory and localization data into quantitative biophysical parameters, enabling rigorous characterization of molecular dynamics, interactions, and states.

---

## See Also

### Core Concepts
- [Architecture Overview](../core-concepts/architecture.md)
- [TR Structure Guide](../core-concepts/tr-structure.md)
- [Data Flow](../core-concepts/data-flow.md)

### Workflows
- [SPT Tracking Workflow](../workflows/spt-tracking.md)
- [Diffusion Analysis](../examples/tracking-diffusion.md)

### Related APIs
- [+smi Namespace](smi-namespace.md) - Main workflow classes
- [+smi_core Namespace](smi-core.md) - Core processing
- `smi.SPT` - Tracking pipeline

### External Resources
- MSD analysis theory and best practices
- Hidden Markov Model tutorials
- Bayesian change point detection literature
