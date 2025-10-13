---
title: "API Reference: +smi_sim Namespace"
category: "api-reference"
level: "beginner"
tags: ["api", "smi_sim", "simulation", "testing", "gaussblobs", "simsmlm", "simspt"]
prerequisites: ["../core-concepts/smd-structure.md", "../core-concepts/smf-structure.md"]
related: ["../how-to/simulate-data.md", "../workflows/smlm-analysis.md", "../workflows/spt-tracking.md"]
summary: "Complete API reference for the +smi_sim namespace covering simulation classes for generating test data, validating algorithms, and understanding SMLM/SPT analysis"
estimated_time: "20 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# API Reference: +smi_sim Namespace

## Purpose

The +smi_sim namespace provides simulation tools for generating synthetic single molecule localization microscopy (SMLM) and single particle tracking (SPT) data. These tools are essential for algorithm validation, testing analysis pipelines, parameter tuning, and understanding how imaging conditions affect results. This reference documents three simulation classes with increasing sophistication: GaussBlobs for simple test patterns, SimSMLM for realistic photophysics, and SimSPT for tracking with diffusion.

## Prerequisites

- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Familiarity with [SMF parameters](../core-concepts/smf-structure.md)
- Basic knowledge of SMLM and SPT concepts
- MATLAB object-oriented programming basics

## Overview

The +smi_sim namespace provides three complementary simulation tools:

**GaussBlobs:** Fast, simple Gaussian blob generation
- Static methods for rapid test data creation
- Random emitter placement with Poisson noise
- No photophysics modeling
- Ideal for quick algorithm validation

**SimSMLM:** Full SMLM photophysics simulation
- Realistic blinking and bleaching behavior
- Structured patterns (Siemen's star, tetramers, circles)
- Labeling efficiency modeling
- Complete data flow from true positions to noisy images

**SimSPT:** Single particle tracking simulation
- Brownian diffusion trajectories
- Photokinetics (blinking, bleaching)
- Oligomerization (dimer formation)
- Motion blur and measurement noise

All tools generate ground truth data for quantitative validation of localization and tracking algorithms.

## Class Reference

### GaussBlobs

**Purpose:** Fast generation of image stacks containing 2D Gaussian blobs with random placement and Poisson noise.

**Key Concept:** GaussBlobs provides the simplest path to test data. Use when you need speed over realism and don't require photophysics modeling.

**Class Definition:**
```matlab
classdef GaussBlobs
```

**Properties:**
None (all methods are static).

**Static Methods:**

#### gaussBlobROIStack

Generate ROI stack with single blob per image.

**Signature:**
```matlab
[Model, Data] = smi_sim.GaussBlobs.gaussBlobROIStack(SZ, SMD, VarianceIm, Covariance, PixType)
```

**Inputs:**
- `SZ`: Box size (pixels), scalar or [Y, X]
- `SMD`: SingleMoleculeData structure with positions and parameters
- `VarianceIm`: Read noise variance image (optional)
- `Covariance`: Include position covariance (optional, default false)
- `PixType`: Pixel data type (optional, default 'single')

**Outputs:**
- `Model`: Noiseless Gaussian blob stack (SZ × SZ × NBlobs)
- `Data`: Stack with Poisson noise added

**Usage:**
```matlab
% Create SMD with blob parameters
SMD = smi_core.SingleMoleculeData.createSMD();
SMD.X = [4.5; 3.2];  % X positions within box
SMD.Y = [3.8; 5.1];  % Y positions
SMD.Photons = [1000; 800];
SMD.Bg = [5; 5];
SMD.PSFSigma = [1.3; 1.3];

% Generate ROI stack
[Model, Data] = smi_sim.GaussBlobs.gaussBlobROIStack(7, SMD);
% Model: 7×7×2, noiseless
% Data: 7×7×2, with Poisson noise
```

#### gaussBlobImage

Generate full-frame images with multiple blobs per frame.

**Signature:**
```matlab
[Model, Data] = smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF, Bg, Density)
```

**Inputs:**
- `SMD`: SingleMoleculeData with blob positions and parameters
- `SMF`: SingleMoleculeFitting structure (for image size, ROI)
- `Bg`: Background level override (optional)
- `Density`: Density flag (optional)

**Outputs:**
- `Model`: Noiseless image stack
- `Data`: Stack with Poisson noise

**Usage:**
```matlab
% Setup SMF for image dimensions
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.DataROI = [1, 1, 128, 128];  % YStart, XStart, YEnd, XEnd

% Create SMD with multiple blobs
SMD = smi_core.SingleMoleculeData.createSMD();
SMD.NFrames = 100;
for nn = 1:100
    N = poissrnd(10);  % ~10 blobs per frame
    SMD.FrameNum = cat(1, SMD.FrameNum, nn*ones(N,1));
    SMD.X = cat(1, SMD.X, 128*rand(N,1));
    SMD.Y = cat(1, SMD.Y, 128*rand(N,1));
    SMD.Photons = cat(1, SMD.Photons, 1000*ones(N,1));
    SMD.Bg = cat(1, SMD.Bg, 5*ones(N,1));
    SMD.PSFSigma = cat(1, SMD.PSFSigma, 1.3*ones(N,1));
end

% Generate image stack
[Model, Data] = smi_sim.GaussBlobs.gaussBlobImage(SMD, SMF);
```

#### genRandomBlobImage

Generate random blob image stack (convenience wrapper).

**Signature:**
```matlab
BlobStack = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, Photons, PSFSigma, Bg)
```

**Inputs:**
- `SZ`: Image size (pixels), scalar or [Y, X] (Default: 256)
- `NFrames`: Number of frames (Default: 1000)
- `Rho`: Density (blobs/pixel) (Default: 0.001)
- `Photons`: Photons per blob (Default: 1000)
- `PSFSigma`: PSF sigma (pixels) (Default: 1.0)
- `Bg`: Background (photons/pixel) (Default: 0)

**Outputs:**
- `BlobStack`: Image stack (SZ × SZ × NFrames) with Poisson noise

**Usage:**
```matlab
% Quick test data generation
SZ = 128;
NFrames = 100;
Rho = 0.01;           % ~160 blobs per 128×128 frame
Photons = 1000;
PSFSigma = 1.3;
Bg = 5;

BlobStack = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, ...
    Photons, PSFSigma, Bg);

% Use with localization
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraGain = 1;  % Already in photons
SMF.Fitting.PSFSigma = 1.3;
LD = smi_core.LocalizeData(BlobStack, SMF);
[SMD, ~] = LD.genLocalizations();

fprintf('Generated %d frames, found %d localizations\n', ...
    NFrames, length(SMD.X));
```

**Typical Parameters:**

| Parameter | Sparse | Moderate | Dense | Use Case |
|-----------|--------|----------|-------|----------|
| Rho | 0.001 | 0.01 | 0.05 | Algorithm stress test |
| Photons | 500 | 1000 | 5000 | Bright/dim emitters |
| Bg | 1 | 5 | 50 | High/low SNR |
| PSFSigma | 1.0 | 1.3 | 1.5 | Diffraction limit variation |

**Performance Note:**

GaussBlobs is the fastest simulation tool. On typical hardware:
- 100 frames (128×128): < 1 second
- 1000 frames (256×256): ~5 seconds

**See Also:**
- [How to Simulate Data](../how-to/simulate-data.md#tool-1-gaussblobs-simple-random-patterns) - Usage examples
- [Basic Localization Example](../examples/basic-localization.md) - Uses GaussBlobs

---

### SimSMLM

**Purpose:** Comprehensive SMLM data generation with realistic photophysics including blinking, bleaching, and structured patterns.

**Key Concept:** SimSMLM models complete photophysics of fluorophores, producing realistic datasets for validating frame connection, testing structured illumination patterns, and understanding labeling efficiency effects. The data flow progresses through stages: SMD_True → SMD_Labeled → SMD_Model → Model/Data.

**Class Definition:**
```matlab
classdef SimSMLM < handle
```

**Properties:**

Imaging parameters:
```matlab
SZ                  % Linear image size (pixels) (Default: 256)
NFrames             % Frames per dataset (Default: 1000)
NDatasets           % Number of datasets (Default: 1)
EmissionRate        % Photons per frame when ON (Default: 1000)
Bg                  % Background (photons/pixel) (Default: 5)
PSFSigma            % PSF sigma (pixels) (Default: 1.3)
```

Pattern parameters:
```matlab
Rho                 % Fluorophore density (fluorophores/pixel) (Default: 30)
ZoomFactor          % Internal zoom factor (Default: 20)
```

Photophysics rates:
```matlab
K_OnToOff           % Turn-off rate (frames^-1) (Default: 1)
K_OffToOn           % Turn-on rate (frames^-1) (Default: 0.0005)
K_OnToBleach        % Bleaching rate (frames^-1) (Default: 0.2)
StartState          % Initial state: 'On' or 'Equib' (Default: 'Equib')
```

Labeling and efficiency:
```matlab
LabelingEfficiency  % Fraction labeled (0-1) (Default: 1)
```

Optimization:
```matlab
SparseFlag          % Use sparse matrices for large datasets (Default: false)
Verbose             % Verbosity level (Default: 1)
```

Data structures (outputs):
```matlab
SMD_True            % True emitter positions (from simStar, etc.)
SMD_Labeled         % After applying labeling efficiency
SMD_Model           % After generating blinks (photophysics)
NOnEvents           % Number of ON events per emitter
```

**Constructor:**
```matlab
obj = smi_sim.SimSMLM()
```

Creates SimSMLM object with default parameters.

**Example:**
```matlab
obj = smi_sim.SimSMLM();
obj.SZ = 128;
obj.NFrames = 1000;
obj.EmissionRate = 1000;
obj.Bg = 10;
obj.K_OnToOff = 1;
obj.K_OffToOn = 0.001;
obj.K_OnToBleach = 0.01;
```

**Key Methods:**

#### simStar

Generate Siemen's star pattern for resolution testing.

**Signature:**
```matlab
obj.simStar(NWings)
```

**Inputs:**
- `NWings`: Number of star wings (Default: 16)

**Outputs:**
Populates `obj.SMD_True` with star pattern coordinates.

**Algorithm:**
1. Generates uniform random emitters in circle
2. Removes emitters in alternating wedges
3. Wing length = SZ/3
4. Number of emitters follows Poisson(π × R² × Rho)

**Usage:**
```matlab
obj = smi_sim.SimSMLM();
obj.SZ = 256;
obj.Rho = 50;  % Total emitters, not density here
obj.NFrames = 1000;
obj.EmissionRate = 1000;
obj.K_OnToOff = 1;
obj.K_OffToOn = 0.001;
obj.K_OnToBleach = 0.01;

% Generate star
NWings = 16;
obj.simStar(NWings);

fprintf('Generated star with %d emitters\n', length(obj.SMD_True.X));

% Visualize true positions
figure;
plot(obj.SMD_True.X, obj.SMD_True.Y, 'k.', 'MarkerSize', 2);
axis equal; axis([0 obj.SZ 0 obj.SZ]);
title(sprintf('Siemens Star: %d Wings', NWings));
```

#### simCircle

Generate circular pattern.

**Signature:**
```matlab
obj.simCircle(Radius, NEmitters)
```

**Inputs:**
- `Radius`: Circle radius (pixels)
- `NEmitters`: Number of emitters on circle

**Outputs:**
Populates `obj.SMD_True`.

**Usage:**
```matlab
obj = smi_sim.SimSMLM();
obj.simCircle(30, 100);  % 100 emitters on radius=30 circle
```

#### simkTets

Generate k-mer patterns (dimers, trimers, tetramers, etc.).

**Signature:**
```matlab
obj.simkTets(k, radius_kTet)
```

**Inputs:**
- `k`: Order of oligomer (2=dimer, 3=trimer, 4=tetramer, etc.)
- `radius_kTet`: Size of k-mer (pixels)

**Outputs:**
Populates `obj.SMD_True` with k-mer patterns distributed across image.

**Usage:**
```matlab
obj = smi_sim.SimSMLM();
obj.SZ = 128;
obj.Rho = 20;  % Number of k-mers
k = 4;         % Tetramers
radius = 0.1;  % 0.1 pixel size

obj.simkTets(k, radius);

fprintf('Generated %d emitters in tetramers\n', length(obj.SMD_True.X));
```

#### kTet (Static)

Generate single k-mer at specified location.

**Signature:**
```matlab
SMD_True = smi_sim.SimSMLM.kTet(k, center, radius, startAngle)
```

**Inputs:**
- `k`: Oligomer order
- `center`: [X, Y] center position (pixels)
- `radius`: k-mer radius (pixels)
- `startAngle`: Rotation angle (radians) (Default: 0)

**Outputs:**
- `SMD_True`: SMD structure with k emitters

**Usage:**
```matlab
% Create tetramer at position [50, 50]
k = 4;
center = [50, 50];
radius = 0.2;
startAngle = pi/4;  % 45 degree rotation

SMD = smi_sim.SimSMLM.kTet(k, center, radius, startAngle);

% Visualize
figure;
plot(SMD.X, SMD.Y, 'ro', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
title('Tetramer');
```

#### applyLabelEffic

Apply labeling efficiency to remove unlabeled emitters.

**Signature:**
```matlab
obj.applyLabelEffic()
```

**Inputs:**
Uses `obj.SMD_True` and `obj.LabelingEfficiency`.

**Outputs:**
Populates `obj.SMD_Labeled` with subset of emitters that are labeled.

**Usage:**
```matlab
obj = smi_sim.SimSMLM();
obj.simStar(16);
fprintf('True emitters: %d\n', length(obj.SMD_True.X));

% Apply 70% labeling efficiency
obj.LabelingEfficiency = 0.7;
obj.applyLabelEffic();

fprintf('Labeled emitters: %d\n', length(obj.SMD_Labeled.X));
fprintf('Labeling rate: %.1f%%\n', ...
    100 * length(obj.SMD_Labeled.X) / length(obj.SMD_True.X));
```

#### genBlinks

Generate blinking time traces based on photophysics rates.

**Signature:**
```matlab
obj.genBlinks(StartState)
```

**Inputs:**
- `StartState`: 'On' or 'Equib' (equilibrium) (Default: 'Equib')

Uses `obj.SMD_Labeled`, photophysics rates, and frame parameters.

**Outputs:**
Populates `obj.SMD_Model` with blink events. Each blink becomes a localization.

**Key Fields Set:**
- `SMD_Model.X`, `SMD_Model.Y`: Position for each blink
- `SMD_Model.Photons`: Photon count for each blink
- `SMD_Model.FrameNum`: Frame where blink occurs
- `SMD_Model.ConnectID`: Links blinks to parent emitter in SMD_Labeled

**Algorithm:**
For each emitter in SMD_Labeled:
1. Start in specified state (On or random equilibrium)
2. Simulate state transitions (On ↔ Off, On → Bleached)
3. Record ON periods as blink events
4. Assign photons and frame numbers to each blink

**Usage:**
```matlab
obj = smi_sim.SimSMLM();
obj.SZ = 128;
obj.NFrames = 500;
obj.simStar(16);
obj.applyLabelEffic();

% Generate blinks starting in equilibrium
obj.genBlinks('Equib');

fprintf('Emitters: %d\n', length(obj.SMD_Labeled.X));
fprintf('Blinks: %d\n', length(obj.SMD_Model.X));
fprintf('Blinks per emitter: %.1f\n', ...
    length(obj.SMD_Model.X) / length(obj.SMD_Labeled.X));

% Access ground truth via ConnectID
true_X = obj.SMD_Labeled.X(obj.SMD_Model.ConnectID);
true_Y = obj.SMD_Labeled.Y(obj.SMD_Model.ConnectID);
% true_X(i) is the true position of blink i
```

#### genImageStack

Generate image stacks from blinking model.

**Signature:**
```matlab
[Model, Data] = obj.genImageStack()
```

**Inputs:**
Uses `obj.SMD_Model` and imaging parameters.

**Outputs:**
- `Model`: Noiseless image stack (SZ × SZ × TotalFrames)
- `Data`: Stack with Poisson noise added

Where TotalFrames = NDatasets × NFrames.

**Usage:**
```matlab
obj = smi_sim.SimSMLM();
obj.SZ = 128;
obj.NFrames = 100;
obj.NDatasets = 1;
obj.simStar(16);
obj.applyLabelEffic();
obj.genBlinks('Equib');

% Generate images
[Model, Data] = obj.genImageStack();

% Visualize
figure;
subplot(1,2,1);
imagesc(sum(Model, 3)); axis image; colormap gray;
title('Model (No Noise)');
colorbar;

subplot(1,2,2);
imagesc(sum(Data, 3)); axis image; colormap gray;
title('Data (Poisson Noise)');
colorbar;

% Use Data for localization
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraGain = 1;  % Already in photons
LD = smi_core.LocalizeData(Data, SMF);
[SMD, ~] = LD.genLocalizations();
```

#### genNoisySMD

Generate SMD with positional and intensity noise (alternative to genImageStack).

**Signature:**
```matlab
SMD_Data = obj.genNoisySMD(SMD_Model)
```

**Inputs:**
- `SMD_Model`: Model SMD structure (typically `obj.SMD_Model`)

**Outputs:**
- `SMD_Data`: SMD with Gaussian position noise and Poisson photon noise

**Usage:**
```matlab
obj = smi_sim.SimSMLM();
obj.simStar(16);
obj.applyLabelEffic();
obj.genBlinks('Equib');

% Generate noisy coordinates directly (no image)
SMD_Data = obj.genNoisySMD(obj.SMD_Model);

% Compare to ground truth
dx = SMD_Data.X - obj.SMD_Labeled.X(obj.SMD_Model.ConnectID);
dy = SMD_Data.Y - obj.SMD_Labeled.Y(obj.SMD_Model.ConnectID);
errors = sqrt(dx.^2 + dy.^2);

fprintf('Mean localization error: %.3f pixels\n', mean(errors));
```

**Complete Workflow Example:**

```matlab
% 1. Create and configure
obj = smi_sim.SimSMLM();
obj.SZ = 128;
obj.NFrames = 1000;
obj.NDatasets = 1;
obj.EmissionRate = 1000;
obj.Bg = 5;
obj.PSFSigma = 1.3;
obj.K_OnToOff = 1;
obj.K_OffToOn = 0.001;
obj.K_OnToBleach = 0.01;
obj.LabelingEfficiency = 1;

% 2. Generate pattern
NWings = 16;
obj.simStar(NWings);
fprintf('True emitters: %d\n', length(obj.SMD_True.X));

% 3. Apply labeling efficiency
obj.applyLabelEffic();
fprintf('Labeled emitters: %d\n', length(obj.SMD_Labeled.X));

% 4. Generate blinks
obj.genBlinks('Equib');
fprintf('Total blinks: %d\n', length(obj.SMD_Model.X));

% 5. Generate image stacks
[Model, Data] = obj.genImageStack();
fprintf('Image stack size: %d × %d × %d\n', size(Data));

% 6. Analyze ground truth
NEmitters = length(obj.SMD_Labeled.X);
NBlinks = length(obj.SMD_Model.X);
fprintf('Average blinks per emitter: %.1f\n', NBlinks / NEmitters);
```

**Tuning Photophysics:**

Long-lived ON states (PAINT):
```matlab
obj.K_OnToOff = 0.5;      % Stays ON longer
obj.K_OffToOn = 0.01;     % Turns ON more often
obj.K_OnToBleach = 0.001; % Bleaches slowly
```

Short blinks (dSTORM):
```matlab
obj.K_OnToOff = 2;        % Blinks off quickly
obj.K_OffToOn = 0.001;    % Rarely turns on
obj.K_OnToBleach = 0.01;  % Moderate bleaching
```

Always ON (photoactivation):
```matlab
obj.K_OnToOff = 0;        % Never turns off
obj.K_OnToBleach = 0.01;  % Only bleaching
obj.StartState = 'On';    % Start in ON state
```

**Performance Note:**

SimSMLM performance depends on number of emitters and frames:
- 100 emitters, 1000 frames: ~2 seconds
- 1000 emitters, 1000 frames: ~20 seconds

For very large simulations (>10,000 emitters), set `obj.SparseFlag = true` to use sparse matrices.

**See Also:**
- [How to Simulate Data](../how-to/simulate-data.md#tool-2-simsmlm-realistic-smlm-with-photophysics)
- [SMLM Workflow](../workflows/smlm-analysis.md)

---

### SimSPT

**Purpose:** Single particle tracking simulation with Brownian diffusion, photokinetics, oligomerization, and realistic measurement effects.

**Key Concept:** SimSPT generates complete tracking datasets including ground truth trajectories with diffusion, blinking/bleaching kinetics, motion blur, and measurement noise. Essential for validating tracking algorithms and understanding how experimental parameters affect tracking performance.

**Class Definition:**
```matlab
classdef SimSPT < handle
```

**Properties:**

Main parameter structure:
```matlab
SimParams           % Structure of all simulation parameters
```

Trajectory structures (outputs, SetAccess = protected):
```matlab
TrajStructSubTrue      % True sub-frame trajectories
TrajStructSubLabeled   % After applying labeling efficiency
TrajStructSubModel     % After photokinetics (blinking, bleaching)
TrajStructModel        % Motion-blurred frame trajectories
TrajStruct             % With measurement noise (final output)
```

Dependent properties (auto-computed):
```matlab
SMD                 % Trajectory data in SMD format
TR                  % Trajectory data in TR format
```

**Constructor:**
```matlab
obj = smi_sim.SimSPT(SimParams)
```

**Inputs:**
- `SimParams`: Parameter structure (optional, uses defaults if empty)

**Example:**
```matlab
% Create with defaults
obj = smi_sim.SimSPT();

% Or with custom parameters
SP = smi_sim.SimSPT.defineDefaultParams();
SP.NFrames = 200;
SP.D = 0.5;
obj = smi_sim.SimSPT(SP);
```

**Static Method: defineDefaultParams**

Create parameter structure with default values.

**Signature:**
```matlab
ParamStruct = smi_sim.SimSPT.defineDefaultParams()
```

**Outputs:**
- `ParamStruct`: Structure with all simulation parameters

**Parameters Reference:**

Simulation geometry:
```matlab
ParticleDensity     % Density (particles/pixel²) (Default: 0.005)
NFrames             % Number of frames (Default: 100)
FrameSize           % Image size [Y, X] (pixels) (Default: [32, 32])
SubframeDensity     % Subframes per frame (Default: 1)
InitialDensityMask  % Binary mask for initial placement (Default: ones(FrameSize))
```

Imaging parameters:
```matlab
PSFSigma            % PSF sigma (pixels) (Default: 1.3)
Intensity           % Photons per frame (Default: 1000)
MinIntensity        % Minimum photons (Default: 50)
Bg                  % Background (photons/pixel) (Default: 5)
```

Diffusion and motion:
```matlab
D                   % Diffusion coefficient(s) (pixels²/frame) (Default: 0.1)
                    % Can be scalar or array for mixed populations
BoundaryCondition   % 'Periodic', 'Reflecting', or 'Free' (Default: 'Periodic')
```

Photokinetics:
```matlab
KOnToBleach         % Bleaching rate (frames^-1) (Default: 0.01)
KOnToOff            % Turn-off rate (frames^-1) (Default: 0.2)
KOffToOn            % Turn-on rate (frames^-1) (Default: 0.8)
PMiss               % Missed detection probability (Default: 0.01)
```

Oligomerization:
```matlab
InteractionDistance % Distance for dimerization (pixels) (Default: 0.5)
InteractionProb     % Probability of dimer formation (Default: 0.5)
KDisconnect         % Dimer dissociation rate (frames^-1) (Default: 0.1)
RestrictToDimers    % Only allow dimers, no higher oligomers (Default: true)
```

Labeling:
```matlab
LabelingEfficiency  % Fraction of particles labeled (Default: 1)
```

**Usage:**
```matlab
% Get defaults
SP = smi_sim.SimSPT.defineDefaultParams();

% Customize
SP.ParticleDensity = 0.01;
SP.NFrames = 200;
SP.FrameSize = [64, 64];
SP.D = 0.3;              % Diffusion coefficient
SP.Intensity = 2000;     % Bright particles
SP.Bg = 5;
SP.KOnToBleach = 0.005;  % Slow bleaching
SP.KOnToOff = 0.1;       % Occasional blinking
SP.KOffToOn = 0.9;       % Comes back quickly

% Create simulation
obj = smi_sim.SimSPT(SP);
```

**Key Method: createSimulation**

Execute complete simulation pipeline.

**Signature:**
```matlab
obj.createSimulation()
```

**Inputs:**
Uses `obj.SimParams`.

**Outputs:**
Populates all trajectory structures. Access results via:
- `obj.SMD` - Localization format
- `obj.TR` - Trajectory format
- `obj.TrajStruct` - Full trajectory structure

**Pipeline Steps:**
1. `simTrajectories`: Generate Brownian diffusion trajectories
2. `applyLabelingEfficiency`: Remove unlabeled particles
3. `simEmitterKinetics`: Apply blinking and bleaching
4. `applyMeasurementModel`: Add motion blur and noise

**Usage:**
```matlab
% Setup parameters
SP = smi_sim.SimSPT.defineDefaultParams();
SP.ParticleDensity = 0.01;
SP.NFrames = 150;
SP.FrameSize = [64, 64];
SP.D = 0.2;
SP.Intensity = 1500;

% Create and run simulation
obj = smi_sim.SimSPT(SP);
obj.createSimulation();

% Access results
SMD = obj.SMD;
TR = obj.TR;

fprintf('Generated %d trajectories\n', length(TR));
fprintf('Total localizations: %d\n', length(SMD.X));

% Visualize trajectories
figure; hold on;
for i = 1:min(20, length(TR))  % Plot first 20
    plot(TR(i).X, TR(i).Y, '-', 'LineWidth', 1.5);
end
axis equal;
xlim([0, SP.FrameSize(2)]);
ylim([0, SP.FrameSize(1)]);
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Simulated Trajectories');
```

**Static Method: simTrajectories**

Generate Brownian diffusion trajectories.

**Signature:**
```matlab
TrajStruct = smi_sim.SimSPT.simTrajectories(SimParams)
```

**Inputs:**
- `SimParams`: Parameter structure

**Outputs:**
- `TrajStruct`: Trajectory structure with coordinates, frames, IDs

**Usage:**
```matlab
SP = smi_sim.SimSPT.defineDefaultParams();
SP.D = [0.1, 0.5];  % Two diffusion populations
TrajStruct = smi_sim.SimSPT.simTrajectories(SP);

fprintf('Generated %d trajectories\n', length(unique(TrajStruct.ID)));
```

**Static Method: simEmitterKinetics**

Apply photokinetics to trajectories.

**Signature:**
```matlab
TrajStruct = smi_sim.SimSPT.simEmitterKinetics(TrajStruct, SimParams)
```

**Inputs:**
- `TrajStruct`: Input trajectory structure
- `SimParams`: Parameters with kinetic rates

**Outputs:**
- `TrajStruct`: Trajectories with blinking/bleaching applied

**Usage:**
```matlab
% After generating trajectories
SP = smi_sim.SimSPT.defineDefaultParams();
TrajStruct = smi_sim.SimSPT.simTrajectories(SP);

% Apply kinetics
TrajStruct = smi_sim.SimSPT.simEmitterKinetics(TrajStruct, SP);

fprintf('Frames before kinetics: %d\n', max(TrajStruct.FrameNum));
fprintf('Frames after kinetics: %d\n', sum(TrajStruct.IsVisible));
```

**Static Method: applyMeasurementModel**

Apply motion blur and measurement noise.

**Signature:**
```matlab
[TrajStruct, TrajStructModel] = smi_sim.SimSPT.applyMeasurementModel(TrajStruct, SimParams)
```

**Inputs:**
- `TrajStruct`: Sub-frame trajectories
- `SimParams`: Measurement parameters

**Outputs:**
- `TrajStruct`: Frame-averaged with noise
- `TrajStructModel`: Frame-averaged without noise

**Usage:**
```matlab
SP = smi_sim.SimSPT.defineDefaultParams();
SP.SubframeDensity = 10;  % 10 subframes for motion blur

% Generate and apply
obj = smi_sim.SimSPT(SP);
obj.createSimulation();

% Access pre-noise model
TrajModel = obj.TrajStructModel;
TrajNoisy = obj.TrajStruct;

% Compare
fprintf('Model positions: %.3f, %.3f\n', TrajModel.X(1), TrajModel.Y(1));
fprintf('Noisy positions: %.3f, %.3f\n', TrajNoisy.X(1), TrajNoisy.Y(1));
```

**Static Method: convertTrajToSMD**

Convert trajectory structure to SMD format.

**Signature:**
```matlab
SMD = smi_sim.SimSPT.convertTrajToSMD(TrajStruct, SimParams)
```

**Inputs:**
- `TrajStruct`: Trajectory structure
- `SimParams`: Parameters for metadata

**Outputs:**
- `SMD`: SingleMoleculeData structure

**Note:** Automatically called via `obj.SMD` dependent property.

**Complete Example with Analysis:**

```matlab
% 1. Configure simulation
SP = smi_sim.SimSPT.defineDefaultParams();
SP.ParticleDensity = 0.01;
SP.NFrames = 200;
SP.FrameSize = [64, 64];
SP.D = 0.3;
SP.Intensity = 2000;
SP.Bg = 5;
SP.PSFSigma = 1.3;
SP.KOnToBleach = 0.005;  % Slow bleaching
SP.KOnToOff = 0.1;       % Occasional blinking
SP.KOffToOn = 0.9;       % Quick recovery
SP.PMiss = 0.01;         % 1% missed detections

% 2. Run simulation
obj = smi_sim.SimSPT(SP);
obj.createSimulation();

% 3. Get results
SMD = obj.SMD;
TR = obj.TR;

fprintf('Simulation complete:\n');
fprintf('  Trajectories: %d\n', length(TR));
fprintf('  Total localizations: %d\n', length(SMD.X));

% 4. Analyze trajectory statistics
Lengths = smi_core.TrackingResults.computeTrajLengths(TR);
fprintf('  Median trajectory length: %d frames\n', median(Lengths));

% 5. Visualize trajectories
figure;
subplot(1,2,1);
hold on;
for i = 1:length(TR)
    plot(TR(i).X, TR(i).Y, '-', 'LineWidth', 1);
end
axis equal;
xlim([0, SP.FrameSize(2)]);
ylim([0, SP.FrameSize(1)]);
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Simulated Trajectories');

subplot(1,2,2);
histogram(Lengths, 20);
xlabel('Trajectory Length (frames)');
ylabel('Count');
title('Trajectory Length Distribution');

% 6. Calculate mean squared displacement
MSD = zeros(length(TR), 1);
for i = 1:length(TR)
    dx = diff(TR(i).X);
    dy = diff(TR(i).Y);
    MSD(i) = mean(dx.^2 + dy.^2);
end

% 7. Estimate diffusion coefficient
D_estimated = MSD / 4;  % For 2D diffusion, MSD = 4*D*dt (dt=1 frame)
fprintf('  True D: %.3f pixels²/frame\n', SP.D);
fprintf('  Estimated D: %.3f pixels²/frame\n', mean(D_estimated));
```

**Boundary Conditions:**

Periodic (particles wrap around edges):
```matlab
SP.BoundaryCondition = 'Periodic';
% Useful for avoiding edge effects
```

Reflecting (particles bounce off edges):
```matlab
SP.BoundaryCondition = 'Reflecting';
% Models physical barriers
```

Free (particles can leave):
```matlab
SP.BoundaryCondition = 'Free';
% Realistic for unconstrained diffusion
```

**Mixed Diffusion Populations:**

```matlab
SP.D = [0.05, 0.2, 0.5];  % Slow, medium, fast

obj = smi_sim.SimSPT(SP);
obj.createSimulation();

% Each trajectory randomly assigned one D value
% Useful for heterogeneous populations
```

**Oligomerization (Dimer Formation):**

```matlab
SP.InteractionDistance = 0.5;  % Within 0.5 pixels
SP.InteractionProb = 0.5;      % 50% chance when close
SP.KDisconnect = 0.1;          % Dissociation rate
SP.RestrictToDimers = true;    % Only pairs, no trimers

obj = smi_sim.SimSPT(SP);
obj.createSimulation();

% Trajectories will show binding/unbinding events
```

**Motion Blur Control:**

```matlab
% No motion blur (fast acquisition)
SP.SubframeDensity = 1;

% Moderate motion blur
SP.SubframeDensity = 10;

% High accuracy motion blur (slower simulation)
SP.SubframeDensity = 20;
```

Higher SubframeDensity = more accurate motion blur but slower simulation.

**Performance Note:**

SimSPT performance scales with:
- Number of particles (ParticleDensity × FrameSize)
- Number of frames
- SubframeDensity (motion blur accuracy)

Typical performance:
- 30 particles, 100 frames, SubframeDensity=1: ~2 seconds
- 100 particles, 200 frames, SubframeDensity=10: ~15 seconds

**See Also:**
- [How to Simulate Data](../how-to/simulate-data.md#tool-3-simspt-single-particle-tracking)
- [SPT Workflow](../workflows/spt-tracking.md)

---

## Common Usage Patterns

### Quick Algorithm Validation

Use GaussBlobs for fast testing:

```matlab
% Generate test data
Data = smi_sim.GaussBlobs.genRandomBlobImage(128, 100, 0.01, 1000, 1.3, 5);

% Test your algorithm
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraGain = 1;
LD = smi_core.LocalizeData(Data, SMF);
[SMD, ~] = LD.genLocalizations();

% Validate
fprintf('Found %d localizations in 100 frames\n', length(SMD.X));
```

### Frame Connection Testing

Use SimSMLM with realistic blinking:

```matlab
obj = smi_sim.SimSMLM();
obj.SZ = 128;
obj.NFrames = 500;
obj.K_OnToOff = 1;
obj.K_OffToOn = 0.01;     % Multiple blinks per emitter
obj.K_OnToBleach = 0.005; % Long-lived emitters
obj.simStar(16);
obj.applyLabelEffic();
obj.genBlinks('Equib');
[~, Data] = obj.genImageStack();

% Test frame connection
SMF = smi_core.SingleMoleculeFitting();
SMF.FrameConnection.On = true;
LD = smi_core.LocalizeData(Data, SMF);
[SMD, ~] = LD.genLocalizations();
FC = smi_core.FrameConnection(SMD, SMF);
[SMDCombined, SMD] = FC.performFrameConnection();

% Validate against ground truth
fprintf('True emitters: %d\n', length(obj.SMD_Labeled.X));
fprintf('Connected emitters: %d\n', length(SMDCombined.X));
```

### Tracking Algorithm Benchmark

Use SimSPT for tracking validation:

```matlab
% Generate ground truth
SP = smi_sim.SimSPT.defineDefaultParams();
SP.ParticleDensity = 0.01;
SP.NFrames = 200;
SP.D = 0.2;
SP.PMiss = 0.05;  % Some missed detections

obj = smi_sim.SimSPT(SP);
obj.createSimulation();
TR_true = obj.TR;

% Run your tracking algorithm on obj.SMD
% Compare results to TR_true
```

### Parameter Sensitivity Study

Test how SNR affects precision:

```matlab
Backgrounds = [1, 5, 10, 50];
Precisions = zeros(size(Backgrounds));

for i = 1:length(Backgrounds)
    % Simulate
    Data = smi_sim.GaussBlobs.genRandomBlobImage(128, 50, 0.005, ...
        1000, 1.3, Backgrounds(i));

    % Localize
    SMF = smi_core.SingleMoleculeFitting();
    SMF.Data.CameraGain = 1;
    LD = smi_core.LocalizeData(Data, SMF, 0);
    [SMD, ~] = LD.genLocalizations();

    % Measure precision
    Precisions(i) = mean(sqrt(SMD.X_SE.^2 + SMD.Y_SE.^2));
end

figure;
plot(Backgrounds, Precisions, 'o-', 'LineWidth', 2);
xlabel('Background (photons/pixel)');
ylabel('Mean Precision (pixels)');
title('Precision vs Background');
```

---

## Comparison: When to Use Each Tool

| Feature | GaussBlobs | SimSMLM | SimSPT |
|---------|------------|---------|--------|
| **Speed** | Fastest (< 1s) | Fast (~2-20s) | Moderate (~2-60s) |
| **Photophysics** | None | Full (blink/bleach) | Kinetics only |
| **Patterns** | Random only | Structured (star, k-mers) | Trajectories |
| **Ground Truth** | Manual | Automatic | Automatic |
| **Motion** | Static | Static | Brownian diffusion |
| **Use Case** | Quick tests | SMLM validation | SPT validation |

**Choose GaussBlobs when:**
- You need test data immediately
- Photophysics don't matter
- Testing localization algorithms

**Choose SimSMLM when:**
- Testing frame connection
- Validating structured patterns
- Studying blinking/bleaching effects
- Need realistic SMLM data

**Choose SimSPT when:**
- Testing tracking algorithms
- Studying diffusion processes
- Validating trajectory analysis
- Need motion blur modeling

---

## Validation Workflows

### Localization Error Analysis

```matlab
% Simulate with GaussBlobs
SZ = 128;
NFrames = 50;
Rho = 0.01;
Photons = 1000;
PSFSigma = 1.3;
Bg = 5;

% Create ground truth
SMD_true = smi_core.SingleMoleculeData.createSMD();
SMD_true.NFrames = NFrames;
for nn = 1:NFrames
    N = poissrnd(Rho * SZ * SZ);
    SMD_true.FrameNum = cat(1, SMD_true.FrameNum, nn*ones(N,1));
    SMD_true.X = cat(1, SMD_true.X, SZ*rand(N,1));
    SMD_true.Y = cat(1, SMD_true.Y, SZ*rand(N,1));
end

% Generate data
Data = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, Photons, PSFSigma, Bg);

% Localize
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraGain = 1;
LD = smi_core.LocalizeData(Data, SMF);
[SMD, ~] = LD.genLocalizations();

% Match to ground truth
max_distance = 2;  % pixels
matched = 0;
errors = [];

for i = 1:length(SMD_true.X)
    same_frame = (SMD.FrameNum == SMD_true.FrameNum(i));
    if sum(same_frame) > 0
        dx = SMD.X(same_frame) - SMD_true.X(i);
        dy = SMD.Y(same_frame) - SMD_true.Y(i);
        dist = sqrt(dx.^2 + dy.^2);
        [min_dist, ~] = min(dist);
        if min_dist < max_distance
            matched = matched + 1;
            errors = [errors; min_dist];
        end
    end
end

fprintf('Detection rate: %.1f%%\n', 100 * matched / length(SMD_true.X));
fprintf('Mean error: %.3f pixels\n', mean(errors));
```

### Photon Recovery Validation

```matlab
% Simulate with SimSMLM
obj = smi_sim.SimSMLM();
obj.SZ = 128;
obj.NFrames = 100;
obj.EmissionRate = 1000;  % True photons
obj.simStar(16);
obj.applyLabelEffic();
obj.genBlinks('On');
[~, Data] = obj.genImageStack();

% Localize
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraGain = 1;
LD = smi_core.LocalizeData(Data, SMF);
[SMD, ~] = LD.genLocalizations();

% Compare
true_photons = obj.EmissionRate;
estimated_photons = SMD.Photons;

figure;
histogram(estimated_photons, 50);
xline(true_photons, 'r--', 'LineWidth', 2);
xlabel('Detected Photons');
ylabel('Count');
title(sprintf('Photon Recovery (True: %d)', true_photons));

fprintf('Mean detected: %.1f photons\n', mean(estimated_photons));
fprintf('Bias: %.1f photons (%.1f%%)\n', ...
    mean(estimated_photons) - true_photons, ...
    100 * (mean(estimated_photons) - true_photons) / true_photons);
```

---

## Performance Optimization

### For Large Simulations

When generating very large datasets:

```matlab
% SimSMLM with many emitters
obj = smi_sim.SimSMLM();
obj.Rho = 1000;  % Many emitters
obj.SparseFlag = true;  % Use sparse matrices
obj.simStar(16);
```

### Memory-Efficient Batch Processing

```matlab
% Generate data in chunks
NFramesTotal = 10000;
ChunkSize = 1000;

for chunk = 1:NFramesTotal/ChunkSize
    obj = smi_sim.SimSMLM();
    obj.NFrames = ChunkSize;
    obj.simStar(16);
    obj.applyLabelEffic();
    obj.genBlinks('Equib');
    [~, Data] = obj.genImageStack();

    % Process chunk
    % ...

    % Save results
    save(sprintf('chunk_%d.mat', chunk), 'Data');
    clear obj Data
end
```

---

## Troubleshooting

### Simulation Too Slow

**Problem:** SimSMLM or SimSPT taking too long

**Solutions:**
```matlab
% Reduce frames
SP.NFrames = 100;  % Instead of 1000

% Reduce particle density
SP.ParticleDensity = 0.005;  % Instead of 0.01

% Reduce motion blur accuracy
SP.SubframeDensity = 1;  % Instead of 10

% For SimSMLM, use sparse matrices
obj.SparseFlag = true;
```

### Unrealistic Data

**Problem:** Simulated data doesn't match experiments

**Solutions:**
```matlab
% Match your experimental parameters
obj.EmissionRate = 800;     % Measured photons/frame
obj.Bg = 15;                % Measured background
obj.PSFSigma = 1.4;         % Measured PSF width
obj.K_OnToOff = 2.5;        % Fitted blink rate
obj.K_OffToOn = 0.002;      % Fitted turn-on rate
obj.K_OnToBleach = 0.015;   % Fitted bleach rate
```

### Ground Truth Mismatch

**Problem:** Can't match localizations to ground truth

**Solution:**
```matlab
% For SimSMLM, use ConnectID
true_X = obj.SMD_Labeled.X(obj.SMD_Model.ConnectID);
true_Y = obj.SMD_Labeled.Y(obj.SMD_Model.ConnectID);
% Now true_X(i) matches SMD_Model.X(i)

% For SimSPT, use TR format
TR = obj.TR;
% Each TR(i) is one trajectory
```

---

## See Also

### How-To Guides
- [How to Simulate Data](../how-to/simulate-data.md) - Detailed usage examples
- [How to Localize Molecules](../how-to/localize-molecules.md) - Using simulated data

### Workflows
- [SMLM Analysis Workflow](../workflows/smlm-analysis.md) - Full analysis pipeline
- [SPT Tracking Workflow](../workflows/spt-tracking.md) - Tracking pipeline

### Related APIs
- [+smi_core Namespace](smi-core.md) - Core processing classes
- [SMD Structure Guide](../core-concepts/smd-structure.md) - Result format
- [TR Structure Guide](../core-concepts/tr-structure.md) - Trajectory format

---

## Summary

The +smi_sim namespace provides three complementary tools for simulation:

**GaussBlobs** offers the fastest path to test data with static methods for generating random Gaussian blobs. Use for quick algorithm validation when photophysics don't matter.

**SimSMLM** provides comprehensive SMLM simulation with realistic blinking, bleaching, and structured patterns. Essential for testing frame connection, understanding photophysics effects, and generating publication-quality test datasets.

**SimSPT** simulates complete tracking experiments with Brownian diffusion, photokinetics, oligomerization, and motion blur. Critical for validating tracking algorithms and understanding how experimental parameters affect trajectory analysis.

All three tools generate ground truth data enabling quantitative validation of smite algorithms. Choose the tool that matches your validation needs, balancing simulation speed against realism requirements.
