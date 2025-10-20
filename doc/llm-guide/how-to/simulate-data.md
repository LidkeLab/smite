---
title: "How to Simulate Data"
category: "how-to"
level: "beginner"
tags: ["simulation", "testing", "validation", "gaussblobs", "simsmlm", "simspt"]
prerequisites: ["../getting-started/installation.md", "../core-concepts/smd-structure.md"]
related: ["../examples/basic-localization.md", "../workflows/smlm-analysis.md", "../workflows/spt-tracking.md"]
summary: "Generate simulated SMLM and SPT data with realistic noise and photophysics for testing and validation"
estimated_time: "15 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# How to Simulate Data

## Purpose

Simulated data is essential for testing localization algorithms, validating tracking pipelines, and understanding how imaging parameters affect analysis results. smite provides three simulation tools with increasing sophistication: GaussBlobs for simple test patterns, SimSMLM for realistic blinking behavior, and SimSPT for single particle tracking with diffusion. This guide shows how to use each tool and generate appropriate simulations for different validation tasks.

## Prerequisites

- smite installed and working
- Basic understanding of [SMD structure](../core-concepts/smd-structure.md)
- Familiarity with SMLM/SPT concepts

## Overview

smite's simulation tools:

1. **GaussBlobs**: Fast generation of Gaussian blobs with Poisson noise
   - Use for: Quick testing, algorithm validation, simple data
   - Speed: Very fast (< 1 second for 100 frames)

2. **SimSMLM**: Full photophysics simulation with blinking and bleaching
   - Use for: Realistic SMLM data, testing frame connection, structured patterns
   - Speed: Fast to moderate (seconds for 1000 frames)

3. **SimSPT**: Tracking simulation with diffusion and photokinetics
   - Use for: SPT algorithm validation, trajectory analysis, complex behaviors
   - Speed: Moderate (seconds to minutes)

All tools generate ground truth data for quantitative validation.

## Tool 1: GaussBlobs (Simple Random Patterns)

### When to Use GaussBlobs

Use GaussBlobs when you need:
- Quick test data for localization algorithms
- Random emitter patterns
- Simple validation without photophysics
- Fast generation of large datasets

### Basic Usage

```matlab
% Generate random blob image
SZ = 128;              % Image size (pixels)
NFrames = 100;         % Number of frames
Rho = 0.01;            % Density (emitters/pixel)
Photons = 1000;        % Photons per emitter
PSFSigma = 1.3;        % PSF sigma (pixels)
Bg = 5;                % Background (photons/pixel)

imageStack = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, ...
    Photons, PSFSigma, Bg);

% imageStack is SZ × SZ × NFrames with Poisson noise
fprintf('Generated %d frames\n', NFrames);
fprintf('Expected emitters per frame: %.0f\n', Rho * SZ * SZ);
```

### Creating Ground Truth

GaussBlobs generates images but doesn't return ground truth positions. Create your own SMD for validation:

```matlab
% Generate ground truth SMD matching the simulation
SMD_true = smi_core.SingleMoleculeData.createSMD();
SMD_true.NFrames = NFrames;
SMD_true.NDatasets = 1;

for nn = 1:NFrames
    N = poissrnd(Rho * SZ * SZ);  % Random number of emitters
    SMD_true.FrameNum = cat(1, SMD_true.FrameNum, nn*ones(N,1));
    SMD_true.X = cat(1, SMD_true.X, SZ*rand(N,1));  % Uniform random positions
    SMD_true.Y = cat(1, SMD_true.Y, SZ*rand(N,1));
    SMD_true.Photons = cat(1, SMD_true.Photons, Photons*ones(N,1));
    SMD_true.Bg = cat(1, SMD_true.Bg, Bg*ones(N,1));
    SMD_true.PSFSigma = cat(1, SMD_true.PSFSigma, PSFSigma*ones(N,1));
end

fprintf('Ground truth: %d total emitters\n', length(SMD_true.X));
```

### Parameter Effects

**Density (Rho)**: Controls emitter overlap
```matlab
% Sparse (no overlap)
Rho = 0.001;  % ~16 emitters per 128×128 frame

% Dense (significant overlap)
Rho = 0.01;   % ~160 emitters per frame

% Very dense (stress test)
Rho = 0.05;   % ~800 emitters per frame
```

**Photons**: Controls signal strength and localization precision
```matlab
% Dim emitters (challenging)
Photons = 500;

% Bright emitters (easy)
Photons = 5000;
```

**Background**: Affects signal-to-noise ratio
```matlab
% Low background (high SNR)
Bg = 1;

% High background (low SNR, more challenging)
Bg = 50;
```

### Complete Example

```matlab
% Simulate and localize with GaussBlobs
SZ = 128;
NFrames = 50;
Rho = 0.01;
Photons = 1000;
PSFSigma = 1.3;
Bg = 5;

% Generate data
imageStack = smi_sim.GaussBlobs.genRandomBlobImage(SZ, NFrames, Rho, ...
    Photons, PSFSigma, Bg);

% Configure localization
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraGain = 1;     % Data already in photons
SMF.Data.CameraOffset = 0;
SMF.Data.PixelSize = 0.1;    % 100 nm pixels
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 300;
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Localize
LD = smi_core.LocalizeData(imageStack, SMF);
SMD = LD.genLocalizations();

fprintf('Localizations found: %d\n', length(SMD.X));
fprintf('Median precision: %.1f nm\n', ...
    median(SMD.X_SE) * SMF.Data.PixelSize * 1000);
```

## Tool 2: SimSMLM (Realistic SMLM with Photophysics)

### When to Use SimSMLM

Use SimSMLM when you need:
- Realistic blinking and bleaching behavior
- Structured patterns (Siemen's star, custom shapes)
- Testing frame connection algorithms
- Labeling efficiency effects
- Ground truth for multi-dataset simulations

### Basic Setup

```matlab
% Create SimSMLM object
obj = smi_sim.SimSMLM();

% Configure imaging parameters
obj.SZ = 128;                  % Image size (pixels)
obj.NFrames = 1000;            % Frames per dataset
obj.NDatasets = 1;             % Number of datasets
obj.EmissionRate = 1000;       % Photons per frame when ON
obj.Bg = 5;                    % Background (photons/pixel)
obj.PSFSigma = 1.3;            % PSF sigma (pixels)

% Configure photophysics
obj.K_OnToOff = 1;             % Turn-off rate (frames^-1)
obj.K_OffToOn = 0.001;         % Turn-on rate (frames^-1)
obj.K_OnToBleach = 0.01;       % Bleaching rate (frames^-1)
obj.StartState = 'Equib';      % Start in equilibrium ('Equib' or 'On')

% Labeling efficiency
obj.LabelingEfficiency = 1;    % Fraction of labeled emitters (0-1)
```

### Generating Structured Patterns

**Siemen's Star** (for testing resolution):
```matlab
obj = smi_sim.SimSMLM();
obj.SZ = 256;
obj.Rho = 50;           % Total emitters
obj.NFrames = 1000;
obj.EmissionRate = 1000;
obj.Bg = 10;
obj.K_OnToOff = 1;
obj.K_OffToOn = 0.001;
obj.K_OnToBleach = 0.01;

% Generate Siemen's star pattern
NWings = 16;
obj.simStar(NWings);

% obj.SMD_True now contains ground truth positions
fprintf('Generated star with %d wings\n', NWings);
fprintf('True emitters: %d\n', length(obj.SMD_True.X));
```

**Custom Patterns** (tetramers, etc.):
```matlab
% Generate k-mer patterns
k = 4;                  % Tetramer
radius = 0.1;           % Size (pixels)
startAngle = 0;         % Orientation

% Create tetramers at specific locations
centers = [50, 50; 100, 100; 150, 150];  % Center positions

SMD_combined = smi_core.SingleMoleculeData.createSMD();
for i = 1:size(centers, 1)
    SMD_temp = smi_sim.SimSMLM.kTet(k, centers(i,:), radius, startAngle);
    SMD_combined.X = [SMD_combined.X; SMD_temp.X];
    SMD_combined.Y = [SMD_combined.Y; SMD_temp.Y];
end

% Use in SimSMLM
obj = smi_sim.SimSMLM();
obj.SMD_True = SMD_combined;
```

### Complete Workflow

```matlab
% 1. Create and configure SimSMLM
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

% 2. Generate pattern (Siemen's star)
NWings = 16;
obj.simStar(NWings);

% 3. Apply labeling efficiency
obj.applyLabelEffic();
fprintf('After labeling: %d emitters\n', length(obj.SMD_Labeled.X));

% 4. Generate blinking traces
obj.genBlinks('Equib');
fprintf('Blinks generated: %d localizations\n', length(obj.SMD_Model.X));

% 5. Generate image stacks
[Model, Data] = obj.genImageStack();
% Model = noiseless images
% Data = with Poisson noise

fprintf('Generated %d frames\n', size(Data, 3));

% 6. Visualize
figure;
subplot(1,2,1);
imagesc(sum(Model, 3)); axis image; colormap gray;
title('Model (No Noise)');

subplot(1,2,2);
imagesc(sum(Data, 3)); axis image; colormap gray;
title('Data (With Noise)');
```

### Understanding SimSMLM Data Flow

SimSMLM generates data in stages, each stored as a property:

1. **SMD_True**: True emitter positions (from simStar, etc.)
2. **SMD_Labeled**: After applying labeling efficiency
3. **SMD_Model**: After generating blinking events (photophysics)
4. **Model/Data**: Image stacks (noiseless/noisy)

Access ground truth at any stage:
```matlab
% Compare detected blinks to true positions
ConnectID = obj.SMD_Model.ConnectID;
true_X = obj.SMD_Labeled.X(ConnectID);
true_Y = obj.SMD_Labeled.Y(ConnectID);

% true_X(i) is the true X position of blink i in SMD_Model
```

### Tuning Photophysics

**Long-lived ON states** (easy to localize):
```matlab
obj.K_OnToOff = 0.5;      % Stays ON longer
obj.K_OffToOn = 0.01;     % Turns ON more often
obj.K_OnToBleach = 0.001; % Bleaches slowly
```

**Short blinks** (challenging, like dSTORM):
```matlab
obj.K_OnToOff = 2;        % Blinks off quickly
obj.K_OffToOn = 0.001;    % Rarely turns on
obj.K_OnToBleach = 0.01;  % Moderate bleaching
```

**No photophysics** (always ON until bleach):
```matlab
obj.K_OnToOff = 0;        % Never turns off
obj.K_OffToOn = Inf;      % (doesn't matter)
obj.K_OnToBleach = 0.01;  % Only bleaching
obj.StartState = 'On';    % Start in ON state
```

## Tool 3: SimSPT (Single Particle Tracking)

### When to Use SimSPT

Use SimSPT when you need:
- Trajectories with realistic diffusion
- Testing tracking algorithms
- Photokinetics (blinking) during tracking
- Oligomerization (dimer formation)
- Motion blur effects
- Ground truth tracks

### Basic Setup

```matlab
% Create SimSPT with default parameters
obj = smi_sim.SimSPT();

% Access default parameters
SP = obj.SimParams;
fprintf('Default diffusion: %.3f pixels^2/frame\n', SP.D);
fprintf('Default frames: %d\n', SP.NFrames);

% Modify parameters
SP.ParticleDensity = 0.01;    % Particles per pixel^2
SP.NFrames = 100;             % Number of frames
SP.FrameSize = [64, 64];      % Image size (pixels)
SP.D = 0.5;                   % Diffusion coefficient (pixels^2/frame)
SP.Intensity = 2000;          % Photons per frame
SP.Bg = 10;                   % Background (photons/pixel)
SP.PSFSigma = 1.3;            % PSF sigma (pixels)

obj.SimParams = SP;
```

### Generating Tracking Data

```matlab
% Create and configure SimSPT
SP = smi_sim.SimSPT.defineDefaultParams();
SP.ParticleDensity = 0.005;   % Density
SP.NFrames = 100;             % Frames
SP.FrameSize = [64, 64];      % Size
SP.D = 0.2;                   % Diffusion coefficient
SP.Intensity = 1500;          % Photons
SP.Bg = 5;                    % Background
SP.PSFSigma = 1.3;            % PSF sigma

% Photokinetics parameters
SP.KOnToBleach = 0.01;        % Bleaching rate (frames^-1)
SP.KOnToOff = 0.2;            % Turn-off rate (frames^-1)
SP.KOffToOn = 0.8;            % Turn-on rate (frames^-1)
SP.PMiss = 0.01;              % Missed detection probability

% Create simulation
obj = smi_sim.SimSPT(SP);
obj.createSimulation();

% Access results
SMD = obj.SMD;                % Localization data
TR = obj.TR;                  % Tracking results

fprintf('Generated %d trajectories\n', max(TR.ConnectID));
fprintf('Total localizations: %d\n', length(SMD.X));
```

### Boundary Conditions

```matlab
SP = smi_sim.SimSPT.defineDefaultParams();

% Periodic boundaries (particles wrap around)
SP.BoundaryCondition = 'Periodic';

% Reflecting boundaries (particles bounce)
SP.BoundaryCondition = 'Reflecting';

% Free boundaries (particles leave and don't return)
SP.BoundaryCondition = 'Free';

obj = smi_sim.SimSPT(SP);
obj.createSimulation();
```

### Multiple Diffusion States

```matlab
SP = smi_sim.SimSPT.defineDefaultParams();
SP.ParticleDensity = 0.01;
SP.NFrames = 200;

% Mix of diffusion coefficients
SP.D = [0.05, 0.2, 0.5];  % Slow, medium, fast

% Each particle randomly assigned one D value
obj = smi_sim.SimSPT(SP);
obj.createSimulation();

% Analyze diffusion
TR = obj.TR;
for i = 1:max(TR.ConnectID)
    idx = TR.ConnectID == i;
    coords = [TR.X(idx), TR.Y(idx)];
    % Calculate MSD, etc.
end
```

### Oligomerization (Dimer Formation)

```matlab
SP = smi_sim.SimSPT.defineDefaultParams();
SP.ParticleDensity = 0.01;
SP.NFrames = 200;
SP.D = 0.2;

% Enable dimerization
SP.InteractionDistance = 0.5;  % Distance for dimer formation (pixels)
SP.InteractionProb = 0.5;      % Probability of forming dimer
SP.KDisconnect = 0.1;          % Dimer dissociation rate (frames^-1)
SP.RestrictToDimers = true;    % Only dimers (no higher oligomers)

obj = smi_sim.SimSPT(SP);
obj.createSimulation();

TR = obj.TR;
fprintf('Simulated oligomerization\n');
```

### Motion Blur

SimSPT automatically simulates motion blur using subframes:

```matlab
SP = smi_sim.SimSPT.defineDefaultParams();
SP.NFrames = 100;
SP.SubframeDensity = 10;  % 10 subframes per frame

% Higher SubframeDensity = more accurate motion blur
% Typical values: 1 (no blur) to 20 (high accuracy)

obj = smi_sim.SimSPT(SP);
obj.createSimulation();

% TrajStructModel contains pre-blur coordinates
% TrajStruct contains motion-blurred coordinates
```

### Complete SPT Example

```matlab
% Simulate diffusing particles with tracking
SP = smi_sim.SimSPT.defineDefaultParams();
SP.ParticleDensity = 0.01;
SP.NFrames = 150;
SP.FrameSize = [64, 64];
SP.D = 0.3;
SP.Intensity = 2000;
SP.Bg = 5;
SP.PSFSigma = 1.3;
SP.KOnToBleach = 0.005;  % Slow bleaching
SP.KOnToOff = 0.1;       % Occasional blinking
SP.KOffToOn = 0.9;       % Comes back quickly

% Create simulation
obj = smi_sim.SimSPT(SP);
obj.createSimulation();

% Get tracking results
TR = obj.TR;
SMD = obj.SMD;

% Visualize trajectories
figure;
hold on;
for i = 1:max(TR.ConnectID)
    idx = TR.ConnectID == i;
    plot(TR.X(idx), TR.Y(idx), '-', 'LineWidth', 1.5);
end
axis equal;
xlim([0, SP.FrameSize(1)]);
ylim([0, SP.FrameSize(2)]);
title(sprintf('%d Trajectories', max(TR.ConnectID)));
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Calculate trajectory statistics
track_lengths = zeros(max(TR.ConnectID), 1);
for i = 1:max(TR.ConnectID)
    track_lengths(i) = sum(TR.ConnectID == i);
end

fprintf('Mean track length: %.1f frames\n', mean(track_lengths));
fprintf('Median track length: %.1f frames\n', median(track_lengths));
```

## Validation with Ground Truth

### Comparing Localizations to Truth

```matlab
% After localization, match to ground truth
max_distance = 2;  % pixels

matched = 0;
distances = [];

for i = 1:length(SMD_true.X)
    % Find nearest localization in same frame
    same_frame = (SMD.FrameNum == SMD_true.FrameNum(i));

    if sum(same_frame) > 0
        dx = SMD.X(same_frame) - SMD_true.X(i);
        dy = SMD.Y(same_frame) - SMD_true.Y(i);
        dist = sqrt(dx.^2 + dy.^2);

        [min_dist, ~] = min(dist);

        if min_dist < max_distance
            matched = matched + 1;
            distances = [distances; min_dist];
        end
    end
end

detection_rate = 100 * matched / length(SMD_true.X);
fprintf('Detection rate: %.1f%%\n', detection_rate);
fprintf('Mean error: %.3f pixels (%.1f nm)\n', ...
    mean(distances), mean(distances) * pixel_size * 1000);
```

### Validating Photon Estimates

```matlab
% Compare estimated photons to simulated photons
true_photons = 1000;  % From simulation

figure;
histogram(SMD.Photons, 50);
xline(true_photons, 'r--', 'LineWidth', 2, 'Label', 'True');
xlabel('Detected Photons');
ylabel('Count');
title('Photon Recovery');

bias = mean(SMD.Photons) - true_photons;
fprintf('Photon bias: %.1f (%.1f%%)\n', bias, 100*bias/true_photons);
```

## Comparison: Which Tool to Use?

| Feature | GaussBlobs | SimSMLM | SimSPT |
|---------|------------|---------|--------|
| **Speed** | Fastest | Fast | Moderate |
| **Complexity** | Simple | Moderate | Complex |
| **Photophysics** | No | Yes (full) | Yes (kinetics) |
| **Patterns** | Random only | Structured | Trajectories |
| **Tracking** | No | No | Yes (diffusion) |
| **Ground truth** | Manual | Automatic | Automatic |
| **Use case** | Quick tests | SMLM validation | SPT validation |

**Choose GaussBlobs when**: You need fast, simple data for localization testing

**Choose SimSMLM when**: You need realistic blinking/bleaching or structured patterns

**Choose SimSPT when**: You're testing tracking algorithms or studying diffusion

## Common Scenarios

### Scenario 1: Algorithm Development

Testing a new localization algorithm:
```matlab
% Use GaussBlobs for speed
imageStack = smi_sim.GaussBlobs.genRandomBlobImage(128, 100, 0.01, 1000, 1.3, 5);
% Test your algorithm on imageStack
```

### Scenario 2: Frame Connection Validation

Testing frame connection with realistic blinking:
```matlab
obj = smi_sim.SimSMLM();
obj.SZ = 128;
obj.NFrames = 500;
obj.K_OnToOff = 1;
obj.K_OffToOn = 0.01;  % Multiple blinks from same emitter
obj.simStar(16);
obj.applyLabelEffic();
obj.genBlinks('Equib');
[~, Data] = obj.genImageStack();
% Use Data with frame connection enabled
```

### Scenario 3: Tracking Algorithm Benchmark

Comparing tracking algorithms:
```matlab
SP = smi_sim.SimSPT.defineDefaultParams();
SP.D = [0.1, 0.5, 1.0];  % Mixed diffusion
SP.NFrames = 200;
SP.PMiss = 0.05;         % Some missed detections
obj = smi_sim.SimSPT(SP);
obj.createSimulation();
TR_true = obj.TR;
% Compare your tracking results to TR_true
```

## Tips and Best Practices

### Choosing Parameters

**For realistic SMLM data**:
- Photons: 500-5000 (depends on fluorophore/laser)
- Background: 5-50 (depends on autofluorescence)
- PSFSigma: 1.0-1.5 pixels (depends on NA and pixel size)
- Density: 0.001-0.01 emitters/pixel (depends on activation)

**For SPT data**:
- D: 0.05-2 pixels^2/frame (depends on particle size and viscosity)
- Intensity: 1000-10000 photons/frame
- Frame rate consideration: Higher D needs higher frame rate

### Performance Optimization

```matlab
% For large simulations, use sparse matrices in SimSMLM
obj = smi_sim.SimSMLM();
obj.SparseFlag = true;  % Slower per operation, but handles more emitters

% For SimSPT, reduce SubframeDensity if speed matters
SP.SubframeDensity = 1;  % Faster, less accurate motion blur
```

### Saving Simulations

```matlab
% Save simulation for reproducibility
save('simulation_params.mat', 'obj', 'SP');

% Save image data
save('simulated_data.mat', 'Data', 'SMD_true', '-v7.3');

% Or save as HDF5 for compatibility
smi_core.LoadData.saveRawData(Data, 'simulated_data.h5');
```

## See Also

- [Basic Localization Example](../examples/basic-localization.md) - Uses GaussBlobs
- [SMLM Workflow](../workflows/smlm-analysis.md) - Analysis pipeline
- [SPT Tracking Workflow](../workflows/spt-tracking.md) - Tracking analysis
- MATLAB/+smi_sim/@GaussBlobs/GaussBlobs.m - Source code
- MATLAB/+smi_sim/@SimSMLM/SimSMLM.m - Source code
- MATLAB/+smi_sim/@SimSPT/SimSPT.m - Source code
