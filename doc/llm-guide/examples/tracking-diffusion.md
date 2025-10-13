---
title: "SPT Tracking and Diffusion Analysis"
category: "examples"
level: "intermediate"
tags: ["example", "spt", "tracking", "diffusion", "msd", "tutorial"]
prerequisites: ["../workflows/spt-tracking.md", "../core-concepts/tr-structure.md"]
related: ["basic-localization.md", "../how-to/analyze-diffusion.md"]
summary: "Complete working example of single particle tracking from simulation through diffusion coefficient estimation"
estimated_time: "15 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# SPT Tracking and Diffusion Analysis

## Purpose

This example demonstrates a complete single particle tracking (SPT) workflow from start to finish. You'll generate simulated diffusing particles with a known diffusion coefficient, track them through time, compute mean square displacement (MSD) curves, estimate the diffusion coefficient from the tracked trajectories, and compare the results to the ground truth. This is the canonical workflow for analyzing particle motion and measuring molecular dynamics.

## Prerequisites

- smite installed and working
- Basic MATLAB knowledge
- Understanding of [SPT workflow](../workflows/spt-tracking.md)
- Understanding of [TR structure](../core-concepts/tr-structure.md)
- 10-15 minutes to run the code

## Overview

This example covers:
1. Simulating diffusing particles with known diffusion coefficient
2. Converting simulation to realistic image data
3. Localizing particles in each frame
4. Tracking particles through time (frame-to-frame + gap closing)
5. Computing MSD curves from trajectories
6. Estimating diffusion coefficients
7. Comparing estimated vs true diffusion
8. Visualizing trajectories and motion statistics

The key difference from basic localization is that particles move between frames and must be tracked to reconstruct their paths (trajectories).

## Complete Working Code

Copy and run this complete example:

```matlab
%% SPT Tracking and Diffusion Analysis Example
% Demonstrates complete workflow from simulation to diffusion estimation

%% Step 1: Generate Simulated Diffusing Particles
fprintf('=== Generating Simulated SPT Data ===\n');

% Set random seed for reproducibility
rng(42);

% Simulation parameters
SimParams = struct();
SimParams.FrameSize = [128, 128];         % Image size (pixels)
SimParams.NFrames = 200;                  % Number of frames
SimParams.ParticleDensity = 0.002;        % Particles per pixel^2 (~50 particles)
SimParams.D = 0.15;                       % True diffusion (pixels^2/frame)
SimParams.KOffToOn = 0.9;                 % Blinking: off -> on probability
SimParams.KOnToOff = 0.05;                % Blinking: on -> off probability
SimParams.Intensity = 1000;               % Photons per particle
SimParams.Bg = 10;                        % Background photons/pixel

% Physical parameters for interpretation
PixelSize = 0.1;      % micrometers (100 nm pixels)
FrameRate = 100;      % Hz (10 ms per frame)

% Convert diffusion to physical units for display
D_physical = SimParams.D * (PixelSize^2) * FrameRate;  % μm^2/s

fprintf('Simulation setup:\n');
fprintf('  Frame size: %d × %d pixels\n', SimParams.FrameSize(1), SimParams.FrameSize(2));
fprintf('  Number of frames: %d\n', SimParams.NFrames);
fprintf('  Particle density: %.4f particles/pixel^2\n', SimParams.ParticleDensity);
fprintf('  Expected particles: ~%d\n', ...
    round(SimParams.ParticleDensity * prod(SimParams.FrameSize)));
fprintf('  True diffusion: %.3f pixels^2/frame = %.3f μm^2/s\n', ...
    SimParams.D, D_physical);
fprintf('  Photons per particle: %d\n', SimParams.Intensity);
fprintf('  Background: %d photons/pixel\n', SimParams.Bg);

% Create SPT simulation
SPTSim = smi_sim.SimSPT;
SPTSim.SimParams.FrameSize = SimParams.FrameSize;
SPTSim.SimParams.NFrames = SimParams.NFrames;
SPTSim.SimParams.ParticleDensity = SimParams.ParticleDensity;
SPTSim.SimParams.D = SimParams.D;
SPTSim.SimParams.KOffToOn = SimParams.KOffToOn;
SPTSim.SimParams.KOnToOff = SimParams.KOnToOff;
SPTSim.SimParams.Intensity = SimParams.Intensity;
SPTSim.SimParams.Bg = SimParams.Bg;

% Generate ground truth trajectories
SPTSim.createSimulation();

% Store ground truth TR
TR_true = SPTSim.TR;
fprintf('Generated %d true trajectories\n', length(TR_true));

%% Step 2: Convert to Image Stack
fprintf('\n=== Converting to Image Stack ===\n');

% Create SMF for image generation
SMF_sim = smi_core.SingleMoleculeFitting();
SMF_sim.Data.DataROI = [1, 1, SPTSim.SMD.YSize, SPTSim.SMD.XSize];
SMF_sim.Fitting.PSFSigma = 1.3;  % PSF width (pixels)

% Generate realistic images with PSF and noise
[~, imageStack] = smi_sim.GaussBlobs.gaussBlobImage(SPTSim.SMD, SMF_sim);

fprintf('Generated image stack: %d × %d × %d frames\n', ...
    size(imageStack, 1), size(imageStack, 2), size(imageStack, 3));
fprintf('Total localizations in ground truth: %d\n', length(SPTSim.SMD.X));

%% Step 3: Configure SMF for Tracking
fprintf('\n=== Configuring Tracking Parameters ===\n');

% Create SMF for analysis
SMF = smi_core.SingleMoleculeFitting();

% Camera parameters
SMF.Data.CameraGain = 1;          % Data already in photons
SMF.Data.CameraOffset = 0;
SMF.Data.PixelSize = PixelSize;   % micrometers
SMF.Data.FrameRate = FrameRate;   % Hz

% Localization parameters
SMF.BoxFinding.BoxSize = 7;              % 7×7 pixel boxes
SMF.BoxFinding.MinPhotons = 300;         % Detection threshold
SMF.BoxFinding.BoxOverlap = 2;           % Allow some overlap

SMF.Fitting.PSFSigma = 1.3;              % Match simulation PSF
SMF.Fitting.FitType = 'XYNB';            % Fit X, Y, photons, background
SMF.Fitting.Iterations = 20;             % MLE iterations

% Thresholding (more permissive for tracking than SMLM)
SMF.Thresholding.On = true;
SMF.Thresholding.AutoThreshLogL = true;  % Auto log-likelihood threshold
SMF.Thresholding.MaxXY_SE = 0.3;         % Max position uncertainty (pixels)
SMF.Thresholding.MinPhotons = 200;       % Min photons after fitting
SMF.Thresholding.MinPValue = 0.01;       % Min fit quality

% Tracking parameters (key for SPT!)
SMF.Tracking.Method = 'CostMatrix';
SMF.Tracking.D = 0.15;                   % Expected diffusion (match true value)
SMF.Tracking.TrajwiseD = true;           % Estimate D per trajectory
SMF.Tracking.MaxDistFF = 5;              % Max frame-to-frame distance (pixels)
SMF.Tracking.MaxDistGC = 10;             % Max gap closing distance (pixels)
SMF.Tracking.MaxFrameGap = 5;            % Max frames to bridge
SMF.Tracking.MinTrackLength = 10;        % Keep trajectories ≥10 frames
SMF.Tracking.K_on = SimParams.KOffToOn;  % Match simulation blinking
SMF.Tracking.K_off = SimParams.KOnToOff;

fprintf('Tracking parameters:\n');
fprintf('  Expected D: %.3f pixels^2/frame\n', SMF.Tracking.D);
fprintf('  Max frame-to-frame distance: %d pixels\n', SMF.Tracking.MaxDistFF);
fprintf('  Max gap closing distance: %d pixels\n', SMF.Tracking.MaxDistGC);
fprintf('  Max frame gap: %d frames\n', SMF.Tracking.MaxFrameGap);
fprintf('  Min trajectory length: %d frames\n', SMF.Tracking.MinTrackLength);

%% Step 4: Localize Particles
fprintf('\n=== Localizing Particles ===\n');

% Create LocalizeData object
LD = smi_core.LocalizeData(imageStack, SMF);
LD.Verbose = 1;  % Show progress

% Run localization
tic;
SMD = LD.genLocalizations();
localization_time = toc;

fprintf('Localization complete in %.2f seconds\n', localization_time);
fprintf('Found %d localizations\n', length(SMD.X));

% Filter by quality thresholds
passed = sum(SMD.ThreshFlag == 0);
fprintf('Passed quality filters: %d (%.1f%%)\n', ...
    passed, 100*passed/length(SMD.X));

%% Step 5: Track Particles
fprintf('\n=== Tracking Particles ===\n');

% Create SPT object
SPT = smi.SPT(SMF, false);  % false = don't open GUI
SPT.GenerateMovies = false;  % Skip movies for speed
SPT.GeneratePlots = false;   % Skip plots for now

% Localization already done, so set SMD directly
SPT.SMD = SMD;

% Run tracking (frame-to-frame + gap closing)
tic;
fprintf('Running frame-to-frame tracking...\n');
SMD = smi.SPT.genTrajFF(SMD, SMF, [], SPT.NonLinkMarker);

fprintf('Running gap closing...\n');
SMD = smi.SPT.genTrajGC(SMD, SMF, [], SPT.NonLinkMarker, SPT.UseSparseMatrices);
tracking_time = toc;

fprintf('Tracking complete in %.2f seconds\n', tracking_time);

% Convert SMD to TR (trajectories)
TR = smi_core.TrackingResults.convertSMDToTR(SMD);

% Filter by minimum length
TR = smi_core.TrackingResults.threshTrajLength(TR, SMF.Tracking.MinTrackLength);

fprintf('Found %d trajectories (after filtering)\n', length(TR));

% Trajectory statistics
traj_lengths = arrayfun(@(x) length(x.X), TR);
fprintf('Trajectory lengths: %d to %d frames\n', min(traj_lengths), max(traj_lengths));
fprintf('Mean trajectory length: %.1f frames\n', mean(traj_lengths));
fprintf('Median trajectory length: %d frames\n', median(traj_lengths));

%% Step 6: Compute MSD Curves
fprintf('\n=== Computing Mean Square Displacements ===\n');

% Configure DiffusionEstimator
DE = smi_stat.DiffusionEstimator();
DE.TR = TR;
DE.SMF = SMF;
DE.FitTarget = 'MSD';                   % Fit MSD curves
DE.FitMethod = 'WeightedLS';            % Weighted least squares
DE.FitIndividualTrajectories = true;    % Fit each trajectory
DE.FrameLagRange = [1, 10];             % Compute MSD for lag 1-10 frames
DE.NFitPoints = 5;                      % Fit first 5 points
DE.NDimensions = 2;                     % 2D tracking
DE.UnitFlag = true;                     % Output in physical units
DE.Verbose = 1;

% Compute MSD and estimate diffusion
tic;
DiffusionStruct = DE.estimateDiffusionConstant();
diffusion_time = toc;

fprintf('Diffusion estimation complete in %.2f seconds\n', diffusion_time);

%% Step 7: Analyze Results
fprintf('\n=== Diffusion Analysis Results ===\n');

% Extract trajectory-wise results
D_trajectory = DiffusionStruct(1).DiffusionConstant;  % μm^2/s
D_trajectory_SE = DiffusionStruct(1).DiffusionConstantSE;

% Extract ensemble results
D_ensemble = DiffusionStruct(2).DiffusionConstant;
D_ensemble_SE = DiffusionStruct(2).DiffusionConstantSE;

% Statistics
D_mean = mean(D_trajectory, 'omitnan');
D_median = median(D_trajectory, 'omitnan');
D_std = std(D_trajectory, 'omitnan');

fprintf('True diffusion coefficient: %.3f μm^2/s\n', D_physical);
fprintf('\nEstimated diffusion (ensemble): %.3f ± %.3f μm^2/s\n', ...
    D_ensemble, D_ensemble_SE);
fprintf('Estimated diffusion (trajectory mean): %.3f ± %.3f μm^2/s\n', ...
    D_mean, D_std);
fprintf('Estimated diffusion (trajectory median): %.3f μm^2/s\n', D_median);

% Compute error
ensemble_error = 100 * abs(D_ensemble - D_physical) / D_physical;
mean_error = 100 * abs(D_mean - D_physical) / D_physical;

fprintf('\nAccuracy:\n');
fprintf('  Ensemble error: %.1f%%\n', ensemble_error);
fprintf('  Mean trajectory error: %.1f%%\n', mean_error);

%% Step 8: Visualizations
fprintf('\n=== Creating Visualizations ===\n');

figure('Name', 'SPT Tracking and Diffusion Analysis', ...
    'Position', [50, 50, 1600, 1000]);

% Panel 1: Raw data (frame 1)
subplot(3,4,1);
imagesc(imageStack(:,:,1));
axis image; colormap(gca, 'gray'); colorbar;
title('Raw Data (Frame 1)');
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Panel 2: Raw data (frame 100)
subplot(3,4,2);
imagesc(imageStack(:,:,100));
axis image; colormap(gca, 'gray'); colorbar;
title('Raw Data (Frame 100)');
xlabel('X (pixels)'); ylabel('Y (pixels)');

% Panel 3: All trajectories
subplot(3,4,3);
hold on;
for i = 1:min(length(TR), 100)  % Plot first 100
    plot(TR(i).X, TR(i).Y, '-', 'LineWidth', 1);
end
axis equal;
xlim([0 SimParams.FrameSize(2)]);
ylim([0 SimParams.FrameSize(1)]);
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('%d Trajectories', min(length(TR), 100)));
grid on;

% Panel 4: Trajectories colored by time
subplot(3,4,4);
hold on;
for i = 1:min(length(TR), 50)  % Plot first 50
    frames = TR(i).FrameNum;
    scatter(TR(i).X, TR(i).Y, 20, frames, 'filled');
end
colorbar;
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Trajectories (colored by time)');
axis equal;
xlim([0 SimParams.FrameSize(2)]);
ylim([0 SimParams.FrameSize(1)]);

% Panel 5: Example trajectory detail
subplot(3,4,5);
if length(TR) >= 10
    idx = 10;
    plot(TR(idx).X, TR(idx).Y, 'b.-', 'MarkerSize', 10, 'LineWidth', 1.5);
    xlabel('X (pixels)'); ylabel('Y (pixels)');
    title(sprintf('Trajectory %d (%d points)', idx, length(TR(idx).X)));
    axis equal; grid on;
end

% Panel 6: Trajectory length distribution
subplot(3,4,6);
histogram(traj_lengths, 30);
xlabel('Trajectory Length (frames)');
ylabel('Count');
title(sprintf('Lengths (median: %d frames)', median(traj_lengths)));
grid on;

% Panel 7: Ensemble MSD curve
subplot(3,4,7);
MSD_ensemble = DE.MSDEnsemble;
tau = MSD_ensemble.FrameLags / FrameRate;  % Convert to seconds
msd_um2 = MSD_ensemble.MSD * (PixelSize^2);  % Convert to μm^2

plot(tau, msd_um2, 'ko-', 'MarkerSize', 8, 'LineWidth', 2);
hold on;

% Fit line: MSD = 4*D*t for 2D diffusion
fit_points = min(DE.NFitPoints, length(tau));
p = polyfit(tau(1:fit_points), msd_um2(1:fit_points), 1);
plot(tau, polyval(p, tau), 'r--', 'LineWidth', 2);

xlabel('Time Lag (s)');
ylabel('MSD (μm^2)');
title(sprintf('Ensemble MSD (D=%.3f μm^2/s)', D_ensemble));
legend('Data', 'Fit', 'Location', 'northwest');
grid on;

% Panel 8: Example single trajectory MSD
subplot(3,4,8);
if length(DE.MSDSingleTraj) >= 1
    msd_single = DE.MSDSingleTraj(1);
    tau_single = msd_single.FrameLags / FrameRate;
    msd_single_um2 = msd_single.MSD * (PixelSize^2);

    plot(tau_single, msd_single_um2, 'bo-', 'MarkerSize', 6, 'LineWidth', 1.5);
    xlabel('Time Lag (s)');
    ylabel('MSD (μm^2)');
    title(sprintf('Single Trajectory MSD (D=%.3f μm^2/s)', D_trajectory(1)));
    grid on;
end

% Panel 9: Diffusion coefficient distribution
subplot(3,4,9);
D_valid = D_trajectory(~isnan(D_trajectory));
histogram(D_valid, 30);
hold on;
xline(D_physical, 'g--', 'LineWidth', 2, 'Label', 'True');
xline(D_mean, 'r--', 'LineWidth', 2, 'Label', 'Mean');
xlabel('Diffusion Coefficient (μm^2/s)');
ylabel('Count');
title(sprintf('D Distribution (mean=%.3f)', D_mean));
legend('Location', 'best');
grid on;

% Panel 10: Diffusion error vs trajectory length
subplot(3,4,10);
D_error = abs(D_trajectory - D_physical) ./ D_physical * 100;  % Percent error
D_valid_idx = ~isnan(D_error);
scatter(traj_lengths(D_valid_idx), D_error(D_valid_idx), 30, 'filled', ...
    'MarkerFaceAlpha', 0.5);
xlabel('Trajectory Length (frames)');
ylabel('Error (%)');
title('Estimation Error vs Length');
grid on;

% Panel 11: Step size distribution
subplot(3,4,11);
all_steps = [];
for i = 1:length(TR)
    dx = diff(TR(i).X);
    dy = diff(TR(i).Y);
    steps = sqrt(dx.^2 + dy.^2);
    all_steps = [all_steps; steps];
end
steps_nm = all_steps * PixelSize * 1000;  % Convert to nm
histogram(steps_nm, 50);
xlabel('Step Size (nm)');
ylabel('Count');
title(sprintf('Steps (mean=%.1f nm)', mean(steps_nm)));
grid on;

% Panel 12: Localization precision
subplot(3,4,12);
precision_nm = SMD.X_SE * PixelSize * 1000;  % Convert to nm
histogram(precision_nm, 40);
xlabel('Localization Precision (nm)');
ylabel('Count');
title(sprintf('Precision (median=%.1f nm)', median(precision_nm)));
grid on;

%% Step 9: Summary Statistics
fprintf('\n=== Summary Statistics ===\n');
fprintf('========================================\n');
fprintf('DATA:\n');
fprintf('  True trajectories: %d\n', length(TR_true));
fprintf('  Frames: %d\n', SimParams.NFrames);
fprintf('  Localizations found: %d\n', length(SMD.X));
fprintf('  Tracked trajectories: %d\n', length(TR));
fprintf('  Mean trajectory length: %.1f frames\n', mean(traj_lengths));
fprintf('  Median trajectory length: %d frames\n', median(traj_lengths));

fprintf('\nLOCALIZATION QUALITY:\n');
fprintf('  Median precision: %.1f nm\n', median(precision_nm));
fprintf('  Mean photons: %.0f\n', mean(SMD.Photons));
fprintf('  Mean background: %.1f photons/pixel\n', mean(SMD.Bg));

fprintf('\nDIFFUSION RESULTS:\n');
fprintf('  True D: %.3f μm^2/s\n', D_physical);
fprintf('  Estimated D (ensemble): %.3f ± %.3f μm^2/s\n', ...
    D_ensemble, D_ensemble_SE);
fprintf('  Estimated D (mean): %.3f ± %.3f μm^2/s\n', D_mean, D_std);
fprintf('  Ensemble error: %.1f%%\n', ensemble_error);
fprintf('  Mean trajectory error: %.1f%%\n', mean_error);

fprintf('\nPERFORMANCE:\n');
fprintf('  Localization time: %.2f s\n', localization_time);
fprintf('  Tracking time: %.2f s\n', tracking_time);
fprintf('  Diffusion estimation time: %.2f s\n', diffusion_time);
fprintf('  Total time: %.2f s\n', ...
    localization_time + tracking_time + diffusion_time);
fprintf('========================================\n');

fprintf('\nExample complete!\n');
```

## What This Example Demonstrates

### Simulation
- Generates realistic diffusing particles with known diffusion coefficient
- Includes blinking behavior (on/off switching)
- Creates photon-realistic images with PSF and noise
- Provides ground truth for validation

### Localization
- Detects particles in each frame independently
- Fits PSF to determine precise positions
- Filters by quality metrics
- Same pipeline as SMLM but with different parameters

### Tracking
- **Frame-to-frame linking**: Connects nearby particles in consecutive frames
- **Gap closing**: Bridges gaps from blinking or missed detections
- **Trajectory filtering**: Removes short, unreliable trajectories
- Uses cost matrix approach based on diffusion model

### MSD Analysis
- Computes trajectory-wise MSD curves
- Computes ensemble MSD from all trajectories
- Fits linear model to extract diffusion coefficient
- Provides uncertainty estimates

### Validation
- Compares estimated vs true diffusion coefficient
- Shows accuracy improves with longer trajectories
- Demonstrates ensemble vs single-trajectory estimates

## Expected Results

When you run this example, you should see:

**Trajectories**: 40-60 trajectories with 10+ frames (depends on blinking and tracking success)

**Trajectory lengths**: 10-100 frames, median ~20 frames

**Diffusion accuracy**: Ensemble estimate within 5-15% of true value

**Localization precision**: ~10-20 nm (depends on photons and background)

**Processing speed**: ~5-10 seconds total (depends on GPU)

## Understanding the Results

### Why Isn't the Estimate Exact?

Several factors cause deviation from the true diffusion coefficient:

1. **Finite trajectory lengths**: Short trajectories have noisy MSD estimates
2. **Localization uncertainty**: Position errors propagate to MSD
3. **Motion blur**: Particle moves during exposure
4. **Model assumptions**: Assumes pure Brownian motion
5. **Sampling**: Limited number of trajectories

Typical accuracy: 5-20% for well-tracked particles with 20+ frame trajectories.

### Ensemble vs Trajectory-Wise Estimates

**Ensemble**:
- Averages MSD across all trajectories
- More stable, lower variance
- Best for homogeneous populations

**Trajectory-wise**:
- Independent estimate per trajectory
- Reveals heterogeneity
- Useful for detecting subpopulations

### MSD Curve Interpretation

For 2D Brownian diffusion:
```
MSD(τ) = 4Dτ + offset
```

Where:
- D = diffusion coefficient
- τ = time lag
- offset = localization uncertainty (should be small)

The slope of MSD vs time gives 4D, so D = slope/4.

## Modifications to Try

### 1. Vary Diffusion Coefficient

```matlab
% Faster diffusion
SimParams.D = 0.3;  % pixels^2/frame
SMF.Tracking.D = 0.3;  % Update tracking parameter!
SMF.Tracking.MaxDistFF = 8;  % Increase search radius

% Slower diffusion
SimParams.D = 0.05;  % pixels^2/frame
SMF.Tracking.D = 0.05;
SMF.Tracking.MaxDistFF = 3;  % Decrease search radius
```

### 2. Change Blinking Behavior

```matlab
% Less blinking (longer continuous tracks)
SimParams.KOffToOn = 0.95;
SimParams.KOnToOff = 0.02;

% More blinking (shorter segments, more gaps)
SimParams.KOffToOn = 0.7;
SimParams.KOnToOff = 0.2;
SMF.Tracking.MaxFrameGap = 10;  % Allow longer gaps
```

### 3. Vary Signal Quality

```matlab
% Brighter particles (better precision)
SimParams.Intensity = 2000;

% Dimmer particles (worse precision)
SimParams.Intensity = 500;
SMF.BoxFinding.MinPhotons = 150;  % Lower threshold

% Higher background (worse SNR)
SimParams.Bg = 30;
```

### 4. Change Particle Density

```matlab
% Higher density (more crowding, harder tracking)
SimParams.ParticleDensity = 0.005;  % ~160 particles
SMF.Tracking.MaxDistFF = 4;  % Be more conservative

% Lower density (easier tracking)
SimParams.ParticleDensity = 0.001;  % ~30 particles
```

### 5. More Frames

```matlab
% Longer acquisition (more trajectories, better statistics)
SimParams.NFrames = 500;
```

### 6. Different MSD Fit Targets

```matlab
% Use CDF of jumps (better for heterogeneous populations)
DE.FitTarget = 'CDFOfJumps';
DE.NComponents = 1;  % Can increase for multiple populations

% Use maximum likelihood (most accurate)
DE.FitTarget = 'LikelihoodOfJumps';
DE.NComponents = 1;
```

### 7. Multiple Diffusive States

```matlab
% Simulate two-state diffusion
SimParams.D = 0.1;  % Base diffusion
% After creating simulation, modify some trajectories:
for i = 1:length(TR_true)/2
    % Make half the trajectories diffuse faster
    TR_true(i).X = TR_true(i).X * 1.5;
    TR_true(i).Y = TR_true(i).Y * 1.5;
end

% Estimate with two components
DE.FitTarget = 'CDFOfJumps';
DE.NComponents = 2;
```

## Troubleshooting

### Issue: Low number of trajectories

**Diagnosis**:
```matlab
fprintf('True trajectories: %d\n', length(TR_true));
fprintf('Found trajectories: %d\n', length(TR));
```

**Solutions**:
- Increase `MaxFrameGap` (allow longer gaps)
- Increase `MaxDistGC` (allow larger gap jumps)
- Lower `MinTrackLength` (accept shorter tracks)
- Check localization quality (need detections in most frames)

### Issue: Poor diffusion estimate accuracy

**Diagnosis**:
```matlab
% Plot MSD curves
figure;
plot(MSD_ensemble.FrameLags, MSD_ensemble.MSD, 'o-');
xlabel('Frame Lag'); ylabel('MSD (pixels^2)');
% Should be linear for first few points
```

**Solutions**:
- Increase trajectory lengths (more frames or lower `MinTrackLength`)
- Check tracking isn't linking wrong particles (plot trajectories)
- Verify `FrameLagRange` and `NFitPoints` are appropriate
- Ensure tracking parameters match true diffusion
- Use `FitTarget = 'LikelihoodOfJumps'` for better accuracy

### Issue: Trajectories broken inappropriately

**Diagnosis**:
```matlab
% Check for multiple short trajectories in same region
starts = arrayfun(@(x) [x.X(1), x.Y(1)], TR, 'UniformOutput', false);
starts = vertcat(starts{:});
figure; scatter(starts(:,1), starts(:,2), 'filled');
% Clustered starts suggest broken tracks
```

**Solutions**:
- Increase `MaxFrameGap` (bridge longer gaps)
- Increase `MaxDistGC` (allow larger gap jumps)
- Check blinking parameters `K_on`, `K_off`
- Verify localization works in all frames

### Issue: Wrong particles linked

**Diagnosis**:
```matlab
% Plot overlapping trajectories
figure; hold on;
for i = 1:20
    plot(TR(i).X, TR(i).Y, '.-', 'LineWidth', 1.5);
end
% Crossing trajectories indicate wrong linking
```

**Solutions**:
- Decrease `MaxDistFF` (be more conservative)
- Decrease `MaxDistGC` (reduce aggressive gap closing)
- Reduce particle density (less crowding)
- Improve localization precision (more photons)

## Advanced Topics

### Heterogeneous Diffusion

Detect multiple diffusive populations:

```matlab
DE.FitTarget = 'CDFOfJumps';
DE.NComponents = 2;  % Two diffusive states
DiffusionStruct = DE.estimateDiffusionConstant();

% Extract results
D1 = DiffusionStruct(2).DiffusionConstant(1);
D2 = DiffusionStruct(2).DiffusionConstant(2);
ratio1 = DiffusionStruct(2).PopulationRatios(1);
ratio2 = DiffusionStruct(2).PopulationRatios(2);

fprintf('State 1: D=%.3f μm^2/s (%.1f%%)\n', D1, ratio1*100);
fprintf('State 2: D=%.3f μm^2/s (%.1f%%)\n', D2, ratio2*100);
```

### Hidden Markov Model Analysis

Detect state switches within trajectories:

```matlab
HMM = smi_stat.HMM();
HMM.NStates = 2;
HMM.PixelSize = PixelSize;
HMM.FrameRate = FrameRate;

% Analyze trajectory
[states, D_states] = HMM.analyzeTrajectory(TR(10), SMF);

% Plot trajectory with states
figure;
subplot(2,1,1);
plot(TR(10).X, TR(10).Y, '.-');
title('Trajectory Path');

subplot(2,1,2);
plot(states, 'o-');
ylabel('State'); xlabel('Time Point');
title(sprintf('States (D1=%.3f, D2=%.3f μm^2/s)', D_states(1), D_states(2)));
```

## Complete Analysis Pipeline

For a publication-quality analysis:

```matlab
% 1. Generate/load data
% ... (from example above)

% 2. Track particles
SPT = smi.SPT(SMF, false);
SPT.performFullAnalysis();

% 3. Load results
load(fullfile(SMF.Data.FileDir, 'Results', 'Results.mat'), 'TR', 'SMD', 'SMF');

% 4. Estimate diffusion
DE = smi_stat.DiffusionEstimator(TR, SMF, 1, true);  % AutoRun=true

% 5. Visualize ensemble MSD
figure;
smi_stat.DiffusionEstimator.plotEnsembleMSD(gca, ...
    DE.MSDEnsemble, DE.DiffusionStruct, 'Brownian', true);

% 6. Save results
DE.SaveDir = fullfile(SMF.Data.FileDir, 'Results');
DE.saveResults();

% 7. Generate trajectory movies
MovieMaker = smi_vis.GenerateMovies();
MovieMaker.TR = TR;
MovieMaker.SMD = SMD;
MovieMaker.SMF = SMF;
MovieMaker.RawData = imageStack;
MovieMaker.gui();  % Interactive movie generation
```

## See Also

- [SPT Tracking Workflow](../workflows/spt-tracking.md) - Complete tracking pipeline
- [TR Structure](../core-concepts/tr-structure.md) - Trajectory data format
- [Basic Localization Example](basic-localization.md) - Localization fundamentals
- MATLAB/examples/Example_SPT.m - Additional SPT examples
- MATLAB/+smi_stat/@DiffusionEstimator/ - Diffusion analysis methods
- MATLAB/+smi/@SPT/ - SPT class documentation
