---
title: "Single Particle Tracking Workflow"
category: "workflows"
level: "intermediate"
tags: ["spt", "tracking", "trajectories", "diffusion"]
prerequisites: ["../core-concepts/smf-structure.md", "../core-concepts/smd-structure.md"]
related: ["smlm-analysis.md", "../how-to/localize-molecules.md"]
summary: "Complete workflow for tracking single particles through time and analyzing their motion"
estimated_time: "20 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Single Particle Tracking Workflow

## Purpose

Single Particle Tracking (SPT) analyzes the motion of individual molecules or particles over time. This workflow guides you through localizing particles in each frame, linking localizations into trajectories, and analyzing the resulting motion. You'll learn how to configure tracking parameters, understand the TR (Tracking Results) structure, and extract meaningful dynamics information.

## Prerequisites

- Understanding of [SMF structure](../core-concepts/smf-structure.md)
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Basic tracking concepts (diffusion, mean square displacement)
- Completion of [SMLM workflow](smlm-analysis.md) helpful but not required

## Overview

SPT differs from SMLM in key ways:
- **SMLM**: Emitters blink, appear in many frames, are stationary
- **SPT**: Particles move continuously, appear in consecutive frames

The SPT workflow:
1. **Localize** particles in each frame (like SMLM)
2. **Track** particles frame-to-frame and through gaps (unique to SPT)
3. **Organize** results as trajectories (TR structure)
4. **Analyze** motion (diffusion, state changes, etc.)

The output is **TR** (Tracking Results): an array where each element is an SMD structure for one trajectory.

## SPT Pipeline Architecture

### Data Flow

```
┌─────────────────┐
│ Raw Data (.mat) │ Movie of moving particles
└────────┬────────┘
         │ LoadData + LocalizeData
         ▼
┌─────────────────┐
│  SMD (all locs) │ Localizations per frame
└────────┬────────┘
         │ Frame-to-Frame Tracking
         ├─► Build cost matrix (distances)
         ├─► Solve linear assignment problem
         └─► Link nearby localizations
         ▼
┌─────────────────┐
│  Initial Tracks │ Short trajectory segments
└────────┬────────┘
         │ Gap Closing
         ├─► Build cost matrix (allow gaps)
         ├─► Solve assignment problem
         └─► Merge trajectory segments
         ▼
┌─────────────────┐
│  TR (raw)       │ Array of trajectories
└────────┬────────┘
         │ Filter by length, quality
         ▼
┌─────────────────┐
│  TR (final)     │ Filtered trajectories
└────────┬────────┘
         │ Analysis + Visualization
         ▼
┌─────────────────┐
│ MSD, Movies     │ Motion analysis results
└─────────────────┘
```

### Control via SMF.Tracking

All tracking behavior is controlled by `SMF.Tracking`:

```matlab
SMF.Tracking.Method = 'CostMatrix';  % Tracking algorithm
SMF.Tracking.D = 0.1;                % Expected diffusion (pixels²/frame)
SMF.Tracking.MaxDistFF = 5;          % Max frame-to-frame distance
SMF.Tracking.MaxDistGC = 10;         % Max gap closing distance
SMF.Tracking.MaxFrameGap = 5;        % Max frames to bridge
SMF.Tracking.MinTrackLength = 10;    % Min trajectory length to keep
```

## Step-by-Step SPT Workflow

### Step 1: Generate or Load Particle Data

**For simulated data:**

```matlab
% Create simulation
SPTSim = smi_sim.SimSPT;
SPTSim.SimParams.FrameSize = [128, 128];
SPTSim.SimParams.ParticleDensity = 0.002;  % particles/pixel²
SPTSim.SimParams.D = 0.1;                  % pixels²/s
SPTSim.SimParams.KOnToOff = 0.05;          % Blinking off rate
SPTSim.SimParams.KOffToOn = 0.9;           % Blinking on rate
SPTSim.SimParams.Intensity = 1000;         % Photons per particle

% Generate
SPTSim.createSimulation();

% Convert to image stack
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.DataROI = [1, 1, 128, 128];  % Required for gaussBlobImage
SMF.Fitting.PSFSigma = 1.3;
[~, sequence] = smi_sim.GaussBlobs.gaussBlobImage(SPTSim.SMD, SMF);

% Save
save('tracking_data.mat', 'sequence');
```

**For real data:**
Prepare .mat or .h5 file with image stack.

### Step 2: Configure SMF for Tracking

```matlab
% Create SMF
SMF = smi_core.SingleMoleculeFitting();

% Data location
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'tracking_data.mat'};
SMF.Data.PixelSize = 0.1;    % micrometers
SMF.Data.FrameRate = 100;    % Hz

% Localization (similar to SMLM)
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 100;  % Lower for tracking (particles can be dim)

% Thresholding (more permissive for tracking)
SMF.Thresholding.AutoThreshLogL = true;
SMF.Thresholding.MaxXY_SE = 0.3;  % Looser than SMLM

% Disable frame connection (SPT does its own linking)
SMF.FrameConnection.On = false;

% Tracking parameters
SMF.Tracking.Method = 'CostMatrix';
SMF.Tracking.D = 0.1;              % Expected diffusion (pixels²/frame)
SMF.Tracking.TrajwiseD = true;     % Estimate D per trajectory
SMF.Tracking.MaxDistFF = 5;        % Max distance frame-to-frame
SMF.Tracking.MaxDistGC = 10;       % Max distance for gap closing
SMF.Tracking.MaxFrameGap = 5;      % Max gap to bridge
SMF.Tracking.MinTrackLength = 10;  % Keep trajectories ≥10 frames

% Blinking parameters
SMF.Tracking.K_on = 0.9;   % Prob of appearing when off
SMF.Tracking.K_off = 0.1;  % Prob of disappearing when on
```

**Choosing diffusion constant D:**

Convert from physical units:
```
D_pixels = D_μm²/s / (pixel_size_μm² × frame_rate_Hz)
```

Example: D = 0.1 μm²/s, pixel = 0.1 μm, rate = 100 Hz
```
D_pixels = 0.1 / (0.01 × 100) = 0.1 pixels²/frame
```

### Step 3: Run Tracking Analysis

**Using smi.SPT class:**

```matlab
% Create SPT object
SPT = smi.SPT(SMF);

% Optionally disable movies/plots for faster processing
SPT.GenerateMovies = false;
SPT.GeneratePlots = true;

% Run complete analysis
SPT.performFullAnalysis();
```

This performs:
1. Localization (via `LocalizeData`)
2. Frame-to-frame tracking
3. Gap closing
4. Filtering by trajectory length
5. Saving results

**Results saved to:**
```
FileDir/Results/
├── [FileName]_Results.mat   % Contains TR, SMD, SMF
└── [FileName]_Plots/         % Diagnostic plots
```

### Step 4: Understanding TR (Tracking Results)

Load results:

```matlab
load('FileDir/Results/tracking_data_Results.mat', 'TR', 'SMD', 'SMF');
```

**TR structure:**
- TR is an **array** where each element is one trajectory
- Each TR(i) is an SMD structure with fields: X, Y, Photons, FrameNum, etc.
- Length of TR = number of trajectories
- Length of TR(i).X = number of localizations in trajectory i

```matlab
% Number of trajectories
N_traj = length(TR);
fprintf('Found %d trajectories\n', N_traj);

% First trajectory
fprintf('Trajectory 1: %d localizations\n', length(TR(1).X));
fprintf('  Frames: %d to %d\n', min(TR(1).FrameNum), max(TR(1).FrameNum));

% Access positions
X_traj1 = TR(1).X;  % X positions (pixels)
Y_traj1 = TR(1).Y;  % Y positions (pixels)
t_traj1 = (TR(1).FrameNum - 1) / SMF.Data.FrameRate;  % Time (seconds)
```

### Step 5: Visualize Trajectories

**Plot all trajectories:**

```matlab
figure; hold on;
for i = 1:length(TR)
    plot(TR(i).X, TR(i).Y, '-', 'LineWidth', 1);
end
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('%d Trajectories', length(TR)));
axis equal;
```

**Plot with time color-coding:**

```matlab
figure; hold on;
for i = 1:min(100, length(TR))  % Plot first 100
    frames = TR(i).FrameNum;
    scatter(TR(i).X, TR(i).Y, 20, frames, 'filled');
end
colorbar; xlabel('X (pixels)'); ylabel('Y (pixels)');
title('Trajectories (colored by time)');
```

**Individual trajectory analysis:**

```matlab
% Pick one trajectory
idx = 10;  % 10th trajectory

% Position vs time
t = (TR(idx).FrameNum - 1) / SMF.Data.FrameRate;
figure;
subplot(2,1,1);
plot(t, TR(idx).X, '.-'); ylabel('X (pixels)'); xlabel('Time (s)');
subplot(2,1,2);
plot(t, TR(idx).Y, '.-'); ylabel('Y (pixels)'); xlabel('Time (s)');
sgtitle(sprintf('Trajectory %d', idx));
```

### Step 6: Analyze Motion

**Mean Square Displacement (MSD):**

```matlab
% Compute MSD for one trajectory
function msd = computeMSD(X, Y)
    N = length(X);
    max_tau = floor(N/4);  % Compute up to 1/4 of trajectory length
    msd = zeros(max_tau, 1);

    for tau = 1:max_tau
        displacements = (X(1+tau:end) - X(1:end-tau)).^2 + ...
                       (Y(1+tau:end) - Y(1:end-tau)).^2;
        msd(tau) = mean(displacements);
    end
end

% Example for trajectory 1
msd = computeMSD(TR(1).X, TR(1).Y);
tau = (1:length(msd)) / SMF.Data.FrameRate;  % Time lags (seconds)

figure;
plot(tau, msd * SMF.Data.PixelSize^2);  % Convert to μm²
xlabel('Time lag (s)'); ylabel('MSD (μm²)');
title('Mean Square Displacement');

% Fit to extract diffusion coefficient
% For 2D diffusion: MSD = 4*D*tau
p = polyfit(tau, msd * SMF.Data.PixelSize^2, 1);
D_fit = p(1) / 4;
fprintf('Estimated D = %.3f μm²/s\n', D_fit);
```

**Using DiffusionEstimator:**

```matlab
DE = smi_stat.DiffusionEstimator();
DE.FrameRate = SMF.Data.FrameRate;
DE.PixelSize = SMF.Data.PixelSize;

% Estimate D for one trajectory
TR_single = TR(1);
[D_est, ~] = DE.estimateDiffusion(TR_single);
fprintf('Diffusion: %.3f μm²/s\n', D_est);

% Batch estimate for all trajectories
D_all = zeros(length(TR), 1);
for i = 1:length(TR)
    [D_all(i), ~] = DE.estimateDiffusion(TR(i));
end

% Plot distribution
histogram(D_all, 30);
xlabel('Diffusion Coefficient (μm²/s)');
ylabel('Count');
title('Diffusion Distribution');
```

**Step size distribution:**

```matlab
% Compute step sizes for all trajectories
all_steps = [];
for i = 1:length(TR)
    dx = diff(TR(i).X);
    dy = diff(TR(i).Y);
    steps = sqrt(dx.^2 + dy.^2);
    all_steps = [all_steps; steps];
end

% Convert to physical units
steps_nm = all_steps * SMF.Data.PixelSize * 1000;

% Plot
histogram(steps_nm, 50);
xlabel('Step Size (nm)'); ylabel('Count');
title('Step Size Distribution');
```

### Step 7: Generate Movies (Optional)

```matlab
% Load raw data
LD = smi_core.LoadData();
[~, sequence, ~] = LD.loadRawData(SMF, 1);

% Create movie maker
MovieMaker = smi_vis.GenerateMovies();
MovieMaker.TR = TR;
MovieMaker.SMD = SMD;
MovieMaker.SMF = SMF;
MovieMaker.RawData = sequence;

% Open interactive GUI
MovieMaker.gui();

% Or generate programmatically
MovieMaker.StartFrame = 1;
MovieMaker.EndFrame = 100;
MovieMaker.ShowTrails = true;
MovieMaker.TrailLength = 10;
movie_frames = MovieMaker.generateMovie();
```

## Tracking Algorithm Details

### Frame-to-Frame Tracking

**Cost matrix approach:**

1. For frames t and t+1, compute distance matrix:
   ```
   D(i,j) = distance from localization i in frame t
                      to localization j in frame t+1
   ```

2. Build cost matrix:
   ```
   C(i,j) = -log(P(link i→j | D, parameters))
   ```
   Based on diffusion model and blinking probabilities.

3. Solve linear assignment problem (LAP):
   Find assignments minimizing total cost subject to:
   - Each localization in frame t assigned to ≤1 in frame t+1
   - Each localization in frame t+1 assigned from ≤1 in frame t

4. Accept links with cost below threshold.

**Parameters affecting this:**
- `MaxDistFF`: Maximum distance to consider
- `D`: Expected diffusion (affects probability model)
- `K_on`, `K_off`: Blinking probabilities

### Gap Closing

After frame-to-frame tracking, handle gaps from blinking or missed detections:

1. Identify trajectory segments (between gaps)

2. Build cost matrix for linking segments:
   ```
   C(i,j) = cost to link end of segment i to start of segment j
   ```
   Considers:
   - Distance between endpoints
   - Time gap
   - Expected motion during gap

3. Solve LAP to merge segments

4. Accept merges with:
   - Gap ≤ `MaxFrameGap` frames
   - Distance ≤ `MaxDistGC` pixels

### Iterative Refinement

If `TrajwiseD = true`:
1. Perform initial tracking with global D
2. Estimate D for each trajectory from MSD
3. Re-track using trajectory-specific D values
4. Iterate until convergence

This improves tracking when particles have heterogeneous diffusion.

## Troubleshooting

### Issue: Too many short trajectories

**Diagnose:**
```matlab
traj_lengths = arrayfun(@(x) length(x.X), TR);
histogram(traj_lengths);
xlabel('Trajectory Length (frames)');
```

**Solutions:**
- Increase `MaxDistFF` (particles moving faster than expected)
- Increase `MaxFrameGap` (more blinking than expected)
- Increase `MaxDistGC` (allow larger gap closing jumps)
- Lower `MinTrackLength` if short trajectories are acceptable
- Check localization quality (poor locs → missed detections)

### Issue: Incorrect trajectory linking

**Diagnose:**
```matlab
% Plot overlapping trajectories (indicates wrong linking)
figure; hold on;
for i = 1:min(50, length(TR))
    plot(TR(i).X, TR(i).Y, '.-', 'LineWidth', 1.5);
end
% Look for trajectories that cross or swap
```

**Solutions:**
- Decrease `MaxDistFF` (too permissive, linking wrong particles)
- Decrease `MaxDistGC` (gap closing too aggressive)
- Adjust `D` to better match true diffusion
- Reduce particle density (too crowded)
- Improve localization precision

### Issue: Trajectories broken inappropriately

**Diagnose:**
```matlab
% Check for many similar short trajectories in same region
% Compute trajectory start/end positions
starts = arrayfun(@(x) [x.X(1), x.Y(1)], TR, 'UniformOutput', false);
starts = vertcat(starts{:});
scatter(starts(:,1), starts(:,2), 'filled');
title('Trajectory Start Positions');
% Clustered starts suggest broken trajectories
```

**Solutions:**
- Increase `MaxFrameGap` (allow longer gaps)
- Increase `MaxDistGC` (allow larger jumps during gaps)
- Check blinking parameters `K_on`, `K_off`
- Ensure localization works on all frames

## Advanced Analysis

### State Detection with HMM

Detect diffusive state changes:

```matlab
HMM = smi_stat.HMM();
HMM.NStates = 2;  % Two diffusive states

% Analyze one trajectory
TR_single = TR(10);
[states, D_states] = HMM.analyzeTrajectory(TR_single, SMF);

% Plot
t = (TR_single.FrameNum - 1) / SMF.Data.FrameRate;
figure;
subplot(3,1,1);
plot(t, TR_single.X, '.-'); ylabel('X');
subplot(3,1,2);
plot(t, TR_single.Y, '.-'); ylabel('Y');
subplot(3,1,3);
plot(t, states, '.-'); ylabel('State');
xlabel('Time (s)');
```

### Change Point Detection

Find abrupt motion changes:

```matlab
CD = smi_stat.ChangeDetection();
[change_points, D_segments] = CD.detectChanges(TR(10));

% Plot with change points marked
figure;
plot(TR(10).X, TR(10).Y, '.-');
hold on;
for cp = change_points'
    plot(TR(10).X(cp), TR(10).Y(cp), 'ro', 'MarkerSize', 10);
end
xlabel('X'); ylabel('Y'); title('Change Points');
```

## Complete Example

```matlab
% ========== Setup ==========
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/tracking';
SMF.Data.FileName = {'particle_movie.mat'};
SMF.Data.PixelSize = 0.1;
SMF.Data.FrameRate = 100;

SMF.Fitting.PSFSigma = 1.3;
SMF.BoxFinding.MinPhotons = 100;
SMF.Thresholding.AutoThreshLogL = true;

SMF.Tracking.D = 0.1;
SMF.Tracking.MaxDistFF = 5;
SMF.Tracking.MaxDistGC = 10;
SMF.Tracking.MaxFrameGap = 5;
SMF.Tracking.MinTrackLength = 10;

% ========== Run Tracking ==========
SPT = smi.SPT(SMF);
SPT.GenerateMovies = false;
SPT.performFullAnalysis();

% ========== Load Results ==========
results_file = fullfile(SMF.Data.FileDir, 'Results', ...
    'particle_movie_Results.mat');
load(results_file, 'TR', 'SMD', 'SMF');

% ========== Analyze ==========
fprintf('=== Tracking Summary ===\n');
fprintf('Total trajectories: %d\n', length(TR));
traj_lengths = arrayfun(@(x) length(x.X), TR);
fprintf('Mean trajectory length: %.1f frames\n', mean(traj_lengths));

% Estimate diffusion
DE = smi_stat.DiffusionEstimator();
DE.FrameRate = SMF.Data.FrameRate;
DE.PixelSize = SMF.Data.PixelSize;
D_all = zeros(length(TR), 1);
for i = 1:length(TR)
    [D_all(i), ~] = DE.estimateDiffusion(TR(i));
end
fprintf('Mean diffusion: %.3f μm²/s\n', mean(D_all));

% ========== Visualize ==========
figure;
subplot(2,2,1);
hold on;
for i = 1:min(100, length(TR))
    plot(TR(i).X, TR(i).Y, '-');
end
xlabel('X'); ylabel('Y'); title('Trajectories');

subplot(2,2,2);
histogram(traj_lengths, 30);
xlabel('Length (frames)'); ylabel('Count');
title('Trajectory Lengths');

subplot(2,2,3);
histogram(D_all, 30);
xlabel('D (μm²/s)'); ylabel('Count');
title('Diffusion Distribution');

subplot(2,2,4);
msd_example = computeMSD(TR(1).X, TR(1).Y);
tau = (1:length(msd_example)) / SMF.Data.FrameRate;
plot(tau, msd_example * SMF.Data.PixelSize^2);
xlabel('Time lag (s)'); ylabel('MSD (μm²)');
title('Example MSD');
```

## See Also

- [SMLM Workflow](smlm-analysis.md) - Localization pipeline
- [SMF Structure](../core-concepts/smf-structure.md) - Tracking parameters
- [SMD Structure](../core-concepts/smd-structure.md) - Localization format
- doc/DataStructures/TR.md - TR structure reference
- MATLAB/+smi/@SPT/README.md - SPT class documentation
- MATLAB/examples/Example_SPT.m - Complete example script
