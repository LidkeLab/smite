---
title: "TR Structure: Tracking Results"
category: "core-concepts"
level: "intermediate"
tags: ["tr", "tracking", "spt", "data-structures", "trajectories"]
prerequisites: ["smd-structure.md", "architecture.md"]
related: ["../workflows/spt-tracking.md", "smd-structure.md"]
summary: "Complete reference to the Tracking Results (TR) structure for single particle tracking data"
estimated_time: "10 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# TR Structure: Tracking Results

## Purpose

The TR (Tracking Results) structure is smite's format for organizing single particle tracking (SPT) data by trajectory. While SMD stores all localizations as a flat list, TR organizes them into separate trajectories, making it much easier to analyze individual particle paths, compute trajectory statistics, and visualize tracks.

## Prerequisites

- Understanding of [SMD structure](smd-structure.md)
- Familiarity with [smite architecture](architecture.md)
- Basic SPT concepts (trajectories, frame connection)

## Overview

TR is fundamentally **an array of SMD structures**, where each array element represents one trajectory:

```matlab
% After tracking analysis
load('Results.mat', 'TR');

% TR is an array - each element is one trajectory
N_trajectories = length(TR);

% Access first trajectory
trajectory_1 = TR(1);

% trajectory_1 is an SMD structure with all localizations for that track
positions_1 = [TR(1).X, TR(1).Y];
frames_1 = TR(1).FrameNum;
```

**Key concept:** TR transforms SMD from "all localizations" to "localizations grouped by trajectory."

### SMD vs TR

| Aspect | SMD | TR |
|--------|-----|-----|
| **Organization** | Flat list of all localizations | Array of trajectories |
| **Indexing** | `SMD.X(i)` = i-th localization | `TR(n).X(i)` = i-th point in n-th trajectory |
| **Best for** | SMLM, visualization, all-data operations | SPT analysis, trajectory statistics |
| **Size** | Single structure | Array of structures |
| **ConnectID** | Links localizations to tracks | All points have same ConnectID |

## Creating and Accessing TR

### From SPT Analysis

TR is created automatically by tracking workflows:

```matlab
% Option 1: Full SPT workflow
SPT = smi.SPT(SMF);
SPT.performFullAnalysis();
load(SPT.SMF.Data.ResultsDir, 'TR', 'SMD');

% Option 2: Frame connection only
SMD = ...; % From localization
FC = smi_core.FrameConnection();
FC.SMF = SMF;
SMD_combined = FC.runFrameConnection(SMD);
TR = smi_core.TrackingResults.convertSMDToTR(SMD_combined);
```

### Loading from Files

```matlab
% Load tracking results
load('FileDir/Results/Results_SPT.mat', 'TR', 'SMD', 'SMF');

% Now TR is available
fprintf('Found %d trajectories\n', length(TR));
```

### Creating Empty TR

```matlab
% Create empty TR structure
TR = smi_core.TrackingResults.createTR();

% TR is just an SMD structure
% Can be used as template for array elements
```

### Converting Between SMD and TR

```matlab
% SMD with ConnectID → TR
SMD = ...; % Must have ConnectID field
TR = smi_core.TrackingResults.convertSMDToTR(SMD);

% TR → SMD (flatten back to list)
SMD = smi_core.TrackingResults.convertTRToSMD(TR);
```

## TR Structure Reference

### Array Organization

```matlab
% TR is an array
TR(1)   % First trajectory
TR(2)   % Second trajectory
TR(end) % Last trajectory

% Each element is an SMD structure
N_points_in_traj_5 = length(TR(5).X);
```

### Fields (Same as SMD)

Each TR(n) element contains all standard SMD fields for that trajectory:

**Position fields:**
- `TR(n).X` - X positions (pixels)
- `TR(n).Y` - Y positions (pixels)
- `TR(n).Z` - Z positions (micrometers, if 3D)
- `TR(n).X_SE` - Position uncertainties
- `TR(n).Y_SE`
- `TR(n).Z_SE`

**Photometry:**
- `TR(n).Photons` - Detected photons
- `TR(n).Bg` - Background per pixel
- `TR(n).PSFSigma` - PSF width

**Temporal:**
- `TR(n).FrameNum` - Frame numbers for each localization
- `TR(n).DatasetNum` - Dataset indices

**Quality:**
- `TR(n).PValue` - Fit quality
- `TR(n).ThreshFlag` - Threshold pass/fail

**Tracking:**
- `TR(n).ConnectID` - All have same value (the trajectory ID)

**Metadata (same across trajectory):**
- `TR(n).NFrames` - Total frames in dataset
- `TR(n).PixelSize` - Pixel size
- `TR(n).FrameRate` - Acquisition rate

## Working with TR

### Basic Trajectory Analysis

```matlab
% Count trajectories
N_traj = length(TR);

% Trajectory lengths
for n = 1:length(TR)
    lengths(n) = length(TR(n).X);
end

% Or use built-in method
lengths = smi_core.TrackingResults.computeTrajLengths(TR);

% Trajectory durations (in time units)
durations = smi_core.TrackingResults.computeTrajDurations(TR);
```

### Filtering Trajectories

```matlab
% Filter by minimum length
MinTrackLength = 10; % At least 10 localizations
TR_filtered = smi_core.TrackingResults.threshTrajLength(TR, MinTrackLength);

% Filter by starting frame
MaxStartFrame = 100; % Must start in first 100 frames
TR_early = smi_core.TrackingResults.windowStartTR(TR, MaxStartFrame, 0);

% Filter by frame range
MinFrame = 50;
MaxFrame = 150;
TR_windowed = smi_core.TrackingResults.windowTR(TR, MinFrame, MaxFrame, 0);
```

### Accessing Trajectory Data

```matlab
% Plot first trajectory
plot(TR(1).X, TR(1).Y, 'o-');
title(sprintf('Trajectory 1 (%d points)', length(TR(1).X)));

% Plot all trajectories
figure;
hold on;
for n = 1:length(TR)
    plot(TR(n).X, TR(n).Y, '-');
end
title(sprintf('%d Trajectories', length(TR)));
```

### Trajectory Statistics

```matlab
% Compute MSD for each trajectory
for n = 1:length(TR)
    % Extract positions
    x = TR(n).X;
    y = TR(n).Y;

    % Compute displacements
    dx = diff(x);
    dy = diff(y);

    % Mean squared displacement
    msd(n) = mean(dx.^2 + dy.^2);
end

% Diffusion analysis
% Use smi_stat.DiffusionEstimator for proper MSD analysis
```

### Combining Trajectories

```matlab
% Concatenate two TR arrays
TR_combined = smi_core.TrackingResults.catTR(TR1, TR2, true);

% Join specific trajectories into one
TrajectoryIDs = [5, 12, 18]; % Which trajectories to join
TR_joined = smi_core.TrackingResults.joinTraj(TR, TrajectoryIDs, 0);
```

## Practical Examples

### Example 1: Trajectory Length Distribution

```matlab
% Compute trajectory lengths
lengths = smi_core.TrackingResults.computeTrajLengths(TR);

% Visualize distribution
figure;
histogram(lengths, 50);
xlabel('Trajectory Length (localizations)');
ylabel('Count');
title(sprintf('%d Trajectories, Median Length: %d', ...
    length(TR), median(lengths)));
grid on;

% Statistics
fprintf('Trajectory Statistics:\n');
fprintf('  Total trajectories: %d\n', length(TR));
fprintf('  Mean length: %.1f localizations\n', mean(lengths));
fprintf('  Median length: %d localizations\n', median(lengths));
fprintf('  Min length: %d, Max length: %d\n', min(lengths), max(lengths));
```

### Example 2: Long Trajectory Analysis

```matlab
% Filter for long trajectories
MinLength = 20; % At least 20 points
TR_long = smi_core.TrackingResults.threshTrajLength(TR, MinLength);

fprintf('Long trajectories: %d (%.1f%%)\n', ...
    length(TR_long), 100*length(TR_long)/length(TR));

% Analyze only long trajectories
figure;
subplot(1,2,1);
hold on;
for n = 1:length(TR_long)
    plot(TR_long(n).X, TR_long(n).Y, '-', 'LineWidth', 1.5);
end
title(sprintf('%d Long Trajectories (≥%d points)', length(TR_long), MinLength));
xlabel('X (pixels)'); ylabel('Y (pixels)');
axis equal;

subplot(1,2,2);
lengths_long = smi_core.TrackingResults.computeTrajLengths(TR_long);
histogram(lengths_long, 20);
xlabel('Length (localizations)');
ylabel('Count');
title('Length Distribution');
```

### Example 3: Convert TR to SMD for Visualization

```matlab
% Convert back to SMD for super-resolution imaging
SMD = smi_core.TrackingResults.convertTRToSMD(TR);

% Generate SR image
SMD.XSize = 256;
SMD.YSize = 256;
SR_image = smi_vis.GenerateImages.gaussianImage(SMD, 10, 0);

figure;
imagesc(SR_image);
colormap hot;
axis image;
title(sprintf('Super-Resolution Image (%d trajectories)', length(TR)));
```

### Example 4: Trajectory Duration Analysis

```matlab
% Compute durations in seconds
durations = smi_core.TrackingResults.computeTrajDurations(TR);

% Visualize
figure;
subplot(2,1,1);
histogram(durations, 50);
xlabel('Trajectory Duration (seconds)');
ylabel('Count');
title('Duration Distribution');

subplot(2,1,2);
lengths = smi_core.TrackingResults.computeTrajLengths(TR);
scatter(lengths, durations, 20, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Trajectory Length (localizations)');
ylabel('Duration (seconds)');
title('Length vs Duration');
grid on;

% Statistics
fprintf('Duration Statistics:\n');
fprintf('  Mean duration: %.2f s\n', mean(durations));
fprintf('  Median duration: %.2f s\n', median(durations));
fprintf('  Total observation time: %.1f s\n', sum(durations));
```

## TR Utility Methods

### Trajectory Statistics

```matlab
% Compute trajectory lengths
lengths = smi_core.TrackingResults.computeTrajLengths(TR);
% Returns: vector of integers, one per trajectory

% Compute trajectory durations
durations = smi_core.TrackingResults.computeTrajDurations(TR);
% Returns: vector of durations in seconds

% Compute trajectory fidelity
fidelity = smi_core.TrackingResults.computeTrajFidelity(TR);
% Returns: fraction of frames with localizations
```

### Trajectory Filtering

```matlab
% Filter by minimum length
TR_filtered = smi_core.TrackingResults.threshTrajLength(TR, MinLength);

% Window by frame range
TR_windowed = smi_core.TrackingResults.windowTR(TR, MinFrame, MaxFrame, Verbose);

% Filter by starting frame
TR_early = smi_core.TrackingResults.windowStartTR(TR, MaxStartFrame, Verbose);
```

### Trajectory Manipulation

```matlab
% Get indices for specific trajectories
IDs = [1, 5, 10];
indices = smi_core.TrackingResults.getTRIndex(TR, IDs);

% Join multiple trajectories
TR_joined = smi_core.TrackingResults.joinTraj(TR, IDs, Verbose);

% Pad trajectories (add padding localizations)
TRPadding = struct('Pre', 2, 'Post', 2); % Add 2 points before and after
TR_padded = smi_core.TrackingResults.padTR(TR, TRPadding);

% Concatenate TR arrays
TR_combined = smi_core.TrackingResults.catTR(TR1, TR2, CheckDims);
```

## Common Patterns

### Iterate Over All Trajectories

```matlab
for n = 1:length(TR)
    % Access trajectory data
    x = TR(n).X;
    y = TR(n).Y;
    frames = TR(n).FrameNum;

    % Do analysis
    % ...
end
```

### Filter and Analyze

```matlab
% 1. Filter trajectories
TR_filtered = smi_core.TrackingResults.threshTrajLength(TR, 10);

% 2. Convert to SMD if needed
SMD = smi_core.TrackingResults.convertTRToSMD(TR_filtered);

% 3. Visualize
SR_image = smi_vis.GenerateImages.gaussianImage(SMD, 10, 0);
```

### Extract Single Trajectory

```matlab
% Get trajectory with specific ConnectID
targetID = 42;
idx = find([TR.ConnectID] == targetID); % Note: This won't work directly

% Better: iterate
for n = 1:length(TR)
    if TR(n).ConnectID(1) == targetID % First element is trajectory ID
        traj = TR(n);
        break;
    end
end
```

## Common Issues

**Issue: TR is empty after tracking**

Solutions:
- Check if frame connection succeeded: `length(SMD.ConnectID)`
- Verify ConnectID field exists before conversion
- Check minimum track length threshold wasn't too strict

**Issue: Can't access TR like SMD**

Remember: TR is an **array of structures**, not a single structure.

```matlab
% ✗ Wrong:
positions = [TR.X, TR.Y]; % This concatenates all trajectories!

% ✓ Correct:
positions = [TR(1).X, TR(1).Y]; % Single trajectory
```

**Issue: Trajectory indices don't match ConnectID**

TR array index (1, 2, 3...) is not the same as ConnectID. Use `getTRIndex` to find specific IDs.

## See Also

- [SMD Structure](smd-structure.md) - Base structure for TR
- [SPT Tracking Workflow](../workflows/spt-tracking.md) - How TR is generated
- [Frame Connection](../workflows/frame-connection.md) - Creates trajectories
- smi_core.TrackingResults - Class with all utility methods
- smi_stat.DiffusionEstimator - Analyze diffusion from TR

## Next Steps

- [SPT Tracking Workflow](../workflows/spt-tracking.md) - Complete tracking analysis
- [Frame Connection](../workflows/frame-connection.md) - How trajectories are built
- Use smi_stat.DiffusionEstimator for MSD and diffusion analysis
- Explore smi_stat.HMM for hidden state analysis of trajectories
