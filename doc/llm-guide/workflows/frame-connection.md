---
title: "Frame Connection for SMLM Blinking Analysis"
category: "workflows"
level: "intermediate"
tags: ["frame-connection", "smlm", "blinking", "precision", "lap-fc"]
prerequisites: ["../core-concepts/smf-structure.md", "../core-concepts/smd-structure.md", "smlm-analysis.md"]
related: ["smlm-analysis.md", "spt-tracking.md", "../how-to/localize-molecules.md"]
summary: "Comprehensive guide to frame connection, which links localizations across frames from the same blinking emitter to improve precision and reduce data size"
estimated_time: "20 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Frame Connection for SMLM Blinking Analysis

## Purpose

Frame connection is a critical post-processing step in SMLM analysis that identifies when the same fluorescent emitter appears in multiple frames. By recognizing that several localizations come from a single blinking molecule, frame connection achieves two key goals: improved localization precision through averaging multiple observations, and reduced data volume by combining redundant measurements. This guide explains what frame connection does, when to use it, how to configure it, and how to interpret the results.

## Prerequisites

- Understanding of [SMF structure](../core-concepts/smf-structure.md)
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Completion of [SMLM workflow](smlm-analysis.md)
- Basic understanding of SMLM blinking behavior
- Familiarity with localization precision concepts

## Overview

### What is Frame Connection?

In SMLM imaging, fluorescent emitters exhibit stochastic blinking: they switch between bright (on) and dark (off) states. A single molecule might appear in frame 10, disappear in frames 11-13, reappear in frames 14-16, and then photobleach. Without frame connection, these appearances would be treated as separate emitters, leading to:

- **Overcounting**: One molecule counted as multiple emitters
- **Reduced precision**: Each localization has individual uncertainty
- **Large datasets**: Redundant measurements increase data size
- **Analysis complications**: Difficult to measure true emitter density

Frame connection solves this by:

1. **Identifying** localizations from the same emitter across frames
2. **Linking** them with a unique ConnectID
3. **Combining** multiple observations into a single high-precision localization
4. **Improving** precision by factors of √N where N is the number of linked localizations

### Frame Connection vs SPT Tracking

It's important to distinguish frame connection (for SMLM) from tracking (for SPT):

| Aspect | Frame Connection (SMLM) | Tracking (SPT) |
|--------|------------------------|----------------|
| **Use case** | Blinking emitters, stationary targets | Moving particles |
| **Assumption** | Emitters are fixed in space | Particles diffuse/move |
| **Goal** | Improve precision, count emitters | Follow trajectories, measure dynamics |
| **Output** | Single combined localization per emitter | Time series of positions per particle |
| **Blinking** | Expected and corrected for | Optional, handled during tracking |
| **Frame gaps** | Small gaps (1-10 frames) from blinking | Larger gaps possible from missed detections |
| **Distance threshold** | Small (< 1-2 pixels) | Larger (depends on diffusion) |

**When to use which:**

- **DNA-PAINT, dSTORM, PALM**: Use frame connection (molecules don't move)
- **Single particle tracking**: Use SPT (particles move continuously)
- **Live-cell SMLM**: Depends on target - use frame connection for static structures, SPT for moving proteins

## How Frame Connection Works

### The Frame Connection Pipeline

Frame connection happens after localization but before drift correction in the SMLM pipeline:

```
Localization (LocalizeData)
    ↓
Raw SMD (all localizations)
    ↓
Frame Connection ← YOU ARE HERE
    ↓
Connected SMD (ConnectID assigned)
    ↓
Combined SMD (one per emitter)
    ↓
Drift Correction
```

### Algorithm Overview

Frame connection uses the **LAP-FC** (Linear Assignment Problem Frame Connection) method, which formulates the problem as an optimization:

**Step 1: Pre-clustering**
- Group localizations that are spatially and temporally close
- Uses `MaxFrameGap` and `NSigmaDev` parameters
- Creates initial clusters to reduce computational complexity

**Step 2: Rate parameter estimation**
- Estimates fluorophore kinetics from pre-clustered data:
  - `KOn`: Probability of transitioning from off to on state
  - `KOff`: Probability of transitioning from on to off state
  - `KBleach`: Probability of permanent photobleaching
  - `PMiss`: Probability of missing a detection

**Step 3: Density estimation**
- Computes local emitter density around each cluster
- Uses `NNearestClusters` parameter
- Accounts for crowding effects on connection probabilities

**Step 4: Cost matrix construction**
- For each cluster, builds cost matrix for linking localizations
- Cost reflects probability that two localizations come from same emitter
- Considers:
  - Spatial distance (smaller = lower cost)
  - Temporal gap (larger = higher cost)
  - Photon counts (similar = lower cost)
  - Local density (higher density = higher cost to connect)
  - Kinetic rates (consistent with blinking = lower cost)

**Step 5: Linear assignment problem**
- Solves LAP to find optimal connections
- Each localization assigned to at most one other
- Minimizes total cost across all connections

**Step 6: Localization combining**
- Groups localizations with same ConnectID
- Combines using precision-weighted averaging:
  ```
  X_combined = Σ(X_i / σ_i²) / Σ(1 / σ_i²)
  σ_combined = √(1 / Σ(1 / σ_i²))
  ```
- Resulting precision improved by √N for N linked localizations

## Configuration in SMF

Frame connection is controlled by `SMF.FrameConnection` sub-structure. Here's a complete parameter guide:

### Essential Parameters

```matlab
% Enable/disable frame connection
SMF.FrameConnection.On = true;  % Enable (default: true)
```

**When to disable:**
- Very sparse data where blinking is rare
- Testing localization parameters
- Comparing with and without frame connection

```matlab
% Method selection
SMF.FrameConnection.Method = 'LAP-FC';  % Default and recommended
```

**Available methods:**
- `'LAP-FC'`: Linear assignment problem (most sophisticated, recommended)
- `'Classical'`: Simple distance/time thresholding (fast but less accurate)
- `'Revised Classical'`: Classical with sigma-based thresholds
- `'Hypothesis Test'`: Statistical hypothesis testing approach

```matlab
% Maximum spatial separation for connections
SMF.FrameConnection.MaxSeparation = 1.0;  % pixels (default: 1.0)
```

**Interpretation:**
- Maximum distance between localizations to consider linking
- Typical values: 0.5-2.0 pixels
- **Too small**: Misses legitimate connections, fragments emitters
- **Too large**: Links different emitters, overcombines
- Rule of thumb: 2-3× expected localization precision

```matlab
% Maximum frame gap to bridge
SMF.FrameConnection.MaxFrameGap = 5;  % frames (default: 5)
```

**Interpretation:**
- Maximum number of frames an emitter can be dark while still connected
- Typical values: 2-10 frames
- **Too small**: Breaks legitimate connections, loses blinking emitters
- **Too large**: Increases false connections, computational cost
- Choose based on known blinking statistics of your fluorophore

### Advanced Parameters

```matlab
% Level of significance for hypothesis testing
SMF.FrameConnection.LoS = 0.01;  % 0-1 (default: 0.01)
```

Used primarily for 'Hypothesis Test' method. Lower values are more conservative.

```matlab
% Pre-clustering threshold
SMF.FrameConnection.NSigmaDev = 5;  % sigma (default: 5)
```

**Interpretation:**
- Spatial threshold for initial clustering in units of localization uncertainty
- Localizations within NSigmaDev standard deviations are pre-clustered
- Higher values create larger pre-clusters (more connections considered)
- Lower values create tighter pre-clusters (faster but may miss connections)

```matlab
% Number of nearest clusters for density estimation
SMF.FrameConnection.NNearestClusters = 2;  % integer (default: 2)
```

**Interpretation:**
- How many neighboring clusters used to estimate local density
- Higher values smooth density estimates
- Lower values make density more local
- Typical range: 1-5

```matlab
% Number of iterations for LAP-FC
SMF.FrameConnection.NIterations = 1;  % integer (default: 1)
```

**Interpretation:**
- LAP-FC can iterate, updating rate/density estimates between iterations
- Usually 1 iteration is sufficient
- More iterations improve results slightly at computational cost

```matlab
% Minimum number of frame connections to retain
SMF.FrameConnection.MinNFrameConns = 1;  % integer (default: 1)
```

**Interpretation:**
- Minimum number of localizations an emitter must have to be kept
- Setting to 2+ filters out single-appearance localizations
- Useful for reducing noise and transient binding events
- Typical values: 1-3 for most SMLM, higher for DNA-PAINT

## Running Frame Connection

### Method 1: Automatic (via smi.SMLM)

Frame connection runs automatically when enabled:

```matlab
% Create and configure SMF
SMF = smi_core.SingleMoleculeFitting();

% Configure data and localization parameters
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = {'dSTORM_data.h5'};
SMF.Data.PixelSize = 0.108;  % 108 nm
SMF.BoxFinding.MinPhotons = 250;
SMF.Fitting.PSFSigma = 1.3;

% Configure frame connection
SMF.FrameConnection.On = true;
SMF.FrameConnection.Method = 'LAP-FC';
SMF.FrameConnection.MaxSeparation = 1.0;
SMF.FrameConnection.MaxFrameGap = 5;
SMF.FrameConnection.MinNFrameConns = 2;  % Keep emitters appearing ≥2 times

% Run full analysis
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();

% Results include frame-connected and combined localizations
SMD = SMLMobj.SMD;
```

### Method 2: Standalone Frame Connection

Apply frame connection to existing localization results:

```matlab
% Load raw localizations
load('Results.mat', 'SMD', 'SMF');

% Create FrameConnection object
FC = smi_core.FrameConnection(SMD, SMF);
FC.Verbose = 1;  % Show progress messages

% Perform frame connection
[SMDCombined, SMD] = FC.performFrameConnection();

% SMD now has ConnectID field added
% SMDCombined contains one localization per unique ConnectID
```

### Method 3: Alternative "Function" Style

Quick one-line frame connection:

```matlab
% Automatically runs frame connection and returns results
[~, SMDCombined, SMD] = smi_core.FrameConnection(SMD, SMF, 1);
```

## Understanding the Results

### The ConnectID Field

Frame connection adds a critical field to the SMD structure:

```matlab
% After frame connection
fprintf('Total localizations: %d\n', length(SMD.X));
fprintf('Unique ConnectIDs: %d\n', length(unique(SMD.ConnectID)));

% Example ConnectID values
% SMD.ConnectID = [1, 1, 2, 2, 2, 3, 4, 4, ...]
%                  └─┬─┘  └──┬──┘  │  └─┬─┘
%                    │       │     │    └─ Emitter 4 (2 localizations)
%                    │       │     └────── Emitter 3 (1 localization)
%                    │       └──────────── Emitter 2 (3 localizations)
%                    └──────────────────── Emitter 1 (2 localizations)
```

**Properties:**
- Integer array same length as SMD.X, SMD.Y, etc.
- Unique value assigned to each distinct emitter
- Localizations sharing a ConnectID come from same emitter
- Exact numeric value is arbitrary (only matching matters)

### Finding Connected Localizations

To find all localizations for a specific emitter:

```matlab
% Get all localizations from emitter with ConnectID = 42
emitter_id = 42;
indices = find(SMD.ConnectID == emitter_id);

% Extract positions
X_emitter = SMD.X(indices);
Y_emitter = SMD.Y(indices);
frames = SMD.FrameNum(indices);
photons = SMD.Photons(indices);

% Visualize this emitter's appearances
figure;
plot(frames, X_emitter, 'o-');
xlabel('Frame Number');
ylabel('X Position (pixels)');
title(sprintf('Emitter %d X Position vs Frame', emitter_id));
```

Or use the helper function:

```matlab
% Find localizations connected to first localization in SMDCombined
indices = smi_core.FrameConnection.findConnected(SMDCombined, SMD, 1);

% indices contains positions in SMD arrays for this emitter
X_connected = SMD.X(indices);
Y_connected = SMD.Y(indices);
```

### The SMDCombined Structure

`SMDCombined` contains the precision-enhanced results:

```matlab
% SMDCombined has fewer entries than SMD
fprintf('Original localizations: %d\n', length(SMD.X));
fprintf('Combined localizations: %d\n', length(SMDCombined.X));

% Compression ratio
compression = length(SMD.X) / length(SMDCombined.X);
fprintf('Compression: %.1f:1\n', compression);
% Typical: 3:1 to 10:1 for well-behaved dSTORM/PAINT
```

**Key differences:**
- Length reduced (one entry per emitter, not per localization)
- Precision improved (standard errors smaller)
- Photons summed (total photons across all appearances)
- FrameNum = last frame where emitter appeared
- New field: `NCombined` = number of localizations combined

```matlab
% Check precision improvement
figure;
subplot(1,2,1);
histogram(SMD.X_SE * SMD.PixelSize * 1000, 50);
xlabel('X Precision (nm)'); ylabel('Count');
title('Before Frame Connection');

subplot(1,2,2);
histogram(SMDCombined.X_SE * SMDCombined.PixelSize * 1000, 50);
xlabel('X Precision (nm)'); ylabel('Count');
title('After Frame Connection');
```

### The NCombined Field

This field shows how many localizations contributed to each combined entry:

```matlab
% Distribution of combination counts
histogram(SMDCombined.NCombined, max(SMDCombined.NCombined));
xlabel('Number of Localizations Combined');
ylabel('Count');
title('Frame Connection Statistics');

% Summary statistics
fprintf('Mean localizations per emitter: %.1f\n', mean(SMDCombined.NCombined));
fprintf('Median: %.1f\n', median(SMDCombined.NCombined));
fprintf('Max: %d\n', max(SMDCombined.NCombined));

% Find "super-blinkers" that appeared many times
super_blinkers = find(SMDCombined.NCombined > 20);
fprintf('Found %d emitters appearing >20 times\n', length(super_blinkers));
```

## Practical Examples

### Example 1: Basic Frame Connection

```matlab
% Configure for dSTORM imaging
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/experiments/2024-01-15';
SMF.Data.FileName = {'Cell1_Actin_dSTORM.h5'};
SMF.Data.PixelSize = 0.108;  % 108 nm pixel size
SMF.Data.FrameRate = 100;    % 100 Hz acquisition

% Localization settings
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 300;
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Frame connection for dSTORM (moderate blinking)
SMF.FrameConnection.On = true;
SMF.FrameConnection.Method = 'LAP-FC';
SMF.FrameConnection.MaxSeparation = 1.0;   % 108 nm threshold
SMF.FrameConnection.MaxFrameGap = 5;       % Bridge 5-frame gaps
SMF.FrameConnection.MinNFrameConns = 2;    % Keep emitters seen ≥2 times

% Run analysis
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();

% Check results
SMD = SMLMobj.SMD;
fprintf('=== Frame Connection Results ===\n');
fprintf('Total localizations: %d\n', length(SMD.X));
unique_emitters = length(unique(SMD.ConnectID));
fprintf('Unique emitters: %d\n', unique_emitters);
fprintf('Compression ratio: %.2f:1\n', length(SMD.X) / unique_emitters);
```

### Example 2: DNA-PAINT Optimization

DNA-PAINT has longer, more frequent blinking:

```matlab
% DNA-PAINT specific settings
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/data/PAINT/2024-02-20';
SMF.Data.FileName = {'Target_PAINT_10k_frames.h5'};
SMF.Data.PixelSize = 0.130;  % 130 nm pixels
SMF.Data.FrameRate = 10;     % 10 Hz (slower acquisition)

% Localization
SMF.BoxFinding.MinPhotons = 500;  % PAINT is typically brighter
SMF.Fitting.PSFSigma = 1.2;

% Frame connection for DNA-PAINT (aggressive)
SMF.FrameConnection.On = true;
SMF.FrameConnection.MaxSeparation = 1.5;   % Slightly larger (195 nm)
SMF.FrameConnection.MaxFrameGap = 15;      % Bridge longer gaps
SMF.FrameConnection.MinNFrameConns = 5;    % Keep only stable binders
SMF.FrameConnection.NIterations = 2;       # Refine connections

% Run
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();
```

### Example 3: Comparing With and Without

See the impact of frame connection:

```matlab
% Run localization once
SMF = smi_core.SingleMoleculeFitting();
% ... configure data and localization parameters ...

% Localize data
SMLMobj = smi.SMLM(SMF);
SMF.FrameConnection.On = false;  % Disable for first run
SMLMobj.SMF = SMF;
SMLMobj.analyzeAll();
SMD_no_FC = SMLMobj.SMD;

% Apply frame connection
SMF.FrameConnection.On = true;
SMF.FrameConnection.MaxSeparation = 1.0;
SMF.FrameConnection.MaxFrameGap = 5;
FC = smi_core.FrameConnection(SMD_no_FC, SMF);
[SMDCombined, SMD_with_FC] = FC.performFrameConnection();

% Compare results
figure('Position', [100, 100, 1200, 400]);

% Original localizations
subplot(1,3,1);
smi_vis.GenerateImages.gaussImage([SMD_no_FC.X, SMD_no_FC.Y], ...
    SMD_no_FC.PixelSize, 20, [SMD_no_FC.YSize, SMD_no_FC.XSize]);
title('Without Frame Connection');
axis image;

% Combined localizations
subplot(1,3,2);
smi_vis.GenerateImages.gaussImage([SMDCombined.X, SMDCombined.Y], ...
    SMDCombined.PixelSize, 20, [SMDCombined.YSize, SMDCombined.XSize]);
title('With Frame Connection');
axis image;

% Precision comparison
subplot(1,3,3);
hold on;
histogram(SMD_no_FC.X_SE * SMD_no_FC.PixelSize * 1000, 50, ...
    'DisplayName', 'Without FC');
histogram(SMDCombined.X_SE * SMDCombined.PixelSize * 1000, 50, ...
    'DisplayName', 'With FC');
xlabel('X Precision (nm)');
ylabel('Count');
legend();
title('Precision Improvement');

% Statistics
fprintf('=== Comparison ===\n');
fprintf('Localizations - No FC: %d\n', length(SMD_no_FC.X));
fprintf('Localizations - With FC: %d\n', length(SMDCombined.X));
fprintf('Reduction: %.1f%%\n', ...
    100 * (1 - length(SMDCombined.X)/length(SMD_no_FC.X)));
fprintf('\nMedian precision - No FC: %.1f nm\n', ...
    median(SMD_no_FC.X_SE * SMD_no_FC.PixelSize * 1000));
fprintf('Median precision - With FC: %.1f nm\n', ...
    median(SMDCombined.X_SE * SMDCombined.PixelSize * 1000));
fprintf('Improvement: %.1fx\n', ...
    median(SMD_no_FC.X_SE) / median(SMDCombined.X_SE));
```

### Example 4: Analyzing Blinking Statistics

Extract information about fluorophore behavior:

```matlab
% Load frame-connected results
load('Results.mat', 'SMD', 'SMDCombined', 'SMF');

% Analyze blinking statistics
N_combined = SMDCombined.NCombined;

% Distribution of appearances
figure;
subplot(2,2,1);
histogram(N_combined, max(N_combined));
xlabel('Appearances per Emitter');
ylabel('Count');
title('Blinking Distribution');

% On-time analysis
subplot(2,2,2);
on_times = N_combined / SMF.Data.FrameRate;  % Convert to seconds
histogram(on_times, 50);
xlabel('Total On-Time (s)');
ylabel('Count');
title('Cumulative On-Time per Emitter');

% Frame gaps (off-times)
subplot(2,2,3);
gaps = [];
for id = unique(SMD.ConnectID)'
    frames = SMD.FrameNum(SMD.ConnectID == id);
    if length(frames) > 1
        gaps = [gaps; diff(sort(frames)) - 1];  % Gap length
    end
end
histogram(gaps, 0:max(gaps));
xlabel('Gap Length (frames)');
ylabel('Count');
title('Off-Time Distribution');

% Photon statistics
subplot(2,2,4);
scatter(SMDCombined.NCombined, SMDCombined.Photons, 10, 'filled', 'MarkerFaceAlpha', 0.3);
xlabel('Number of Appearances');
ylabel('Total Photons');
title('Photon Budget vs Blinking');

% Summary statistics
fprintf('=== Blinking Statistics ===\n');
fprintf('Mean appearances: %.1f\n', mean(N_combined));
fprintf('Single-appearance emitters: %d (%.1f%%)\n', ...
    sum(N_combined == 1), 100 * sum(N_combined == 1) / length(N_combined));
fprintf('Mean gap length: %.1f frames (%.0f ms)\n', ...
    mean(gaps), mean(gaps) * 1000 / SMF.Data.FrameRate);
fprintf('Mean total photons: %.0f\n', mean(SMDCombined.Photons));
```

## Troubleshooting

### Issue: Too Many Single-Appearance Localizations

**Diagnose:**
```matlab
single_appearances = sum(SMDCombined.NCombined == 1);
total = length(SMDCombined.NCombined);
fprintf('Single appearances: %d / %d (%.1f%%)\n', ...
    single_appearances, total, 100 * single_appearances / total);
% Problematic if >50% for dSTORM, >30% for PAINT
```

**Solutions:**
- **Increase MaxSeparation**: Connections too restrictive spatially
- **Increase MaxFrameGap**: Missing connections due to long dark times
- **Check localization quality**: Poor precision prevents accurate connection
- **Reduce MinNFrameConns**: Currently filtering them out post-connection
- **Verify fluorophore**: Some dyes genuinely blink once then bleach

### Issue: Over-connection (Compression Too High)

**Diagnose:**
```matlab
compression = length(SMD.X) / length(unique(SMD.ConnectID));
fprintf('Compression ratio: %.1f:1\n', compression);
% Problematic if >20:1 for dSTORM, >50:1 for PAINT
```

**Solutions:**
- **Decrease MaxSeparation**: Connecting distinct nearby emitters
- **Decrease MaxFrameGap**: Bridging too-long gaps, merging different emitters
- **Check for drift**: Uncorrected drift causes false connections
- **Reduce NSigmaDev**: Pre-clustering too aggressively
- **Verify localization**: Very poor precision causes spurious connections

### Issue: Inconsistent Positions Within Connected Groups

**Diagnose:**
```matlab
% Check position variation for highly-connected emitters
for id = unique(SMD.ConnectID)'
    indices = find(SMD.ConnectID == id);
    if length(indices) > 10  % Look at emitters with many appearances
        X_std = std(SMD.X(indices)) * SMD.PixelSize * 1000;  % nm
        Y_std = std(SMD.Y(indices)) * SMD.PixelSize * 1000;  % nm
        if X_std > 50 || Y_std > 50  % More than 50nm variation
            fprintf('ConnectID %d: X_std=%.1f nm, Y_std=%.1f nm\n', ...
                id, X_std, Y_std);
        end
    end
end
```

**Solutions:**
- **Decrease MaxSeparation**: Linking drifting or moving structures
- **Check sample stability**: Physical drift during acquisition
- **Apply drift correction first**: Run pipeline in correct order
- **Verify sample is static**: Frame connection assumes immobile emitters
- **Check for molecular motion**: Use SPT instead if targets move

### Issue: Frame Connection Too Slow

**Diagnose:**
Frame connection scales with localization density and MaxFrameGap.

**Solutions:**
- **Reduce MaxFrameGap**: Largest impact on computational cost
- **Lower NSigmaDev**: Smaller pre-clusters process faster
- **Reduce NIterations**: Usually 1 iteration sufficient
- **Use 'Classical' method**: Faster but less sophisticated
- **Filter before frame connection**: Apply stricter thresholding first
- **Process smaller ROIs**: Divide and conquer approach

### Issue: Poor Precision Improvement

**Diagnose:**
```matlab
improvement = median(SMD.X_SE) / median(SMDCombined.X_SE);
fprintf('Precision improvement: %.2fx\n', improvement);
% Expected: 1.5x-3x for typical blinking
```

**Solutions:**
- **Not enough connections**: See "Too Many Single-Appearance" above
- **Check NCombined**: If mean < 3, limited improvement possible
- **Verify localization precision**: If precision already excellent, little room to improve
- **Increase MinNFrameConns**: Filters weakly-connected emitters
- **Check photon counts**: Low photons limit fundamental precision

## Best Practices

### Parameter Selection Guidelines

**Conservative (high precision, fewer connections):**
```matlab
SMF.FrameConnection.MaxSeparation = 0.75;  % Tight spatial
SMF.FrameConnection.MaxFrameGap = 3;       # Short gaps
SMF.FrameConnection.MinNFrameConns = 3;    % Require multiple appearances
```

**Aggressive (more connections, risk over-combination):**
```matlab
SMF.FrameConnection.MaxSeparation = 1.5;   % Generous spatial
SMF.FrameConnection.MaxFrameGap = 10;      % Long gaps
SMF.FrameConnection.MinNFrameConns = 1;    # Keep all
```

**Standard dSTORM starting point:**
```matlab
SMF.FrameConnection.MaxSeparation = 1.0;
SMF.FrameConnection.MaxFrameGap = 5;
SMF.FrameConnection.MinNFrameConns = 2;
```

**DNA-PAINT starting point:**
```matlab
SMF.FrameConnection.MaxSeparation = 1.5;
SMF.FrameConnection.MaxFrameGap = 15;
SMF.FrameConnection.MinNFrameConns = 5;
```

### Optimization Workflow

1. **Start with defaults**: Run initial analysis with standard parameters
2. **Check compression ratio**: Should be 3:1 to 10:1 for dSTORM
3. **Examine NCombined distribution**: Most emitters should have >2 appearances
4. **Verify precision improvement**: Expect 1.5x-3x improvement
5. **Visualize super-resolution images**: Check for artifacts or over-smoothing
6. **Adjust parameters**: Iterate based on observations
7. **Document final settings**: Record parameters for reproducibility

### Quality Control Checks

Always verify frame connection worked correctly:

```matlab
% QC Script
fprintf('=== Frame Connection QC ===\n');

% 1. Compression ratio
compression = length(SMD.X) / length(unique(SMD.ConnectID));
fprintf('Compression: %.1f:1 ', compression);
if compression < 2
    fprintf('[LOW - consider more permissive settings]\n');
elseif compression > 20
    fprintf('[HIGH - consider stricter settings]\n');
else
    fprintf('[OK]\n');
end

% 2. Single appearances
single_pct = 100 * sum(SMDCombined.NCombined == 1) / length(SMDCombined.NCombined);
fprintf('Single appearances: %.1f%% ', single_pct);
if single_pct > 50
    fprintf('[HIGH - check MaxFrameGap and MaxSeparation]\n');
else
    fprintf('[OK]\n');
end

% 3. Precision improvement
improvement = median(SMD.X_SE) / median(SMDCombined.X_SE);
fprintf('Precision improvement: %.2fx ', improvement);
if improvement < 1.3
    fprintf('[LOW - limited benefit from frame connection]\n');
else
    fprintf('[OK]\n');
end

% 4. Mean appearances
mean_appearances = mean(SMDCombined.NCombined);
fprintf('Mean appearances: %.1f ', mean_appearances);
if mean_appearances < 2
    fprintf('[LOW - limited blinking]\n');
else
    fprintf('[OK]\n');
end
```

## Advanced Topics

### Custom Frame Connection Methods

While LAP-FC is recommended, other methods exist:

```matlab
% Classical method (fast, simple thresholding)
SMF.FrameConnection.Method = 'Classical';
% Uses only MaxSeparation and MaxFrameGap
% No rate estimation or density calculation
% Faster for large datasets

% Revised Classical (classical with sigma-based thresholds)
SMF.FrameConnection.Method = 'Revised Classical';
% Similar to Classical but thresholds scale with localization precision
% Good middle ground between speed and sophistication
```

### Integration with Drift Correction

Frame connection happens before drift correction in the pipeline:

```matlab
% Correct order in smi.SMLM:
% 1. Localization
% 2. Frame Connection
% 3. Drift Correction (intra-dataset)
% 4. Concatenate datasets
% 5. Drift Correction (inter-dataset)
```

This order is important because:
- Drift correction needs stable reference points (frame-connected emitters)
- Pre-connection drift would artificially increase separation between localizations
- Post-connection drift correction gives better alignment

### Accessing Internal Parameters

The LAP-FC algorithm estimates kinetic parameters:

```matlab
FC = smi_core.FrameConnection(SMD, SMF);
[SMDCombined, SMD] = FC.performFrameConnection();

% Access estimated parameters
params = FC.InternalParams;
fprintf('Estimated KOn: %.3f (prob on→off per frame)\n', params.KOn);
fprintf('Estimated KOff: %.3f (prob off→on per frame)\n', params.KOff);
fprintf('Estimated KBleach: %.3f (prob bleach per frame)\n', params.KBleach);
fprintf('Estimated PMiss: %.3f (prob miss detection)\n', params.PMiss);
fprintf('Estimated number of emitters: %d\n', params.NEmitters);

% These can provide insight into fluorophore behavior
```

## Citations and Further Reading

Frame connection in smite uses the LAP-FC algorithm described in:

**Primary Citation:**
David J. Schodt and Keith A. Lidke, "Spatiotemporal Clustering of Repeated Super-Resolution Localizations via Linear Assignment Problem", *Frontiers in Bioinformatics*, 2021
https://doi.org/10.3389/fbinf.2021.724325

This paper provides:
- Complete mathematical formulation
- Comparison with other methods
- Validation on simulated and experimental data
- Performance benchmarks

## See Also

- [SMLM Analysis Workflow](smlm-analysis.md) - Complete pipeline including frame connection
- [SPT Tracking Workflow](spt-tracking.md) - Tracking moving particles (different from frame connection)
- [SMF Structure](../core-concepts/smf-structure.md) - All SMF.FrameConnection parameters
- [SMD Structure](../core-concepts/smd-structure.md) - Understanding ConnectID field
- MATLAB/+smi_core/@FrameConnection/README.md - Technical class documentation
- MATLAB/+smi_core/@FrameConnection/FrameConnection.m - Source code
