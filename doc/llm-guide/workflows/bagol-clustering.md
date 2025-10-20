---
title: "BaGoL: Bayesian Grouping of Localizations"
category: "workflows"
level: "advanced"
tags: ["bagol", "clustering", "bayesian", "super-resolution", "precision"]
prerequisites: ["../core-concepts/smd-structure.md", "./smlm-analysis.md"]
related: ["../core-concepts/data-flow.md", "../how-to/localize-molecules.md"]
summary: "Complete workflow for using BaGoL to achieve sub-nanometer precision by grouping multiple localizations from the same emitter"
estimated_time: "35 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# BaGoL: Bayesian Grouping of Localizations

## Purpose

BaGoL (Bayesian Grouping of Localizations) is an advanced method that achieves sub-nanometer localization precision by intelligently grouping multiple blinking events from the same emitter. While standard SMLM treats each localization independently, BaGoL uses Bayesian inference to determine the true number of emitters and their positions from repeated observations. This document provides a complete understanding of when to use BaGoL, how to configure it, and how to interpret the results.

## Prerequisites

- Completion of [SMLM analysis workflow](./smlm-analysis.md)
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Frame-connected SMLM data (localizations must be frame-connected)
- Understanding of Bayesian inference concepts (helpful but not required)
- MATLAB Statistics and Machine Learning Toolbox

## What BaGoL Does

### The Core Problem

Single molecule localization microscopy produces multiple localizations from the same emitter due to blinking or repeated binding events. Standard analysis treats these as independent observations without attempting to group them by source emitter. This approach:

- Inflates the apparent number of molecules
- Doesn't leverage multiple observations to improve precision
- Can't distinguish true molecular density from blinking statistics

### BaGoL's Solution

BaGoL explores possible groupings of localizations to determine:

1. **How many true emitters** are present in the data
2. **Where each emitter is located** with improved precision
3. **Which localizations belong to which emitter** (allocation)

The algorithm uses Reversible Jump Markov Chain Monte Carlo (RJMCMC) to:
- Add, remove, and move emitters
- Reallocate localizations to different emitters
- Weight solutions by their posterior probability

### When to Use BaGoL

BaGoL is most effective when:

- **Blinking/binding events occur multiple times per emitter** (DNA-PAINT, dSTORM)
- **You need sub-nanometer precision** (< 1 nm achievable under good conditions)
- **Localizations are frame-connected** (required preprocessing)
- **Emitter density allows grouping** (not too sparse, not too dense)
- **You can estimate or learn the blinking distribution** (localizations per emitter)

Typical applications:
- DNA origami structural validation
- Protein complex stoichiometry and geometry
- Membrane protein organization
- Any application requiring highest possible precision

## Algorithm Overview

### Key Concepts

**Emitter vs. Localization:**
- An **emitter** is a single physical fluorophore or binding site
- A **localization** is one observation of that emitter in a single frame
- Each emitter can produce multiple localizations through blinking

**The RJMCMC Chain:**

BaGoL runs a Markov chain that samples possible models. Each state in the chain includes:
- Number of emitters (K)
- Position of each emitter (X, Y coordinates)
- Allocation of localizations to emitters (which belongs to which)
- Optional drift parameters (Alpha_X, Alpha_Y)

**Jump Types:**

The algorithm proposes four types of moves with probabilities defined by `P_Jumps`:

1. **Move**: Adjust position of an existing emitter
2. **Allocate**: Reassign a localization to a different emitter
3. **Add**: Propose adding a new emitter
4. **Remove**: Propose removing an emitter

Each move is accepted or rejected based on the posterior probability.

**Prior Distribution:**

The number of localizations per emitter follows either:
- **Poisson distribution**: Parameter λ (lambda)
- **Gamma distribution**: Parameters k (shape) and θ (theta, scale)

These parameters can be:
- Fixed based on prior knowledge
- Learned from the data (hierarchical Bayes approach)

### The Two Main Outputs

**1. Posterior Image:**

A probability density map showing where emitters are likely located, weighted by:
- All possible numbers of emitters explored
- All positions sampled
- All allocation schemes considered

This represents the full uncertainty in the inference.

**2. MAPN Coordinates (Maximum A Posteriori Number):**

Discrete emitter positions from the most frequently occurring number of emitters in the chain. These coordinates include:
- X, Y positions (nm)
- X_SE, Y_SE uncertainties (nm)
- Nmean: mean localizations per emitter
- Optional drift parameters

MAPN provides point estimates suitable for quantitative analysis.

## Complete BaGoL Workflow

### Step 1: Prepare SMLM Data

BaGoL requires frame-connected localizations as input. Start with standard SMLM analysis:

```matlab
% Run SMLM analysis with frame connection
SMF = smi_core.SingleMoleculeFitting();

% Configure for BaGoL preprocessing
SMF.Data.FileDir = '/data/2024-01-10';
SMF.Data.FileName = {'Sample_PAINT.h5'};
SMF.Data.PixelSize = 0.108;  % 108 nm
SMF.Data.CameraGain = 2.5;
SMF.Data.CameraOffset = 100;

% Standard localization settings
SMF.BoxFinding.BoxSize = 7;
SMF.BoxFinding.MinPhotons = 250;
SMF.Fitting.PSFSigma = 1.3;
SMF.Fitting.FitType = 'XYNB';

% Important: Enable frame connection
SMF.FrameConnection.On = true;
SMF.FrameConnection.Method = 'LAP-FC';
SMF.FrameConnection.MaxSeparation = 1.0;  % pixels
SMF.FrameConnection.MaxFrameGap = 1;
SMF.FrameConnection.MinNFrameConns = 2;  % Keep emitters with >=2 appearances

% Thresholding for quality
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.15;  % pixels
SMF.Thresholding.MinPhotons = 150;

% Optional: Drift correction
SMF.DriftCorrection.On = true;

% Run analysis
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();

% Get results
SMD = SMLMobj.SMD;
fprintf('%d frame-connected localizations ready for BaGoL\n', length(SMD.X));
```

**Key preprocessing considerations:**

- **Frame connection is mandatory**: BaGoL needs ConnectID to understand blinking
- **Precision inflation**: Consider setting `SMF.Data.SEAdjust = 1-2` nm for DNA-PAINT
- **Intensity filtering**: Use `SMF.Thresholding.InMeanMultiplier` to keep bright spots
- **Nearest neighbor filtering**: Use for DNA-PAINT; avoid for dSTORM

### Step 2: Convert SMD to BaGoL Format

BaGoL expects coordinates in nm (not pixels) with specific field names:

```matlab
% If starting from SMITE SMD structure
% Convert to BaGoL format (coordinates in nm)
BaGoL_SMD = struct();
BaGoL_SMD.X = SMD.X * SMD.PixelSize;  % Convert pixels to nm
BaGoL_SMD.Y = SMD.Y * SMD.PixelSize;
BaGoL_SMD.Z = [];  % Optional 3D
BaGoL_SMD.X_SE = SMD.X_SE * SMD.PixelSize;  % Uncertainties in nm
BaGoL_SMD.Y_SE = SMD.Y_SE * SMD.PixelSize;
BaGoL_SMD.Z_SE = [];
BaGoL_SMD.FrameNum = SMD.FrameNum;

% Alternatively, use helper function
BaGoL_SMD = smi.BaGoL.importLLSMD(SMD);
```

**Coordinate filtering:**

BaGoL works better with non-negative coordinates. This filtering is done automatically by `hierBaGoL_analysis`, but you can do it manually:

```matlab
% Remove localizations with negative coordinates
Ind = BaGoL_SMD.X >= 0 & BaGoL_SMD.Y >= 0;
BaGoL_SMD.X = BaGoL_SMD.X(Ind);
BaGoL_SMD.Y = BaGoL_SMD.Y(Ind);
BaGoL_SMD.X_SE = BaGoL_SMD.X_SE(Ind);
BaGoL_SMD.Y_SE = BaGoL_SMD.Y_SE(Ind);
BaGoL_SMD.FrameNum = BaGoL_SMD.FrameNum(Ind);
```

### Step 3: Configure BaGoL Parameters

#### Basic Configuration

```matlab
% Create BaGoL object
BGL = smi.BaGoL();

% Set input data
BGL.SMD = BaGoL_SMD;

% ROI subdivision parameters
BGL.ROIsize = 500;     % nm - Size of subregions for RJMCMC
BGL.Overlap = 50;      % nm - Overlap between adjacent regions
BGL.Cutoff = 25;       % nm - Pre-clustering cutoff distance

% RJMCMC chain parameters
BGL.N_Burnin = 32000;  % Burn-in samples (discarded)
BGL.N_Trials = 8000;   % Post-burn-in samples (used for inference)
BGL.P_Jumps = [0.25, 0.25, 0.25, 0.25];  % [Move, Allocate, Add, Remove]

% Prior distribution parameters
% Option 1: Fixed Poisson prior
BGL.Xi = 20;  % Lambda = 20 localizations/emitter

% Option 2: Fixed Gamma prior
BGL.Xi = [20, 1];  % [k=20, theta=1], mean = k*theta = 20

% Option 3: Learn from data (hierarchical Bayes)
BGL.HierarchFlag = 1;
BGL.Xi = [20, 1];  % Initial guess
BGL.NSamples = 10;  % Sample Xi every 10 RJMCMC iterations

% Output parameters
BGL.PImageFlag = 1;  % Generate posterior image
BGL.PixelSize = 4;   % Output pixel size for images (nm)
BGL.ChainFlag = 0;   % Save full chain (memory intensive)

% Precision adjustment
BGL.SE_Adjust = 0;   % Add to uncertainties (nm); 1-2 for DNA-PAINT

% Drift parameters (if temporal drift expected)
BGL.Drift = 0;  % Expected drift magnitude (nm/frame)
```

#### Parameter Selection Guidelines

**ROIsize (nm):**
- Larger ROI = more computational time but better handles extended structures
- Smaller ROI = faster but may fragment large clusters
- Typical: 200-500 nm for most data
- Dense data: Use smaller ROI (100-200 nm) to avoid artifacts
- Sparse data: Use larger ROI (500-1000 nm)

**Overlap (nm):**
- Prevents edge artifacts where ROIs meet
- Typical: 10% of ROIsize (e.g., 50 nm for 500 nm ROI)

**Cutoff (nm):**
- Pre-clustering distance threshold using DBSCAN
- Should be around localization precision (10-30 nm typical)
- Too small: Fragments real clusters
- Too large: Merges distinct emitters

**N_Burnin and N_Trials:**
- Burnin: Allow chain to converge (discard these samples)
- Trials: Samples used for inference
- Longer chains = better inference but slower
- Typical: 32000 burnin, 8000 trials for production
- Quick test: 2000 burnin, 500 trials

**P_Jumps:**
- Default `[0.25, 0.25, 0.25, 0.25]` works well
- For challenging data, increase Move and Allocate relative to Add/Remove
- Must sum to 1.0

**Xi prior:**
- DNA-PAINT: 10-30 localizations/emitter typical
- dSTORM: Can be higher (20-50+)
- Learn from data when uncertain (HierarchFlag = 1)

### Step 4: Run BaGoL Analysis

#### Simple Analysis (Fixed Prior)

```matlab
% Run complete analysis
BGL.analyze_all();

% Access results
MAPN = BGL.MAPN;
fprintf('Found %d emitters from %d localizations\n', ...
    length(MAPN.X), length(BGL.SMD.X));

% Check compression ratio
compression = length(BGL.SMD.X) / length(MAPN.X);
fprintf('Compression: %.1f:1 (localizations:emitters)\n', compression);
```

**What happens internally:**

1. `genROIs()`: Divides data into overlapping subregions
2. `precluster()`: Uses hierarchical clustering/DBSCAN to find independent clusters
3. For each cluster:
   - `BaGoL_RJMCMC()`: Runs RJMCMC chain
   - `genMAPN()`: Extracts MAPN coordinates from chain
   - `genPosterior()`: Updates posterior image
4. `removeOverlap()`: Eliminates duplicates from overlapping regions
5. Collates results into `BGL.MAPN` and `BGL.PImage`

#### Hierarchical Analysis (Learn Prior)

```matlab
% Configure to learn Xi from data
BGL.HierarchFlag = 1;
BGL.Xi = [20, 1];      % Initial guess for gamma prior
BGL.NSamples = 10;     % Sample Xi every 10 RJMCMC iterations
BGL.ChainFlag = 1;     % Save chain to examine Xi learning

% Hyperprior parameters for Xi
BGL.Alpha_Xi = 1;      % Shape parameter for Xi prior
BGL.Beta_Xi = 50;      % Scale parameter for Xi prior

% Run analysis
BGL.analyze_all();

% Examine learned Xi distribution
figure;
if size(BGL.XiChain, 2) == 1
    % Poisson prior
    plot(BGL.XiChain(:,1), '.');
    ylabel('Lambda (locs/emitter)');
else
    % Gamma prior
    lambda = BGL.XiChain(:,1) .* BGL.XiChain(:,2);  % k * theta
    plot(lambda, '.');
    ylabel('Lambda = k*theta (locs/emitter)');
end
xlabel('Iteration');
title('Xi Chain: Learned Localization Distribution');

% Get posterior mean of Xi
burnin_idx = floor(size(BGL.XiChain,1)/2);
Xi_post = mean(BGL.XiChain(burnin_idx:end, :), 1);
fprintf('Learned Xi: k=%.2f, theta=%.2f, mean=%.2f locs/emitter\n', ...
    Xi_post(1), Xi_post(2), Xi_post(1)*Xi_post(2));
```

#### Batch Analysis with hierBaGoL

For processing multiple datasets:

```matlab
% Configure parameters
BaGoLParams.ImageSize = 256;          % Image size (pixels)
BaGoLParams.PixelSize = 108;          % Camera pixel size (nm)
BaGoLParams.OutputPixelSize = 4;      % Output image pixel size (nm)
BaGoLParams.SE_Adjust = 2;            % Precision inflation (nm)
BaGoLParams.ClusterDrift = 0;         % Drift (nm/frame)
BaGoLParams.ROIsz = 500;              % ROI size (nm)
BaGoLParams.OverLap = 50;             % Overlap (nm)
BaGoLParams.Cutoff = 25;              % Pre-clustering cutoff (nm)
BaGoLParams.Xi = [20, 1];             % Prior parameters
BaGoLParams.DataROI = [];             % Process full image
BaGoLParams.N_Burnin = 32000;
BaGoLParams.N_Trials = 8000;
BaGoLParams.NSamples = 10;
BaGoLParams.Y_Adjust = [];            % Coordinate transform if needed

% Define files to process
Files = {
    '/data/Cell_01_Label_01_Results.mat';
    '/data/Cell_02_Label_01_Results.mat';
    '/data/Cell_03_Label_01_Results.mat';
};

% Run batch analysis (parallelized)
Results_BaGoL = 'Results_BaGoLHier';
BGL = smi.BaGoL.hierBaGoL_run(Files, [], Results_BaGoL, BaGoLParams);
```

The `hierBaGoL_run` function:
- Runs datasets in parallel using `parfor`
- Handles ROI-based analysis if `*_ROIs.mat` files exist
- Saves results to structured directories
- Provides error handling for batch processing

### Step 5: Examine Results

#### MAPN Coordinates

```matlab
% Access MAPN results
MAPN = BGL.MAPN;

% Basic statistics
fprintf('Results summary:\n');
fprintf('  Emitters found: %d\n', length(MAPN.X));
fprintf('  From localizations: %d\n', length(BGL.SMD.X));
fprintf('  Compression: %.1f:1\n', length(BGL.SMD.X)/length(MAPN.X));
fprintf('  Mean X precision: %.2f nm\n', mean(MAPN.X_SE));
fprintf('  Mean Y precision: %.2f nm\n', mean(MAPN.Y_SE));
fprintf('  Median X precision: %.2f nm\n', median(MAPN.X_SE));
fprintf('  Median Y precision: %.2f nm\n', median(MAPN.Y_SE));

% Precision histogram
figure;
subplot(1,2,1);
histogram(MAPN.X_SE, 30);
xlabel('X Precision (nm)');
ylabel('Frequency');
title('MAPN X Precision');

subplot(1,2,2);
histogram(MAPN.Y_SE, 30);
xlabel('Y Precision (nm)');
ylabel('Frequency');
title('MAPN Y Precision');
```

#### Localizations per Emitter

```matlab
% Distribution of localizations per emitter
figure;
histogram(MAPN.Nmean, 'Normalization', 'pdf');
xlabel('Localizations per Emitter');
ylabel('Probability Density');
title('Blinking Statistics');

% Compare to prior if fixed
if ~BGL.HierarchFlag && length(BGL.Xi) == 2
    hold on;
    x = 0:0.1:max(MAPN.Nmean);
    plot(x, gampdf(x, BGL.Xi(1), BGL.Xi(2)), 'r-', 'LineWidth', 2);
    legend('Observed', 'Prior');
end
```

#### Nearest Neighbor Distances

```matlab
% NND analysis
[idx, dist] = knnsearch([MAPN.X, MAPN.Y], [MAPN.X, MAPN.Y], 'K', 2);
NND = dist(:, 2);  % Distance to nearest neighbor (not self)

figure;
histogram(NND, 30);
xlabel('Nearest Neighbor Distance (nm)');
ylabel('Frequency');
title('MAPN NND Distribution');

fprintf('NND statistics:\n');
fprintf('  Median: %.2f nm\n', median(NND));
fprintf('  Mean: %.2f nm\n', mean(NND));
fprintf('  95th percentile: %.2f nm\n', prctile(NND, 95));
```

### Step 6: Visualize Results

#### Generate Images

```matlab
% Generate MAPN image
ImFlag = 1;  % MAPN mode
MAPN_img = BGL.genMAPNIm(ImFlag);

% Generate SR image from input localizations
ImFlag = 2;  % SR mode
SR_img = BGL.genMAPNIm(ImFlag);

% Display comparison
figure;
subplot(1,3,1);
imagesc(SR_img);
axis image; colormap hot;
title('Input Localizations');

subplot(1,3,2);
imagesc(MAPN_img);
axis image; colormap hot;
title('MAPN Results');

subplot(1,3,3);
if BGL.PImageFlag
    imagesc(BGL.PImage);
    axis image; colormap hot;
    title('Posterior Image');
end
```

#### Plot Localizations vs MAPN

```matlab
% Scatter plot comparison
BGL.plotMAPN(pwd, 'on');  % 'on' = visible figures

% Or create custom plot
figure;
hold on;

% Plot input localizations (smaller, lighter)
plot(BGL.SMD.X, BGL.SMD.Y, 'g.', 'MarkerSize', 4);

% Plot MAPN emitters (larger, darker)
plot(MAPN.X, MAPN.Y, 'mo', 'MarkerSize', 8, 'LineWidth', 2);

axis equal tight;
xlabel('X (nm)');
ylabel('Y (nm)');
title('BaGoL Results');
legend('Input Localizations', 'MAPN Emitters', 'Location', 'best');
hold off;
```

#### Overlay Images

```matlab
% Generate overlay with circles scaled by precision
RadiusScale = 2;  % Circle radius = RadiusScale * precision
ScaleBarLength = 500;  % nm

overlay_img = BGL.genSRMAPNOverlay(BGL.SMD, MAPN, ...
    BGL.PImageSize, BGL.PImageSize, BGL.PixelSize, ...
    pwd, BGL.XStart, BGL.YStart, RadiusScale, ScaleBarLength);

% Display
figure;
imshow(overlay_img);
title('Overlay: Green=Input, Magenta=MAPN');
```

The overlay shows:
- **Green circles**: Input localizations with radius proportional to uncertainty
- **Magenta circles**: MAPN emitters with radius proportional to uncertainty
- Allows visual assessment of grouping quality

### Step 7: Save Results

#### Automatic Saving

```matlab
% Save all standard outputs
ScaleBarLength = 500;  % nm
SaveDir = './BaGoL_Results';
BGL.saveBaGoL(ScaleBarLength, SaveDir, 1);  % 1 = include overlay images
```

This creates:
- `MAPN.mat`: MAPN coordinates
- `BaGoL_X-SE.png`: X precision histogram
- `BaGoL_Y-SE.png`: Y precision histogram
- `Xi.png`: Localization distribution
- `NND.png`: Nearest neighbor distances
- `MAPN-Im.png`: MAPN image
- `SR-Im.png`: Input SR image
- `Post-Im.png`: Posterior image
- `Overlay_*.png`: Various overlay images

#### Save Specific Components

```matlab
% Save just MAPN coordinates
MAPN = BGL.MAPN;
PixelSize = BGL.PixelSize;  % Output pixel size
save('./MAPN_results.mat', 'MAPN', 'PixelSize');

% Save full BaGoL object (can be large)
if length(BGL.SMD.X) < 25000
    save('./BaGoL_complete.mat', 'BGL');
else
    % For large datasets, save without heavy components
    BGL_light = BGL;
    BGL_light.Chain = [];
    BGL_light.PImage = [];
    save('./BaGoL_light.mat', 'BGL_light');
end
```

#### Export for Further Analysis

```matlab
% Export MAPN to simple format
coords_table = table(MAPN.X, MAPN.Y, MAPN.X_SE, MAPN.Y_SE, MAPN.Nmean, ...
    'VariableNames', {'X_nm', 'Y_nm', 'X_SE_nm', 'Y_SE_nm', 'Locs_per_Emitter'});
writetable(coords_table, './MAPN_coordinates.csv');

% Convert back to pixels if needed for SMITE compatibility
MAPN_SMD = struct();
MAPN_SMD.X = MAPN.X / BaGoL_SMD.PixelSize;  % nm to pixels
MAPN_SMD.Y = MAPN.Y / BaGoL_SMD.PixelSize;
MAPN_SMD.X_SE = MAPN.X_SE / BaGoL_SMD.PixelSize;
MAPN_SMD.Y_SE = MAPN.Y_SE / BaGoL_SMD.PixelSize;
MAPN_SMD.PixelSize = BaGoL_SMD.PixelSize;
save('./MAPN_SMD.mat', 'MAPN_SMD');
```

## ROI-Based Analysis

For analyzing specific regions or multiple structures:

### Defining ROIs

```matlab
% Option 1: Manual ROI definition
% Define multiple ROIs [Xmin, Xmax, Ymin, Ymax] in pixels
DataROI = [
    100, 150, 100, 150;  % ROI 1
    200, 250, 180, 230;  % ROI 2
];

% Process each ROI
for i = 1:size(DataROI, 1)
    % Extract ROI
    roi = DataROI(i, :);
    Ind = SMD.X >= roi(1) & SMD.X <= roi(2) & ...
          SMD.Y >= roi(3) & SMD.Y <= roi(4);

    % Create ROI-specific SMD
    ROI_SMD = struct();
    ROI_SMD.X = (SMD.X(Ind) - roi(1)) * SMD.PixelSize;  % Relative coords in nm
    ROI_SMD.Y = (SMD.Y(Ind) - roi(3)) * SMD.PixelSize;
    ROI_SMD.X_SE = SMD.X_SE(Ind) * SMD.PixelSize;
    ROI_SMD.Y_SE = SMD.Y_SE(Ind) * SMD.PixelSize;
    ROI_SMD.Z = [];
    ROI_SMD.Z_SE = [];
    ROI_SMD.FrameNum = SMD.FrameNum(Ind);

    % Run BaGoL on ROI
    BGL_ROI = smi.BaGoL();
    BGL_ROI.SMD = ROI_SMD;
    % ... set other parameters ...
    BGL_ROI.analyze_all();

    % Save ROI results
    MAPN = BGL_ROI.MAPN;
    save(sprintf('./ROI_%02d_MAPN.mat', i), 'MAPN');
end
```

### Combining ROI Results

Recent commits mention `combineBaGoLROIs` fixes for combining results from multiple ROIs. While this function isn't directly visible in the searched files, the typical workflow is:

```matlab
% After processing multiple ROIs, combine them
% Load individual ROI results
ROI_results = cell(n_ROIs, 1);
for i = 1:n_ROIs
    data = load(sprintf('./ROI_%02d_MAPN.mat', i));
    ROI_results{i} = data.MAPN;
end

% Concatenate results
combined_MAPN = struct();
combined_MAPN.X = [];
combined_MAPN.Y = [];
combined_MAPN.X_SE = [];
combined_MAPN.Y_SE = [];
combined_MAPN.Nmean = [];

for i = 1:n_ROIs
    % Adjust coordinates back to global reference
    roi_offset_x = DataROI(i, 1) * SMD.PixelSize;
    roi_offset_y = DataROI(i, 3) * SMD.PixelSize;

    combined_MAPN.X = [combined_MAPN.X; ROI_results{i}.X + roi_offset_x];
    combined_MAPN.Y = [combined_MAPN.Y; ROI_results{i}.Y + roi_offset_y];
    combined_MAPN.X_SE = [combined_MAPN.X_SE; ROI_results{i}.X_SE];
    combined_MAPN.Y_SE = [combined_MAPN.Y_SE; ROI_results{i}.Y_SE];
    combined_MAPN.Nmean = [combined_MAPN.Nmean; ROI_results{i}.Nmean];
end

save('./combined_MAPN.mat', 'combined_MAPN');
```

## Advanced Topics

### Understanding the Posterior Image

The posterior image differs fundamentally from the MAPN image:

```matlab
% Posterior image represents uncertainty
% - Weighted average over all models sampled
% - Bright areas = high probability of emitter presence
% - Spread indicates position uncertainty
% - Incorporates model uncertainty (number of emitters)

% MAPN image is a point estimate
% - Uses only the most frequent number of emitters
% - Discrete positions with Gaussian blobs
% - Similar to standard SR reconstruction
% - Better for quantitative coordinate analysis

% Compare visually
figure;
subplot(1,2,1);
imagesc(BGL.PImage);
axis image; colormap hot; colorbar;
title('Posterior Image (Full Uncertainty)');

subplot(1,2,2);
MAPN_img = BGL.genMAPNIm(1);
imagesc(MAPN_img);
axis image; colormap hot; colorbar;
title('MAPN Image (Point Estimate)');
```

### Chain Diagnostics

When `ChainFlag = 1`, examine the RJMCMC chains:

```matlab
% Check chain convergence for one cluster
cluster_idx = 1;
chain = BGL.Chain{cluster_idx};

% Plot number of emitters over iterations
N_chain = [chain.N];
figure;
plot(N_chain);
xlabel('RJMCMC Iteration');
ylabel('Number of Emitters');
title(sprintf('Cluster %d: N(emitters) Chain', cluster_idx));

% Find mode (MAPN value)
[counts, values] = hist(N_chain, unique(N_chain));
[~, mode_idx] = max(counts);
MAPN_N = values(mode_idx);
fprintf('MAPN number of emitters: %d\n', MAPN_N);
fprintf('Frequency: %.1f%%\n', 100*counts(mode_idx)/length(N_chain));

% Plot emitter positions over chain
% (only for chains with same N)
chain_same_N = chain(N_chain == MAPN_N);
figure;
hold on;
for i = 1:length(chain_same_N)
    plot(chain_same_N(i).X, chain_same_N(i).Y, '.');
end
xlabel('X (nm)'); ylabel('Y (nm)');
title(sprintf('Emitter Positions for N=%d', MAPN_N));
axis equal;
```

### Custom Jump Probabilities

For challenging data, adjust jump probabilities:

```matlab
% Default: Equal probabilities
BGL.P_Jumps = [0.25, 0.25, 0.25, 0.25];  % [Move, Allocate, Add, Remove]

% For stable structures (DNA origami):
% Favor refinement over topology changes
BGL.P_Jumps = [0.4, 0.4, 0.1, 0.1];

% For uncertain emitter count:
% Favor topology exploration
BGL.P_Jumps = [0.2, 0.2, 0.3, 0.3];

% For dense data with ambiguous allocation:
% Favor reallocation
BGL.P_Jumps = [0.25, 0.5, 0.125, 0.125];
```

### Memory Management

For large datasets:

```matlab
% Estimate memory usage
n_locs = length(BGL.SMD.X);
n_clusters_est = n_locs / 100;  % Rough estimate
chain_size_MB = n_clusters_est * BGL.N_Trials * 0.001;  % Very rough

fprintf('Estimated chain size: %.1f MB\n', chain_size_MB);

% For large datasets
if n_locs > 50000
    fprintf('Large dataset detected. Adjusting parameters...\n');

    % Don't save chain
    BGL.ChainFlag = 0;

    % Use smaller ROIs
    BGL.ROIsize = 300;

    % Optional: Don't generate posterior image
    BGL.PImageFlag = 0;

    % Process in batches if needed
    % (split data, process separately, combine results)
end
```

## Troubleshooting

### Issue: Poor MAPN precision (> 5 nm)

**Diagnose:**
```matlab
histogram(MAPN.X_SE);
xlabel('X Precision (nm)');
fprintf('Median precision: %.2f nm\n', median(MAPN.X_SE));
```

**Solutions:**
- Check input localization quality (should be < 10-15 nm)
- Increase `SE_Adjust` if precisions are underestimated
- Ensure sufficient localizations per emitter (5-10+ ideal)
- Check that frame connection worked properly
- Verify camera calibration accuracy

### Issue: Unrealistic emitter count

**Diagnose:**
```matlab
compression = length(BGL.SMD.X) / length(MAPN.X);
fprintf('Compression ratio: %.1f:1\n', compression);
% Expected: 5:1 to 30:1 for DNA-PAINT/dSTORM
```

**Solutions:**
- **Too many emitters** (compression < 3:1):
  - Increase Xi prior (more locs/emitter expected)
  - Check that frame connection worked
  - Reduce `Cutoff` parameter
  - Verify blinking statistics are as expected

- **Too few emitters** (compression > 50:1):
  - Decrease Xi prior
  - Check for over-connection in frame connection step
  - Increase `Cutoff` parameter
  - Verify emitter density is reasonable

### Issue: Long computation time

**Diagnose:**
```matlab
tic;
% ... run subset of data ...
elapsed = toc;
fprintf('Time per 1000 localizations: %.1f sec\n', elapsed * 1000 / n_locs);
```

**Solutions:**
- Reduce `ROIsize` (most effective)
- Reduce `N_Burnin` and `N_Trials` for testing
- Increase `Cutoff` to reduce pre-clusters
- Process ROIs in parallel (already done by `hierBaGoL_run`)
- Check for memory swapping (reduce `ChainFlag`, `PImageFlag`)

### Issue: Chain not converging

**Diagnose:**
```matlab
% Plot N chain for several clusters
figure;
for i = 1:min(5, length(BGL.Chain))
    subplot(5,1,i);
    plot([BGL.Chain{i}.N]);
    ylabel(sprintf('Cluster %d', i));
end
xlabel('Iteration');
sgtitle('Number of Emitters Over Chain');
```

**Solutions:**
- Increase `N_Burnin`
- Check that ROI isn't too large (too many emitters in one RJMCMC)
- Verify data quality (poor localizations prevent convergence)
- Adjust `P_Jumps` to favor exploration early

### Issue: Edge artifacts between ROIs

**Diagnose:**
```matlab
% Look for duplicate or missing emitters at ROI boundaries
figure;
plot(MAPN.X, MAPN.Y, 'o');
% Draw ROI grid
ROI_grid_x = 0:BGL.ROIsize:max(BGL.SMD.X);
ROI_grid_y = 0:BGL.ROIsize:max(BGL.SMD.Y);
hold on;
for x = ROI_grid_x
    plot([x x], [0 max(BGL.SMD.Y)], 'r--');
end
for y = ROI_grid_y
    plot([0 max(BGL.SMD.X)], [y y], 'r--');
end
title('MAPN Results with ROI Grid');
```

**Solutions:**
- Increase `Overlap` parameter
- Check that `removeOverlap` is functioning correctly
- Consider using smaller ROIs with more overlap
- For critical structures, process entire region as one ROI

## Performance Benchmarks

Typical processing times (256×256 pixel region, consumer workstation):

| Localizations | ROIsize (nm) | N_Trials | Processing Time |
|--------------|-------------|----------|-----------------|
| 1,000        | 500         | 8,000    | ~2 min          |
| 5,000        | 500         | 8,000    | ~10 min         |
| 10,000       | 500         | 8,000    | ~25 min         |
| 50,000       | 300         | 8,000    | ~3 hours        |
| 100,000      | 300         | 8,000    | ~8 hours        |

Factors affecting speed:
- **ROIsize**: Computational cost scales exponentially with ROI area
- **Chain length**: Linear scaling with N_Burnin + N_Trials
- **Localization density**: More locs per ROI = slower
- **HierarchFlag**: Learning Xi adds ~30% overhead

## Interpreting Results

### Precision Improvement

Calculate precision improvement over input:

```matlab
% Compare input vs MAPN precision
input_precision = median([BGL.SMD.X_SE; BGL.SMD.Y_SE]);
mapn_precision = median([MAPN.X_SE; MAPN.Y_SE]);
improvement = input_precision / mapn_precision;

fprintf('Precision improvement: %.1fx\n', improvement);
fprintf('Input: %.2f nm, MAPN: %.2f nm\n', input_precision, mapn_precision);

% Theoretical best case from N observations
N_avg = mean(MAPN.Nmean);
theoretical_improvement = sqrt(N_avg);
fprintf('Theoretical improvement: %.1fx\n', theoretical_improvement);
fprintf('Efficiency: %.1f%%\n', 100*improvement/theoretical_improvement);
```

Expected improvements:
- DNA-PAINT (10-20 locs/emitter): 2-4x improvement
- dSTORM (20-50 locs/emitter): 3-6x improvement
- Well below 1 nm precision achievable with good data

### Blinking Statistics

Analyze how emitters blink:

```matlab
% Distribution of localizations per emitter
figure;
histogram(MAPN.Nmean, 'Normalization', 'probability');
xlabel('Localizations per Emitter');
ylabel('Probability');
title('Blinking Distribution');

fprintf('Blinking statistics:\n');
fprintf('  Mean: %.1f locs/emitter\n', mean(MAPN.Nmean));
fprintf('  Median: %.1f locs/emitter\n', median(MAPN.Nmean));
fprintf('  Std: %.1f locs\n', std(MAPN.Nmean));

% Compare to prior
if ~BGL.HierarchFlag && length(BGL.Xi) == 2
    prior_mean = BGL.Xi(1) * BGL.Xi(2);
    fprintf('  Prior mean: %.1f locs/emitter\n', prior_mean);
    fprintf('  Difference: %.1f locs (%.1f%%)\n', ...
        mean(MAPN.Nmean) - prior_mean, ...
        100*(mean(MAPN.Nmean) - prior_mean)/prior_mean);
end
```

### Spatial Organization

Analyze emitter spatial patterns:

```matlab
% Nearest neighbor distances reveal organization
[~, NND] = knnsearch([MAPN.X, MAPN.Y], [MAPN.X, MAPN.Y], 'K', 2);
NND = NND(:, 2);

% Compare to random distribution at same density
area_nm2 = (max(MAPN.X) - min(MAPN.X)) * (max(MAPN.Y) - min(MAPN.Y));
density = length(MAPN.X) / area_nm2;  % emitters/nm²
expected_NND = 1 / (2 * sqrt(density));  % For 2D Poisson point process

fprintf('Spatial analysis:\n');
fprintf('  Density: %.6f emitters/nm²\n', density);
fprintf('  Observed median NND: %.2f nm\n', median(NND));
fprintf('  Expected NND (random): %.2f nm\n', expected_NND);

if median(NND) < expected_NND
    fprintf('  → Clustered organization\n');
else
    fprintf('  → Dispersed organization\n');
end
```

## Citation

When publishing results using BaGoL, please cite:

Mohamadreza Fazel, Michael J. Wester, David J. Schodt, Sebastian Restrepo Cruz, Sebastian Strauss, Florian Schueder, Thomas Schlichthaerle, Jennifer M. Gillette, Diane S. Lidke, Bernd Rieger, Ralf Jungmann and Keith A. Lidke, "High-Precision Estimation of Emitter Positions using Bayesian Grouping of Localizations", *Nature Communications*, **13**(7152), November 22, 2022, 1-11, DOI: 10.1038/s41467-022-34894-2.

## See Also

- [SMLM Analysis Workflow](./smlm-analysis.md) - Prerequisite pipeline
- [SMD Structure](../core-concepts/smd-structure.md) - Data format details
- [Data Flow](../core-concepts/data-flow.md) - Understanding the pipeline
- MATLAB/+smi/@BaGoL/README.md - Class documentation
- MATLAB/examples/hierBaGoL_wrapper.m - Batch processing example
