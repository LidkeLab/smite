---
title: "API Reference: +smi_cluster Namespace"
category: "api-reference"
level: "intermediate"
tags: ["api", "smi_cluster", "clustering", "dbscan", "voronoi", "hset", "pair-correlation", "spatial-analysis"]
prerequisites: ["../core-concepts/architecture.md", "../core-concepts/smd-structure.md", "smi-core.md"]
related: ["../workflows/bagol-clustering.md", "../examples/clustering-analysis.md", "../how-to/visualize-results.md"]
summary: "Complete API reference for the +smi_cluster namespace covering spatial clustering algorithms, statistical analysis, and pair correlation methods"
estimated_time: "30 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# API Reference: +smi_cluster Namespace

## Purpose

The +smi_cluster namespace provides comprehensive tools for spatial clustering analysis of single molecule localization data. This reference documents clustering algorithms (DBSCAN, Voronoi, H-SET, Hierarchical), statistical analysis methods, pair correlation functions, and cluster characterization tools. These methods are essential for quantifying protein organization, membrane domains, and spatial patterns in SMLM data.

## Prerequisites

- Understanding of [smite architecture](../core-concepts/architecture.md)
- Familiarity with [SMD structure](../core-concepts/smd-structure.md)
- Knowledge of [+smi_core namespace](smi-core.md)
- Basic spatial statistics concepts
- MATLAB class usage

## Overview

The +smi_cluster namespace provides:

**Clustering Algorithms:**
- `Clustering` - Main class implementing DBSCAN, Voronoi, Hierarchical, H-SET
- Density-based, tessellation-based, and precision-weighted methods
- Comprehensive cluster statistics and characterization

**Statistical Analysis:**
- `StatisticsClustering` - Hopkins statistic, Ripley's K function
- Tests for spatial randomness vs clustering
- Bivariate statistics for multi-channel analysis

**Pair Correlation:**
- `PairCorrelation` - Auto- and cross-correlation functions
- Quantifies spatial relationships and length scales
- Identifies clustering and co-localization patterns

**Analysis Interfaces:**
- `ClusterInterface` - ROI-based analysis workflows
- `PairAnalysis` - Multi-channel correlation analysis
- Integration with BaGoL clustering results

## Class Reference

### Main Clustering Class

#### Clustering

**Purpose:** Implements multiple spatial clustering algorithms and computes comprehensive cluster statistics.

**Key Concept:** Different clustering algorithms have different strengths. DBSCAN excels at density-based detection, Voronoi handles irregular shapes, H-SET uses localization precision, and Hierarchical provides multi-scale analysis. Choose the algorithm that matches your biological question.

**Class Definition:**
```matlab
classdef Clustering < handle
```

**Properties:**
```matlab
% Visualization and output
Font_props      % Font properties for plots
Fig_ext         % Figure format ('png', 'pdf', etc.)
ResultsDir      % Output directory for results
ShrinkFactor    % Boundary compactness (0-1)
Verbose         % Verbosity level (0-3)
Xlim, Ylim      % Optional axis limits

% Statistics computation
DoSigmaActual   % Compute cluster spread (memory intensive)

% Voronoi-specific
Alpha           % Density ratio threshold (default: 2)
Valgorithm      % Algorithm mode (1, 2, or 3)
Plotting        % Generate Voronoi plots
PtIDs           % Label points in plots

% H-SET specific
Method          % Collapse method ('hierarchal_singlelabel')
LoS             % Level of significance (default: 0.01)
PixelSize       % nm per pixel (default: 100)
PlotFigures     % Generate cluster plots
Sigma_Reg       % Registration error [x, y] (nm)
Timing          % Display timing information
```

**Constructor:**
```matlab
c = smi_cluster.Clustering(SMF)
```

Creates clustering object with parameters from SMF (optional).

**Main Methods:**

`cluster(algorithm, SMD, E, minPts)` - Perform clustering
```matlab
% DBSCAN clustering
[nC, C, centers, ptsI] = c.cluster('DBSCAN', SMD, E, minPts);

% Voronoi clustering
[nC, C, centers, ptsI] = c.cluster('Voronoi', SMD, [], minPts);

% Hierarchical clustering
[nC, C, centers, ptsI] = c.cluster('Hierarchical', SMD, E, minPts);

% H-SET clustering (requires SMD.X_SE, SMD.Y_SE)
[nC, C, centers, ptsI] = c.cluster('H-SET', SMD, [], minPts);
```

**Inputs:**
- `algorithm`: Clustering method (see below)
- `SMD`: Coordinates as SMD structure or N×2/N×3 array (nm)
- `E`: Epsilon/cutoff distance (nm), not used for Voronoi/H-SET
- `minPts`: Minimum points per cluster (default: 3)

**Outputs:**
- `nC`: Number of clusters found
- `C`: Cell array of indices for each cluster {nC × 1}
- `centers`: Cluster center coordinates [n_dim × nC]
- `ptsI`: Indices of unclustered (isolated) points

**Supported Algorithms:**

**DBSCAN** - Density-Based Spatial Clustering
```matlab
% Standard DBSCAN with specified epsilon
[nC, C, centers, ptsI] = c.cluster('DBSCAN', SMD, 50, 3);
% E = 50 nm: maximum distance within cluster
% minPts = 3: minimum cluster size

% DBSCAN with automatic epsilon
[nC, C, centers, ptsI] = c.cluster('DBSCAN_Daszykowski_noE', SMD, [], 3);
% Algorithm estimates optimal E from data
```

**When to use DBSCAN:**
- Compact, well-separated clusters
- Known typical cluster density
- Noise tolerance required
- Arbitrary cluster shapes (not just circles)

**Key parameters:**
- E: Too small = everything isolated; too large = everything clusters
- minPts: Noise filtering; typically 3-5 for 2D data

**Voronoi** - Tessellation-Based Clustering
```matlab
% Configure Voronoi
c.Alpha = 2.0;        % Density threshold
c.Valgorithm = 2;     % Consider neighbors
c.Plotting = false;   % No intermediate plots

% Perform clustering
[nC, C, centers, ptsI] = c.cluster('Voronoi', SMD, [], 3);
```

**When to use Voronoi:**
- Irregular or complex cluster shapes
- Variable cluster densities
- Unknown cluster sizes
- Membrane domains or diffuse structures

**Key parameters:**
- Alpha: Higher = more conservative (only very dense regions)
- Valgorithm: 1 (cell only), 2 (cell + neighbors, recommended), 3 (median)

**Hierarchical** - Linkage-Based Clustering
```matlab
% Hierarchical clustering
[nC, C, centers, ptsI] = c.cluster('Hierarchical', SMD, 50, 3);
% Uses MATLAB's linkage algorithm
```

**When to use Hierarchical:**
- Multi-scale cluster structure
- Known hierarchical organization
- Need dendrogram visualization
- Exploratory analysis

**H-SET** - Precision-Weighted Clustering
```matlab
% Requires localization uncertainties
% SMD.X_SE and SMD.Y_SE must be populated

c.LoS = 0.01;              % Level of significance
c.Sigma_Reg = [10, 10];    % Registration error (nm)
c.Method = 'hierarchal_singlelabel';

% Perform H-SET clustering
[nC, C, centers, ptsI] = c.cluster('H-SET', SMD, [], 3);
```

**When to use H-SET:**
- High-precision data with uncertainty estimates
- Need to resolve nearby emitters
- Precision varies across localizations
- Statistical rigor required

**Key parameters:**
- LoS: Significance level for combining points (typically 0.01)
- Sigma_Reg: Account for systematic errors (typically 10 nm)

`clusterStats(SMD, C, centers)` - Compute comprehensive statistics
```matlab
results = c.clusterStats(SMD, C, centers);
```

**Returns extensive cluster characterization:**
```matlab
% Basic counts
results.nC              % Number of clusters
results.n_points        % Total points
results.n_clustered     % Points in clusters
results.n_isolated      % Isolated points
results.n_pts           % Points per cluster [1 × nC]

% Cluster classification
results.numclust        % [singlets, doublets, multiplets]
results.singlet_fraction % Fraction as singlets

% Size metrics
results.clust_width     % Max vertex-to-vertex [1 × nC]
results.areas           % Cluster areas (nm²) [1 × nC]
results.equiv_radii     % Equivalent radius (nm) [1 × nC]
results.sigma_actual    % Intra-cluster spread [1 × nC]

% Density metrics
results.n_pts_per_area  % Density (points/nm²)

% Shape metrics (2D only)
results.compactness     % 4π×Area/Perimeter² [1 for circle]
results.circularity     % Shape regularity
results.convexity       % Perimeter ratio [1 for convex]
results.solidity        % Area ratio [1 for solid]

% Boundaries
results.indices_hull    % Hull boundary indices per cluster
results.indConvex       % Convex hull indices
results.indCompact      % Compact boundary indices
results.perimeters      % Perimeter lengths (nm)

% Inter-cluster distances
results.min_c2c_dists   % Center-to-center [1 × nC]
results.min_e2e_dists   % Edge-to-edge [1 × nC]
results.min_c2c_dist    % Overall minimum c2c
results.min_e2e_dist    % Overall minimum e2e

% Intra-cluster distances
results.nn_within_clust % NN distances within clusters
```

`plotClusters(SMD, C, centers, ptsI, title_str)` - Visualize clusters
```matlab
% Plot with color-coded clusters
h = c.plotClusters(SMD, C, centers, ptsI, 'My Clusters');

% Returns figure handle
set(h, 'Visible', 'on');  % Make visible
saveas(h, 'clusters.png');
```

**Visualization features:**
- Each cluster gets unique color
- Isolated points shown in gray
- Cluster centers marked
- Boundaries drawn (convex or compact)
- Scalable to thousands of clusters

`plotClusters3(SMD, C, centers, ptsI, title_str)` - 3D visualization
```matlab
% For 3D data
h = c.plotClusters3(SMD, C, centers, ptsI, '3D Clusters');
```

**Static Methods:**

`dbscan_Daszykowski(XY, k, Eps)` - DBSCAN implementation
```matlab
[class, type, Eps] = smi_cluster.Clustering.dbscan_Daszykowski(XY, k, Eps);
% Direct DBSCAN call
% class: cluster assignment per point
% type: point type (core, border, noise)
% Eps: used/estimated epsilon
```

`hierarchal(XY, E, minPts)` - Hierarchical clustering
```matlab
[C, ptsI] = smi_cluster.Clustering.hierarchal(XY, E, minPts);
```

`voronoi_Levet(XY, alpha, epsilon, minPts, algorithm)` - Voronoi clustering
```matlab
[area, rho, nC, C] = smi_cluster.Clustering.voronoi_Levet(...
    XY, alpha, epsilon, minPts, algorithm);
```

`nn_distances(xy)` - Nearest neighbor distances
```matlab
min_dists = smi_cluster.Clustering.nn_distances(xy);
% Returns NN distance for each point
```

`edge2edge(hull1, hull2)` - Edge-to-edge distance
```matlab
dist = smi_cluster.Clustering.edge2edge(hull1, hull2);
% Minimum distance between two cluster hulls
```

`vertex2vertex(hull)` - Maximum hull span
```matlab
max_dist = smi_cluster.Clustering.vertex2vertex(hull);
% Maximum distance between vertices
```

**Usage Examples:**

Basic clustering workflow:
```matlab
% Setup
SMF = smi_core.SingleMoleculeFitting();
c = smi_cluster.Clustering(SMF);
c.PixelSize = 100;  % nm per pixel
c.Verbose = 1;

% Load data
load('Results.mat', 'SMD');

% DBSCAN clustering
E = 50;      % nm
minPts = 3;
[nC, C, centers, ptsI] = c.cluster('DBSCAN', SMD, E, minPts);

fprintf('Found %d clusters\n', nC);
fprintf('Isolated points: %d\n', length(ptsI));

% Compute statistics
results = c.clusterStats(SMD, C, centers);

fprintf('Clustered: %d/%d (%.1f%%)\n', ...
    results.n_clustered, results.n_points, ...
    100 * results.n_clustered / results.n_points);

% Visualize
c.plotClusters(SMD, C, centers, ptsI, 'My Clusters');
```

Comparing algorithms:
```matlab
% Same data, multiple algorithms
algorithms = {'DBSCAN', 'Voronoi', 'Hierarchical'};
E = 50;
minPts = 3;

figure;
for i = 1:length(algorithms)
    alg = algorithms{i};

    % Configure algorithm-specific parameters
    if strcmp(alg, 'Voronoi')
        c.Alpha = 2.0;
        [nC, C, centers, ptsI] = c.cluster(alg, SMD, [], minPts);
    else
        [nC, C, centers, ptsI] = c.cluster(alg, SMD, E, minPts);
    end

    % Plot
    subplot(1, 3, i);
    h = c.plotClusters(SMD, C, centers, ptsI, alg);
    set(h, 'Visible', 'on');
    title(sprintf('%s: %d clusters', alg, nC));
end
```

H-SET precision-weighted clustering:
```matlab
% Requires localization uncertainties
% SMD.X_SE, SMD.Y_SE must be present

c = smi_cluster.Clustering(SMF);
c.LoS = 0.01;              % Significance level
c.Sigma_Reg = [10, 10];    % Registration error
c.Timing = true;

% Perform H-SET clustering
[nC, C, centers, ptsI] = c.cluster('H-SET', SMD, [], 3);

% Compute statistics
results = c.clusterStats(SMD, C, centers);

% Compare cluster sizes to DBSCAN
fprintf('H-SET found %d clusters\n', nC);
fprintf('Mean cluster size: %.1f points\n', mean(results.n_pts));
```

Parameter optimization:
```matlab
% Scan epsilon values to find optimal
E_range = 20:10:100;
nC_array = zeros(size(E_range));
clustered_frac = zeros(size(E_range));

for i = 1:length(E_range)
    [nC, C, ~, ~] = c.cluster('DBSCAN', SMD, E_range(i), 3);
    nC_array(i) = nC;

    results = c.clusterStats(SMD, C, []);
    clustered_frac(i) = results.n_clustered / results.n_points;
end

% Plot optimization curves
figure;
subplot(1,2,1);
plot(E_range, nC_array, 'o-', 'LineWidth', 2);
xlabel('Epsilon (nm)');
ylabel('Number of Clusters');
title('Cluster Count vs Epsilon');

subplot(1,2,2);
plot(E_range, 100*clustered_frac, 'o-', 'LineWidth', 2);
xlabel('Epsilon (nm)');
ylabel('Clustered Points (%)');
title('Clustering Fraction vs Epsilon');
```

**See Also:**
- [Clustering Analysis Example](../examples/clustering-analysis.md)
- [BaGoL Workflow](../workflows/bagol-clustering.md)

**Citations:**

DBSCAN:
- M. Ester, H.-P. Kriegel, J. Sander, X. Xu, "A Density-Based Algorithm for Discovering Clusters in Large Spatial Databases with Noise", KDD-96, 1996.
- M. Daszykowski, B. Walczak, D. L. Massart, "Looking for Natural Patterns in Data. Part 1: Density Based Approach", Chemom. Intell. Lab. Syst. 56 (2001) 83-92.

Voronoi:
- F. Levet, E. Hosy, A. Kechkar, C. Butler, A. Beghin, D. Choquet, J.-B. Sibarita, "SR-Tesseler: a method to segment and quantify localization-based super-resolution microscopy data", Nature Methods 12 (2015) 1065-1071.

H-SET:
- J. Lin, M. J. Wester, M. S. Graus, K. A. Lidke, A. K. Neumann, "Nanoscopic cell wall architecture of an immunogenic ligand in Candida albicans during antifungal drug treatment", Molecular Biology of the Cell 27(6) (2016) 1002-1014.

---

### Statistical Analysis

#### StatisticsClustering

**Purpose:** Statistical tests for spatial randomness, clustering, and regularity using established spatial statistics methods.

**Key Concept:** Quantitative tests distinguish true spatial organization from random patterns. Hopkins statistic tests for clustering vs randomness, while Ripley's K function characterizes clustering at multiple length scales.

**Class Definition:**
```matlab
classdef StatisticsClustering < handle
```

**Properties:**
```matlab
% General
ROI             % ROI bounds [x_min, x_max, y_min, y_max] (nm)
Font_props      % Plot font properties
Line_props      % Plot line properties
Fig_ext         % Figure format
BaseName        % Descriptive name for files
ResultsDir      % Output directory
Xlim, Ylim      % Optional axis limits
PixelSize       % Camera pixel size (nm)
Verbose         % Verbosity level

% Plotting control
PlotDo          % Which plots to generate
LinLog          % Plot scaling ('plot', 'semilogx', etc.)
LegendTitle     % Optional legend title
ShowMM          % Show mean/median markers
CSV             % Generate CSV output

% Test parameters
HopTestPts      % Number of test points (default: 10)
HopTests        % Number of Hopkins tests (default: 1000)
Rate            % Sampling rate (default: 20)
Dendro_cutoff   % Dendrogram cutoff (nm, default: 50)
Ripley_cutoff   % Ripley distance cutoff (nm, default: 200)
Confidence      % Confidence interval (default: 2.576 = 99%)
Nsims           % Number of simulations (default: 20)
```

**Constructor:**
```matlab
SC = smi_cluster.StatisticsClustering(SMF)
```

**Main Methods:**

`hopkins(labels, SMD)` - Hopkins statistic for spatial randomness
```matlab
results = SC.hopkins({'MyData'}, SMD);
```

**Hopkins Statistic Interpretation:**
- H < 0.5: Data is clustered
- H ≈ 0.5: Data is randomly distributed
- H > 0.5: Data is regularly spaced/dispersed

**Returns:**
```matlab
results.H{1}         % Hopkins statistic value
results.HopStats{1}  % Individual test statistics
results.A{1}         % Area used
results.n{1}         % Number of points
```

**Algorithm:** Compares distances from random points to data vs distances from random points to random points. Clustering makes data-distances shorter.

`hopkins_ROIcombined(labels, n_ROIs, RoI)` - Hopkins across multiple ROIs
```matlab
results = SC.hopkins_ROIcombined({'Condition'}, n_ROIs, RoI);
% Combines Hopkins tests across ROIs for better statistics
```

`ripley(labels, SMD, desc)` - Ripley's K function
```matlab
results = SC.ripley({'MyData'}, SMD, 'Analysis');
```

**Ripley's K Interpretation:**
- K(r) > πr²: Clustering at distance r
- K(r) = πr²: Random distribution at distance r
- K(r) < πr²: Dispersion at distance r

**Returns:**
```matlab
results.K{1}         % K function vs distance
results.L{1}         % L(r) = sqrt(K(r)/π) - r (linearized)
results.H{1}         % H function
results.r{1}         % Distance values (nm)
results.edges{1}     % Confidence interval edges
results.A{1}         % Area
```

**Algorithm:** For each distance r, counts points within r of each point. Normalized by area and density. L(r) > 0 indicates clustering.

`ripley_ROIcombined(labels, n_ROIs, RoI, desc)` - Combined Ripley's K
```matlab
results = SC.ripley_ROIcombined({'Data'}, n_ROIs, RoI, 'Combined');
% Averages Ripley's K across multiple ROIs
```

`bivariateRipley(labels, SMD1, SMD2, desc)` - Cross-channel Ripley
```matlab
results = SC.bivariateRipley({'Ch1_vs_Ch2'}, SMD1, SMD2, 'Cross');
```

**Bivariate Ripley Interpretation:**
- Tests spatial correlation between two patterns
- L₁₂(r) > 0: Channel 1 and 2 are co-localized at distance r
- L₁₂(r) ≈ 0: Independent distributions
- L₁₂(r) < 0: Anti-correlation/exclusion

**Returns:**
```matlab
results.K_12{1}      % Cross K function
results.L_12{1}      % Cross L function
results.H_12{1}      % Cross H function
results.r{1}         % Distance values
results.edges{1}     % Confidence intervals
```

`bivariateRipley_ROIcombined(...)` - Combined bivariate Ripley
```matlab
results = SC.bivariateRipley_ROIcombined(labels, n_ROIs, RoI, desc);
```

`pairwiseDist(labels, SMD, desc)` - Pairwise distance distribution
```matlab
results = SC.pairwiseDist({'MyData'}, SMD, 'Analysis');
```

**Returns:**
```matlab
results.all_dists{1}   % All pairwise distances
results.nn_dists{1}    % Nearest neighbor distances
results.mean_nn{1}     % Mean NN distance
results.median_nn{1}   % Median NN distance
```

`pairwiseMutualDist(labels, SMD1, SMD2, desc)` - Cross-channel distances
```matlab
results = SC.pairwiseMutualDist({'Ch1_to_Ch2'}, SMD1, SMD2, 'Cross');
```

**Returns distances from points in one channel to nearest point in other channel.**

`plotCombined(data, labels, title_str)` - Multi-format plotting
```matlab
SC.PlotDo = 'fnpcCsSxb';  % All plot types
SC.plotCombined(data, labels, 'My Analysis');
```

**Plot types:**
- f: Frequency histogram
- n: Normalized histogram
- p: PDF (probability density)
- c: CDF (cumulative distribution)
- C: CDF with log scale options
- s: Scatter plot (for 2D data)
- S: PlotSpread (categorical data)
- x: Box plot
- b: Bar plot

**Static Methods:**

`histogram(A, ab, n)` - Enhanced histogram
```matlab
[X, V] = smi_cluster.StatisticsClustering.histogram(data, [min, max], nbins);
```

`hopkinstat(P, A, B, mm)` - Hopkins statistic computation
```matlab
H = smi_cluster.StatisticsClustering.hopkinstat(P, A, B, mm);
% P: number of test points
% A, B: ROI bounds
% mm: data points
```

`hopkinstat3(P, A, B, C, mm)` - 3D Hopkins statistic
```matlab
H = smi_cluster.StatisticsClustering.hopkinstat3(P, A, B, C, mm);
```

**Usage Examples:**

Hopkins test for clustering:
```matlab
% Setup
SMF = smi_core.SingleMoleculeFitting();
SC = smi_cluster.StatisticsClustering(SMF);
SC.HopTests = 1000;        % Number of tests
SC.HopTestPts = 10;        % Points per test
SC.ROI = [0, 10000, 0, 10000];  % nm

% Load data
load('Results.mat', 'SMD');

% Perform Hopkins test
results = SC.hopkins({'MyProtein'}, SMD);

fprintf('Hopkins statistic: %.3f\n', results.H{1});
if results.H{1} < 0.5
    fprintf('Data is CLUSTERED (H < 0.5)\n');
elseif results.H{1} > 0.55
    fprintf('Data is DISPERSED (H > 0.55)\n');
else
    fprintf('Data is RANDOM (H ≈ 0.5)\n');
end

% Plot histogram of test statistics
figure;
histogram(results.HopStats{1}, 30);
xlabel('Hopkins Statistic');
ylabel('Frequency');
title(sprintf('Hopkins Test (mean = %.3f)', results.H{1}));
xline(0.5, 'r--', 'Random', 'LineWidth', 2);
```

Ripley's K analysis:
```matlab
SC = smi_cluster.StatisticsClustering(SMF);
SC.Ripley_cutoff = 500;    % Analyze up to 500 nm
SC.Rate = 10;              % Sampling rate (nm)
SC.Nsims = 20;             % Simulations for confidence
SC.BaseName = 'RipleyAnalysis';

% Perform Ripley analysis
results = SC.ripley({'MyData'}, SMD, 'Analysis');

% Plot L(r) function
figure;
plot(results.r{1}, results.L{1}, 'b-', 'LineWidth', 2);
hold on;
plot(results.r{1}, results.edges{1}(1,:), 'r--', 'LineWidth', 1);
plot(results.r{1}, results.edges{1}(2,:), 'r--', 'LineWidth', 1);
yline(0, 'k--', 'Random');
xlabel('Distance r (nm)');
ylabel('L(r) - r');
title('Ripley L Function');
legend('Data', 'Confidence bounds', '', 'Random', 'Location', 'best');
grid on;

% Identify clustering length scales
clustered = results.L{1} > results.edges{1}(2,:);
if any(clustered)
    cluster_scales = results.r{1}(clustered);
    fprintf('Clustering detected at: %.0f - %.0f nm\n', ...
        min(cluster_scales), max(cluster_scales));
end
```

Bivariate analysis (two channels):
```matlab
% Load two-channel data
load('Results_Ch1.mat', 'SMD');
SMD1 = SMD;
load('Results_Ch2.mat', 'SMD');
SMD2 = SMD;

% Bivariate Ripley
SC = smi_cluster.StatisticsClustering(SMF);
SC.Ripley_cutoff = 300;
SC.Confidence = 2.576;  % 99% confidence

results = SC.bivariateRipley({'Ch1_vs_Ch2'}, SMD1, SMD2, 'Colocalization');

% Plot cross-correlation
figure;
plot(results.r{1}, results.L_12{1}, 'b-', 'LineWidth', 2);
hold on;
plot(results.r{1}, results.edges{1}(1,:), 'r--', 'LineWidth', 1);
plot(results.r{1}, results.edges{1}(2,:), 'r--', 'LineWidth', 1);
yline(0, 'k--', 'Independent');
xlabel('Distance r (nm)');
ylabel('L_{12}(r)');
title('Bivariate Ripley: Channel Cross-Correlation');
legend('Cross L', 'Confidence', '', 'Independent', 'Location', 'best');

% Interpret
if any(results.L_12{1} > results.edges{1}(2,:))
    fprintf('Significant co-localization detected\n');
else
    fprintf('Channels are independent\n');
end
```

Nearest neighbor analysis:
```matlab
SC = smi_cluster.StatisticsClustering(SMF);

% Compute pairwise distances
results = SC.pairwiseDist({'MyData'}, SMD, 'NN_Analysis');

% Analyze NN distances
fprintf('Mean NN distance: %.1f nm\n', results.mean_nn{1});
fprintf('Median NN distance: %.1f nm\n', results.median_nn{1});

% Plot NN distribution
figure;
histogram(results.nn_dists{1}, 50);
xlabel('Nearest Neighbor Distance (nm)');
ylabel('Count');
title('NN Distance Distribution');
xline(results.mean_nn{1}, 'r--', 'Mean', 'LineWidth', 2);
xline(results.median_nn{1}, 'g--', 'Median', 'LineWidth', 2);

% Compare to random expectation
A_ROI = (SC.ROI(2) - SC.ROI(1)) * (SC.ROI(4) - SC.ROI(3));
density = numel(SMD.X) / A_ROI;
expected_nn = 1 / (2 * sqrt(density));
xline(expected_nn, 'k--', 'Random', 'LineWidth', 2);

fprintf('Expected NN (random): %.1f nm\n', expected_nn);
fprintf('Observed/Expected ratio: %.2f\n', results.mean_nn{1} / expected_nn);
```

Multi-ROI combined analysis:
```matlab
% Prepare ROI structure (see smi_helpers.ROITools)
n_ROIs = 5;
RoI = cell(1, n_ROIs);
for i = 1:n_ROIs
    RoI{i}.ROI = [0, 5000, 0, 5000];
    RoI{i}.X = {SMD_list{i}.X * PixelSize};
    RoI{i}.Y = {SMD_list{i}.Y * PixelSize};
end

% Combined Hopkins test
results_hop = SC.hopkins_ROIcombined({'Combined'}, n_ROIs, RoI);
fprintf('Combined Hopkins: %.3f ± %.3f\n', ...
    mean(results_hop.H{1}), std(results_hop.H{1}));

% Combined Ripley
results_rip = SC.ripley_ROIcombined({'Combined'}, n_ROIs, RoI, 'Analysis');
% Plots average Ripley function with error bars
```

**See Also:**
- [Clustering Analysis Example](../examples/clustering-analysis.md)
- [Statistical Methods Documentation](../how-to/statistical-analysis.md)

**Citations:**

Hopkins Statistic:
- B. Hopkins, J. G. Skellam, "A New Method for Determining the Type of Distribution of Plant Individuals", Annals of Botany 18 (1954) 213-227.

Ripley's K Function:
- B. D. Ripley, "Modelling Spatial Patterns", Journal of the Royal Statistical Society B 39(2) (1977) 172-212.
- A. J. Baddeley, E. Rubak, R. Turner, "Spatial Point Patterns: Methodology and Applications with R", Chapman and Hall/CRC, 2015.

---

### Pair Correlation Analysis

#### PairCorrelation

**Purpose:** Computes pair correlation functions (auto- and cross-correlation) to quantify spatial relationships and co-localization between localizations.

**Key Concept:** The pair correlation function g(r) measures how density varies with distance. g(r) > 1 indicates clustering, g(r) = 1 is random, and g(r) < 1 is dispersion. This provides quantitative length scales for spatial organization.

**Class Definition:**
```matlab
classdef PairCorrelation < handle
```

**Properties:**
```matlab
% ROI and data
ROI             % ROI bounds [x_min, x_max, y_min, y_max] (nm)
BaseName        % Descriptive name for files

% Visualization
Font_props      % Font properties
Fig_ext         % Figure format (empty = .fig only, else .fig + ext)
Lines           % Plot lines vs points (default: true)
ResultsDir      % Output directory

% Analysis parameters
HistBinSize     % Internal pixel size (nm, default: 100)
PixelSize       % Camera pixel size (nm, default: 100)
Rmax_axis       % Plot limit (nm, -1 = auto)
RmaxFitFactor   % Fitting limit factor (≤1, default: 1)
Veatch_fit      % Fit model for pair_correlation_Veatch
                % 'exponential_and_gaussian' (default)
                % 'exponential_and_cosine'
                % 'exponential'

Verbose         % Verbosity level (0-2)
```

**Constructor:**
```matlab
PC = smi_cluster.PairCorrelation(SMF)
```

**Main Methods:**

`pair_correlation(Data1, Data2)` - Compute g(r)
```matlab
% Auto-correlation (single dataset)
results = PC.pair_correlation(SMD);
results = PC.pair_correlation(XY);  % XY = N×2 array (nm)

% Cross-correlation (two datasets)
results = PC.pair_correlation(SMD1, SMD2);
results = PC.pair_correlation(XY1, XY2);
```

**Inputs:**
- SMD structure(s) with X, Y fields (pixels)
- N×2 coordinate arrays (nm)

**Returns:**
```matlab
results.r            % Distance values (nm)
results.g            % g(r) values
results.dg           % Uncertainty in g(r)
results.rho          % Overall density (1/nm²)
results.mask         % Analysis mask
results.rmax         % Maximum radius used
results.G            % Unnormalized correlation
results.estimates    % Gaussian fit parameters
results.errors       % Fit uncertainties
results.model        % Fit model object
```

**Interpretation:**
- g(r) > 1: Points are clustered at distance r
- g(r) = 1: Random distribution at distance r
- g(r) < 1: Points are dispersed/excluded at distance r
- Peak position: Characteristic cluster spacing
- Peak height: Strength of clustering
- Peak width: Cluster size distribution

`pair_correlation_Veatch(Data1, Data2, corr_type, fit_model)` - Alternative method
```matlab
% Uses Sarah Veatch's original implementation
results = PC.pair_correlation_Veatch(SMD1, SMD2, 'cross', 'exponential_and_gaussian');
results = PC.pair_correlation_Veatch(SMD1, [], 'auto', 'exponential');
```

**corr_type:**
- 'auto': Auto-correlation
- 'cross': Cross-correlation

**fit_model:**
- 'exponential_and_gaussian': Cluster + background
- 'exponential_and_cosine': Periodic organization
- 'exponential': Simple decay

`pair_correlation_ROIcombined(n_channels, n_ROIs, RoI, channel)` - Combined analysis
```matlab
% Combine pair correlation across multiple ROIs
results = PC.pair_correlation_ROIcombined(2, n_ROIs, RoI);
% More robust statistics, averages g(r) across ROIs
```

**Inputs:**
- `n_channels`: 1 (auto) or 2 (cross)
- `n_ROIs`: Number of ROIs to combine
- `RoI`: ROI structure array
- `channel`: Which channel for auto-correlation (optional)

**Static Methods:**

`get_autocorr(I1, mask, rmax, flag)` - Low-level auto-correlation
```matlab
[G, r, g, dg, mask, rmax] = ...
    smi_cluster.PairCorrelation.get_autocorr(Image, mask, rmax, flag);
```

`get_crosscorr(I1, I2, mask, rmax, flag)` - Low-level cross-correlation
```matlab
[C, r, c, dc, mask, rmax] = ...
    smi_cluster.PairCorrelation.get_crosscorr(Image1, Image2, mask, rmax, flag);
```

`get_corr(n_ROIs, rmax, II1, II2)` - Combined correlation
```matlab
[C, r, c, dc, rmax] = ...
    smi_cluster.PairCorrelation.get_corr(n_ROIs, rmax, Images1, Images2);
```

`pc_GaussFit(r, g_r, rmax, rho)` - Fit Gaussian model
```matlab
[estimates, errors, model] = ...
    smi_cluster.PairCorrelation.pc_GaussFit(r, g, rmax, density);
% Fits 2D Gaussian to g(r) peak
```

**Usage Examples:**

Basic auto-correlation:
```matlab
% Setup
SMF = smi_core.SingleMoleculeFitting();
PC = smi_cluster.PairCorrelation(SMF);
PC.BaseName = 'MyProtein';
PC.PixelSize = 100;        % nm per pixel
PC.HistBinSize = 100;      % nm internal bins
PC.Rmax_axis = 500;        % Plot to 500 nm

% Set ROI
PC.ROI = [0, 10000, 0, 10000];  % nm

% Load data and convert to nm
load('Results.mat', 'SMD');
XY = [SMD.X * PC.PixelSize, SMD.Y * PC.PixelSize];

% Compute auto-correlation
results = PC.pair_correlation(XY);

% Display results
fprintf('Overall density: %.2e points/nm²\n', results.rho);
fprintf('Analysis radius: %.1f nm\n', results.rmax);

% Plot
figure;
if PC.Lines
    plot(results.r, results.g, 'b-', 'LineWidth', 2);
else
    plot(results.r, results.g, 'b.');
end
hold on;
yline(1, 'k--', 'Random', 'LineWidth', 1.5);
xlabel('Distance r (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('g(r)', 'FontSize', 12, 'FontWeight', 'bold');
title('Pair Correlation Function');
xlim([0, PC.Rmax_axis]);
legend('g(r)', 'Random', 'Location', 'best');
grid on;

% Interpret
peak_idx = find(results.g == max(results.g), 1);
if results.g(peak_idx) > 1.5
    fprintf('Strong clustering detected\n');
    fprintf('Peak at r = %.1f nm\n', results.r(peak_idx));
    fprintf('Peak height g(r) = %.2f\n', results.g(peak_idx));
end
```

Cross-correlation (two channels):
```matlab
% Load two-channel data
load('Results_Ch1.mat', 'SMD');
XY1 = [SMD.X * PC.PixelSize, SMD.Y * PC.PixelSize];
load('Results_Ch2.mat', 'SMD');
XY2 = [SMD.X * PC.PixelSize, SMD.Y * PC.PixelSize];

% Setup pair correlation
PC = smi_cluster.PairCorrelation(SMF);
PC.BaseName = 'Ch1_vs_Ch2';
PC.ROI = [0, 10000, 0, 10000];

% Cross-correlation
results_cross = PC.pair_correlation(XY1, XY2);

% Plot
figure;
plot(results_cross.r, results_cross.g, 'r-', 'LineWidth', 2);
hold on;
yline(1, 'k--', 'Independent');
xlabel('Distance r (nm)');
ylabel('g(r)');
title('Cross-Correlation: Channel 1 vs Channel 2');

% Interpretation
if max(results_cross.g) > 1.2
    fprintf('Channels are co-localized\n');
elseif max(results_cross.g) < 0.8
    fprintf('Channels are anti-correlated\n');
else
    fprintf('Channels are independent\n');
end
```

Comparing clustered vs random data:
```matlab
% Generate test data
SMD_Clustered = generateClusteredData();
SMD_Random = generateRandomData();

PC = smi_cluster.PairCorrelation(SMF);
PC.ROI = [0, 10000, 0, 10000];

% Auto-correlation of each
XY_Clust = [SMD_Clustered.X * 100, SMD_Clustered.Y * 100];
XY_Rand = [SMD_Random.X * 100, SMD_Random.Y * 100];

PC.BaseName = 'Clustered';
results_clust = PC.pair_correlation(XY_Clust);

PC.BaseName = 'Random';
results_rand = PC.pair_correlation(XY_Rand);

% Plot comparison
figure;
subplot(1,2,1);
plot(results_clust.r, results_clust.g, 'b-', 'LineWidth', 2);
hold on;
yline(1, 'k--');
xlabel('r (nm)');
ylabel('g(r)');
title('Clustered Data');
xlim([0, 500]);

subplot(1,2,2);
plot(results_rand.r, results_rand.g, 'r-', 'LineWidth', 2);
hold on;
yline(1, 'k--');
xlabel('r (nm)');
ylabel('g(r)');
title('Random Data (Control)');
xlim([0, 500]);
```

Gaussian fitting to extract cluster parameters:
```matlab
% Compute pair correlation
results = PC.pair_correlation(XY);

% Fit Gaussian model
[estimates, errors, model] = PC.pc_GaussFit(...
    results.r, results.g, results.rmax, results.rho);

% Extract parameters
sigma_fit = estimates.sigma;      % Cluster size
amplitude = estimates.amplitude;  % Cluster strength

fprintf('Fitted cluster size: %.1f ± %.1f nm\n', ...
    sigma_fit, errors.sigma);
fprintf('Fitted amplitude: %.2f ± %.2f\n', ...
    amplitude, errors.amplitude);

% Plot with fit
figure;
plot(results.r, results.g, 'b-', 'LineWidth', 2);
hold on;
plot(model);
xlabel('Distance r (nm)');
ylabel('g(r)');
legend('Data', 'Gaussian fit', 'Location', 'best');
title(sprintf('Cluster size: %.1f nm', sigma_fit));
```

Multi-ROI combined analysis:
```matlab
% Prepare ROI structure
n_ROIs = 10;
RoI = cell(1, n_ROIs);
for i = 1:n_ROIs
    % Load or define each ROI
    load(sprintf('ROI_%d.mat', i), 'SMD');
    XY = [SMD.X * PixelSize, SMD.Y * PixelSize];

    RoI{i}.ROI = [0, 5000, 0, 5000];
    RoI{i}.X = {XY(:,1)};
    RoI{i}.Y = {XY(:,2)};
end

% Combined auto-correlation
PC = smi_cluster.PairCorrelation(SMF);
PC.BaseName = 'Combined';
results = PC.pair_correlation_ROIcombined(1, n_ROIs, RoI);

% Plot averaged result
figure;
plot(results.r, results.c, 'b-', 'LineWidth', 2);
hold on;
% Error bars from combining ROIs
errorbar(results.r, results.c, results.dc, 'b.');
yline(1, 'k--');
xlabel('Distance r (nm)');
ylabel('g(r)');
title(sprintf('Combined Analysis (%d ROIs)', n_ROIs));
```

Analyzing length scales:
```matlab
results = PC.pair_correlation(XY);

% Find peaks in g(r)
[peaks, locs] = findpeaks(results.g, results.r, ...
    'MinPeakHeight', 1.2, 'MinPeakDistance', 50);

fprintf('Clustering length scales:\n');
for i = 1:length(peaks)
    fprintf('  Peak %d: r = %.1f nm, g(r) = %.2f\n', ...
        i, locs(i), peaks(i));
end

% Find where g(r) returns to 1 (cluster size)
above_one = results.g > 1.1;
if any(above_one)
    cluster_extent = max(results.r(above_one));
    fprintf('Clustering extends to: %.1f nm\n', cluster_extent);
end
```

**See Also:**
- [Clustering Analysis Example](../examples/clustering-analysis.md)
- [Multi-Channel Example](../examples/multi-channel.md)

**Citations:**

Pair Correlation:
- S. L. Veatch, B. B. Machta, S. A. Shelby, E. N. Chiang, D. A. Holowka, B. A. Baird, "Correlation functions quantify super-resolution images and estimate apparent clustering due to over-counting", PLoS ONE 7(2) (2012) e31457.

---

### Analysis Interface Classes

#### ClusterInterface

**Purpose:** High-level interface functions for ROI-based clustering workflows, particularly for analyzing results from BaGoL clustering.

**Key Concept:** Provides workflows for defining ROIs, combining results across cells/conditions, and generating statistics. Especially useful for batch processing and standardized analysis pipelines.

**Class Definition:**
```matlab
classdef ClusterInterface < handle
```

**Static Methods:**

`defineROIs(pathname, files, Pixel2nm, RT, oneROI, ROI_sizes)` - Interactive ROI selection
```matlab
smi_cluster.ClusterInterface.defineROIs(...
    '/data', {'Cell1.mat', 'Cell2.mat'}, 100, 'rect', false, [1, 1]);
```

**Inputs:**
- `pathname`: Data directory
- `files`: Cell array of file names
- `Pixel2nm`: Conversion factor
- `RT`: ROI type ('rect', 'poly', 'circle')
- `oneROI`: Use single ROI for all files
- `ROI_sizes`: [width, height] in pixels

`defineBaGoLROIs(pathnameR, filesR, pathnameB, filesB, MAPNfile, OriginLLvsUL)` - BaGoL ROI definition
```matlab
smi_cluster.ClusterInterface.defineBaGoLROIs(...
    '/SR', SR_files, '/BaGoL', BaGoL_files, 'MAPN.mat', 'LL');
```

Defines ROIs overlaying super-resolution and BaGoL results.

`combineBaGoLROIs(pathnameR, filesR, pathnameB, filesB, MAPNfile, keep_numbering, GaussIm)` - Combine BaGoL ROIs
```matlab
smi_cluster.ClusterInterface.combineBaGoLROIs(...
    '/SR', SR_files, '/BaGoL', BaGoL_files, 'MAPN.mat', true, []);
```

Combines ROIs from BaGoL analysis across multiple cells.

`combineResults(Files, analysis_dir, out_file)` - Merge analysis results
```matlab
smi_cluster.ClusterInterface.combineResults(...
    file_list, '/results', 'combined_output.mat');
```

`filterROIs(pathname, files, filter)` - Filter ROIs by criteria
```matlab
[n_ROIs, RoI] = smi_cluster.ClusterInterface.filterROIs(...
    '/data', files, struct('MinPoints', 100));
```

`singleCondition(...)` - Single condition analysis workflow
```matlab
smi_cluster.ClusterInterface.singleCondition(...
    pathname, files, algorithm_range, E_range, minPts_range, ...
    Pixel2nm, base_name, A_ROI, doHopkins, doSigmaActual, Alpha, SC);
```

Complete workflow: clustering, statistics, visualization for single condition.

`plotROIDriver(...)` - Batch ROI visualization
```matlab
smi_cluster.ClusterInterface.plotROIDriver(...
    PixelSize, options, start_datadir, SaveDir, IncludeCell);
```

`plotROI(...)` - Individual ROI visualization
```matlab
smi_cluster.ClusterInterface.plotROI(...
    opt, pathnameC, filesC, pathnameB, filesB, PixelSize, SaveDir);
```

`combinedStatistics1/2(...)` - Combined statistical analysis
```matlab
smi_cluster.ClusterInterface.combinedStatistics1(...
    SC, pathname, files, base_name, A_ROI, doHopkins);
```

**Usage Example:**

```matlab
% Define analysis parameters
pathname = '/data/experiment1';
files = {'Cell1.mat', 'Cell2.mat', 'Cell3.mat'};
Pixel2nm = 100;

% Interactive ROI definition
smi_cluster.ClusterInterface.defineROIs(...
    pathname, files, Pixel2nm, 'rect', false, [50, 50]);

% Load ROI structure (saved by defineROIs)
load(fullfile(pathname, 'RoI.mat'), 'RoI');
n_ROIs = length(RoI);

% Single condition analysis
SC = smi_cluster.StatisticsClustering();
algorithm_range = {'DBSCAN'};
E_range = [30, 50, 70];
minPts_range = [3];
A_ROI = 5000^2;  % nm^2

smi_cluster.ClusterInterface.singleCondition(...
    pathname, files, algorithm_range, E_range, minPts_range, ...
    Pixel2nm, 'Experiment1', A_ROI, true, true, 2.0, SC);
```

---

#### PairAnalysis

**Purpose:** High-level interface for two-channel pair analysis workflows including clustering, statistics, and correlation.

**Static Methods:**

`defineROIs2(Files1, Files2, ...)` - Two-channel ROI definition
```matlab
n_ROIs = smi_cluster.PairAnalysis.defineROIs2(...
    Files1, Files2, Pixel2nm, Color, GaussIm, ...
    ROI_sizes, ResultsDir, RegistrationNeeded, RegistrationViaDC);
```

`doAnalysis(n_ROIs, RoI, ...)` - Complete pair analysis
```matlab
[results_pcc, resultsRC_pcc, results_ss, results_c, results_cs, ...
 results_ls, results_o1, results_o2] = ...
    smi_cluster.PairAnalysis.doAnalysis(...
        n_ROIs, RoI, ROI_sizes, desc, particles, results_dir, ...
        options, PixelSize, HistBinSize, RmaxAxisLimit, E, minPts, ...
        PlotNonOverlap, Color);
```

Returns results for:
- pcc: Pair cross-correlation
- ss: Simple statistics
- c: Clustering
- cs: Cluster separation
- ls: Localization separation
- o: Overlap analysis

`doClustering(n_ROIs, RoI, desc, results_dir, plotting, PixelSize, E, minPts)` - Clustering analysis
```matlab
results_c = smi_cluster.PairAnalysis.doClustering(...
    n_ROIs, RoI, 'MyAnalysis', '/results', true, 100, 50, 3);
```

`doPairCorr(n_ROIs, RoI, ...)` - Pair correlation
```matlab
[results_pcc, resultsRC_pcc] = ...
    smi_cluster.PairAnalysis.doPairCorr(...
        n_ROIs, RoI, ROI_sizes, desc, results_dir, combined, ...
        plotting, HistBinSize, RmaxAxisLimit);
```

`doBiStats(n_ROIs, RoI, desc, particles, results_dir, combined)` - Bivariate statistics
```matlab
results_bi = smi_cluster.PairAnalysis.doBiStats(...
    n_ROIs, RoI, 'Analysis', particles, '/results', true);
```

`doOverlap(n_ROIs, RoI, results_c, l12, ...)` - Overlap analysis
```matlab
results_o = smi_cluster.PairAnalysis.doOverlap(...
    n_ROIs, RoI, results_c, l12, desc, particles, ...
    results_dir, PlotNonOverlap, Color, plotting);
```

`doSimpleStats(n_ROIs, RoI, PixelSize, desc, particles, results_dir, combined)` - Basic statistics
```matlab
results_ss = smi_cluster.PairAnalysis.doSimpleStats(...
    n_ROIs, RoI, PixelSize, 'Stats', particles, '/results', true);
```

**Usage Example:**

```matlab
% Setup two-channel analysis
Files1 = {'Cell1_Ch1.mat', 'Cell2_Ch1.mat'};
Files2 = {'Cell1_Ch2.mat', 'Cell2_Ch2.mat'};
Pixel2nm = 100;
Color = [1, 0, 0; 0, 1, 0];  % Red, Green

% Define ROIs
n_ROIs = smi_cluster.PairAnalysis.defineROIs2(...
    Files1, Files2, Pixel2nm, Color, [], ...
    [50, 50], '/results', true, false);

% Load ROI structure
load('/results/RoI.mat', 'RoI');

% Complete analysis
options = struct('DoAll', true);
[results_pcc, ~, results_ss, results_c, ~, ~, results_o] = ...
    smi_cluster.PairAnalysis.doAnalysis(...
        n_ROIs, RoI, [50, 50], 'TwoChannel', ...
        {'Protein1', 'Protein2'}, '/results', options, ...
        100, 100, 500, 50, 3, false, Color);

% Interpret results
fprintf('Pair correlation peak: %.2f\n', max(results_pcc.g{1}));
fprintf('Overlap fraction: %.1f%%\n', 100*results_o.overlap_fraction);
fprintf('Channel 1 clusters: %d\n', results_c.nC{1});
fprintf('Channel 2 clusters: %d\n', results_c.nC{2});
```

---

## Common Usage Patterns

### Basic Clustering Workflow

```matlab
% 1. Setup
SMF = smi_core.SingleMoleculeFitting();
c = smi_cluster.Clustering(SMF);
c.PixelSize = 100;  % nm per pixel

% 2. Load data
load('Results.mat', 'SMD');

% 3. Cluster
E = 50;      % nm
minPts = 3;
[nC, C, centers, ptsI] = c.cluster('DBSCAN', SMD, E, minPts);

% 4. Statistics
results = c.clusterStats(SMD, C, centers);

% 5. Visualize
c.plotClusters(SMD, C, centers, ptsI, 'My Clusters');

% 6. Report
fprintf('Found %d clusters\n', nC);
fprintf('Clustered: %.1f%%\n', ...
    100 * results.n_clustered / results.n_points);
fprintf('Mean size: %.1f points\n', mean(results.n_pts));
```

### Multi-Algorithm Comparison

```matlab
algorithms = {'DBSCAN', 'Voronoi', 'Hierarchical'};
results_all = cell(1, length(algorithms));

for i = 1:length(algorithms)
    alg = algorithms{i};

    if strcmp(alg, 'Voronoi')
        c.Alpha = 2.0;
        [nC, C, centers, ptsI] = c.cluster(alg, SMD, [], 3);
    else
        [nC, C, centers, ptsI] = c.cluster(alg, SMD, 50, 3);
    end

    results_all{i} = c.clusterStats(SMD, C, centers);

    fprintf('%s: %d clusters, %.1f%% clustered\n', ...
        alg, nC, ...
        100 * results_all{i}.n_clustered / results_all{i}.n_points);
end
```

### Statistical Validation

```matlab
% Test if data is clustered
SC = smi_cluster.StatisticsClustering(SMF);
SC.ROI = [0, 10000, 0, 10000];

% Hopkins test
results_hop = SC.hopkins({'MyData'}, SMD);
fprintf('Hopkins: %.3f ', results_hop.H{1});
if results_hop.H{1} < 0.5
    fprintf('(CLUSTERED)\n');
else
    fprintf('(not clustered)\n');
end

% Ripley's K
results_rip = SC.ripley({'MyData'}, SMD, 'Analysis');
fprintf('Max L(r): %.2f at r = %.0f nm\n', ...
    max(results_rip.L{1}), ...
    results_rip.r{1}(results_rip.L{1} == max(results_rip.L{1})));
```

### Two-Channel Colocalization

```matlab
% Load channels
load('Ch1.mat', 'SMD'); SMD1 = SMD;
load('Ch2.mat', 'SMD'); SMD2 = SMD;

% Pair correlation
PC = smi_cluster.PairCorrelation(SMF);
PC.ROI = [0, 10000, 0, 10000];

XY1 = [SMD1.X * 100, SMD1.Y * 100];
XY2 = [SMD2.X * 100, SMD2.Y * 100];

results = PC.pair_correlation(XY1, XY2);

% Bivariate Ripley
SC = smi_cluster.StatisticsClustering(SMF);
results_bi = SC.bivariateRipley({'Ch1_vs_Ch2'}, SMD1, SMD2, 'Coloc');

% Report
if max(results.g) > 1.5
    fprintf('Strong co-localization detected\n');
    fprintf('Peak g(r) = %.2f at r = %.0f nm\n', ...
        max(results.g), results.r(results.g == max(results.g)));
end
```

---

## Performance Considerations

### Large Datasets

For datasets with >100,000 localizations:

```matlab
% Disable expensive calculations
c.DoSigmaActual = false;  % Skip sigma_actual computation

% Use spatial subsampling for initial exploration
subsample_idx = randsample(numel(SMD.X), 10000);
SMD_sub = smi_core.SingleMoleculeData.isolateSubSMD(SMD, subsample_idx);

% Quick clustering test
[nC, C, centers, ptsI] = c.cluster('DBSCAN', SMD_sub, 50, 3);

% Once parameters optimized, run on full dataset
[nC_full, C_full, centers_full, ptsI_full] = ...
    c.cluster('DBSCAN', SMD, 50, 3);
```

### Memory Management

```matlab
% For very dense ROIs, skip memory-intensive calculations
c.DoSigmaActual = false;  % Avoids pdist() on large clusters

% Process ROIs individually
for i = 1:n_ROIs
    SMD_roi = extractROI(SMD, RoI{i});
    [nC, C, centers, ptsI] = c.cluster('DBSCAN', SMD_roi, 50, 3);
    results{i} = c.clusterStats(SMD_roi, C, centers);

    % Clear variables to free memory
    clear SMD_roi C centers
end
```

---

## Troubleshooting

### Clustering Issues

**No clusters found**
```matlab
% Solutions:
% 1. Check units (pixels vs nm)
% 2. Increase epsilon
% 3. Decrease minPts
% 4. Visualize data first

% Debug
fprintf('Data range: X [%.1f, %.1f], Y [%.1f, %.1f] nm\n', ...
    min(SMD.X*c.PixelSize), max(SMD.X*c.PixelSize), ...
    min(SMD.Y*c.PixelSize), max(SMD.Y*c.PixelSize));
```

**Everything clusters together**
```matlab
% Solutions:
% 1. Decrease epsilon
% 2. Increase minPts
% 3. Check for duplicates

% Remove duplicates
[~, unique_idx] = unique([SMD.X, SMD.Y], 'rows');
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, unique_idx);
```

### Pair Correlation Issues

**Noisy g(r)**
```matlab
% Increase bin size for smoothing
PC.HistBinSize = 200;  % nm (larger bins)

% Or use more data
PC.ROI = []; % Use full field of view
```

**g(r) always ~1**
```matlab
% Check if ROI is too small
A_ROI = (PC.ROI(2) - PC.ROI(1)) * (PC.ROI(4) - PC.ROI(3));
fprintf('ROI area: %.2e nm^2\n', A_ROI);
% Need area >> cluster size^2

% Check localization density
density = numel(SMD.X) / A_ROI;
fprintf('Density: %.2e points/nm^2\n', density);
% Need sufficient density
```

---

## See Also

### Core Concepts
- [Architecture Overview](../core-concepts/architecture.md)
- [SMD Structure Guide](../core-concepts/smd-structure.md)
- [Data Flow](../core-concepts/data-flow.md)

### Workflows
- [BaGoL Clustering](../workflows/bagol-clustering.md)
- [ROI Analysis](../workflows/roi-analysis.md)

### Examples
- [Clustering Analysis Example](../examples/clustering-analysis.md)
- [Multi-Channel Analysis](../examples/multi-channel.md)

### How-To Guides
- [Visualize Results](../how-to/visualize-results.md)
- [Work with ROIs](../how-to/work-with-rois.md)

### Other API References
- [+smi_core Namespace](smi-core.md)
- [+smi Namespace](smi-namespace.md)

---

## Summary

The +smi_cluster namespace provides comprehensive spatial analysis for SMLM data:

**Clustering Algorithms:** DBSCAN for density-based detection, Voronoi for irregular shapes, H-SET for precision weighting, and Hierarchical for multi-scale analysis. Each has specific strengths for different biological questions.

**Statistical Methods:** Hopkins statistic tests clustering vs randomness, Ripley's K characterizes multi-scale organization, and bivariate methods quantify co-localization between channels.

**Pair Correlation:** Quantifies spatial relationships through g(r), providing characteristic length scales and clustering strength. Essential for membrane organization and protein complex analysis.

**Comprehensive Statistics:** Cluster characterization includes size metrics (points, area, radius), density metrics (points/area), shape metrics (compactness, circularity, convexity, solidity), and separation metrics (center-to-center, edge-to-edge).

**Analysis Interfaces:** High-level workflows for ROI-based batch processing, particularly useful for BaGoL clustering and multi-channel co-localization studies.

Choose clustering algorithms based on your data characteristics and biological question. Use statistical validation (Hopkins, Ripley) alongside clustering. Compare results across algorithms for robust conclusions. Pair correlation provides complementary information about spatial organization length scales.

**Word count: ~2450 words**
