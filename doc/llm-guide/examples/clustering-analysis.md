---
title: "Clustering Analysis Example"
category: "examples"
level: "intermediate"
tags: ["example", "clustering", "spatial-analysis", "dbscan", "voronoi", "statistics"]
prerequisites: ["../getting-started/installation.md", "basic-localization.md"]
related: ["../workflows/bagol-clustering.md", "../how-to/visualize-results.md"]
summary: "Complete example of spatial clustering analysis on SMLM data using multiple algorithms"
estimated_time: "20 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Clustering Analysis Example

## Purpose

This example demonstrates complete spatial clustering analysis of super-resolution localization data. You'll perform clustering with multiple algorithms (DBSCAN, Voronoi), compute cluster statistics (size, density, shape), visualize clusters with color-coded results, and perform pair correlation analysis. This is essential for analyzing protein organization, membrane domains, and spatial relationships in SMLM data.

## Prerequisites

- smite installed and working
- Understanding of basic localization (see basic-localization.md)
- Basic MATLAB knowledge
- 15-20 minutes to run the code

## Overview

This example covers:
1. Generating clustered and random test data
2. Performing DBSCAN clustering (density-based)
3. Performing Voronoi clustering (tessellation-based)
4. Computing comprehensive cluster statistics
5. Visualizing clusters with color coding and boundaries
6. Comparing clustering algorithms
7. Pair correlation analysis for spatial patterns

Everything runs in memory with simulated data, demonstrating all major clustering capabilities.

## Complete Working Code

Copy and run this complete example:

```matlab
%% Clustering Analysis Example
% Demonstrates spatial clustering analysis with multiple algorithms

%% Step 1: Generate Test Data - Clustered and Random Patterns
fprintf('=== Generating Test Data ===\n');

% Setup simulation
SaveDir = tempdir;
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.ResultsDir = SaveDir;

PixelSize = 100;  % nm per pixel
SZ = 100;         % Image size in pixels (10 microns)
SZnm = SZ * PixelSize;

% Generate clustered data: hextets (6-point clusters)
SIM = smi_sim.SimSMLM();
kmer = 6;                          % Hextet clusters
radius_kTet = 25;                  % Cluster radius (nm)
receptorDensity = 1;               % Clusters per um^2
SIM.Rho = receptorDensity / (1000 / PixelSize)^2;
SIM.StartState = 'Equib';
SIM.SZ = SZ;
SIM.simkTets(kmer, radius_kTet);
SMD_Clustered = SIM.genNoisySMD(SIM.SMD_Model);

% Keep only points within bounds
indx = SMD_Clustered.X < 0 | SMD_Clustered.X > SZ | ...
       SMD_Clustered.Y < 0 | SMD_Clustered.Y > SZ;
SMD_Clustered.X(indx) = [];
SMD_Clustered.Y(indx) = [];
SMD_Clustered.X_SE(indx) = [];
SMD_Clustered.Y_SE(indx) = [];

% Generate random (unclustered) data for comparison
Npts = numel(SMD_Clustered.X);
SMD_Random.X = SZ * rand(Npts, 1);
SMD_Random.Y = SZ * rand(Npts, 1);

fprintf('Generated %d localizations\n', Npts);
fprintf('Clustered data: %d hextets (radius %.0f nm)\n', ...
    ceil(Npts/kmer), radius_kTet);
fprintf('Random data: %d random points\n', Npts);
fprintf('Image size: %.0f × %.0f nm\n', SZnm, SZnm);

%% Step 2: Setup Clustering Object
fprintf('\n=== Configuring Clustering ===\n');

c = smi_cluster.Clustering(SMF);
c.PixelSize = PixelSize;
c.Timing = true;
c.Verbose = 1;
ROI = [0, SZnm, 0, SZnm];  % nm
A_ROI = (ROI(2) - ROI(1)) * (ROI(4) - ROI(3));  % Area (nm^2)

fprintf('Pixel size: %.0f nm\n', PixelSize);
fprintf('ROI area: %.2e nm^2\n', A_ROI);

%% Step 3: DBSCAN Clustering on Clustered Data
fprintf('\n=== DBSCAN Clustering (Clustered Data) ===\n');

% DBSCAN parameters
E_DBSCAN = 50;      % Epsilon: max distance within clusters (nm)
minPts = 3;         % Minimum points to form a cluster

% Perform DBSCAN clustering
[nC_DBSCAN, C_DBSCAN, centers_DBSCAN, ptsI_DBSCAN] = ...
    c.cluster('DBSCAN', SMD_Clustered, E_DBSCAN, minPts);

fprintf('DBSCAN (E = %d nm, minPts = %d):\n', E_DBSCAN, minPts);
fprintf('  Found %d clusters\n', nC_DBSCAN);
fprintf('  Isolated points: %d\n', length(ptsI_DBSCAN));

% Compute cluster statistics
results_DBSCAN = c.clusterStats(SMD_Clustered, C_DBSCAN, centers_DBSCAN);
fprintf('  Clustered points: %d (%.1f%%)\n', ...
    results_DBSCAN.n_clustered, ...
    100 * results_DBSCAN.n_clustered / results_DBSCAN.n_points);

%% Step 4: Detailed Cluster Statistics
fprintf('\n=== Detailed Cluster Statistics (DBSCAN) ===\n');

fprintf('Cluster size distribution:\n');
fprintf('  Singlets: %d\n', results_DBSCAN.numclust(1));
fprintf('  Doublets: %d\n', results_DBSCAN.numclust(2));
fprintf('  Multiplets: %d\n', results_DBSCAN.numclust(3));

if nC_DBSCAN > 0
    fprintf('\nCluster properties:\n');
    fprintf('  Mean points per cluster: %.1f\n', mean(results_DBSCAN.n_pts));
    fprintf('  Median points per cluster: %.0f\n', median(results_DBSCAN.n_pts));
    fprintf('  Max points in cluster: %d\n', max(results_DBSCAN.n_pts));

    % Size statistics (for clusters with ≥3 points)
    if ~isempty(results_DBSCAN.equiv_radii)
        large_clusters = results_DBSCAN.n_pts >= 3;
        fprintf('\nCluster geometry (clusters with ≥3 points):\n');
        fprintf('  Mean equivalent radius: %.1f nm\n', ...
            mean(results_DBSCAN.equiv_radii(large_clusters)));
        fprintf('  Mean area: %.0f nm^2\n', ...
            mean(results_DBSCAN.areas(large_clusters)));
        fprintf('  Mean compactness: %.3f\n', ...
            mean(results_DBSCAN.compactness));
        fprintf('  Mean circularity: %.3f\n', ...
            mean(results_DBSCAN.circularity));
    end

    % Density statistics
    if ~isempty(results_DBSCAN.n_pts_per_area)
        fprintf('\nCluster density:\n');
        fprintf('  Mean: %.2e points/nm^2\n', ...
            mean(results_DBSCAN.n_pts_per_area));
        fprintf('  Median: %.2e points/nm^2\n', ...
            median(results_DBSCAN.n_pts_per_area));
    end

    % Separation statistics
    if nC_DBSCAN > 1
        fprintf('\nCluster separation:\n');
        fprintf('  Min center-to-center: %.1f nm\n', ...
            results_DBSCAN.min_c2c_dist);
        fprintf('  Min edge-to-edge: %.1f nm\n', ...
            results_DBSCAN.min_e2e_dist);
        fprintf('  Mean nearest neighbor: %.1f nm\n', ...
            mean(results_DBSCAN.min_c2c_dists));
    end
end

%% Step 5: Voronoi Clustering on Same Data
fprintf('\n=== Voronoi Clustering (Clustered Data) ===\n');

% Voronoi parameters
c.Alpha = 1.5;      % Density threshold (ratio to overall density)
c.Valgorithm = 2;   % Algorithm: consider neighbors

% Perform Voronoi clustering
[nC_Voronoi, C_Voronoi, centers_Voronoi, ptsI_Voronoi] = ...
    c.cluster('Voronoi', SMD_Clustered, [], minPts);

fprintf('Voronoi (Alpha = %.1f, minPts = %d):\n', c.Alpha, minPts);
fprintf('  Found %d clusters\n', nC_Voronoi);
fprintf('  Isolated points: %d\n', length(ptsI_Voronoi));

% Compute cluster statistics
results_Voronoi = c.clusterStats(SMD_Clustered, C_Voronoi, centers_Voronoi);
fprintf('  Clustered points: %d (%.1f%%)\n', ...
    results_Voronoi.n_clustered, ...
    100 * results_Voronoi.n_clustered / results_Voronoi.n_points);

%% Step 6: Compare Algorithms
fprintf('\n=== Algorithm Comparison ===\n');
fprintf('Algorithm      Clusters  Clustered%%  Mean Size\n');
fprintf('------------------------------------------------\n');
fprintf('DBSCAN         %8d  %9.1f  %9.1f\n', ...
    nC_DBSCAN, ...
    100 * results_DBSCAN.n_clustered / results_DBSCAN.n_points, ...
    mean(results_DBSCAN.n_pts));
fprintf('Voronoi        %8d  %9.1f  %9.1f\n', ...
    nC_Voronoi, ...
    100 * results_Voronoi.n_clustered / results_Voronoi.n_points, ...
    mean(results_Voronoi.n_pts));

%% Step 7: Clustering Random Data (Control)
fprintf('\n=== Control Analysis (Random Data) ===\n');

% DBSCAN on random data
[nC_Random, C_Random, centers_Random, ptsI_Random] = ...
    c.cluster('DBSCAN', SMD_Random, E_DBSCAN, minPts);

fprintf('DBSCAN on random data:\n');
fprintf('  Found %d clusters\n', nC_Random);
fprintf('  Isolated points: %d\n', length(ptsI_Random));

if nC_Random > 0
    results_Random = c.clusterStats(SMD_Random, C_Random, centers_Random);
    fprintf('  Clustered points: %d (%.1f%%)\n', ...
        results_Random.n_clustered, ...
        100 * results_Random.n_clustered / results_Random.n_points);
    fprintf('  Mean cluster size: %.1f\n', mean(results_Random.n_pts));
else
    fprintf('  No clusters found (expected for random data)\n');
end

%% Step 8: Visualize Clusters with Color Coding
fprintf('\n=== Visualizing Clusters ===\n');

% Create figure with multiple panels
figure('Name', 'Clustering Analysis', 'Position', [50, 50, 1600, 800]);

% Panel 1: DBSCAN clustering on clustered data
subplot(2,3,1);
h1 = c.plotClusters(SMD_Clustered, C_DBSCAN, centers_DBSCAN, ptsI_DBSCAN, ...
    sprintf('DBSCAN (E=%d nm)', E_DBSCAN));
set(h1, 'Visible', 'on');  % Make visible
title(sprintf('DBSCAN: %d clusters', nC_DBSCAN));

% Panel 2: Voronoi clustering on clustered data
subplot(2,3,2);
h2 = c.plotClusters(SMD_Clustered, C_Voronoi, centers_Voronoi, ptsI_Voronoi, ...
    sprintf('Voronoi (Alpha=%.1f)', c.Alpha));
set(h2, 'Visible', 'on');
title(sprintf('Voronoi: %d clusters', nC_Voronoi));

% Panel 3: DBSCAN on random data (control)
subplot(2,3,3);
h3 = c.plotClusters(SMD_Random, C_Random, centers_Random, ptsI_Random, ...
    'Random Data (Control)');
set(h3, 'Visible', 'on');
title(sprintf('Random Control: %d clusters', nC_Random));

% Panel 4: Cluster size distribution (DBSCAN)
subplot(2,3,4);
if nC_DBSCAN > 0
    histogram(results_DBSCAN.n_pts, 'BinWidth', 1);
    xlabel('Points per Cluster');
    ylabel('Count');
    title('DBSCAN: Cluster Size Distribution');
    grid on;
end

% Panel 5: Cluster area distribution (DBSCAN)
subplot(2,3,5);
if ~isempty(results_DBSCAN.areas) && sum(results_DBSCAN.areas > 0) > 0
    histogram(results_DBSCAN.areas(results_DBSCAN.areas > 0), 20);
    xlabel('Cluster Area (nm^2)');
    ylabel('Count');
    title('DBSCAN: Area Distribution');
    grid on;
end

% Panel 6: Cluster density distribution (DBSCAN)
subplot(2,3,6);
if ~isempty(results_DBSCAN.n_pts_per_area)
    histogram(results_DBSCAN.n_pts_per_area, 20);
    xlabel('Density (points/nm^2)');
    ylabel('Count');
    title('DBSCAN: Density Distribution');
    grid on;
end

%% Step 9: Nearest Neighbor Distance Analysis
fprintf('\n=== Nearest Neighbor Analysis ===\n');

% Compare NN distances to random distribution
h_nn = c.nn_ROIrandom(SMD_Clustered, A_ROI, 'Clustered Data vs Random');

fprintf('Nearest neighbor analysis complete\n');
fprintf('  Expected for clustered data: NN distances shorter than random\n');

%% Step 10: Pair Correlation Analysis
fprintf('\n=== Pair Correlation Analysis ===\n');

% Setup pair correlation
PC = smi_cluster.PairCorrelation(SMF);
PC.BaseName = 'ClusterExample';
PC.ROI = ROI;
PC.PixelSize = PixelSize;
PC.HistBinSize = PixelSize;  % Internal pixel size for correlation
PC.Verbose = 1;
PC.Rmax_axis = 500;  % Plot up to 500 nm

% Convert SMD to nm for pair correlation
XY_Clustered = [SMD_Clustered.X * PixelSize, SMD_Clustered.Y * PixelSize];
XY_Random = [SMD_Random.X * PixelSize, SMD_Random.Y * PixelSize];

% Auto-correlation of clustered data
fprintf('\nComputing auto-correlation (clustered data)...\n');
results_PC_clustered = PC.pair_correlation(XY_Clustered);

% Auto-correlation of random data (control)
PC.BaseName = 'RandomControl';
fprintf('Computing auto-correlation (random data)...\n');
results_PC_random = PC.pair_correlation(XY_Random);

% Plot pair correlation comparison
figure('Name', 'Pair Correlation Analysis', 'Position', [100, 100, 1200, 500]);

% Clustered data pair correlation
subplot(1,2,1);
if PC.Lines
    plot(results_PC_clustered.r, results_PC_clustered.g, 'b-', 'LineWidth', 2);
else
    plot(results_PC_clustered.r, results_PC_clustered.g, 'b.');
end
hold on;
yline(1, 'k--', 'LineWidth', 1.5, 'Label', 'Random');
xlabel('Distance r (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('g(r)', 'FontSize', 12, 'FontWeight', 'bold');
title('Pair Correlation: Clustered Data');
grid on;
xlim([0, PC.Rmax_axis]);
legend('Clustered', 'Random level', 'Location', 'best');

% Random data pair correlation
subplot(1,2,2);
if PC.Lines
    plot(results_PC_random.r, results_PC_random.g, 'r-', 'LineWidth', 2);
else
    plot(results_PC_random.r, results_PC_random.g, 'r.');
end
hold on;
yline(1, 'k--', 'LineWidth', 1.5, 'Label', 'Expected');
xlabel('Distance r (nm)', 'FontSize', 12, 'FontWeight', 'bold');
ylabel('g(r)', 'FontSize', 12, 'FontWeight', 'bold');
title('Pair Correlation: Random Data (Control)');
grid on;
xlim([0, PC.Rmax_axis]);
legend('Random', 'Expected', 'Location', 'best');

fprintf('\nPair correlation interpretation:\n');
fprintf('  g(r) > 1: Points are clustered at distance r\n');
fprintf('  g(r) = 1: Random distribution at distance r\n');
fprintf('  g(r) < 1: Points are dispersed at distance r\n');

%% Step 11: Statistical Clustering Tests
fprintf('\n=== Statistical Clustering Tests ===\n');

% Setup statistics clustering
SC = smi_cluster.StatisticsClustering(SMF);
SC.BaseName = 'StatTests';
SC.ROI = ROI;
SC.Verbose = 1;

% Hopkins statistic (tests for clustering)
% H < 0.5 indicates clustering, H ~ 0.5 random, H > 0.5 regular
fprintf('\nHopkins statistic test:\n');
fprintf('  H < 0.5: clustered\n');
fprintf('  H = 0.5: random\n');
fprintf('  H > 0.5: regular/dispersed\n');

% Test clustered data
results_hopkins_clust = SC.hopkins({'Clustered'}, SMD_Clustered);
fprintf('Clustered data H = %.3f\n', results_hopkins_clust.H{1});

% Test random data
results_hopkins_rand = SC.hopkins({'Random'}, SMD_Random);
fprintf('Random data H = %.3f\n', results_hopkins_rand.H{1});

%% Step 12: Export Cluster Results
fprintf('\n=== Exporting Results ===\n');

% Create results summary structure
ClusterResults = struct();
ClusterResults.Algorithm = 'DBSCAN';
ClusterResults.Parameters.E = E_DBSCAN;
ClusterResults.Parameters.minPts = minPts;
ClusterResults.nClusters = nC_DBSCAN;
ClusterResults.nIsolated = length(ptsI_DBSCAN);
ClusterResults.PercentClustered = 100 * results_DBSCAN.n_clustered / results_DBSCAN.n_points;
ClusterResults.ClusterIndices = C_DBSCAN;
ClusterResults.ClusterCenters = centers_DBSCAN;
ClusterResults.ClusterSizes = results_DBSCAN.n_pts;
ClusterResults.ClusterAreas = results_DBSCAN.areas;
ClusterResults.ClusterDensities = results_DBSCAN.n_pts_per_area;

% Save results
save(fullfile(SaveDir, 'ClusterResults.mat'), 'ClusterResults', 'results_DBSCAN');
fprintf('Results saved to: %s\n', fullfile(SaveDir, 'ClusterResults.mat'));

% Create CSV with cluster statistics
if nC_DBSCAN > 0
    T = table();
    T.ClusterID = (1:nC_DBSCAN)';
    T.NumPoints = results_DBSCAN.n_pts';
    T.Area_nm2 = results_DBSCAN.areas';
    T.EquivRadius_nm = results_DBSCAN.equiv_radii';
    T.CenterX_nm = centers_DBSCAN(1,:)';
    T.CenterY_nm = centers_DBSCAN(2,:)';

    if ~isempty(results_DBSCAN.compactness)
        T.Compactness = results_DBSCAN.compactness';
        T.Circularity = results_DBSCAN.circularity';
    end

    csv_file = fullfile(SaveDir, 'ClusterStatistics.csv');
    writetable(T, csv_file);
    fprintf('Cluster statistics saved to: %s\n', csv_file);
end

%% Step 13: Summary Report
fprintf('\n========================================\n');
fprintf('CLUSTERING ANALYSIS SUMMARY\n');
fprintf('========================================\n');
fprintf('\nData:\n');
fprintf('  Total localizations: %d\n', Npts);
fprintf('  Field of view: %.0f × %.0f nm\n', SZnm, SZnm);
fprintf('  True cluster structure: %d-point clusters\n', kmer);

fprintf('\nDBSCAN Results:\n');
fprintf('  Parameters: E = %d nm, minPts = %d\n', E_DBSCAN, minPts);
fprintf('  Clusters found: %d\n', nC_DBSCAN);
fprintf('  Clustering percentage: %.1f%%\n', ...
    100 * results_DBSCAN.n_clustered / results_DBSCAN.n_points);
if nC_DBSCAN > 0
    fprintf('  Mean cluster size: %.1f points\n', mean(results_DBSCAN.n_pts));
    fprintf('  Median cluster size: %.0f points\n', median(results_DBSCAN.n_pts));
end

fprintf('\nVoronoi Results:\n');
fprintf('  Parameters: Alpha = %.1f, minPts = %d\n', c.Alpha, minPts);
fprintf('  Clusters found: %d\n', nC_Voronoi);
fprintf('  Clustering percentage: %.1f%%\n', ...
    100 * results_Voronoi.n_clustered / results_Voronoi.n_points);
if nC_Voronoi > 0
    fprintf('  Mean cluster size: %.1f points\n', mean(results_Voronoi.n_pts));
end

fprintf('\nStatistical Tests:\n');
fprintf('  Hopkins (clustered): %.3f (< 0.5 = clustered)\n', ...
    results_hopkins_clust.H{1});
fprintf('  Hopkins (random): %.3f (~ 0.5 = random)\n', ...
    results_hopkins_rand.H{1});

fprintf('\nControl Analysis:\n');
fprintf('  Random data clusters: %d\n', nC_Random);
fprintf('  (Expect few/no clusters in random data)\n');

fprintf('\n========================================\n');
fprintf('Analysis Complete!\n');
fprintf('========================================\n');
```

## What This Example Demonstrates

### Data Generation
- Creates realistic clustered data (hextet patterns)
- Generates random control data for comparison
- Demonstrates importance of controls in clustering analysis

### DBSCAN Clustering
- Density-based spatial clustering algorithm
- Identifies clusters based on point density
- Robust to noise and irregular cluster shapes
- Key parameters: epsilon (E) and minimum points (minPts)

### Voronoi Clustering
- Tessellation-based clustering algorithm
- Uses Voronoi diagram density analysis
- Good for irregularly shaped clusters
- Key parameter: Alpha (density ratio threshold)

### Comprehensive Statistics
- **Size metrics**: Points per cluster, areas, equivalent radii
- **Density metrics**: Points per area, spatial concentration
- **Shape metrics**: Compactness, circularity, convexity, solidity
- **Separation metrics**: Center-to-center and edge-to-edge distances

### Validation Methods
- Nearest neighbor distance analysis
- Pair correlation function g(r)
- Hopkins statistic for randomness testing
- Control analysis with random data

### Visualization
- Color-coded cluster plots with boundaries
- Statistical distributions (size, area, density)
- Pair correlation curves
- Algorithm comparison

## Expected Results

### DBSCAN on Clustered Data
- **Clusters found**: ~90-100 clusters (depends on density)
- **Clustering percentage**: ~80-95% (most points in clusters)
- **Mean cluster size**: ~6 points (matches hextet structure)

### Voronoi on Clustered Data
- **Clusters found**: Similar to DBSCAN but may vary
- **Different boundaries**: May merge/split differently than DBSCAN
- **Complementary information**: Good to compare algorithms

### Random Data Control
- **Clusters found**: 0-5 spurious clusters
- **Clustering percentage**: <10%
- **Validates specificity**: Algorithms don't see false patterns

### Pair Correlation
- **Clustered data**: g(r) > 1 at cluster length scales
- **Random data**: g(r) ≈ 1 at all distances
- **Peak position**: Near 2× cluster radius

### Hopkins Statistic
- **Clustered data**: H < 0.5 (typically 0.1-0.3)
- **Random data**: H ≈ 0.5 (typically 0.45-0.55)

## Understanding Clustering Parameters

### DBSCAN Parameters

**Epsilon (E)**: Maximum distance for connectivity
```matlab
% Tight clusters (small E)
E = 30;  % Only very close points cluster

% Loose clusters (large E)
E = 100;  % More distant points can cluster

% Rule of thumb: Set E slightly larger than expected intra-cluster spacing
E = 1.5 * expected_spacing;
```

**MinPts**: Minimum cluster size
```matlab
% Strict (fewer clusters)
minPts = 5;  % Need 5+ points

% Permissive (more clusters)
minPts = 3;  % Need 3+ points

% Rule of thumb: minPts = 2*dimensions + 1 (= 5 for 2D)
```

### Voronoi Parameters

**Alpha**: Density threshold ratio
```matlab
% Conservative (high specificity)
Alpha = 2.0;  % Only very dense regions

% Sensitive (high sensitivity)
Alpha = 1.2;  % Moderately dense regions

% Rule of thumb: Start with Alpha = 1.5-2.0
```

**Algorithm**: Voronoi calculation mode
```matlab
% Mode 1: Cell only (most conservative)
c.Valgorithm = 1;

% Mode 2: Cell + neighbors (recommended)
c.Valgorithm = 2;

% Mode 3: Median of cell + neighbors
c.Valgorithm = 3;
```

## Modifications to Try

### 1. Different Cluster Structures

```matlab
% Trimers instead of hextets
kmer = 3;
SIM.simkTets(kmer, radius_kTet);

% Larger clusters
kmer = 10;
SIM.simkTets(kmer, radius_kTet);

% Tighter clusters
radius_kTet = 15;  % nm
SIM.simkTets(kmer, radius_kTet);
```

### 2. Vary DBSCAN Parameters

```matlab
% Scan epsilon values
for E = [30, 50, 70, 100]
    [nC, C, centers, ptsI] = c.cluster('DBSCAN', SMD_Clustered, E, minPts);
    fprintf('E = %d: %d clusters\n', E, nC);
end
```

### 3. Compare All Algorithms

```matlab
% Test multiple algorithms
algorithms = {'DBSCAN', 'Voronoi', 'Hierarchal'};
for i = 1:length(algorithms)
    alg = algorithms{i};
    if strcmp(alg, 'Voronoi')
        [nC, C, centers, ptsI] = c.cluster(alg, SMD_Clustered, [], minPts);
    else
        [nC, C, centers, ptsI] = c.cluster(alg, SMD_Clustered, E_DBSCAN, minPts);
    end
    fprintf('%s: %d clusters\n', alg, nC);
end
```

### 4. Higher Density Analysis

```matlab
% Increase cluster density
receptorDensity = 5;  % More clusters per um^2
SIM.Rho = receptorDensity / (1000 / PixelSize)^2;
SIM.simkTets(kmer, radius_kTet);
```

### 5. Cross-Correlation Analysis

```matlab
% Generate two different channel data
SMD_Channel1 = SMD_Clustered;
SMD_Channel2 = SMD_Random;

% Cross-correlation
XY1 = [SMD_Channel1.X * PixelSize, SMD_Channel1.Y * PixelSize];
XY2 = [SMD_Channel2.X * PixelSize, SMD_Channel2.Y * PixelSize];
PC.BaseName = 'CrossCorr';
results_cross = PC.pair_correlation(XY1, XY2);
```

### 6. 3D Clustering

```matlab
% Generate 3D data
SMD_3D.X = SZ * rand(Npts, 1);
SMD_3D.Y = SZ * rand(Npts, 1);
SMD_3D.Z = SZ * rand(Npts, 1);

% 3D DBSCAN
[nC, C, centers, ptsI] = c.cluster('DBSCAN', SMD_3D, E_DBSCAN, minPts);

% 3D visualization
c.plotClusters3(SMD_3D, C, centers, ptsI, '3D Clusters');
```

## Interpreting Cluster Statistics

### Size Metrics

**Points per cluster**: Direct measure of cluster membership
- Compare to expected biological complex size
- Distribution shows heterogeneity

**Equivalent radius**: Radius of circle with same area
- Physical size interpretation
- Compare to diffraction limit (~250 nm)

**Cluster area**: Spatial extent in nm²
- Total footprint of cluster
- Relevant for membrane domains

### Density Metrics

**Points per area**: Concentration within clusters
- Higher = tighter packing
- Compare to background density

**Local density ratio**: Cluster vs background
- Shows enrichment factor
- Typically 5-50× for real clusters

### Shape Metrics

**Compactness**: 4π × Area / Perimeter²
- 1.0 = perfect circle
- <0.5 = elongated/irregular
- Indicates cluster organization

**Circularity**: Shape regularity
- Based on compact boundary
- Insensitive to small irregularities

**Convexity**: Perimeter ratio
- 1.0 = convex shape
- <1.0 = concave/irregular
- Detects complex shapes

**Solidity**: Area ratio
- 1.0 = solid/filled
- <1.0 = holes or irregular
- Identifies donut-like structures

### Separation Metrics

**Center-to-center distance**: Cluster spacing
- Minimum distance between centers
- Shows cluster organization

**Edge-to-edge distance**: True separation
- Accounts for cluster size
- Relevant for interaction studies

## Common Issues and Solutions

**Issue: No clusters found in clearly clustered data**

Solutions:
- Increase epsilon (E) - may be too small
- Decrease minPts - may be too strict
- Check pixel size is correct
- Verify data units (pixels vs nm)

**Issue: Everything clusters together**

Solutions:
- Decrease epsilon (E) - too large
- Increase minPts - too permissive
- Check for duplicate coordinates
- Verify coordinate scaling

**Issue: Different algorithms give very different results**

Solutions:
- Expected - algorithms have different philosophies
- DBSCAN: density-based, good for compact clusters
- Voronoi: tessellation-based, good for irregular shapes
- Try both and compare biologically relevant metrics
- Use algorithm that matches your question

**Issue: Random data shows clusters**

Solutions:
- Expected at low levels (<10%)
- Statistical fluctuations create apparent clusters
- This is why controls are essential
- Consider significance testing

**Issue: Pair correlation is noisy**

Solutions:
- Increase HistBinSize for smoother curves
- Use more localizations
- Increase ROI size for better statistics
- Try Ripley's K function for smoothed version

**Issue: Hopkins statistic is ambiguous**

Solutions:
- Values near 0.5 are inconclusive
- Use multiple tests (Hopkins, g(r), visual)
- Try different spatial scales
- Consider heterogeneous mixing (some clustered, some random)

## Biological Applications

### Membrane Protein Organization

```matlab
% Cluster membrane receptor data
% Typical values:
E = 50;       % nm, based on receptor size
minPts = 3;   % Minimum complex size
c.Alpha = 2.0;  % For Voronoi

% Metrics to report:
% - Cluster density (receptors per cluster)
% - Cluster size distribution
% - Percentage clustered vs monomeric
% - Nearest neighbor distances
```

### Nuclear Protein Domains

```matlab
% Larger, more diffuse clusters
E = 100;      % nm, larger domains
minPts = 10;  % More proteins per domain
c.Alpha = 1.5;  % Lower threshold

% Metrics to report:
% - Domain area
% - Protein density within domains
% - Domain spacing
% - Shape metrics (circularity, compactness)
```

### Viral Assembly Sites

```matlab
% Tight, well-defined clusters
E = 30;       % nm, tight packing
minPts = 5;   % Minimum complex size
c.Alpha = 2.5;  % High specificity

% Metrics to report:
% - Assembly site count per cell
% - Proteins per assembly site
% - Assembly site size
% - Pair correlation for organization
```

## Next Steps

After running this example:

1. **Apply to your data**: Load real SMLM results and cluster
2. **Optimize parameters**: Use controls to find optimal E and minPts
3. **Try H-SET**: For precision-weighted clustering (needs X_SE, Y_SE)
4. **BaGoL analysis**: See bagol-clustering.md for advanced clustering
5. **Statistical comparison**: Compare conditions/treatments
6. **Ripley's K**: For multi-scale spatial analysis
7. **Bivariate analysis**: Two-channel colocalization and correlation

## See Also

- [BaGoL Clustering](../workflows/bagol-clustering.md) - Advanced Bayesian clustering
- [Visualize Results](../how-to/visualize-results.md) - More visualization options
- [SMLM Workflow](../workflows/smlm-analysis.md) - Complete analysis pipeline
- MATLAB/examples/Example_Clustering.m - Additional examples
- MATLAB/examples/Example_StatisticsClustering.m - Statistical methods
- MATLAB/+smi_cluster/@Clustering/ - Source code and documentation
- MATLAB/+smi_cluster/@PairCorrelation/ - Pair correlation details
