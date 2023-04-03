### +smi_cluster/Clustering

See MATLAB/examples/Example_Clustering for an example using cluster, which
dispatches on "algorithm" to the appropriate routines.

- *DBSCAN_Daszykowski* is the recommended DBSCAN algorithm as it is both fast
  and stable under coordinate reordering.

  Martin Ester, Hans-Peter Kriegel and J\"org Sander and Xiaowei Xu, ``A
  Density-Based Algorithm for Discovering Clusters in Large Spatial Databases
  with Noise'', in _Proceedings of 2nd International Conference on Knowledge
  Discovery and Data Mining (KDD-96)_ edited by Evangelos Simoudis, Jiawei Han
  and Usama M. Fayyad, AAAI Press, 1996, 226--231, ISBN:1-57735-004-9,
  DOI:10.1.1.71.1980.

  M. Daszykowski, B. Walczak, D. L. Massart, Looking for Natural Patterns in
  Data. Part 1: Density Based Approach, Chemom. Intell. Lab. Syst. 56 (2001)
  83-92. 

- *Hierarchical* is Matlab's hierarchal clustering algorithm with some small
  additions.

- *Voronoi* is Florian Levet et al's Voronoi based algorithm.

  Florian Levet, Eric Hosy, Adel Kechkar, Corey Butler, Anne Beghin, Daniel
  Choquet and Jean-Baptiste Sibarita, ``SR-Tesseler: a method to segment and
  quantify localization-based super-resolution microscopy data'', _Nature
  Methods_, Volume 12, Number 11, 2015, 1065--1071 (DOI:10.1038/NMETH.3579).

- *H-SET* is the clustering implied by H-SET, in which the nodes that would be
  combined in normal H-SET are taken to be clusters here.

  Jia Lin, Michael J. Wester, Matthew S. Graus, Keith~A. Lidke and Aaron K.
  Neumann, ``Nanoscopic cell wall architecture of an immunogenic ligand in
  _Candida albicans_ during antifungal drug treatment'', _Molecular
  Biology of the Cell_, Volume 27, Number 6, March 15, 2016, 1002--1014
  (DOI: 10.1091/mbc.E15-06-0355, PMID: 26792838, PMCID: PMC4791122).

properties:
```
   % --- Generic properties.
   Font_props = {'FontSize', 15, 'FontWeight', 'bold'};
   Fig_ext = 'png';
   ResultsDir = '.'; % Directory to store results.
   ShrinkFactor = 0.5; % 'boundary' compactness factor in the range 0 - 1
   Verbose = 0; % verbosity level
   Xlim = [];   % x-axis limits if defined
   Ylim = [];   % y-axis limits if defined

   % Properties used by clusterStats.
   DoSigmaActual = true; % this calc. can run out of memory for very dense ROIs
                         % (see clusterStats), so make it optional

   % Properties used by voronoi_Levet.
   Alpha      = 2;
      % Ratio of local density / overall density for a point's Voronoi
      % region to be considered sufficiently dense for clustering purposes
   Valgorithm = 2;       % Voronoi algorithm to apply:
      % 1   [0]  calculations consider Voronoi regions only
      % 2   [1]  calculations consider Voronoi regions and their adjacent
      %          neighbors
      % 3   [1M] consider the median density of each cell and its neighbors
   Plotting   = false;   % Produce Voronoi plots
   PtIDs      = false;   % Label the points in the plots

   % --- H-SET properties.
   % H-SET collapse method: 'hierarchal_singlelabel' or 'trivial'.
   Method = 'hierarchal_singlelabel';
   LoS = 0.01;           % level of significance
   PixelSize = 100;      % conversion from pixels to nm
   PlotFigures = true;   % plot various cluster related figures
   Sigma_Reg = [10, 10]; % registration error in x, y (nm)
   Timing = true;        % produce timings for clustering
```
methods:
- **cluster**:
- **clusterSR**:
  Combine multiple clustered points into single localizations when appropriate
  via a top-down descent through a hierarchal dendrogram relationship between
  points
- **clusterStats**:
  produces large variety of statistics on computed clusters:
```
     nC                number of clusters
     C                 cell array of the indices of the points in each cluster
     SMD               (x, y) coordinates in the format
                        (1) SMD structure: SMD.X, SMD.Y with optional SMD.X_SE,
                            SMD.Y_SE (only needed for H-SET clustering) (pixel)
                        (2) N x 2 or N x 3 array of coordinates (nm)
     n_points          total number of points
     n_clustered       number of points in clusters
     n_isolated        number of points not in clusters
     n_pts             number of points per cluster
     numclust(1,2,3)   number of singlet, double, multiple clusters, where
                       singlet clusters include isolated points (see
                       SRcluster.m for an equivalent definition)
     singlet_faction   numclust(1) / sum(numclust)
     sigma_actual      actual (computed) sigma of each cluster, that is, the
                       standard deviation of the intracluster distances
     indices_hull      cell array of boundary hull indices relative to XY per
                       cluster
     areas             area of each cluster
     equiv_radii       equivalent radius of each cluster
     n_pts_per_area    number of points per area for clusters containing 3 or
                       more points
     perimeters        perimeter of each cluster
     compactness       4 pi area / perimeter^2 of each cluster
     min_c2c_dists     minimum center-to-center distances for each cluster
                       with respect to all the others
     min_e2e_dists     minimum edge-to-edge distances for each cluster convex
                       hull with respect to all the others
     min_c2c_dist      min(min_c2c_dists)
     min_e2e_dist      min(min_e2e_dists)
     nn_within_clust   nearest neighbor distances between points within
                       clusters only
```
- **cluster_HSET**:
  Perform the clustering implied by H-SET, in which the nodes that would be
  combined in normal H-SET are taken to be clusters here
- **cluster_voronoi**:
  Taking the density indices i_rho that identify points to be clustered, 
  generate the clusters C (their number given by nC)
- **dbscan_Daszykowski**:
  Clustering the data with Density-Based Scan Algorithm with Noise (DBSCAN);
  this seems to be the same or very similar to MATLAB's dbscan introduced in
  R2019a
- **edge2edge**:
  the minimum edge-to-edge distance between hull 1 and hull 2
- **hierarchal**:
  Form clusters such that any point in a cluster is within E of some other
  point in the same cluster
- **hierarchalSingleLabel**:
  Combine multiple clustered points into single labels when appropriate via a
  top-down descent through a hierarchal dendrogram relationship between points
- **my_ismemberBuiltinTypes**:
  Extracted and simplified from MATLAB's ismember.m for simplified usage
- **nn_ROIcombined**:
  Plots the mean particle nearest neighbor distances for a series of ROIs
- **nn_ROIrandom**:
  plots the PDF of nearest neighbor distances (NND) for points in a ROI vs a
  theoretical curve based on the same point density
- **nn_distances**:
  minimum nearest neighbor distances from each point in xy to the other points
- **plotClusters**:
  Plot and label the 2D clusters
- **plotClusters3**:
  Plot and label the 3D clusters
- **plotClustersSE**:
  Plot and label the 2D clusters, producing circles with radii proportional to
  the standard error (SE)
- **plot_voronoi**:
  plots the Voronoi diagram corresponding to (X, Y), coloring cells
  according to the density rho
- **plot_voronoi3**:
  plots the Voronoi diagram corresponding to (X, Y, Z), coloring
  the cells according to the density rho
- **singleLabelTest**:
  Tests if a cluster of points came from point source 
- **voronoi_Levet**:
  implements Voronoi diagram based clustering
