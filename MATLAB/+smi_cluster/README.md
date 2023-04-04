### +smi_cluster

+smi_cluster is the namespace for the post-processing clustering classes of
***smite***:
- [@Clustering](@Clustering/README.md):
  clustering algorithms (DBSCAN, Voronoi, H-SET)
- [@ClusterInterface](@ClusterInterface/README.md): 
  interface functions for Clustering class
- [@PairAnalysis](@PairAnalysis/README.md): 
  interface functions for PairCorrelation class
- [@PairCorrelation](@PairCorrelation/README.md): 
  pair correlation (auto- and cross-) on ROIs
- [@StatisticsClustering](@StatisticsClustering/README.md): 
  clustering statistical analyses

In addition, there are functions:
- **[clusterSTDist](clusterSTDist.m)**:
  performs pre-clustering on localizations in SMD based
  on their spatiotemporal separations with cutoffs MaxFrameGap and MaxDist
- **[clusterSTSigma](clusterSTSigma.m)**:
  performs pre-clustering on localizations in SMD based on their spatiotemporal
  separations with cutoffs MaxFrameGap, NSigmaDev, MaxNN
