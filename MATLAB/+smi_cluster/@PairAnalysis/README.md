### +smi_cluster/@PairAnalysis

A collection of functions that interface with the main smi.PairCorrelation
routines.
These are used in the MATLAB/examples simplePairCorrelation.m, which performs a
series of steps (one per section, mostly optional) to analyze pairs of 1-color
SR files (so 2-color) produced
by smi.SMLM by computing clusters and various statistics.

---

methods:
- **[defineROIs2](defineROIs2.m)**:
  Select ROIs simultaneously for label 1 and label 2 over all images
- **[doAnalysis](doAnalysis.m)**:
  Dispatches analyses to various helper functions depending options provided
- **[doBiStats](doBiStats.m)**:
  Pairwise mutual distances and bivariate Ripley's statistics for each ROI
- **[doClusterSep2](doClusterSep2.m)**:
  Find the nearest neighbor of each label 2 cluster to each label 1 cluster
  using center-to-center distances
- **[doClustering](doClustering.m)**:
  Clustering for each label in each ROI given epsilon (E) and minPts
- **[doLocSep2](doLocSep2.m)**:
  Find the nearest neighbor of each label 2 localization to each label 1
  localization
- **[doOverlap](doOverlap.m)**:
  Find overlaps between label 1/2 clusters and label 2/1 localizations
- **[doPairCorr](doPairCorr.m)**:
  Pair cross-correlation for each ROI and combined ROI
- **[doPlot2](doPlot2.m)**:
  Plot 2D ROIs
- **[doSimpleStats](doSimpleStats.m)**:
  Pearson's correlation and Manders' split coefficients for each ROI
- **[overlayBaGoLROIs](overlayBaGoLROIs.m)**:
  overlay label 1 onto label 2 Gaussian images produced from the input ROI file
  names into appropriately named files
