### +smi_cluster/@ClusterInterface

A collection of functions that interface with the main smi.Clustering routines.
These are used in the MATLAB/examples simpleROIcluster.m, which performs a
series of steps (one per section, mostly optional) to analyze SR files produced
by smi.SMLM by computing clusters and various statistics.

---

methods:
- **[combineBaGoLROIs](combineBaGoLROIs.m)**:
  This function will combine multiple BaGoL processed ROIs from multiple
  (biological) cells into appropriately named combined \_ROIs.mat files
- **[combineResults](combineResults.m)**:
  Multiple instantiations of a single condition potentially located in multiple
  directories are collected into a combined instance in a specified directory.
- **[combinedStatistics1](combinedStatistics1.m)**:
  For multiple conditions, experimental (e.g., resting vs. activated) and/or
  analytical (e.g., various DBSCAN parameter combinations), produce plots
  that compare the results for a variety of studies (e.g., Hopkins' statistic,
  number of clusters per ROI, circular equivalent cluster radii, etc.)
- **[combinedStatistics2](combinedStatistics2.m)**:
  This routine is very similar to combinedStatistics1, but allows the user to
  have fine control over line colors and types.
- **[defineBaGoLROIs](defineBaGoLROIs.m)**:
  Often, ROIs are defined from SR data, then need to be transferred to BaGoL
  (MAPN) processing of that SR data.  The BaGoL coordinates will replace the
  SR coordinates in the original ROI files.
- **[defineROIs](defineROIs.m)**:
  Choose ROIs of a fixed size over a series of images.  These are typically
  used for cluster analysis
- **[filterROIs](filterROIs.m)**:
  Filter out ROIs that meet certain criteria on their contents
- **[genMAPNIm1](genMAPNIm1.m)**:
  Produces a Gaussian blob image from either SMD or MAPN
- **[genSRMAPNOverlay1](genSRMAPNOverlay1.m)**:
  generates a multicolor overlay containing circles with radii proportional to
  the average localization precision (the generalized mean between the X and Y
  precisions) for each localization from an SMD and a BaGoL MAPN structure.
  Additional features to smi.BaGoL.genSRMAPNOVerlay have been added, in
  particular, SMD is magenta and MAPN green, the latter emphasized with a
  series of concentric rings for each localization
- **[plotROI](plotROI.m)**:
  Plot dot, Gaussian or circle images of SMD/MAPN coordinates per ROI per cell
- **[plotROI1](plotROI1.m)**:
  As above, but make a single image per cell of all ROIs
- **[plotROIDriver](plotROIDriver.m)**:
  Driver for plotROI and plotROI1
- **[singleCondition](singleCondition.m)**:
  Perform cluster analysis for comparison of experimental conditions
