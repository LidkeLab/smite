### +smi_core/@FrameConnection

FrameConnection performs frame-connection on data in an SMD structure

This class contains methods to perform frame-connection and to do
associated tasks.  More specifically, the main usage of this class
is to combine a time series of localizations arising from the same
emitter into a single localization with precision greater than any
one of the localizations in the time series.

EXAMPLE USAGE:
```
   FC = smi_core.FrameConnection(SMD, SMF);
   [SMDCombined, SMD] = FC.performFrameConnection();
   Alternatively, you can use this class as a "function" as follows:
   [~, SMDCombined, SMD] = smi_core.FrameConnection(SMD, SMF, 1);
```

REQUIRES:
- smi_c_FrameConnection.mex
- Optimization Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

CITATION:
   David J. Schodt and Keith A. Lidke, "Spatiotemporal Clustering of
   Repeated Super-Resolution Localizations via Linear Assignment
   Problem", *Frontiers in Bioinformatics*, 2021 
   [https://doi.org/10.3389/fbinf.2021.724325](https://doi.org/10.3389/fbinf.2021.724325)

---

properties:
```
   SMF % see SingleMoleculeFitting class
   SMD % see SingleMoleculeData class
   Verbose = 1; % (Default = 1) Verbosity level of main workflow
```

---

methods:
- **[classicalFC](classicalFC.m)**:
  connects localizations in 'SMD' by simple thresholds
- **[combineLocalizations](combineLocalizations.m)**:
  combines localizations in SMD with the same ConnectID
- **[createCostMatrix](createCostMatrix.m)**:
  creates the cost matrix for frame connection
- **[estimateLocalDensity](estimateLocalDensity.m)**:
  estimates the local density for clustered localizations
- **[estimateRateParameters](estimateRateParameters.m)**:
  estimates rates from clustered localizations
- **[findConnected](findConnected.m)**:
  finds localizations that were frame connected
- **[hypothesisTestFC](hypothesisTestFC.m)**:
  connects localizations in 'SMD' by the hypothesis test method
- **[lapFC](lapFC.m)**:
  connects localizations in 'SMD' by solving a LAP
- **[linkClusters](linkClusters.m)**:
  updates cluster IDs based on the linkages in Link12
- **[organizeClusterData](organizeClusterData.m)**:
  organizes data according to SMD.ConnectID
- **[performFrameConnection](performFrameConnection.m)**:
  is the "run" method of the FrameConnection class
- **[revisedClassicalFC](revisedClassicalFC.m)**:
  connects localizations in 'SMD' by simple thresholds.  The "revision" with
  respect to classicalFC() is that the spatial thresholds are defined in 
  terms of the position standard errors of the localizations. 
- **[unitTest](unitTest.m)**:
  tests vital functionality of the class smi_core.FrameConnection
