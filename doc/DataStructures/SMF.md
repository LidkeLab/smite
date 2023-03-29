SingleMoleculeFitting: A class defining the Single Molecule Fitting structure

The SMF structure is a structure of structures that collectively contain
all parameters required to go from raw data to an SMD results structure.
The SMF structure is an input of many smi methods. It
intended to be extensible to enable new analysis tools and methods.
The SMF class implements tools for working with SMF structures,
but the data structure itself is not an object of the class.

Parameters of sub-structures are explained in more detail in
the classes and methods that use them.  An incomplete list of classes
that use each sub-structure is listed in {}.

The SMF structure has the following sub-structures and fields:

```
SMF:  Fields that are structures:

Data:             {LoadData}
  FileName:       File name (cell array of char array)
  FileDir:        File directory (char array)
  ResultsDir:     Results directory (char array)(Default='FileDir/Results')
  AnalysisID:     ID tagged onto saved results (char array)(Default='')
  FileType:       Type of data specified by FileName. If using a custom
                  extension, you must set this field manually to the true
                  underlying file type (e.g., if using a .mat file saved
                  as exFile.spt, set obj.Data.FileType = 'mat')
                  (char array)(Default set to extension of FileName{1})
  DataVariable:   Name of variable saved in FileName which contains the
                  raw data. (char array)(Default='sequence')
  DatasetList:    List of datasets of the raw data to be analyzed.
                  (array of int32)(Default=int32([]))
  DatasetMods:    Cell array containing datasets to be used/excluded from
                  analysis (Mods <-> modifiers). This is meant to be the
                  user-facing lists which define DatasetList, meaning that
                  this is what would be set in the GUI. DatasetMods{1} will
                  contain an array of the "inclusion" dataset numbers and
                  DatasetMods{2} will contain an array of the "exclusion"
                  datasets. DatasetList will be set elsewhere (e.g.,
                  smi_core.LoadData) to include the set
                     intersect(intersect(1:NDatasets, DatasetMods{1}), ...
                     setdiff(1:NDatasets, DatasetMods{2}))
                  unless DatasetMods{1} is empty, in which case the first
                  parantheses term is dropped. For example, if
                  NDatasets = 20, and you only want to analyze datasets 1:5,
                  you can set DatasetMods{1} = 1:5. If you further decide to
                  exclude datsaets 2 and 4, you could set
                  DatasetMods{2} = [2, 4].
                  (cell array of int32 arrays)(Default={[]; []})
  CameraType:     'EMCCD','SCMOS' (Default='EMCCD')
  CameraGain:     Camera Gain, scalar or image (Default=1)
  CameraOffset:   Camera Offset, scalar or image (Default=0)
  CameraNoise:    Camera readnoise, scalar or image (Default=0)
  CalibrationFilePath: Path to the camera calibration file (Default='')
  RegistrationFilePath: Path to channel registration file (Default='')
  DataROI:        Region of interest of data file to be used (Default=[])
  FrameRate:      Data Collection Frame Rate (1/s)
  PixelSize:      Camera back-projected pixel size (micrometers)
  SEAdjust:       Standard error inflation per localization (Pixels)(Default=0)

BoxFinding:       {FindROI}
  BoxSize:        Linear box size for fitting (Pixels)(Default=7)
  BoxOverlap:     Overlap of boxes allowed (Pixels)(Default=2)
  MinPhotons:     Minimum number of photons from emitter (Default=200)

Fitting           {GaussMLE}
  PSFSigma:   Initial or fixed Sigma of 2D Gaussian PSF Model (Pixels)
              (Default=1)
  FitType:    See fit class for options  (Default='XYNB')
  NParams:    Number of fitting parameters (auto-set based on FitType)
  Iterations: Newton Raphson iterations (Default=20)
  ZFitStruct: Structure for astigmatic fitting:
      Ax:         Astigmatism fit parameter (see GaussMLE)
      Ay:         Astigmatism fit parameter (see GaussMLE)
      Bx:         Astigmatism fit parameter (see GaussMLE)
      By:         Astigmatism fit parameter (see GaussMLE)
      Gamma:      Astigmatism fit parameter (see GaussMLE)
      D:          Astigmatism fit parameter (see GaussMLE)

Thresholding      {ThresholdFits,SRA}
  On              Perform thresholding? (Default=true)
  MaxXY_SE:       Maximum allowed precision in x,y (Pixels)(Default=.2)
  MaxZ_SE:        Maximum allowed precision in z (Microns)(Default=.5)
  MinPValue:      Minimum accepted p-value from fit (Default=.01)
  AutoThreshLogL: Automatically threshold on LogL and ignore MinPValue
                  (Default = false)
  AutoThreshPrctile: Extrema percentile thrown out when computing LogL
                  auto-threshold (Default = 1e-4)
  MinPSFSigma:    Minimum PSF Sigma from fit (Pixels)(Default=.5);
  MaxPSFSigma:    Maximum PSF Sigma from fit (Pixels)(Default=2);
  MinPhotons:     Minimum accepted photons from fit (Default=100)
  MaxBg:          Maximum background accepted from fit (Default=Inf)
  InMeanMultiplier:   Determines maximum intensity accepted (Default=Inf)
  NNMedianMultiplier: Nearest neighbor acceptance region (Default=3)
  MinNumNeighbors:    Minimum number of neighbors in above (Default=0)

FrameConnection:  {FrameConnect,SRA}
  On              Perform frame connection? (Default=true)
  Method:         Frame connection method being used (Default='LAP-FC')
  MaxSeparation:  Maximum separation for connection (Pixels)(Default=1)
  LoS:            Minimum accepted p-value for connection (Default=.01)
  MaxFrameGap:    Maximum frame gap for connection (Frames)(Default=5)
  NSigmaDev:      SE multiplier for pre-cluster distance threshold (Default=5)
  NNearestClusters: Number of clusters used in density estimates (Default=2)
  NIterations:    Number of iterative FC attempts when Method=lap-fc
                  (Default=1)
  MinNFrameConns  Minimum accepted number of frame connections (Default=1)

DriftCorrection   {DriftCorrection,SRA}
 On               Perform drift correction? (Default=true)
 Method:          Drift correction method being used (Default='DC-KNN')
 BFRegistration   Was brightfield registration performed? (Default=true)
 L_intra          Intra-dataset threshold (Pixel)(Default=1)
 L_inter          Inter-dataset threshold (Pixel)(Default=2)
 PixelSizeZUnit   X/Y pixel size (3D drift correction) (um)(Default=0.1)
 PDegree          Degree intra-dataset fitting poly for drift rate (Default=1)

Tracking          {SPT}
  Method:         Type of method used for tracking (Default='CostMatrix')
  D:              Diffusion Constant (Pixels^2/Frame) (Default=0.01)
  TrajwiseD:      Use traj.-wise value for D (logical)(Default=true)
  K_on:           Off to On Rate (Frame^-1) (Default=.9)
  K_off:          On to Off Rate (Frame^-1) (Default=.1)
  MaxDistFF:      Maximum distance gap for frame-to-frame connection (Pixels)
                  (Default=5)
  MaxDistGC:      Maximum distance gap for Gap Closing (Pixels) (Default=10)
  MaxFrameGap:    Maximum frame gap for Gap Closing (Pixels) (Default=10)
  MinTrackLength: Minimum track length of trajectory (Frames) (Default=3)
  NIterMax:  Max. number of iterative tracking attempts (Integer)(Default=5)
  NIterMaxBatch:  Max. number of batch tracking iterations (Integer)
                  (Default = 5)
  MaxRelativeChange: Max. relative param. change to end iterations
                  (Default = 1e-5)
  MaxZScoreDist:  Max. abs(z-score) x/y jump size (Default=inf)
  MaxZScorePhotons: Max. abs(z-score) for photon diffs. (Default=inf)
  TryLowPValueLocs: Try to incorporate low p-val. locs. (Default=false)
```
