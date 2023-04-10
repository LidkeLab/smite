### +smi/@Publish

Publish contains methods useful for batch processing of SR data.
  This class contains a collection of analysis/visualization methods
  useful for batch-processing of super-resolution data, particularly
  for data with multiple labels stored in independent files (e.g.,
  sequential super-resolution of two labels on multiple cells).

NOTE: This class is designed around .h5 files containing raw data
      stored in directories separating distinct cells and labels,
      with the directory names following the scheme
      Cell*\Label*\Data*.h5

REQUIRES:
- MATLAB 2018a or later (for Publish.genAlignStats())
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
- Curve Fitting Toolbox

---

```
properties:
    % Structure of parameters (see smi_core.SingleMoleculeFitting)
    SMF

    % Directory containing the Cell*\Label*\Data*.h5 sub-directories.
    CoverslipDir

    % Base directory for saving (Default set in performFullAnalysis())
    SaveBaseDir

    % Log file for errors (Default set in performFullAnalysis())
    LogFilePath

    % Label(s) to be analyzed (Default = [], analyze all labels)
    LabelID = [];

    % Cell(s) to be analyzed (Default = [], analyze all cells)
    CellList = [];

    % Zoom factor for output SR images (Default = 20)
    SRImageZoom = 20;

    % Zoom factor for circle images (Default = 50)
    SRCircleImageZoom = 50;

    % Flag to indicate SR results should be generated (Default = true)
    GenerateSR = true;

    % Flag to generate various imaging stats (Default = true)
    GenerateImagingStats = true;

    % Flag to generate overlay info. between channels (Default = false)
    GenerateOverlayStats = false;

    % Flag to perform analysis on bleaching results (Default = false)
    AnalyzeBleaching = false;

    % Apply brightfield drift correction (Default = false)
    % NOTE: If SMF.DriftCorrection.On=true, brightfield DC is still
    %       applied just before the post-processing DC.
    UseBrightfieldDC = false

    % Max. brightfield shift used to define overlay masks (pixels)
    % NOTE: This is defined in terms of brightfield pixels, e.g., units
    %       of obj.SMF.Data.PixelSize.
    MaxBrightfieldShift = inf

    % Shift localizations based on brightfield results (Default = false)
    ShiftToReg = false;

    % Verbosity of the main analysis workflow. (Default = 1)
    Verbose = 1;

    % Structure containing several analysis results.
    ResultsStruct = struct([]);
```

---

methods:
- **[computeBFShifts](computeBFShifts.m)**:
  computes subROI shifts from brightfield data
- **[computeRegCorrection](computeRegCorrection.m)**:
  computes the registration corrections made
- **[defineShiftMask](defineShiftMask.m)**:
  generates a mask from the array of shifts LocalImShifts
- **[estimateLocalCoordShifts](estimateLocalCoordShifts.m)**:
  estimates local shifts between two SMDs
- **[estimateLocalImShifts](estimateLocalImShifts.m)**:
  estimates local shifts between two images
- **[genAlignMovies](genAlignMovies.m)**:
  generates movies related to brightfield registration
- **[genAlignResults](genAlignResults.m)**:
  produces figures/info related to brightfield registration
- **[genAlignStats](genAlignStats.m)**:
  generates interesting plots from AlignRegStruct
- **[genAlignXCorr](genAlignXCorr.m)**:
  generates xcorr curve data for AlignRegStruct
- **[genBFMask](genBFMask.m)**:
  masks overlays based on brightfield shifts
- **[genOverlayResults](genOverlayResults.m)**:
  generates some two-color overlay results
- **[genSROverlays](genSROverlays.m)**:
  generates various types of  SR overlay images
- **[makeOverlayPlots](makeOverlayPlots.m)**:
  makes interesting plots from two color overlays
- **[performFullAnalysis](performFullAnalysis.m)**:
  is the main run method for the smi.Publish class
- **[plotXYRegError](plotXYRegError.m)**:
  plots the x,y registration error in a scatter plot
- **[processCell](processCell.m)**:
  will process data corresponding to CellName
- **[processLabel](processLabel.m)**:
  processes the SR data for the specified label
- **[shiftToBestReg](shiftToBestReg.m)**:
  shifts coordinates based on best alignment results
