### +smi/@SMLM

Single Molecule Localization Microscopy Analysis

This is a high-level class that provides complete analysis of SMLM data.
fullAnalysis/testFit performs an analysis on all/a selected dataset(s).
These routines require a Single Molecule Fitting (SMF) structure describing
the data and defining the analysis parameters.  This structure is either
produced earlier and fed into SMLM or it is created interactively when SMLM
is invoked.  The output is a Results.mat file holding the SMF and Single
Molecule Data (SMD) structures, the latter containing the processed data.  In
addition, various plots describing the processed data are created and placed
in a directory under Results identifying the dataset.  This identification is
derived from the original dataset's name, optionally with an analysis ID
appended.  See generatePlots (below) for more details on the plots produced.

---

```
properties:
    SMDPreThresh      % Keeps track of why localizations were filtered out
    SMD               % SMD structure with final analysis results
    SMF               % Single Molecule Fitting structure
    FitFramePreFC     % Fits per frame pre-frame connection
    PlotDo = []       % Plots to generate (all by default);see generatePlots comments
    SRImageZoom  = 20 % magnification factor for SR     images generated
    SRCircImZoom = 25 % magnification factor for circle images generated
    Verbose = 1       % Verbosity level
    VerboseTest = 3   % Verbosity level for testFits
    FullvsTest        % Logical value set by fullAnalysis or testFit to tell
                      % saveResults to make the proper call to generatePlots
    CalledByGUI=false % Keeps track of how fitting is called
```

---

analyzeAll loops over a list of datasets and creates an SMD.
If DatasetList not provided, use obj.SMD.Data.DatasetList .

analyzeAll flow:

```
analyzeAll:
   for n = DatasetList      (iterate over each dataset)
      analyzeDataset:
         LoadData           (load raw data [ADU] from the camera)
         DataToPhotons      (gain and offset corrections to convert to photons)
         LocalizeData       (produce localizations and put into SMD structure)
            Threshold       (threshold localizations generated)
         SingleMoleculeData (concatenate SMD into structure "SMDPreThresh")
         FrameConnection    (frame connection)
         DriftCorrection    (intra-dataset drift correction)
      SingleMoleculeData    (concatenate SMD into structure "SMD")
   end
   DriftCorrection (inter-dataset drift correction)
   Threshold       (produce statistics for rejected localizations)
```

---

generatePlots creates all histograms and plots for an SMD structure.

```
INPUT:
   obj          SMLM object
      obj.SMD      Single Molecule Data structure
      obj.SMF      Single Molecule Fitting structure
      obj.SRImageZoom    magnification factor for SR     images
      obj.SRCircImZoom   magnification factor for circle images
   PlotSaveDir1 Directory in which to save especially useful (priority 1)
                plots, like GaussIm
   PlotSaveDir2 Directory in which to save all the other (priority 2) plots
                (typically, a subdirectory of PlotSaveDir1)
   AnalysisID   Analysis ID, if non-empty, to add to the filenames generated
   ShowPlots:   Flag for showing plots on the screen (Default = false)
   PlotDo:      Plots to make chosen from the following list:
                "Photons"    intensity (estimated photons) histogram
                "Bg"         background intensity histogram
                "PSFSigma"   sigma of 2D Gaussian PSF model histogram
                "PValue"     P-value of fit histogram
                "X_SE"       standard error in estimated X position histogram
                "Y_SE"       standard error in estimated Y position histogram
                "Z_SE"       standard error in estimated Z position histogram
                "NCombined"  number of connection localizations histogram
                "DriftX"     cumulative x-drift
                "DriftY"     cumulative y-drift
                "DriftZ"     cumulative z-drift
                "CumDrift"   estimated 2D or 3D cumulative drift
                "Drift"      estimated 2D or 3D absolute drift
                "FitFrame"   number of fits per frame
                "DriftIm"    2D drift image from SR data
                "GaussIm"    2D Gaussian blob image from SR data
                "HistIm"     2D histogram image from SR data
                "CircleIm"   2D Circle image from SR data
                "CircleImDrift" 2D circle image color coded by time
                (Default is to make all plots)
                For example, PlotDo = ["PValue", "FitFrame", "DriftIm"]
                NOTE: plots will only be produced if there is corresponding
                     data in the SMD structure!

OUTPUT:
   The figures are saved in .png format in PlotSaveDir1/2.
```
