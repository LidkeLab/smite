## Publish

The Publish class batch-processes SR data assuming the data-containing
.h5 files follow a standard naming convention
(obj.CoverslipDir/Cell\*/Label\*/Data\*.h5).

[+smi/@Publish](../../MATLAB/+smi/@Publish/README.md) summarizes the
properties and methods of the Publish class.  See also tutorial
information on [SMLM](SMLM.md) and information on the
[SMF](DataStructures/SMF.md) data structure and its properties.

Below is a more annotated version of the script
[Example_Publish.m](../MATLAB/examples/Example_Publish.m).
[Example_Publish_generic.m](../MATLAB/examples/Example_Publish_generic.m)
is a similar example, but includes code at the end to perform additional
analyses using the same SMF parameters.

```
% This script will call smi.Publish to generate misc. results for an
% experiment on the sequential microscope.

%% Define the analysis parameters.
% Define the 'CoverslipDir'.
% NOTE: 'CoverslipDir' is the top-level directory which contains the
%        sub-directories 'CoverslipDir'/Cell*/Label*, which themselves
%        contain data in .h5 files 'CoverslipDir'/Cell*/Label*/Data*.h5
%        [MUST BE PROVIDED]
CoverslipDir = '/Experiment1';

% Prepare the SMF structure.
SMF = smi_core.SingleMoleculeFitting;
   % CameraType:     'EMCCD','SCMOS' (Default='EMCCD')
SMF.Data.CameraType = 'SCMOS';
   % Path to the camera calibration file (Default='')   [MUST BE PROVIDED]
SMF.Data.CalibrationFilePath = ...
   '/mnt/nas/lidkelab/sCMOS Calibrations/SequentialSR/GainCalibration_medianGain_2022_05_26.mat';
   % Camera back-projected pixel size (micrometers)
SMF.Data.PixelSize = 0.0954; % microns
   % Linear box size for fitting (Pixels)(Default=7)
SMF.BoxFinding.BoxSize = 8; % pixels
   % See fit class for options  (Default='XYNB')
SMF.Fitting.FitType = 'XYNBS';
   % Initial or fixed Sigma of 2D Gaussian PSF Model (Pixels) (Default=1)
SMF.Fitting.PSFSigma = 1.3; % pixels
   % Maximum allowed precision in x,y (Pixels)(Default=.2)
SMF.Thresholding.MaxXY_SE = 0.15; % pixels
   % Minimum accepted photons from fit (Default=100)
SMF.Thresholding.MinPhotons = 200;
   % Frame connection method being used (Default='LAP-FC')
SMF.FrameConnection.Method = 'LAP-FC';

% Alternatively, you may wish to prepare the SMF using the GUI.
% SMF.gui()

%% Prepare the smi.Publish class and run the standard analysis.
% The smi.Publish class requires the SMF (defined above) as well as the
% 'CoverslipDir'.  All other class properties specify which analyses to do.
    % SMF is a structure of parameters (see smi_core.SingleMoleculeFitting)
Publish = smi.Publish(SMF);
    % Directory containing the Cell*\Label*\Data*.h5 sub-directories.
Publish.CoverslipDir = CoverslipDir;
    % Verbosity of the main analysis workflow. (Default = 1)
Publish.Verbose = 1;
    % Flag to indicate SR results should be generated (Default = true)
Publish.GenerateSR = 1;
    % Flag to generate various imaging stats (Default = true)
Publish.GenerateImagingStats = 1;
    % Flag to generate overlay info. between channels (Default = false)
Publish.GenerateOverlayStats = 0;
    % Shift localizations based on brightfield results (Default = false)
Publish.ShiftToReg = 0;% can be useful for color overlay data, use with caution!

% Define trust regions, so that anything that seems to have a shift above the
% value below will be masked out.  Note that PixelSize is in microns/pixel, so
% the multiplying factor (0.2) is in units of microns.
    % Max. brightfield shift used to define overlay masks (pixels)
    % NOTE: This is defined in terms of brightfield pixels, e.g., units
    %       of obj.SMF.Data.PixelSize.
%Publish.MaxBrightfieldShift = 0.2 / SMF.Data.PixelSize; % pixels

% smi.Publish contains several useful methods, however we'll almost always
% just call performFullAnalysis().
Publish.performFullAnalysis();
```
Supposing the CoverslipDir is '/Experiment1', and data for 2 cells are
organized as follows:
```
/Experiment1/
  Cell_01/
    Label_01/
      Data_2022-8-11-17-10-2.h5
      Data_2022-8-12-9-1-51_bleaching.h5
    Label_02/
      Data_2022-8-12-11-32-50.h5
  Cell_02/
    Label_01/
      Data_2022-8-11-18-45-22.h5
      Data_2022-8-12-9-12-58_bleaching.h5
    Label_02/
      Data_2022-8-12-11-49-56.h5
```
Running the above script will produce files like:
```
/Experiment1/
  Results/
    Cell01_CircleOverlay_GM.png   (label 1 is green; label 2 is magenta)
    Cell01_GaussianOverlay_GM.png
    Cell_01/
      Label_01/
        Data_2022-8-11-17-10-2_GaussImage.png
        Data_2022-8-11-17-10-2/
          AlignRegDiffImageMovie.mp4
          AlignRegErrorPerDataset.png
          AlignRegErrorSignal.png
          AlignRegHistorySum.png
          AlignRegOverlayMovie.mp4
          AlignRegXCorrMaxima.png
          Data_2022-8-11-17-10-2_Background_Hist.png
          Data_2022-8-11-17-10-2_CircleImageDrift.png
          Data_2022-8-11-17-10-2_CircleImage.png
          Data_2022-8-11-17-10-2_Connected_emitters_Hist.png
          Data_2022-8-11-17-10-2_CumDriftCorrection.png
          Data_2022-8-11-17-10-2_DriftCorrection.png
          Data_2022-8-11-17-10-2_DriftImage.png
          Data_2022-8-11-17-10-2_FitsPerFrame.png
          Data_2022-8-11-17-10-2_GaussImage.png
          Data_2022-8-11-17-10-2_Intensity_Hist.png
          Data_2022-8-11-17-10-2_PSFSigma_Hist.png
          Data_2022-8-11-17-10-2_P_value_Hist.png
          Data_2022-8-11-17-10-2_Results.mat
          Data_2022-8-11-17-10-2_X_std_error_Hist.png
          Data_2022-8-11-17-10-2_Y_std_error_Hist.png
          DiffImageHistogram.png
          SSIM.png
          XCorrPlots/
            XCorrSequence1.fig
            XCorrSequence1.png
            XCorrSequence2.fig
            XCorrSequence2.png
            ...
      Label_02/
        Data_2022-8-12-11-32-50_GaussImage.png
        Data_2022-8-12-11-32-50/
          ...
    Cell02_CircleOverlay_GM.png
    Cell02_GaussianOverlay_GM.png
    Cell_02/
      ...
    ResultsStructs/
      Cell_01_Label_01_Data_2022-8-11-17-10-2_Results.mat
      Cell_01_Label_02_Data_2022-8-12-11-32-50_Results.mat
      Cell_02_Label_01_Data_2022-8-11-18-45-22_Results.mat
      Cell_02_Label_02_Data_2022-8-12-11-49-56_Results.mat

```
- AlignReg\* are various summary plots dealing with brightfield registration
  when using a sequential microscope (see also XCorrPlots below).
- Circle images use circles with radii proportional to the standard error to
  represent localizations.
- Drift images are color coded to show the movement over time (frame number)
  of localizations.
- Gaussian images use Gaussian blobs to represent localizations.
- Overlays are produced for 2-color (2-label) images.
- PSFSigma is the point spread function sigma or standard error.
- P_value refers to P-values from emitter fitting.
- Results,mat files contain the SMF and SMD structures generated for the
  given analysis.
- The single color Gaussian images and the Results files are duplicated for
  convenience in managing directories containing multiple analyses.
- XCorrPlots contains plots related to the brightfield registration process of
  the acquistion when using a sequential microscope, in which information about
  the cross-correlation process used to correct for sample drift during
  brightfield registration is displayed.
