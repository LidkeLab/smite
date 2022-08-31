% This script will call smi.Publish to generate misc. results for an
% experiment on the sequential microscope.

%% Define the analysis parameters.
% Define the 'CoverslipDir'.
% NOTE: 'CoverslipDir' is the top-level directory which contains the
%        sub-directories 'CoverslipDir'\Cell*\Label*, which themselves
%        contain data in .h5 files 'CoverslipDir'\Cell*\Label*\Data*.h5
CoverslipDir = 'Y:\dschodt\analysis_test_directory\SMA_Publish_testing\21_5_4_HeLa_alpha_beta_tubulin';

% Prepare the SMF structure.
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = 'Y:\sCMOS Calibrations\Sequential SR\GainCalibration_medianGain_2022_05_26_14_51_08.mat';
SMF.Data.PixelSize = 0.0954; % microns
SMF.BoxFinding.BoxSize = 8; % pixels
SMF.Fitting.FitType = 'XYNBS';
SMF.Fitting.PSFSigma = 1.3; % pixels
SMF.Thresholding.MaxXY_SE = 0.15; % pixels
SMF.Thresholding.MinPhotons = 200;
SMF.FrameConnection.Method = 'LAP-FC';

% Alternatively, you may wish to prepare the SMF using the GUI.
% SMF.gui()

%% Prepare the smi.Publish class and run the standard analysis.
% The smi.Publish class requires the SMF (defined above) as well as the
% 'CoverslipDir'.  All other class properties specify which analyses to do.
Publish = smi.Publish(SMF);
Publish.CoverslipDir = CoverslipDir;
Publish.Verbose = 1;
Publish.GenerateSR = 1;
Publish.GenerateImagingStats = 1;
Publish.GenerateOverlayStats = 0;
Publish.ShiftToReg = 0; % can be useful for color overlay data, use with caution!

% Define trust regions, so that anything that seems to have a shift above the
% value below will be masked out.  Note that PixelSize is in microns/pixel, so
% the multiplying factor (0.2) is in units of microns.
%Publish.MaxBrightfieldShift = 0.2 / SMF.Data.PixelSize; % pixels

% smi.Publish contains several useful methods, however we'll almost always
% just call performFullAnalysis().
Publish.performFullAnalysis();
