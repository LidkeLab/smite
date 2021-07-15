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
SMF.Data.CalibrationFilePath = 'Y:\sCMOS Calibrations\Sequential SR\GainCalibration-2015-12-10-17-19-23.mat';
SMF.Data.PixelSize = 0.0954; % microns
SMF.BoxFinding.BoxSize = 8; % pixels
SMF.Fitting.FitType = 'XYNB';
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
Publish.GenerateOverlayStats = 1;

% smi.Publish contains several useful methods, however we'll almost always
% just call performFullAnalysis().
Publish.performFullAnalysis();