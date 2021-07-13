% This script will call smi.Publish to generate misc. results for an
% experiment on the sequential microscope.  

% Define the analysis parameters.
% Directory containing the Cell*\Label*\Data*.h5 sub-directories:
CoverslipDir = 'Y:\dschodt\analysis_test_directory\SMA_Publish_testing\21_5_4_HeLa_alpha_beta_tubulin';
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = ['Y:\sCMOS Calibrations\Sequential SR\', ...
    'GainCalibration-2015-12-10-17-19-23.mat'];
SMF.Data.PixelSize = 0.0954; % microns
SMF.BoxFinding.BoxSize = 8; % pixels
SMF.Fitting.FitType = 'XYNB';
SMF.Fitting.PSFSigma = 1.3; % pixels
SMF.Thresholding.MaxXY_SE = 0.15; % pixels
SMF.Thresholding.MinPhotons = 200;
SMF.FrameConnection.Method = 'LAP-FC';

% Prepare the smi.Publish class and run the standard analysis.
Publish = smi.Publish(SMF);
Publish.Verbose = 1;
Publish.GenerateSR = 1;
Publish.GenerateImagingStats = 1;
Publish.GenerateOverlayStats = 0;
Publish.CoverslipDir = CoverslipDir;
Publish.performFullAnalysis();
    
