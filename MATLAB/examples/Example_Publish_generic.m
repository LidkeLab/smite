% A generic Publish script to analyze multiple directories of data taken under
% the same conditions.
%
% This script will call smi.Publish to generate misc. results for an
% experiment on the sequential microscope.

%% Define the analysis parameters.
% Define the 'CoverslipDir'.
% NOTE: 'CoverslipDir' is the top-level directory which contains the
%        sub-directories 'CoverslipDir'\Cell*\Label*, which themselves
%        contain data in .h5 files 'CoverslipDir'\Cell*\Label*\Data*.h5
CoverslipDir = '/mnt/nas/cellpath/Actin Asters/22-3-8/HeLa_GFP-GPI_Actin';
SCMOSDir = '/mnt/nas/lidkelab/sCMOS Calibrations/Sequential SR/GainCalibration-2015-12-10-17-19-23.mat';

% Prepare the SMF structure.
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.AnalysisID = ''; % ### Set as desired ###
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = SCMOSDir;
SMF.Data.PixelSize = 0.0954; % microns
SMF.Fitting.PSFSigma = 1.3; % pixels
SMF.Fitting.FitType = 'XYNBS';
SMF.FrameConnection.Method = 'LAP-FC';
SMF.FrameConnection.On = true;
SMF.DriftCorrection.On = true;
SMF.Thresholding.MaxPSFSigma = 1.5; % pixels
SMF.Thresholding.MinPSFSigma = 0.9; % pixels
SMF.Thresholding.MinPhotons = 200; % dSTORM
%SMF.Thresholding.MinPhotons = 500; % DNA-PAINT
SMF.Thresholding.MaxXY_SE = 0.15; % pixels

% Additional filtering
%SMF.Data.SEAdjust = 3 / (1000 * SMF.Data.PixelSize); % nm converted to pixels
%SMF.FrameConnection.MinNFrameConns = 2;
%SMF.Thresholding.InMeanMultiplier = 2;
%SMF.Thresholding.NNMedianMultiplier = 3; % DNA-PAINT
%SMF.Thresholding.MinNumNeighbors = 2; % DNA-PAINT

%% Prepare the smi.Publish class and run the standard analysis.
% The smi.Publish class requires the SMF (defined above) as well as the
% 'CoverslipDir'.  All other class properties specify which analyses to do.
Publish = smi.Publish(SMF);
Publish.CoverslipDir = CoverslipDir;
Publish.Verbose = 1;
Publish.GenerateSR = 1;
Publish.GenerateImagingStats = 0;
Publish.GenerateOverlayStats = 0;
Publish.ShiftToReg = 0;% can be useful for color overlay data, use with caution!
Publish.SRImageZoom = 20;

% Define trust regions, so that anything that seems to have a shift above the
% value below will be masked out.  Note that PixelSize is in microns/pixel, so
% the multiplying factor (0.2) is in units of microns.
%Publish.MaxBrightfieldShift = 0.2 / SMF.Data.PixelSize; % pixels

% smi.Publish contains several useful methods, however we'll almost always
% just call performFullAnalysis().
Publish.performFullAnalysis()

% ----------

% Additional directories to Publish can be done by the following lines.

CoverslipDir = '/mnt/nas/cellpath/Actin Asters/22-3-8-b/HeLa_GFP-GPI_Actin';
Publish = smi.Publish(SMF);
Publish.CoverslipDir = CoverslipDir;
Publish.Verbose = 1;
Publish.GenerateSR = 1;
Publish.GenerateImagingStats = 0;
Publish.GenerateOverlayStats = 0;
Publish.ShiftToReg = 0;% can be useful for color overlay data, use with caution!
Publish.SRImageZoom = 20;
Publish.performFullAnalysis()

CoverslipDir = '/mnt/nas/cellpath/Actin Asters/22-3-8-b/HeLa_GFP-GPI_Actinb';
Publish = smi.Publish(SMF);
Publish.CoverslipDir = CoverslipDir;
Publish.Verbose = 1;
Publish.GenerateSR = 1;
Publish.GenerateImagingStats = 0;
Publish.GenerateOverlayStats = 0;
Publish.ShiftToReg = 0;% can be useful for color overlay data, use with caution!
Publish.SRImageZoom = 20;
Publish.performFullAnalysis()

CoverslipDir = '/mnt/nas/cellpath/Actin Asters/22-3-8-c/HeLa_GFP-GPI_Actin';
Publish = smi.Publish(SMF);
Publish.CoverslipDir = CoverslipDir;
Publish.Verbose = 1;
Publish.GenerateSR = 1;
Publish.GenerateImagingStats = 0;
Publish.GenerateOverlayStats = 0;
Publish.ShiftToReg = 0;% can be useful for color overlay data, use with caution!
Publish.SRImageZoom = 20;
Publish.performFullAnalysis()
