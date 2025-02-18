% This script will call smi.Publish to generate misc. results for an
% experiment on the sequential microscope.

%% Define the analysis parameters.
% Define the 'CoverslipDir'.
% NOTE: 'CoverslipDir' is the top-level directory which contains the
%        sub-directories 'CoverslipDir'\Cell*\Label*, which themselves
%        contain data in .h5 files 'CoverslipDir'\Cell*\Label*\Data*.h5
CoverslipDir = 'NEEDS_TO_BE_SET!';

%% Prepare the SMF structure.
% Load an existing SMF structure and modify from this starting point.
load('NEEDS_TO_BE_SET!');
SMF = smi_core.SingleMoleculeFitting.reloadSMF(SMF);
SMF.Data.PixelSize = 0.0954; % microns
SMF.Data.SEAdjust = 7 / (1000 * SMF.Data.PixelSize); % pixels

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

% smi.Publish contains several useful methods, however we'll almost always
% just call performFullAnalysis().
Publish.performFullAnalysis();
