% A real example used to develop SMLM.  See also smi.SMLM.unitTest

% If true, save plots produced into the ResultsDir defined below for a full
% analysis or a test fit with SMLMobj.Verbose < 5.
Saving = false;

% --- 2D ---

% Create an SMF (Single Molecule Fitting) structure with default values, some
% of which will be overridden below.  See smi_core.SingleMoleculeFitting for
% details on the possible parameters that can be set.
SMF = smi_core.SingleMoleculeFitting();

   % Results directory (char array)(Default='FileDir/Results')
SMF.Data.ResultsDir = tempdir;

if Saving
   if ~exist(SMF.Data.ResultsDir, 'dir')
      mkdir(SMF.Data.ResultsDir);
   end
end

% See
%    https://digitalrepository.unm.edu/physics_data/3/#attach_additional_files
%    (DOI: 10.25827/cs2a-dh13)
% for some example sequence files that can be used for testing.

% Set some Single Molecule Fitting (SMF) parameters:

   % File directory (char array)
SMF.Data.FileDir           = ...
   '\\rayleigh.phys.unm.edu\cell-path\Genmab\Data\10082020\Wien133_LQT_CD52_HexElect1\Cell_01\Label_01';
   % File name (cell array of char array)
SMF.Data.FileName         = 'Data_2020-10-8-17-58-54.h5';
   % ID tagged onto saved results (char array)(Default='')
SMF.Data.AnalysisID       = 'ID';
   % 'EMCCD','SCMOS' (Default='EMCCD')
SMF.Data.CameraType       = 'EMCCD';
   % Camera Gain, scalar or image (Default=1)
SMF.Data.CameraGain       = 1;
   % Camera Offset, scalar or image (Default=0)
SMF.Data.CameraOffset     = 0;
   % Perform thresholding? (Default=true)
SMF.Thresholding.On       = true;
   % Maximum allowed precision in x,y (Pixels)(Default=.2)
SMF.Thresholding.MaxXY_SE = 0.1;
   % Perform frame connection? (Default=true)
SMF.FrameConnection.On    = true;
   % Perform drift correction? (Default=true)
SMF.DriftCorrection.On    = true;

% Create an SMLM object using the values in the SMF structure.
SMLMobj = smi.SMLM(SMF);
% Do a test fit, displaying all the results to the screen (if Verbose >= 5,
% otherwise saving the results in ResultsDir/TestFit).
SMLMobj.Verbose = 5;
SMLMobj.testFit(1);
% Do a full analysis, saving results in ResultsDir.
%SMLMobj.fullAnalysis();
