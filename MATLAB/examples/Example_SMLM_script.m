% A real example used to develop SMLM.  See also smi.SMLM.unitTest

% If true, save plots produced into the ResultsDir defined below for a full
% analysis or a test fit with SMLMobj.VerboseTest < 5.
Saving = false;

% --- 2D ---

% Create an SMF (Single Molecule Fitting) structure with default values, some
% of which will be overridden below.  See smi_core.SingleMoleculeFitting for
% details on the possible parameters that can be set.
SMF = smi_core.SingleMoleculeFitting();

   % Results directory (char array)(Default='FileDir/Results')
SMF.Data.ResultsDir = tempdir;

   % Also, SMF.Data.FileDir and SMF.Data.Filename below need to be set with
   % the input data location.

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
SMF.Data.FileDir           = 'NEEDS_TO_BE_SET!';
   % File name (cell array of char array)
SMF.Data.FileName         = 'NEEDS_TO_BE_SET!';
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
SMLMobj.Verbose = 1;
% Do a test fit, displaying all the results to the screen (if VerboseTest >= 5,
% otherwise saving the results in ResultsDir/TestFit).  Note that calling
% testFit from the smi.gui is equivalent to what we are doing here.
SMLMobj.VerboseTest = 5;
SMLMobj.testFit(1);
% Do a full analysis, saving results in ResultsDir.
SMLMobj.fullAnalysis();
