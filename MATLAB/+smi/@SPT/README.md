### +smi/@SPT

SPT contains methods useful for single-particle tracking analysis.
  This class contains a collection of analysis/visualization methods
  useful for the analysis of single-particle tracking data.

REQUIRES:
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

```
properties:
    % Structure of parameters (see smi_core.SingleMoleculeFitting)
    SMF
    
    % Directory containing channel reg. transforms (char array)
    % NOTE: This property is only used in obj.batchTrack().
    TransformDir = '';
    
    % Pattern to match for transform files in TransformDir.
    % (see obj.batchTrack() for usage)
    % NOTE: This is used in the MATLAB built-in method dir(), which
    %       allows for a wildcard (*) in the name, but not a full
    %       regexp.
    TransformPattern = 'RegistrationTransform*.mat';
    
    % Indicate we should search for files to track. (Default = true)
    % NOTES: This is only used in obj.batchTrack().  If true,
    %        obj.batchTrack() will search for the files matching
    %        obj.FilePattern in obj.SMF.Data.FileDir.  If false, the
    %        files in obj.SMF.Data.FileName will be used instead.
    FindFiles = true;
    
    % Pattern to match for file names in obj.SMF.FileDir.
    % (see obj.batchTrack() for usage)
    % NOTE: This is used in the MATLAB built-in method dir(), which
    %       allows for a wildcard (*) in the name, but not a full
    %       regexp.
    FilePattern = '*.mat';
    
    % Diffusion estimator class for when UseTrackByTrackD is set.
    % NOTE: This is used here so that the user can change properties of
    %       the DiffusionEstimator class as needed when using
    %       obj.SMF.Tracking.TrajwiseD
    DiffusionEstimator
    
    % Structure of parameters used when generating movies.
    % NOTE: This only gets used when calling obj.saveResults() when
    %       obj.GenerateMovies = true.
    MovieParams
    
    % Marker to ignore entries in cost matrices (Default = -1)
    % NonlinkMarker can't be inf or NaN.
    NonLinkMarker = -1;
    
    % Flag to indicate movies should be made (Default = true)
    GenerateMovies = true;
    
    % Flag to indicate plots should be made (Default = true)
    GeneratePlots = true;
    
    % Flag to indicate test run (Default = false)
    IsTestRun = false;
    
    % Flag to make some outputs in physical units (Default = false)
    UnitFlag = false;
    
    % Flag to indicate sparse matrix usage (Default = true)
    % For now, this only applies to gap-closing.
    UseSparseMatrices = true;
    
    % Verbosity of the main analysis workflow. (Default = 1)
    Verbose = 1;
```
