classdef SPT < handle
    % SPT contains methods useful for single-particle tracking analysis.
    %   This class contains a collection of analysis/visualization methods
    %   useful for the analysis of single-particle tracking data.
    
    
    properties
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
    end
    
    properties (SetAccess = protected)
        % Single Molecule Data structure (see smi_core.SingleMoleculeData)
        % NOTE: SMD contains the localizations in SMDPreThresh for which
        %       ThreshFlag was equal to 0 (i.e., the localizations which
        %       we wish to keep)
        SMD
        
        % Pre-threshold SMD structure (see smi_core.SingleMoleculeData)
        SMDPreThresh
        
        % Tracking Results structure (see smi_core.TrackingResults)
        TR
        
        % Tracking Results structure before channel registration.
        TRPreCR = struct([]);
        
        % 'SMD' before channel registration is applied.
        SMDPreCR = struct([]);
        
        % Pre-threshold pre-channel registration 'SMD'.
        SMDPreThreshPreCR = struct([]);
        
        % Scaled raw data.
        ScaledData;
    end
    
    properties (Hidden)
        % Diffusion coefficients for each localization in trajectories.
        % This array is organized as a two-column array as [D, D_SE]
        % NOTE: This is only used when UseTrackByTrackD is set to true.
        %       I've made this hidden because the user shouldn't really be
        %       using these values to do anything.  If they're needed, the
        %       user should produce them in the diffusion estimator class,
        %       or access them in the appropriate properties of
        %       obj.DiffusionEstimator.
        DiffusionCoefficients = [];
        
        % Copy of the SMF structure.
        % This is used for a few random tests/things like
        % obj.TryLowPValueLocs which, when enabled, requires us to modify
        % the SMF provided by the user.  We will want to revert to the
        % user's original SMF if that's done.
        SMFCopy
        
        % Regular expression for filename timestamps.
        % NOTE: This is used in obj.batchTrack()
        TimeStampRegExp = '\d{4,}-\d{1,2}-\d{1,2}-\d{1,2}-\d{1,2}-\d{1,2}';
        
        % Delimiter in the timestamp regular expression.
        % NOTE: This is used in obj.batchTrack()
        TimeStampDelimiter = '-';
    end
    
    methods
        
        function obj = SPT(SMF, StartGUI)
            %SPT is the class constructor for SPT.
            %
            % INPUTS:
            %   SMF: Single Molecule Fitting structure.
            %        (Default = smi_core.SingleMoleculeFitting with some
            %        minor changes to parameters).
            %   StartGUI: Flag to decide if GUI should be opened
            %             automatically. (Default = true)
            
            % Set defaults if needed.
            if (~exist('SMF', 'var') || isempty(SMF))
                SMF = smi_core.SingleMoleculeFitting;
            end
            if (~exist('StartGUI', 'var') || isempty(StartGUI))
                StartGUI = true;
            end
            
            % Store the SMF as a class property.
            % NOTE: There's some strangeness related to property
            %       instantiation in the SMF which causes our input value
            %       of ResultsDir to be ovewritten by a default.  Until I
            %       fix this, I'm adding the (hopefully) temporary fix
            %       below.
            obj.SMF = SMF;
            obj.SMF.Data.ResultsDir = SMF.Data.ResultsDir;
            
            % Create an instance of the diffusion estimator class.
            obj.DiffusionEstimator = smi_stat.DiffusionEstimator;
            obj.DiffusionEstimator.FitIndividualTrajectories = true;
            obj.DiffusionEstimator.UnitFlag = false;
            
            % Set some default movie parameters.
            obj.MovieParams = smi_vis.GenerateMovies.prepDefaults();
            
            % Start the GUI if needed.
            if StartGUI
                obj.gui();
            end
            
        end
        
        function set.SMF(obj, SMFInput)
            %set method for the property SMF.
            % We want to ensure some fields of the SMF are always turned
            % off for tracking. Also, we want to just make a copy of the
            % SMF instead of keeping the original reference to an SMF
            % class instance.
            SMFInput.DriftCorrection.On = false;
            SMFInput.FrameConnection.On = false;
            obj.SMF = smi_core.SingleMoleculeFitting.reloadSMF(SMFInput);
        end
        
        [TR, SMD, SMDPreThresh] = performFullAnalysis(obj);
        [TR, SMD, SMDPreThresh, FileList, TransformList] = batchTrack(obj);
        autoTrack(obj)
        generateTrajectories(obj)
        updateTrackingParams(obj, SMD, TR)
        saveResults(obj)
        gui(obj)
        
    end
    
    methods(Static)
        [Success] = unitTestFFGC()
        [SMD] = genTrajFF(SMD, SMF, DiffusionCoefficients, NonLinkMarker);
        [SMD] = genTrajGC(SMD, SMF, DiffusionCoefficients, ...
            NonLinkMarker, UseSparseMatrices);
        [CostMatrix] = createCostMatrixFF(SMD, SMF, ...
            DiffusionCoefficients, FrameNumber, NonLinkMarker);
        [CostMatrix, StartEndIndices] = createCostMatrixGC(SMD, SMF, ...
            DiffusionCoefficients, NonLinkMarker, CreateSparseMatrix);
        [Assign12, Cost12] = solveLAP(CostMatrix, NonlinkMarker);
        [SMD] = connectTrajFF(SMD, Link12, FrameNumber);
        [SMD] = connectTrajGC(SMD, Link12);
        [KOn, KOff] = estimateRateParameters(SMD);
        [RhoOff, RhoOn] = estimateDensities(SMD, SMF);
        [DiffusionStruct] = ...
            estimateDiffCoeffs(TR, DiffusionEstimator, DReset)
    end
    
    
end