classdef SPT < handle
    % SPT contains methods useful for single-particle tracking analysis.
    %   This class contains a collection of analysis/visualization methods
    %   useful for the analysis of single-particle tracking data.
    
    
    properties
        % Structure of parameters (see smi_core.SingleMoleculeFitting)
        SMF
        
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
        NonlinkMarker = -1;
        
        % Flag to indicate movies should be made (Default = true)
        GenerateMovies = true;
        
        % Flag to indicate plots should be made (Default = true)
        GeneratePlots = true;
        
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
        
        % History of parameters used with iterative tracking.
        % (see obj.autoTrack() for usage and organization)
        ParamsHistory = cell(0);
    end
    
    properties (Hidden)
        % Diffusion constants for each localization in trajectories.
        % This array is organized as a two-column array as [D, D_SE]
        % NOTE: This is only used when UseTrackByTrackD is set to true.
        %       I've made this hidden because the user shouldn't really be
        %       using these values to do anything.  If they're needed, the
        %       user should produce them in the diffusion estimator class,
        %       or access them in the appropriate properties of
        %       obj.DiffusionEstimator.
        DiffusionConstant = [];
        
        % Copy of the SMF structure.
        % This is used for a few random tests/things like 
        % obj.TryLowPValueLocs which, when enabled, requires us to modify
        % the SMF.  We will need to revert to the users original SMF if
        % that's done.
        SMFCopy
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
            obj.SMF = SMF;
            
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
        
        [TR, SMD] = performFullAnalysis(obj);
        autoTrack(obj)
        generateTrajectories(obj)
        saveResults(obj)
        gui(obj)
        
    end
    
    methods(Static)
        [Success] = unitTest();
        [Success] = unitTestFFGC()
        [CostMatrix] = createCostMatrixFF(SMD, SMF, ...
            DiffusionConstants, FrameNumber, NonLinkMarker);
        [CostMatrix, StartEndIndices] = createCostMatrixGC(SMD, SMF, ...
            DiffusionConstants, NonLinkMarker, CreateSparseMatrix);
        [Assign12, Cost12] = solveLAP(CostMatrix, NonlinkMarker);
        [SMD] = connectTrajFF(SMD, Link12, FrameNumber);
        [SMD] = connectTrajGC(SMD, Link12);
        [ConnectID] = validifyConnectID(ConnectID);
        [KOn, KOff] = estimateRateParameters(SMD);
    end
    
    
end