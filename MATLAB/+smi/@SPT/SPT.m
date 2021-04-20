classdef SPT < handle
    % SPT contains methods useful for single-particle tracking analysis.
    %   This class contains a collection of analysis/visualization methods
    %   useful for the analysis of single-particle tracking data.
    
    
    properties
        % ROI of the raw data used to generate SMD. (Pixels)
        % Default will be set in obj.generateTrajectories() if needed.
        % Organized as [YStart, XStart, YEnd, XEnd]
        DataROI
        
        % Structure of parameters (see smi_core.SingleMoleculeFitting)
        SMF
        
        % Indicate SMF.Tracking.Rho_off can be overwritten (Default = true)
        % NOTE: As of this writing, this only makes an appearance in
        %       obj.generateTrajectories(). If you aren't using that
        %       method, the dark emitter density Rho_off will be taken from
        %       SMF.Tracking.Rho_off.
        EstimateRhoFromData = true;
        
        % Marker to ignore entries in cost matrices (Default = -1)
        % NonlinkMarker can't be inf or NaN. 
        NonlinkMarker = -1;
        
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
        
        % Pre gap-closing SMD structure (see smi_core.SingleMoleculeData)
        SMDPreGapClosing
        
        % Pre-threshold SMD structure (see smi_core.SingleMoleculeData)
        SMDPreThresh
       
        % Tracking Results structure (see smi_core.TrackingResults)
        TR
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
            
            % Store the SMF as a class property, ensuring that some vital
            % "flags" are turned off first.
            SMF.DriftCorrection.On = false;
            SMF.FrameConnection.On = false;
            obj.SMF = SMF;
            
            % Start the GUI if needed.
            if StartGUI
                obj.gui();
            end
            
        end
        
        performFullAnalysis(obj)
        [TR, SMD] = generateTrajectories(obj);
        saveResults(obj)
        gui(obj)
        
    end
    
    methods(Static)
        [Success] = unitTest();
        [Success] = unitTestFFGC()
        [CostMatrix] = createCostMatrixFF(SMD, SMF, FrameNumber, ...
            NonLinkMarker);
        [CostMatrix] = createCostMatrixGC(SMD, SMF, ...
            NonLinkMarker, CreateSparseMatrix);
        [Assign12, Cost12] = solveLAP(CostMatrix, NonlinkMarker);
        [SMD] = connectTrajFF(SMD, Link12, FrameNumber);
        [SMD] = connectTrajGC(SMD, Link12);
    end
    
    
end