classdef SPT < handle
    % SPT contains methods useful for single-particle tracking analysis.
    %   This class contains a collection of analysis/visualization methods
    %   useful for the analysis of single-particle tracking data.
    
    
    properties
        % Structure of parameters (see smi_core.SingleMoleculeFitting)
        SMF
        
        % Verbosity of the main analysis workflow. (Default = 1)
        Verbose = 1;
    end
    
    properties (Access = protected)
        % Single Molecule Data structure (see smi_core.SingleMoleculeData)
        % NOTE: SMD contains the localizations in SMDPreThresh for which
        %       ThreshFlag was equal to 0 (i.e., the localizations which
        %       we wish to keep)
        SMD
        
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
            %        (Default = smi_core.SingleMoleculeFitting)
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
            
            % Start the GUI if needed.
            if StartGUI
                obj.gui();
            end
            
        end
        
        performFullAnalysis(obj)
        saveResults(obj)
        gui(obj)
        
    end
    
    methods(Static)
        [Success] = unitTest();
        [Success] = unitTestFFGC()
        [CostMatrix] = createCostMatrixFF(SMD, SMF, FrameNumber, ...
            RhoOff, NonLinkMarker);
        [CostMatrix] = createCostMatrixGC(SMD, SMF, ...
            RhoOff, NonLinkMarker, CreateSparseMatrix);
        [Assign12, Cost12] = solveLAP(CostMatrix, NonlinkMarker);
        [SMD] = connectTrajFF(SMD, Link12, FrameNumber);
        [SMD] = connectTrajGC(SMD, Link12);
    end
    
    
end