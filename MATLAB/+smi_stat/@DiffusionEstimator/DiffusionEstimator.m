classdef DiffusionEstimator < handle
    %DiffusionEstimator contains methods useful for diffusion estimation.
    % This class contains several methods used to estimate diffusion
    % constants from single-particle tracking data.
    
    properties        
        % Fit method for fitting MSD results (char array/string)
        FitMethod = 'weightedLS';
        
        % Max. frame lag of the MSD (scalar, integer)
        MaxFrameLag = inf;
        
        % Directory in which results will be saved by saveResults().
        SaveDir = pwd();
        
        % Name of results file. Default defined in obj.saveDir().
        SaveName
        
        % Tracking results structure.
        TR
        
        % Boolean flag to indicate units of outputs (boolean)(Default = 1)
        % 1 (true) will make the outputs of estimateDiffusionConstant()
        %   micrometers and seconds.
        % 0 (false) will make the outputs of estimateDiffusionConstant()
        %   pixels and frames.
        % NOTE: Most methods of this class will use pixels and frames
        %       regardless of obj.UnitFlag. This property will only affect
        %       the "user-facing" wrapper methods, such as
        %       estimateDiffusionConstant() and saveResults()
        UnitFlag = true;
        
        % Verbosity level for estimateDiffusionConstant() (Default = 0)
        %   Verbose 0: no Command Window outputs
        %   Verbose 1: General progress printed to Command Window
        %   Verbose 2: Display some intermediate results
        %   Verbose 3: Debugging mode, extensive outputs
        Verbose = 0;
    end
    
    properties (SetAccess = protected)
        % Structure array containing diffusion estimates.
        DiffusionStruct
        
        % Structure array containing trajectory-wise MSDs.
        MSDSingleTraj
        
        % Structure array containing the ensemble MSD.
        MSDEnsemble
    end
    
    methods
        function [obj, DiffusionStruct] = ...
                DiffusionEstimator(TR, MaxFrameLag, FitMethod, ...
                Verbose, AutoRun)
            %DiffusionEstimator is the class constructor.
            % Several optional inputs can be provided to directly set class
            % properties.  'AutoRun' is a boolean flag which specifies 
            % whether or not this constructor should call
            % obj.estimateDiffusionConstant() if all requisite class
            % properties were provided.
            
            % Set defaults if needed.
            if (~exist('AutoRun', 'var') || isempty(AutoRun))
                AutoRun = false;
            end
            
            % Set class properties based on the inputs.
            AllFieldsSet = true;
            if (exist('TR', 'var') && ~isempty(TR))
                obj.TR = TR;
            else
                AllFieldsSet = false;
            end
            if (exist('MaxFrameLag', 'var') && ~isempty(MaxFrameLag))
                obj.MaxFrameLag = MaxFrameLag;
            end
            if (exist('FitMethod', 'var') && ~isempty(FitMethod))
                obj.FitMethod = FitMethod;
            end
            if (exist('Verbose', 'var') && ~isempty(Verbose))
                obj.Verbose = Verbose;
            end
            
            % Run obj.estimateDiffusionConstant() if requested.
            if (AutoRun && AllFieldsSet)
                [DiffusionStruct] = obj.estimateDiffusionConstant();
            end
        end
        
        [DiffusionStruct] = estimateDiffusionConstant(obj, ...
            SaveFlag);
        saveResults(obj)
        
    end
    
    methods (Static)
        [FitParams, FitParamsSE] = fitMSD(MSDStruct, FitMethod, Verbose);
        [MSDSingleTraj, MSDEnsemble] = ...
            computeMSD(TR, MaxFrameLag, Verbose);
    end
    
    methods (Static, Hidden)
        % These methods are 'Hidden' because they are primarily intended
        % for use within other more user-centric methods (i.e., we don't
        % want to distract the user with these options, but if they need
        % them they are still accessible).
        
        [MSDSingleTraj] = computeSingleTrajMSD(TR, Verbose);
        
    end
    
    
end