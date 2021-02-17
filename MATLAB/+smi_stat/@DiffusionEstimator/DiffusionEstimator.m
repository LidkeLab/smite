classdef DiffusionEstimator
    %DiffusionEstimator contains methods useful for diffusion estimation.
    % This class contains several methods used to estimate diffusion
    % constants from single-particle tracking data.
    
    properties
        % Tracking results structure.
        TR
        
        % Max. frame lag of the MSD (scalar, integer)
        MaxFrameLag
        
        % Fit method for fitting MSD results (char array/string)
        FitMethod = 'weightedLS';
    end
    
    properties (SetAccess = protected)
        % Structure array containing trajectory-wise MSDs.
        MSDSingleTraj
        
        % Structure array containing the ensemble MSD.
        MSDEnsemble
    end
    
    methods
        function [obj, DiffusionStruct] = ...
                DiffusionEstimator(TR, MaxFrameLag, AutoRun)
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
            else
                AllFieldsSet = false;
            end
            if (exist('FitMethod', 'var') && ~isempty(FitMethod))
                obj.FitMethod = FitMethod;
            end
            
            % Run obj.estimateDiffusionConstant() if requested.
            if (AutoRun && AllFieldsSet)
                [DiffusionStruct] = obj.estimateDiffusionConstant();
            end
        end
        
        [DiffusionConstant, DiffusionConstantSE] = ...
            estimateDiffusionConstant(obj);
        
    end
    
    methods (Static)
        [FitParams, FitParamsSE] = fitMSD(MSDStruct, Method);
        [MSDSingleTraj, MSDEnsemble] = computeMSD(TR, MaxFrameLag)
    end
    
    methods (Static, Hidden)
        % These methods are 'Hidden' because they are primarily intended
        % for use within other more user-centric methods (i.e., we don't
        % want to distract the user with these options, but if they need
        % them they are still accessible).
        
        [MSDSingleTraj] = computeSingleTrajMSD(TR);
        
    end
    
    
end