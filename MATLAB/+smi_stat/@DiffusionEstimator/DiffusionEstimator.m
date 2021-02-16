classdef DiffusionEstimator
    %DiffusionEstimator contains methods useful for diffusion estimation.
    % This class contains several methods used to estimate diffusion
    % constants from single-particle tracking data.
    
    properties
        % Tracking results structure.
        TR
    end
    
    methods
        function obj = DiffusionEstimator()
            %DiffusionEstimator is the class constructor
        end
    end
    
    methods (Static)
        [DiffusionConstant, DiffusionConstantSE] = fitMSD(TR);
        [MSDStruct, TR] = computeMSD(TR, MaxLag);
    end
    
    methods (Static, Hidden)
        % These methods are 'Hidden' because they are primarily intended
        % for use within other more user-centric methods (i.e., we don't
        % want to distract the user with these options, but if they need
        % them they are still accessible).
        
        [MSD, NCount, SquaredDisplacement] = computeSingleTrajMSD(TR);
        [MSD] = computeEnsembleMSD(TR, MaxLag);
    end
    
    
end