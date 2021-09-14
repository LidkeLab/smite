classdef SimSPT < handle
    %SimSPT contains methods useful for simulating SPT data.
    % SimSPT is a collection of methods that can be used to simulate
    % single-particle tracking (SPT) data, including realistic raw data and
    % ground-truth trajectories.
    
    % Created by:
    %   David J. Schodt (Lidke Lab, 2021)
    
    
    properties
        % Structure of parameters used in the simulation.
        SimParams = struct();
    end
    
    properties (SetAccess = 'protected')
        % True locations of the underlying diffusing targets.
        % NOTE: This structure has fields in units of subframes!
        SMDTrue
        
        % obj.SMDTrue after applying labeling efficiencies.
        % NOTE: This structure has fields in units of subframes!
        SMDLabeled
        
        % obj.SMDLabeled after applying blinking/missed locs./frame avg.ing
        SMDModel
        
        % obj.SMDModel after noisy measurement.
        SMD
    end
    
    methods
        function obj = SimSPT(SimParams)
            %SimSPT class constructor.
            
            % Set defaults.
            if (~exist('SimParams', 'var') || isempty(SimParams))
                obj.SimParams = smi_sim.SimSPT.defineDefaultParams();
            else
                obj.SimParams = smi_helpers.padStruct(SimParams, ...
                    smi_sim.SimSPT.defineDefaultParams());
            end
        end
        
        [SMD, SMDModel, SMDLabeled, SMDTrue, OligomerStruct] = ...
            simulateTrajectories(obj);
        
    end
    
    methods (Static)
        [SimParams] = defineDefaultParams();
    end
    
    
end