classdef SimSPT < handle
    %SimSPT contains methods useful for simulating SPT data.
    % SimSPT is a collection of methods that can be used to simulate
    % single-particle tracking (SPT) data, including realistic raw data and
    % ground-truth trajectories.
    
    % Created by:
    %   David J. Schodt (Lidke Lab, 2021) with class structure modeled on
    %       smi_sim.SimSMLM
    
    
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
        
        % obj.SMDLabeled after applying photokinetics model.
        SMDModel
        
        % obj.SMDModel after frame averaging and noisy measurement.
        SMD
    end
    
    methods
        function obj = SimSPT(SimParams)
            %SimSPT class constructor.
            
            % Set defaults.
            if (~exist('SimParams', 'var') || isempty(SimParams))
                obj.SimParams = smi_sim.SimSPT.defineDefaultParams();
            else
                obj.SimParams = SimParams;
            end
        end
        
        function set.SimParams(obj, SimParams)
            % This method ensures obj.SimParams has all parameters.
            obj.SimParams = smi_helpers.padStruct(SimParams, ...
                smi_sim.SimSPT.defineDefaultParams());
        end
        
        [SMD, SMDModel, SMDLabeled, SMDTrue, OligomerStruct] = ...
            createSimulation(obj);
        
    end
    
    methods (Static)
        [SimParams] = defineDefaultParams();
        [Coordinates, MaskedCoordinates] = ...
            applyCoordMask(Coordinates, Mask, FrameSize);
        [TrajectoryStruct] = simTrajectories(SimParams);
        [TrajectoryStruct] = applyLabelingEfficiency(TrajectoryStruct, ...
            SimParams);
        [TrajectoryStruct] = simEmitterKinetics(TrajectoryStruct, ...
            SimParams);
        [SMD] = applyMeasurementModel(TrajectoryStruct, SimParams);
    end
    
    methods (Static, Hidden)
        [TrajectoryStruct] = simTrajBrownian(InitialPositions, SimParams);
        [TrajectoryStruct] = simOligoTrajBrownian(InitialPositions, ...
            SimParams);
        [Trajectories] = smi_sim.SimSPT.enforcePeriodicBoundary(...
            Trajectories, PeriodicityMapT);
    end
    
    
end