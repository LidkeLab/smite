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
        % obj.TrajStructModel after applying a measurement noise model.
        TrajStruct
        
        % obj.TrajStructSubModel after motion blurring, before noising.
        TrajStructModel
        
        % obj.TrajStructSubLabeled after applying photokinetics model.
        % NOTE: This structure has fields in units of subframes!
        TrajStructSubModel
        
        % obj.TrajStructSubTrue after applying labeling efficiencies.
        % NOTE: This structure has fields in units of subframes!
        TrajStructSubLabeled
        
        % True locations of the underlying diffusing targets.
        % NOTE: This structure has fields in units of subframes!
        TrajStructSubTrue
    end
    
    properties (Dependent)
        % obj.TrajStruct converted to the more useable SMD format.
        SMD
        
        % SMD converted to a Tracking Results structure.
        TR
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
        
        function [SMD] = get.SMD(obj)
            % This method converts obj.TrajStruct to an SMD.
            SMD = obj.convertTrajToSMD(obj.TrajStruct, obj.SimParams);
        end
        
        function [TR] = get.TR(obj)
            % This method converts obj.SMD to a TR.
            TR = smi_core.TrackingResults.convertSMDToTR(obj.SMD);
        end
        
        createSimulation(obj);
        
    end
    
    methods (Static)
        [TrajStruct] = simTrajectories(SimParams);
        [TrajStruct, KeepInd] = applyLabelingEfficiency(...
            TrajStruct, LabelingEfficiency);
        [TrajStruct] = simEmitterKinetics(TrajStruct, SimParams);
        [TrajStruct, TrajStructModel] = applyMeasurementModel(...
            TrajStruct, SimParams);
        [SimParams] = defineDefaultParams();
        [Coordinates, MaskedCoordinates] = ...
            applyCoordMask(Coordinates, Mask, FrameSize);
        [SMD] = convertTrajToSMD(TrajStruct, SimParams);
        [Files, SimParams, DataParams] = makeExampleSim(...
            SimParams, DataParams, SaveDir);
    end
    
    methods (Static, Hidden)
        [TrajStruct] = simTrajBrownian(InitialPositions, SimParams);
        [Trajectories, PeriodicityMap, ConnectionMap] = ...
            applyBoundaryCondition(...
            Trajectories, BoundaryCondition, FrameSize, ConnectionMap);
        [Trajectories, ConnectionMapT, IsOn, TrajMap] = ...
            enforcePeriodicBoundary(...
            Trajectories, PeriodicityMapT, ConnectionMapT);
        [Trajectories, ConnectionMap] = simOligomers(...
            Trajectories, TrajectoryUpdates, NTraj, ...
            ConnectionMap, InteractionDistance, InteractionProb, ...
            KDisconnect, RestrictToDimers);
    end
    
    
end