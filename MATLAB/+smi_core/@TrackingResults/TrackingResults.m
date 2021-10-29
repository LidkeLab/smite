classdef TrackingResults
    % TrackingResults A class defining the Tracking Results structure
    %
    % This datatype is one of the primary results structures in the smite
    % enviroment. The Tracking Results (TR) structure is an input/output of
    % many methods in smite which are related to single-particle tracking
    % (SPT) analysis.  TR structures are organized as follows: each
    % structure element corresponds to a single trajectory, i.e., TR(n)
    % contains all relevant properties of the n-th trajectory. The TR
    % structure is intended to carry all necessary information from an SMD
    % structure (see smi_core.SingleMoleculeData) but organized in a more
    % user-friendly manner for SPT data.
    %
    % The TR structure is just an array of SMD structures, with each array
    % element being an SMD structure corresponding to a single trajectory.
    %
    % SEE ALSO:
    %   smi_core.SingleMoleculeData, smi_core.SingleMoleculeFitting
    
    properties
        
    end
    
    methods (Static)
        function [TR] = createTR()
            %createTR creates an empty Tracking Results (TR) structure.
            TR = smi_core.SingleMoleculeData.createSMD();
        end
        
        [Durations] = computeTrajDurations(TR);
        [TR] = convertSMDToTR(SMD);
        [TRIndex] = getTRIndex(TR, TrajectoryIDs);
        [TR] = joinTraj(TR, TrajectoryIDs, Verbose);
        [TR] = threshTrajLength(TR, MinTrackLength);
        [TR] = padTR(TR, TRPadding);
        [TR] = catTR(TR1, TR2, CheckDims);
        
    end
end