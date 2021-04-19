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
    % The structure contains (at least) the following properties:
    %   TrajectoryID: An integer identifier unique to each trajectory.
    %   X: X position of the trajectory over time. (pixels)
    %   Y: Y position of the trajectory over time. (pixels)
    %   X_SE: The standard error of the X position estimate. (pixels)
    %   Y_SE: The standard error of the Y position estimate. (pixels)
    %   FrameNum: The frame in which localizations of the trajectory
    %             appear. (frames)
    %   Photons: The integrated photons found in each localization of
    %            the trajectory. (photons)
    %   Bg: The background found for each localization of the trajectory.
    %       (photons)
    %   IndTD: An index back into the associated TD structure, i.e.,
    %          TR(nn).X(ii) can be found in TD.X(TR(nn).IndTD(ii)).
    %   FrameRate: The frame rate used to acquire this data.
    %              (frames / second)
    %   PixelSize: The pixel size of the sensor used to acquire this data
    %              BACK PROJECTED to the sample plane. (micrometers)
    %   DataROI: ROI of original raw data in which the trajectories were
    %            found, formatted as [YStart, XStart, YEnd, XEnd]. (pixels)
    %   XSize: The number of pixels along the X dimension of the ROI.
    %          (pixels)
    %   YSize: The number of pixels along the Y dimensino of the ROI.
    %          (pixels)
    %   FileDir: Directory containing the raw data corresponding to this TR
    %   FileName: Name of the file corresponding to this TR.
    %
    % SEE ALSO:
    %   smi_core.SingleMoleculeData, smi_core.SingleMoleculeFitting
    
    properties
        
    end
    
    methods (Static)
        function [TR] = createTR()
            %createTR creates an empty Tracking Results (TR) structure.
            TR.TrajectoryID = [];
            TR.X = [];
            TR.Y = [];
            TR.X_SE = [];
            TR.Y_SE = [];
            TR.FrameNum = [];
            TR.Photons = [];
            TR.Bg = [];
            TR.IndTD = [];
            TR.FrameRate = [];
            TR.PixelSize = [];
            TR.DataROI = [];
            TR.XSize = [];
            TR.YSize = [];
            TR.FileDir = '';
            TR.FileName = '';
        end
        
        [TR] = convertSMDToTR(SMD, FileInfoStruct);
        [TRIndex] = getTRIndex(TR, TrajectoryIDs);
        [TR] = joinTraj(TR, TrajectoryIDs, Verbose);
        
    end
end