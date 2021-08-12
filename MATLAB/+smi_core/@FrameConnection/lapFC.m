function [SMD, InternalParams] = lapFC(SMD, SMF, Verbose, InternalParams)
%lapFC connects localizations in 'SMD' by solving a LAP.
% This method solves the frame-connection problem by formulating it as a
% linear assignment problem, where the costs of connection/no connection
% are influenced by kinetic rates and local densities.
% NOTE: This method will add an additional field to SMD called
%       "ConnectID".  SMD.ConnectID is an integer array indicating
%       which localizations were connected during the frame connection
%       process.  For example, if
%       (SMD.ConnectID(nn) == SMD.ConnectID(mm)), the localizations
%       in SMD identified by the indices nn and mm were connected during
%       frame connection.  The exact value of the field "ConnectID" is
%       itself arbitrary and carries no meaning further than associating
%       localizations. This field is directly related to
%       SMDCombined.ConnectID as follows:
%           For a given ConnectID, say nn, the indices in arrays of SMD
%           that were combined to generate a field in SMDCombined can be
%           found as IndicesSMD = find(SMD.ConnectID == nn) (alternatively,
%           IndicesSMD = smi_core.FrameConnection.findConnected(...
%               SMDCombined, SMD, nn) )
%
% INPUTS:
%   SMD: SingleMoleculeData structure with the localizations that we wish
%        to frame-connect.
%   SMF: SingleMoleculeFitting structure defining relevant parameters.
%   Verbose: Integer specifying the verbosity level. (Default = 1)
%   InternalParams: Structure of parameters which are used in this method
%                   but are not present in SMF.  Typically, this input
%                   should not be provided (can be entered as []), however
%                   it is allowed as an input for testing purposes (e.g.,
%                   testing the algorithm performance with perfect inputs).
%                   (Default = [] so that these parameters are estimated
%                   internally from the data in SMD).
%
% OUTPUTS:
%   SMD: SMD but with the field 'ConnectID' populated.
%   InternalParams: Structure of parameters which are computed internally
%                   in this method (e.g., the rates KOn, KOff, and
%                   KBleach).  The field 'NEmitters' is provided for
%                   convenience (unless 'InternalParams' is provided as an
%                   input) but is not actually used anywhere in the
%                   algorithm.
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('SMF', 'var') || isempty(SMF))
    SMF = smi_core.SingleMoleculeFitting;
end
if (~exist('InternalParams', 'var') || isempty(InternalParams))
    InternalParams = struct();
    EstimateParams = true;
else
    % This input was provided, so we'll extract some parameters for later
    % use.
    KOn = InternalParams.KOn;
    KOff = InternalParams.KOff;
    KBleach = InternalParams.KBleach;
    PMiss = InternalParams.PMiss; 
    InitialDensity = InternalParams.InitialDensity;
    EstimateParams = false;
end
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 1;
end    

% Perform an initial pre-clustering of the data.
if (Verbose > 2)
    fprintf('\tFrameConnection.lapFC(): Pre-clustering localizations...\n')
end
SMDPreClustered = smi_core.FrameConnection.preClusterCoords(SMD, SMF);
SMD.ConnectID = SMDPreClustered.ConnectID;

% Recursively perform the frame-connection, using updated density and
% fluorophore transition rate estimates from the previous iteration.
for rr = 1:SMF.FrameConnection.NRecursions
    if (Verbose > 2)
        fprintf(['\tFrameConnection.lapFC(): ', ...
            'LAP-FC recursion %i of %i...\n'], ...
            rr, SMF.FrameConnection.NRecursions)
    end
    % Estimate the number of emitters present at the start of the
    % experiment.  For now, I'll just use the total number of clusters,
    % even though this breaks several assumptions made in the rate
    % parameter estimation step.
    if (Verbose > 1)
        fprintf(['\tFrameConnection.lapFC(): ', ...
            'Re-organizing precluster data...\n'])
    end
    ClusterData = smi_core.FrameConnection.organizeClusterData(SMD);
    
    % Estimate some needed parameters from the data.
    if EstimateParams
        % Estimate the fluorophore transition rate parameters.
        if (Verbose > 1)
            fprintf(['\tFrameConnection.lapFC(): ', ...
                'Estimating rate parameters from preclusters...\n'])
        end
        [KOn, KOff, KBleach, PMiss, NEmitters] = ...
            smi_core.FrameConnection.estimateRateParameters(...
            ClusterData, Verbose);
        InternalParams.NEmitters = NEmitters;
       
        % Estimate the local density around each cluster.
        if (Verbose > 1)
            fprintf(['\tFrameConnection.lapFC(): ', ...
                'Estimating local emitter densities...\n'])
        end
        InitialDensity = smi_core.FrameConnection.estimateLocalDensity(...
            ClusterData, SMF.FrameConnection.NNearestClusters, ...
            KOn, KOff, KBleach, PMiss);
    end

    % Grab a few arrays out of obj and SMD (to improve speed within the
    % for loop below).
    ConnectID = SMD.ConnectID;
    MaxFrameGap = SMF.FrameConnection.MaxFrameGap;
    NFramesMax = SMD.NFrames;
    if isempty(NFramesMax)
        NFramesMax = max(SMD.FrameNum);
    end

    % Loop through all clusters, construct and minimize the cost matrix, and
    % then cluster the localizations.
    if (Verbose > 1)
        fprintf(['\tFrameConnection.lapFC(): ', ...
            'Looping over clusters and solving the LAP...\n'])
    end
    UniqueIDs = unique(ConnectID, 'sorted');
    NIDsInitial = numel(UniqueIDs);
    MaxConnectID = NIDsInitial;
    ClustersToAnalyze = find(cellfun(@(X) size(X, 1) > 1, ClusterData));
    for nn = ClustersToAnalyze.'
        % Construct the cost matrix.
        CostMatrix = smi_core.FrameConnection.createCostMatrix(...
            ClusterData{nn}, ...
            KOn, KOff, KBleach, PMiss, InitialDensity(nn), ...
            MaxFrameGap, NFramesMax, -1);
        
        % Solve the linear assignment problem.
        Link12 = smi.SPT.solveLAP(CostMatrix, -1);
        
        % Update the connect IDs to indicate the linkage.
        [ConnectID, MaxConnectID] = ...
            smi_core.FrameConnection.linkClusters(...
            ConnectID, MaxConnectID, ...
            ClusterData{nn}(:, 8), Link12);
    end
    SMD.ConnectID = smi.SPT.validifyConnectID(ConnectID);
end
InternalParams.KOn = KOn;
InternalParams.KOff = KOff;
InternalParams.KBleach = KBleach;
InternalParams.PMiss = PMiss;
InternalParams.InitialDensity = InitialDensity;


end