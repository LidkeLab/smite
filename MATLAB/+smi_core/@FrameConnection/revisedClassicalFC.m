function [SMD] = revisedClassicalFC(SMD, SMF, Verbose)
%revisedClassicalFC connects localizations in 'SMD' by simple thresholds.
% This method solves the frame-connection problem by connecting
% localizations within hard spatiotemporal thresholds.  The "revision" with
% respect to classicalFC() is that the spatial thresholds are defined in
% terms of the position standard errors of the localizations.
%
% NOTE: This method is nearly identical to
%       FrameConnection.preClusterCoords(), with the only difference being
%       that localizations within the same frame are prohibited from being
%       connected to one another.
%
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
%
% OUTPUTS:
%   SMD: SMD but with the field 'ConnectID' populated.
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('SMF', 'var') || isempty(SMF))
    SMF = smi_core.SingleMoleculeFitting;
end
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 1;
end    

% Gather/revise/reorganize some arrays for further use.
NLocalizations = numel(SMD.FrameNum);
[DatasetNum, SortIndices] = sort(SMD.DatasetNum);
X = SMD.X(SortIndices);
Y = SMD.Y(SortIndices);
X_SE = SMD.X_SE(SortIndices);
Y_SE = SMD.Y_SE(SortIndices);
FrameNum = SMD.FrameNum(SortIndices);
MeanXYSE = mean([X_SE, Y_SE], 2);
MaxFrameGap = SMF.FrameConnection.MaxFrameGap;
NSigmaDev = SMF.FrameConnection.NSigmaDev;

% Initialize each localization as a new cluster.
ConnectID = (1:NLocalizations).';

% Loop through datasets and perform the pre-clustering.
[NLocPerDataset, DatasetArray] = groupcounts(DatasetNum);
CumulativeDatasetLocs = [0; cumsum(NLocPerDataset)];
MaxID = NLocalizations;
for ii = 1:numel(DatasetArray)
    % Provide a Command Window update if needed.
    if (Verbose > 2)
        fprintf(['\tFrameConnection.revisedClassicalFC(): ', ...
            'Performing frame connection for dataset %i...\n'], ii)
    end
    
    % Isolate some arrays for the current dataset (CDS = current dataset)
    CurrentDatasetInd = (1:NLocPerDataset(ii)) + CumulativeDatasetLocs(ii);
    [FrameNumCDs, SortIndicesFN] = sort(FrameNum(CurrentDatasetInd));
    XCDs = X(CurrentDatasetInd(SortIndicesFN));
    YCDs = Y(CurrentDatasetInd(SortIndicesFN));
    MeanXYSECDs = MeanXYSE(CurrentDatasetInd(SortIndicesFN));
    ConnectIDCDs = ConnectID(CurrentDatasetInd(SortIndicesFN));
    
    % Loop through frames and add localizations to clusters.
    IsClustered = zeros(NLocPerDataset(ii), 1, 'logical');
    [NLocPerFrame, FrameArray] = groupcounts(FrameNumCDs);
    CumulativeLocs = [0; cumsum(NLocPerFrame)];
    for ff = 1:numel(FrameArray)
        % Determine which localizations should be considered for clustering.
        % NOTE: Even though we don't want clusters with multiple 
        %       localizations in the same frame for the final results, we 
        %       don't want to exclude those until later (since inclusion of
        %       WRONG localizations now can exclude CORRECT localizations 
        %       if we restrict same frame localizations).
        CurrentFrameInd = (1:NLocPerFrame(ff)) + CumulativeLocs(ff);
        CandidateFrameInd = ...
            find((FrameNumCDs >= (FrameArray(ff)-MaxFrameGap)) ...
            & (FrameNumCDs<FrameArray(ff)));
        if isempty(CandidateFrameInd)
            MaxID = MaxID + 1;
            ConnectID(CurrentFrameInd) = (1:NLocPerFrame(ff)).' + MaxID;
            continue
        end
        
        % Determine the nearest neighbor to the current localizations in
        % all candidate frames (noting that we're allowing comparisons to
        % the current frame as well).
        [NNIndices, NNDistances] = knnsearch(...
            [XCDs(CandidateFrameInd), YCDs(CandidateFrameInd)], ...
            [XCDs(CurrentFrameInd), YCDs(CurrentFrameInd)], ...
            'k', 2);
        NNIndices = NNIndices(:, 2:end);
        NNDistances = NNDistances(:, 2:end);
        
        % Place the CurrentFrameInd localizations into clusters.
        ValidNNInd = find(NNDistances ...
            <= (NSigmaDev*MeanXYSECDs(CurrentFrameInd)));
        if isempty(ValidNNInd)
            continue
        end
        for nn = ValidNNInd.'
            % Place this localization into the same cluster as its nearest
            % neighbor.
            NNIndex = CandidateFrameInd(NNIndices(nn));
            if IsClustered(NNIndex)
                ConnectIDCDs(CurrentFrameInd(nn)) = ConnectIDCDs(NNIndex);
                IsClustered(CurrentFrameInd(nn)) = true;
            else
                MaxID = MaxID + 1;
                ConnectIDCDs([CurrentFrameInd(nn), NNIndex]) = MaxID;
                IsClustered([CurrentFrameInd(nn), NNIndex]) = true;
            end
        end
    end
    ConnectID(CurrentDatasetInd(SortIndicesFN)) = ConnectIDCDs;
end
SMD.ConnectID(SortIndices, 1) = smi.SPT.validifyConnectID(ConnectID);


end