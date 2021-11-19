function [ConnectID] = clusterSTDist(SMD, MaxFrameGap, MaxDist)
%clusterSTSigma performs pre-clustering on localizations in SMD.
% This method clusters localizations in SMD based on their spatiotemporal
% separations.  Localizations within MaxDist of one another which 
% appear within MaxFrameGap frames will be assigned to the same cluster.  
% The assignment is designated by a shared integer value of the output 
% ConnectID.
%
% NOTE: This function was originally written in the context of
%       FrameConnection.lapFC().
%
% INPUTS:
%   SMD: SingleMoleculeData structure with the localizations that we wish
%        to frame-connect.
%   MaxFrameGap: Maximum frame gap allowed between cluster members.
%   MaxDist: Maximum distance allowed between members of the same cluster.
%
% OUTPUTS:
%   ConnectID: Set of integers defining links between localizations in SMD,
%              with indexing matching the indices of SMD localizations.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Gather/revise/reorganize some arrays for further use.
[DatasetNum, SortIndices] = sort(SMD.DatasetNum);
X = SMD.X(SortIndices);
Y = SMD.Y(SortIndices);
X_SE = SMD.X_SE(SortIndices);
Y_SE = SMD.Y_SE(SortIndices);
FrameNum = SMD.FrameNum(SortIndices);
MeanXYSE = mean([X_SE, Y_SE], 2);

% Initialize each localization as a new cluster.
NLocalizations = numel(SMD.FrameNum);
ConnectID = (1:NLocalizations).';

% Loop through datasets and perform the pre-clustering.
[NLocPerDataset, DatasetArray] = groupcounts(DatasetNum);
CumulativeDatasetLocs = [0; cumsum(NLocPerDataset)];
MaxID = NLocalizations;
for ii = 1:numel(DatasetArray)
    % Isolate some arrays for the current dataset (CDS = current dataset)
    CurrentDatasetInd = (1:NLocPerDataset(ii)) + CumulativeDatasetLocs(ii);
    [FrameNumCDs, SortIndicesFN] = sort(FrameNum(CurrentDatasetInd));
    CurrentDatasetInd = CurrentDatasetInd(SortIndicesFN);
    XCDs = X(CurrentDatasetInd);
    YCDs = Y(CurrentDatasetInd);
    MeanXYSECDs = MeanXYSE(CurrentDatasetInd);
    ConnectIDCDs = ConnectID(CurrentDatasetInd);
    
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
            & (FrameNumCDs<=FrameArray(ff)));
        if isempty(CandidateFrameInd)
            MaxID = MaxID + 1;
            ConnectID(CurrentFrameInd) = (1:NLocPerFrame(ff)).' + MaxID;
            MaxID = MaxID + NLocPerFrame(ff);
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
        ValidNNInd = find(NNDistances <= MaxDist);
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
    ConnectID(CurrentDatasetInd) = ConnectIDCDs;
end
ConnectID(SortIndices, 1) = ConnectID;
ConnectID = smi_helpers.compressToRange(ConnectID);


end