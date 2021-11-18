function [SMD] = classicalFC(SMD, SMF, Verbose)
%classicalFC connects localizations in 'SMD' by simple thresholds.
% This method solves the frame-connection problem by connecting
% localizations within hard spatiotemporal thresholds.
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
FrameNum = SMD.FrameNum(SortIndices);
MaxFrameGap = SMF.FrameConnection.MaxFrameGap;
MaxSeparation = SMF.FrameConnection.MaxSeparation;

% Initialize each localization as a new cluster.
ConnectID = (1:NLocalizations).';

% Loop through datasets and perform frame-connection.
[NLocPerDataset, DatasetArray] = groupcounts(DatasetNum);
CumulativeDatasetLocs = [0; cumsum(NLocPerDataset)];
MaxID = NLocalizations;
for ii = 1:numel(DatasetArray)
    % Provide a Command Window update if needed.
    if (Verbose > 2)
        fprintf(['\tFrameConnection.classicalFC(): ', ...
            'Performing frame connection for dataset %i...\n'], ii)
    end
    
    % Isolate some arrays for the current dataset (CDS = current dataset)
    CurrentDatasetInd = (1:NLocPerDataset(ii)) + CumulativeDatasetLocs(ii);
    [FrameNumCDs, SortIndicesFN] = sort(FrameNum(CurrentDatasetInd));
    XCDs = X(CurrentDatasetInd(SortIndicesFN));
    YCDs = Y(CurrentDatasetInd(SortIndicesFN));
    ConnectIDCDs = ConnectID(CurrentDatasetInd(SortIndicesFN));
    
    % Loop through frames and add localizations to clusters.
    IsClustered = zeros(NLocPerDataset(ii), 1, 'logical');
    [NLocPerFrame, FrameArray] = groupcounts(FrameNumCDs);
    CumulativeLocs = [0; cumsum(NLocPerFrame)];
    for ff = 1:numel(FrameArray)
        % Determine which localizations should be considered for clustering
        CurrentFrameInd = (1:NLocPerFrame(ff)) + CumulativeLocs(ff);
        CandidateFrameInd = ...
            find((FrameNumCDs >= (FrameArray(ff)-MaxFrameGap)) ...
            & (FrameNumCDs<FrameArray(ff)));
        if isempty(CandidateFrameInd)
            ConnectID(CurrentFrameInd) = (1:NLocPerFrame(ff)).' + MaxID;
            MaxID = MaxID + NLocPerFrame(ff);
            continue
        end
        
        % Determine the nearest neighbor to the current localizations in
        % all candidate frames.
        [NNIndices, NNDistances] = knnsearch(...
            [XCDs(CandidateFrameInd), YCDs(CandidateFrameInd)], ...
            [XCDs(CurrentFrameInd), YCDs(CurrentFrameInd)], ...
            'k', 1);
        
        % Place the CurrentFrameInd localizations into clusters.
        ValidNNInd = find(NNDistances <= MaxSeparation);
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
SMD.ConnectID(SortIndices, 1) = smi_helpers.compressToRange(ConnectID);


end