function [MSDStruct] = computeCDFOfJumps(MSDStruct, FrameLagRange)
%computeCDFOfMSD computes the CDF (CPD) of the jumps in 'MSDStruct'.
% This method will compute the empirical cumulative distribution function
% of the displacements provided in 'MSDStruct'.
%
% INPUTS:
%   MSDStruct: A structure array of MSD data as output from computeMSD()
%              (see computeMSD() outputs for more details).
%   FrameLagRange: Range of frame lags included in the CDF computation. 
%                  (2 element array, [min. frame lag, max. frame lag])
%                  (Default = [2, 2], so only jumps across 2 frames are 
%                  included)
%
% OUTPUTS:
%   MSDStruct: A structure array of MSD data as output from computeMSD()
%              (see computeMSD() for more details), with a new field(s)
%              added containing the CDF of the MSD.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Set defaults/revise inputs if needed.
if (~exist('FrameLagRange', 'var') || isempty(FrameLagRange))
    FrameLagRange = [2, 2];
else
    FrameLagRange = sort(FrameLagRange);
end

% Loop through the provided MSDStruct and compute the CDF.
for ii = 1:numel(MSDStruct)
    if ~isempty(MSDStruct(ii).MSD)
        % Sort arrays as needed.
        [MSDStruct(ii).SortedSquaredDisp, SortIndices] = ...
            sort(MSDStruct(ii).SquaredDisplacement);
        MSDStruct(ii).FrameLagsAll = ...
            MSDStruct(ii).FrameLagsAll(SortIndices);
        
        % Isolate the data corresponding to frame lags in 'FrameLagRange'.
        KeepBool = ((MSDStruct(ii).FrameLagsAll>=FrameLagRange(1)) ...
            & (MSDStruct(ii).FrameLagsAll<=FrameLagRange(2)));
        MSDStruct(ii).SortedSquaredDisp = ...
            MSDStruct(ii).SortedSquaredDisp(KeepBool);
        MSDStruct(ii).FrameLagsAll = ...
            MSDStruct(ii).FrameLagsAll(KeepBool);
        
        % Compute the CDF of the jumps.
        NJumps = numel(MSDStruct(ii).SortedSquaredDisp);
        MSDStruct(ii).CDFOfJumps = cumsum(ones(NJumps, 1)) / NJumps;
    end
end


end