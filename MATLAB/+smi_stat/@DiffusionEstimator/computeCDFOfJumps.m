function [MSDStruct] = computeCDFOfJumps(MSDStruct)
%computeCDFOfMSD computes the CDF (CPD) of the jumps in 'MSDStruct'.
% This method will compute the empirical cumulative distribution function
% of the displacements provided in 'MSDStruct'.
%
% INPUTS:
%   MSDStruct: A structure array of MSD data as output from computeMSD()
%              (see computeMSD() outputs for more details).
%
% OUTPUTS:
%   MSDStruct: A structure array of MSD data as output from computeMSD()
%              (see computeMSD() for more details), with a new field(s)
%              added containing the CDF of the MSD.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Loop through the provided MSDStruct and compute the CDF.
for ii = 1:numel(MSDStruct)
    if ~isempty(MSDStruct(ii).MSD)
        % Sort arrays as needed.
        [MSDStruct(ii).SortedSquaredDisp, SortIndices] = ...
            sort(MSDStruct(ii).SquaredDisplacement);
        MSDStruct(ii).FrameLagsAll = ...
            MSDStruct(ii).FrameLagsAll(SortIndices);
        
        % Compute the CDF of the jumps.
        NJumps = numel(MSDStruct(ii).SortedSquaredDisp);
        MSDStruct(ii).CDFOfJumps = cumsum(ones(NJumps, 1)) / NJumps;
    end
end


end