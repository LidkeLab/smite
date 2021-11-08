function [SMD] = isolateSubSMD(SMD, SubIndices)
%isolateSubSMD isolates a subset of SMD defined by SubIndices.
% This method will isolate a subset of SMD defined by the indices in
% SubIndices.  This is done automatically by indexing all SMD arrays which
% share the same size as SMD.FrameNum.
%
% WARNING: This method can fail if there are "constant" SMD fields that
%          happen to share the same array size as SMD.FrameNum.
%
% INPUTS:
%   SMD: A Single Molecule Data structure.
%   SubIndices: Set of indices/booleans corresponding to the desired 
%               entries of SMD. (integer array/logical array)
%
% OUTPUTS:
%   SMD: A Single Molecule Data structure whose vector fields correspond to
%        the input SMD.*(SubIndices), with other fields left unchanged.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Determine the length of vector fields and ensure consistency with the
% provided SubIndices.
ArrayLength = numel(SMD.FrameNum);
if islogical(SubIndices)
    assert(numel(SubIndices) == ArrayLength, ...
        ['If ''SubIndices'' is a boolean array, it must be the same ', ...
        'length as SMD.FrameNum'])
else
    assert(max(SubIndices) <= ArrayLength, ...
        ['Provided ''SubIndices'' contains indices greater than the ', ...
        'length of SMD.FrameNum'])
end

% Loop through all SMD fields and update the output SMD as necessary.
SMDIn = SMD;
SMDFields = fieldnames(SMD);
for ff = 1:numel(SMDFields)
    if (numel(SMDIn.(SMDFields{ff})) == ArrayLength)
        SMD.(SMDFields{ff}) = SMDIn.(SMDFields{ff})(SubIndices);
    end
end


end