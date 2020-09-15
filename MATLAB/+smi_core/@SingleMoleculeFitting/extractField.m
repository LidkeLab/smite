function [FieldValue] = extractField(SMF, FieldName)
%ExtractField extracts the field FieldName from the given SMF structure.
% This method will search through an SMF structure and return the value
% stored in the field given by FieldName.
%
% NOTE: If there are identically named fields in SMF, this method will
%       return the first instance found, and thus special care should be
%       taken for those fields (i.e., this method might not return what you
%       want it to in that case!).
%
% EXAMPLE USAGE:
%   PSFSigma = smi_core.SingleMoleculeFitting.extractField(SMF, 'PSFSigma')
%
% INPUTS:
%   SMF: Single Molecule Fitting structure (or similarly organized struct)
%   FieldName: Name of field to be extracted from SMF (Char. array)
%
% OUTPUTS:
%   FieldValue: The value stored in SMF in the field named FieldName.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Initialize the output FieldValue to an empty (in case we never find it
% inside of SMF).
FieldValue = [];

% Generate a list of the fields in the SMF (these should all be
% sub-structures at this stage, unless the SMF structure comes from an
% older version of the code).
TopLevelFields = fieldnames(SMF);

% Search for FieldValue at this level of the SMF.
IsDesiredField = strcmp(TopLevelFields, FieldName);
if any(IsDesiredField)
    FieldValue = SMF.(TopLevelFields{IsDesiredField});
    return
end

% Search for the desired field at the top level.
for ff = 1:numel(TopLevelFields)
    if isstruct(SMF.(TopLevelFields{ff}))
        % If this is a sub-structure, we need to iteratively call this
        % method to keep digging deeper.
        FieldValue = smi_core.SingleMoleculeFitting.extractField(...
            SMF.(TopLevelFields{ff}), FieldName);
        if ~isempty(FieldValue)
            % The desired field was found so we should return to caller.
            return
        end
    end
end


end