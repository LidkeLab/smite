function [SMD] = catSMD(SMD1, SMD2)
%catSMD concatenates two SMD structures into one.
% This method will concatenate two SMD structures, SMD1 and SMD2, into a
% single output structure, SMD.  Vector fields are concatenated directly,
% i.e., SMD.X = [SMD1.X; SMD2.X] .  Scalar fields are checked for
% consistency (i.e., ensured they are equal) before being stored in the
% output SMD.  Scalar fields which are not consistent are stored in the
% output SMD concatenated in a cell array.
%
% NOTE: There is a hard-coded set of fields below called StackableFields,
%       which are the fields that can be safely concatenated as 
%       SMD.Example = [SMD1.Example; SMD2.Example]. If these are changed in 
%       SingleMoleculeData.m (or if they are used differently than usual in
%       some context), this set of fields may need to be revised. 
%       Similarly, there is a set of "special" fields defined by 
%       SpecialFields which are treated differently towards the end of this
%       method (their treatment will be different on a case-by-case basis).
%
% INPUTS:
%   SMD1: A Single Molecule Data structure.
%   SMD2: A Single Molecule Data structure.
%
% OUTPUTS:
%   SMD: A Single Molecule Data structure whose vector fields are the
%        concatenated vector fields from SMD1 and SMD2.  The scalar fields
%        will be kept as scalar fields when equal in both SMD1 and SMD2,
%        and kept as cell arrays when they are not.

% Created by:
%   David J. Schodt (Lidke lab, 2020)


% Ensure neither of the input SMD's are empty.  If one of them is empty,
% output the other one without changes.  If both are empty, throw an error.
if (isempty(SMD1) && isempty(SMD2))
    error('Inputs SMD1 and SMD2 cannot both be empty');
elseif isempty(SMD1)
    SMD = SMD2;
    return
elseif isempty(SMD2)
    SMD = SMD1;
    return
end

% Define which SMD fields can be stacked directly (see note at top of code)
% NOTE: We could just check which fields are numeric vectors with 
%       (numel() > 1), but that convenience limits us to SMD structures 
%       with more than one localization.
StackableFields = {'XBoxCorner', 'YBoxCorner', ...
    'X', 'Y', 'Z', ...
    'X_SE', 'Y_SE', 'Z_SE', ...
    'Photons', 'Photons_SE', ...
    'Bg', 'Bg_SE', ...
    'PSFSigma', 'PSFSigmaX', 'PSFSigmaY', ...
    'PSFSigma_SE', 'PSFSigmaX_SE', 'PSFSigmaY_SE', ...
    'PValue', 'LogLikelihood', 'ThreshFlag', ...
    'NCombined', 'IndSMD', 'DatasetNum', 'FrameNum', 'ConnectID'};

% Define a set of "special" fields in SMD: these fields are treated
% differently on a case-by-case basis towards the end of this method.
% NOTE: I mostly added this to bypass the warnings that happen in the main
%       for loop below. In fact, not all of the "special" fields are even
%       in this list, e.g., PSFSigma is also treated differently in some
%       cases, but not included in this list.
SpecialFields = {'DriftX', 'DriftY', 'DriftZ', 'NDatasets', ...
    'ConnectID', 'IndSMD'};

% Create a list of fields present in each of SMD1 and SMD2.
FieldNames1 = fieldnames(SMD1);
FieldNames2 = fieldnames(SMD2);

% Check if SMD1 and SMD2 have the same fields, warning the user if this
% isn't the case.
UniqueFields1 = FieldNames1(~ismember(FieldNames1, FieldNames2));
UniqueFields2 = FieldNames2(~ismember(FieldNames2, FieldNames1));
if ~isempty(UniqueFields1)
    PrintFriendlyFields1 = cell2mat(cellfun(@(X) [X, ' '], ...
        UniqueFields1, 'UniformOutput', false).');
    warning(['smi_core.SingleMoleculeData.catSMD() - ', ...
        'The following SMD1 fields aren''t present in SMD2: ', ...
        PrintFriendlyFields1])
end
if ~isempty(UniqueFields2)
    PrintFriendlyFields2 = cell2mat(cellfun(@(X) [X, ' '], ...
        UniqueFields2, 'UniformOutput', false).');
    warning(['smi_core.SingleMoleculeData.catSMD() - ', ...
        'The following SMD2 fields aren''t present in SMD1: ', ...
        PrintFriendlyFields2])
end

% Place the unique fields (fields that are only in one of SMD1 or SMD2) in
% the SMD structure directly.
SMD = smi_core.SingleMoleculeData.createSMD();
for uu = 1:numel(UniqueFields1)
    SMD.(UniqueFields1{uu}) = SMD1.(UniqueFields1{uu});
end
for uu = 1:numel(UniqueFields2)
    SMD.(UniqueFields2{uu}) = SMD2.(UniqueFields2{uu});
end

% Loop through fields present in SMD1 and SMD2 and begin constructing the
% concatenated output SMD.
% NOTE: There are a several "special" fields (which aren't defined in
%       VectorFields above) that will be defined illogically/incorrectly in
%       this loop.  Those fields will be revised later on.
AllFieldNames = unique([FieldNames1; FieldNames2]);
SharedFields = AllFieldNames(ismember(AllFieldNames, FieldNames1) ...
    & ismember(AllFieldNames, FieldNames2));
IsStackable = ismember(SharedFields, StackableFields);
IsSpecialField = ismember(SharedFields, SpecialFields);
for ff = 1:numel(SharedFields)
    % Place the current field in the output SMD, concatenating when
    % appropriate.
    CurrentField = SharedFields{ff};
    if IsStackable(ff)
        % Stackable fields should be concatenated directly.
        SMD.(CurrentField) = [SMD1.(CurrentField); SMD2.(CurrentField)];
    elseif isequal(SMD1.(CurrentField), SMD2.(CurrentField))
        % We want to ensure non-stackable fields are consistent with one
        % another before setting the output field.
        SMD.(CurrentField) = SMD1.(CurrentField);
    elseif ~IsSpecialField(ff)
        % This appears to be a non-stackable field, but the field is not
        % equal in SMD1 and SMD2; we will trust the user and concatenate
        % these two as though it were a vector field.
        warning(['smi_core.SingleMoleculeData.catSMD() - ', ...
            'SMD field ''%s'' is different in each of ', ...
            'SMD1 and SMD2, attempting to concatenate.'], CurrentField)
        try
            SMD.(CurrentField) = [SMD1.(CurrentField); ...
                SMD2.(CurrentField)];
        catch Exception
            warning(['smi_core.SingleMoleculeData.catSMD() - ', ...
                'Concatenation of field ''%s'' failed, storing as ', ...
                'cell array.\n%s, %s'], CurrentField, ...
                Exception.identifier, Exception.message)
            SMD.(CurrentField) = {SMD1.(CurrentField); ...
                SMD2.(CurrentField)};
        end
    end
end

% Revise fields that may or may not be vector fields (e.g., PSFSigma is
% either scalar or vector depending on the fit that was used).
if (numel(unique(SMD.PSFSigma)) == 1)
    SMD.PSFSigma = unique(SMD.PSFSigma);
    SMD.PSFSigma_SE = unique(SMD.PSFSigma_SE);
end
if (numel(unique(SMD.DatasetNum)) == 1)
    SMD.DatasetNum = unique(SMD.DatasetNum);
end

% Revise the drift related fields, which should be treated differently.
% NOTE: The dimensions of drift will be [NFrames, NDatasets], so if one of 
SMD.DriftX = [SMD1.DriftX, SMD2.DriftX];
SMD.DriftY = [SMD1.DriftY, SMD2.DriftY];
SMD.DriftZ = [SMD1.DriftZ, SMD2.DriftZ];

% Define the output SMD.NDatasets in a manner consistent with
% SMD.DatasetNum.
SMD.NDatasets = numel(unique(SMD.DatasetNum));
if (SMD.NDatasets ~= (SMD1.NDatasets+SMD2.NDatasets))
    warning(['smi_core.SingleMoleculeData.catSMD() - ', ...
        'Concatenated SMD field ''NDatasets'' may not match ', ...
        'expectations: SMD.NDatasets set to ', ...
        'numel(unique(SMD.DatasetNum))'])
end

% Define the outputs SMD.ConnectID and SMD.IndSMD ensuring consistency
% (their values represent indices which lose meaning after concatenation).
% NOTE: While this might make more sense to put at the top of this code, I 
%       wanted to do it this way (with 'ConnectID' and 'IndSMD' added to
%       'SpecialFields') so that it's more clear to readers that something
%       special is being done with these fields.
MaxConnectIDSMD1 = max(SMD1.ConnectID);
MaxIndSMD1 = max(cell2mat(SMD1.IndSMD));
SMD2.ConnectID = SMD2.ConnectID + MaxConnectIDSMD1;
SMD2.IndSMD = cellfun(@(X) X + MaxIndSMD1, SMD2.IndSMD, ...
    'UniformOutput', false);
SMD.ConnectID = [SMD1.ConnectID; SMD2.ConnectID];
SMD.IndSMD = [SMD1.IndSMD; SMD2.IndSMD];


end