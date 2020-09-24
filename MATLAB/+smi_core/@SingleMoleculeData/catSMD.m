function [SMD] = catSMD(SMD1, SMD2)
%catSMD concatenates two SMD structures into one.
% This method will concatenate two SMD structures, SMD1 and SMD2, into a
% single output structure, SMD.  Vector fields are concatenated directly,
% i.e., SMD.X = [SMD1.X; SMD2.X] .  Scalar fields are checked for
% consistency (i.e., ensured they are equal) before being stored in the
% output SMD.  Scalar fields which are not consistent are stored in the
% output SMD concatenated in a cell array.
%
% NOTE: There is a hard-coded set of vector fields below called
%       VectorFields!!!  If these are changed in SingleMoleculeData.m (or
%       if they are used differently than usual in some context), this
%       set of fields may need to be revised.
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


% Define which SMD fields are vector fields vs. scalar fields.
% NOTE: We could just check which fields are numeric vectors with 
%       (numel() > 1), but that convenience limits us to SMD structures 
%       with more than one localization.
VectorFields = {'XBoxCorner', 'YBoxCorner', ...
    'X', 'Y', 'Z', ...
    'X_SE', 'Y_SE', 'Z_SE', ...
    'Photons', 'Photons_SE', ...
    'Bg', 'Bg_SE', ...
    'PSFSigma', 'PSFSigmaX', 'PSFSigmaY', ...
    'PSFSigma_SE', 'PSFSigmaX_SE', 'PSFSigmaY_SE', ...
    'PValue', 'LogL', 'ThreshFlag', ...
    'DatasetNum', 'FrameNum', 'ConnectID', ...
    'DriftX', 'DriftY', 'DriftZ'};

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
AllFieldNames = unique([FieldNames1; FieldNames2]);
SharedFields = AllFieldNames(ismember(AllFieldNames, FieldNames1) ...
    & ismember(AllFieldNames, FieldNames2));
IsVectorField = ismember(SharedFields, VectorFields);
for ff = 1:numel(SharedFields)
    % Place the current field in the output SMD, concatenating when
    % appropriate.
    if IsVectorField(ff)
        % Vector fields should be concatenated directly.
        SMD.(SharedFields{ff}) = [SMD1.(SharedFields{ff}); ...
            SMD2.(SharedFields{ff})];
    elseif isequal(SMD1.(SharedFields{ff}), SMD2.(SharedFields{ff}))
        % We want to ensure non-vector fields are consistent with one
        % another before setting the output field.
        SMD.(SharedFields{ff}) = SMD1.(SharedFields{ff});
    else
        % This appears to be a non-vector field, but the field is not equal
        % in SMD1 and SMD2; we will trust the user and concatenate these
        % two as though it were a vector field.
        warning(['smi_core.SingleMoleculeData.catSMD() - ', ...
            'SMD field ''%s'' is different in each of ', ...
            'SMD1 and SMD2, attempting to concatenate.'], ...
            SharedFields{ff})
        try
            SMD.(SharedFields{ff}) = ...
                [SMD1.(SharedFields{ff}); ...
                SMD2.(SharedFields{ff})];
        catch Exception
            warning(['smi_core.SingleMoleculeData.catSMD() - ', ...
                'Concatenation of field ''%s'' failed, storing as ', ...
                'cell array.\n%s, %s'], SharedFields{ff}, ...
                Exception.identifier, Exception.message)
            SMD.(SharedFields{ff}) = {SMD1.(SharedFields{ff}); ...
                SMD2.(SharedFields{ff})};
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


end