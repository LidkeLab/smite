function [SMFPadded] = padSMF(SMF, SMFPadding)
%padSMF adds fields in SMFPadding to SMF that weren't already present.
% This method will search for fields present in the union of the SMF
% fields and SMFPadding fields and set them in the output SMFPadded.  When
% fields are in both SMF and SMFPadding, the output SMFPadded will contain
% the value present in SMFPadding.
%
% INPUTS:
%   SMF: Single Molecule Fitting structure whose existing entries are
%        retained in the output SMFPadded.
%   SMFPadding: An SMF structure with additional fields not present in SMF
%               (e.g., SMF might be an incomplete SMF structure and
%               SMFPadding might be an SMF with all default values).
%               (Default = smi_core.SingleMoleculeFitting.createSMF())
%
% OUTPUTS:
%   SMFPadded: An SMF structure containing all fields present in the union
%              of SMF and SMFPadding, with values of the intersecting
%              fields always being taken from the input SMF.  Additionally,
%              this output is guaranteed to have all fields defined in
%              smi_core.SingleMoleculeFitting.createSMF(), even if they
%              weren't in either of SMF or SMFPadding.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set defaults if needed.
if (~exist('SMFPadding', 'var') || isempty(SMFPadding))
    SMFPadding = smi_core.SingleMoleculeFitting.createSMF();
end

% Create a default SMF structure to initialize the output.
SMFPadded = smi_core.SingleMoleculeFitting.createSMF();

% Merge the input SMF and SMFPadding structs, treating the SMF as the
% "primary" struct (i.e., the merged field values are only taken from
% SMFPadding when they aren't present in SMF).
[SMFInputMerged] = mergeStructs(SMF, SMFPadding);

% Merge the complete (but default valued) SMFPadded with the merged input
% structures to generate the desired output SMFPadded.
[SMFPadded] = mergeStructs(SMFInputMerged, SMFPadded);


    function [MergedStruct] = mergeStructs(PrimaryStruct, SecondaryStruct)
        % This function will merge PrimaryStruct and SecondaryStruct, with
        % field values coming from SecondaryStruct only if not present in
        % PrimaryStruct.
        
        % Create a list of fields in the input structs.
        FieldsPrimary = fieldnames(PrimaryStruct);
        FieldsSecondary = fieldnames(SecondaryStruct);
        
        % Loop through the fields in SecondaryStruct and add them to
        % PrimaryStruct if needed.
        for ii = 1:numel(FieldsSecondary)
            CurrentSubField = FieldsSecondary{ii};
            if isstruct(SecondaryStruct.(CurrentSubField))
                % If this field is itself a structure, we need to dig
                % deeper into this structure to extract fields, unless it's
                % unique to PrimaryStruct (in which case we just leave it
                % as is).
                if ismember(CurrentSubField, FieldsPrimary)
                    MergedStruct.(CurrentSubField) = mergeStructs(...
                        PrimaryStruct.(CurrentSubField), ...
                        SecondaryStruct.(CurrentSubField));
                else
                    MergedStruct.(CurrentSubField) = ...
                        SecondaryStruct.(CurrentSubField);
                end
            elseif ismember(CurrentSubField, FieldsPrimary)
                % If this field exists in FieldsPrimary, we'll keep it in
                % the output structure.
                MergedStruct.(CurrentSubField) = ...
                    PrimaryStruct.(CurrentSubField);
            else
                MergedStruct.(CurrentSubField) = ...
                    SecondaryStruct.(CurrentSubField);
            end
        end
    end


end