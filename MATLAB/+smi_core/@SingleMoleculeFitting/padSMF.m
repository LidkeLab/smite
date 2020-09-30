function [SMFPadded, PaddedFields] = padSMF(SMF, SMFPadding, ...
    DisplayMessages)
%padSMF adds fields in SMFPadding to SMF that weren't already present.
% This method will search for fields present in the union of the SMF
% fields and SMFPadding fields and set them in the output SMFPadded.  When
% fields are in both SMF and SMFPadding, the output SMFPadded will contain
% the value present in SMFPadding.
%
% EXAMPLE USAGE:
%   SMF.Data.CameraGain = 1.23;
%   SMFPadding.Data.CameraGain = 10;
%   SMFPadding.Fitting.FitType = 'XYNBS';
%   [SMFPadded] = smi_core.SingleMoleculeFitting.padSMF(SMF, SMFPadding);
%       SMFPadded.Data.CameraGain==1.23, but SMFPadded.Fitting.FitType will
%       be set to 'XYNBS'.
%   [SMFPadded] = smi_core.SingleMoleculeFitting.padSMF(SMF);
%       SMFPadded will be a complete SMF structure with all default values
%       present except Data.CameraGain, which will be set to
%       SMF.Data.CameraGain (this can be used to ensure an incomplete SMF
%       is padded to contain all fields defined in the
%       SingleMoleculeFitting class).
%
% INPUTS:
%   SMF: Single Molecule Fitting structure whose existing entries are
%        retained in the output SMFPadded.
%   SMFPadding: An SMF structure with additional fields not present in SMF
%               (e.g., SMF might be an incomplete SMF structure and
%               SMFPadding might be an SMF with all default values).
%               (Default = smi_core.SingleMoleculeFitting)
%   DisplayMessages: A flag to specify whether or not a message should be
%                    displayed in the Command Window when a field gets
%                    padded. (Default = 0)
%
% OUTPUTS:
%   SMFPadded: An SMF structure containing all fields present in the union
%              of SMF and SMFPadding, with values of the intersecting
%              fields always being taken from the input SMF.  Additionally,
%              this output is guaranteed to have all fields defined in
%              smi_core.SingleMoleculeFitting, even if they weren't in
%              either of SMF or SMFPadding.
%   PaddedFields: A structure with similar organization to an SMF structure
%                 whose fields are the fields in SMFPadded that were not
%                 present in the input SMF (i.e., these are all of the
%                 fields that were added to SMF to generate SMFPadded).

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set defaults if needed.
if (~exist('SMFPadding', 'var') || isempty(SMFPadding))
    SMFPadding = smi_core.SingleMoleculeFitting;
end
if (~exist('DisplayMessages', 'var') || isempty(DisplayMessages))
    DisplayMessages = 0;
end
if ((nargin>1) && ~isempty(inputname(2)))
    SecondInputName = inputname(2);
else
    SecondInputName = 'default SMF';
end

% Create a default SMF structure to initialize the output.
SMFDefault = smi_core.SingleMoleculeFitting;

% Merge the input SMF and SMFPadding structs, treating the SMF as the
% "primary" struct (i.e., the merged field values are only taken from
% SMFPadding when they aren't present in SMF).
[SMFInputMerged, PaddingStruct] = mergeStructs(SMF, SMFPadding, ...
    inputname(1), SecondInputName, DisplayMessages);

% Merge the complete (but default valued) SMFDefault with the merged input
% structures to generate the desired output SMFPadded.
[SMFPadded, PaddingStructDefaults] = ...
    mergeStructs(SMFInputMerged, SMFDefault, ...
    inputname(1), 'default SMF', DisplayMessages);

% Merge the padding structures to generate the output PaddedFields.
[PaddedFields] = mergeStructs(PaddingStruct, PaddingStructDefaults);

    function [MergedStruct, PaddingStruct] = mergeStructs(...
            PrimaryStruct, SecondaryStruct, ...
            PrimaryStructName, SecondaryStructName, DisplayMessages)
        % This function will merge PrimaryStruct and SecondaryStruct, with
        % field values coming from SecondaryStruct only if not present in
        % PrimaryStruct.  The additional output PaddingStruct will be a
        % structure of similar organization to
        % PrimaryStruct/SecondaryStruct whose fields are those which were
        % padded to PrimaryStruct to generate MergedStruct. If the optional
        % flag DisplayMessages==1, a warning will be displayed in the
        % Command Window each time a field is padded.
        
        % Set default inputs if needed.
        if (~exist('DisplayMessages', 'var') || isempty(DisplayMessages))
            DisplayMessages = 0;
        end
        if (~exist('PrimaryStructName', 'var') ...
                || isempty(PrimaryStructName))
            PrimaryStructName = 'PrimaryStruct';
        end
        if (~exist('SecondaryStructName', 'var') ...
                || isempty(SecondaryStructName))
            SecondaryStructName = 'SecondaryStruct';
        end
        
        % Create a list of fields in the input structs.
        FieldsPrimary = fieldnames(PrimaryStruct);
        FieldsSecondary = fieldnames(SecondaryStruct);
        
        % Loop through the fields in SecondaryStruct and add them to
        % PrimaryStruct if needed.
        PaddingStruct = struct();
        MergedStruct = PrimaryStruct;
        for ii = 1:numel(FieldsSecondary)
            CurrentSubField = FieldsSecondary{ii};
            AddToOutput = 0;
            if isstruct(SecondaryStruct.(CurrentSubField))
                if ismember(CurrentSubField, FieldsPrimary)
                    % This sub-structure is present in PrimaryStruct, but
                    % we still want to ensure that it's complete, so we'll
                    % call the mergeStructs() function recursively.
                    PrimarySubStructName = sprintf('%s.%s', ...
                        PrimaryStructName, CurrentSubField);
                    SecondarySubStructName = sprintf('%s.%s', ...
                        SecondaryStructName, CurrentSubField);
                    [MergedStruct.(CurrentSubField), SubPadStruct] = ...
                        mergeStructs(PrimaryStruct.(CurrentSubField), ...
                        SecondaryStruct.(CurrentSubField), ...
                        PrimarySubStructName, SecondarySubStructName, ...
                        DisplayMessages);
                    if ~isempty(fieldnames(SubPadStruct))
                        % NOTE: When we're only padding part of the
                        %       sub-structure, we don't want to set the
                        %       AddToOutput flag to 1.
                        PaddingStruct.(CurrentSubField) = SubPadStruct;
                    end
                else
                    % In this case, the entire sub-structure didn't exist
                    % in PrimaryStruct, so we'll add it in its entirety.
                    AddToOutput = 1;
                end
            elseif ~ismember(CurrentSubField, FieldsPrimary)
                % This field is unique to SecondaryStruct so we need to add
                % it to the output structure.
                AddToOutput = 1;
            end
            
            % Update the output structure (note that partial sub-structures
            % were already padded in the if/else block(s) above).
            if AddToOutput
                MergedStruct.(CurrentSubField) = ...
                    SecondaryStruct.(CurrentSubField);
                PaddingStruct.(CurrentSubField) = ...
                    SecondaryStruct.(CurrentSubField);
                if DisplayMessages
                    fprintf('Field ''%s'' in %s padded from %s\n', ...
                        CurrentSubField, ...
                        PrimaryStructName, ...
                        SecondaryStructName)
                end
            end
        end
    end


end