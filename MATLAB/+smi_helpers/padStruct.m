function [Struct] = padStruct(Struct, DefaultStruct)
%padStruct pads the input 'Struct' with defaults in 'DefaultStruct'.
% This method merges the two structures 'Struct' and 'DefaultStruct', with
% values in 'Struct' taking precedent unless not available.
%
% INPUTS:
%   Struct: Structure array which is to be padded.
%   DefaultStruct: Structure array of "default" fields that will be
%                  added to 'Struct' unless already present in 'Struct', in
%                  which case the "default" values are ignored.
%
% OUTPUTS:
%   Struct: Input 'Struct' padded with values from 'DefaultStruct'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Loop through the fields of the default fields and add them to
% 'Struct' as needed.
DefaultFieldNames = fieldnames(DefaultStruct);
InputFieldNames = fieldnames(Struct);
for pp = 1:numel(DefaultFieldNames)
    if ~any(ismember(DefaultFieldNames{pp}, InputFieldNames))
        % The field DefaultParameterNames{pp} is not present in the
        % Struct structure and so the default must be added.
        Struct.(DefaultFieldNames{pp}) = ...
            DefaultStruct.(DefaultFieldNames{pp});
    elseif iscell(DefaultStruct.(DefaultFieldNames{pp}))
        % Type casting can't be done on cells as above, so we'll just place
        % whatever is present in a cell array.
        if ~iscell(Struct.(DefaultFieldNames{pp}))
            Struct.(DefaultFieldNames{pp}) = ...
                {Struct.(DefaultFieldNames{pp})};
        end
    end
end


end