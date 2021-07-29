function [TR] = padTR(TR, TRPadding)
%padTR pads the input 'TR' with fields present in 'TRPadding'.
% This method ensures that the output 'TR' contains all of the fields
% present in 'TRPadding'.  This method is intended to fix some forward and
% backward compatability issues, where old/new TR structures have a
% differing number of fields but still need to be concatenated. 
%
% INPUTS:
%   TR: TR structure which is to be padded.
%   TRPadding: TR structure containing the "complete" set of fields, some
%              of which may need to be added to 'TR'.
%
% OUTPUTS:
%   TR: Input 'TR' padded with fields from 'TRPadding'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Loop through the fields of the padding fields and add them to 'TR' as 
% needed.
PaddingFieldNames = fieldnames(TRPadding);
InputFieldNames = fieldnames(TR);
for pp = 1:numel(PaddingFieldNames)
    if ~any(ismember(PaddingFieldNames{pp}, InputFieldNames))
        % The field PaddingFieldNames{pp} is not present so the default 
        % must be added.
        if iscell(TRPadding(1).(PaddingFieldNames{pp}))
            [TR.(PaddingFieldNames{pp})] = deal({});
        elseif isstruct(TRPadding(1).(PaddingFieldNames{pp}))
            [TR.(PaddingFieldNames{pp})] = deal(struct([]));
        else
            [TR.(PaddingFieldNames{pp})] = deal([]);
        end
    end
end


end