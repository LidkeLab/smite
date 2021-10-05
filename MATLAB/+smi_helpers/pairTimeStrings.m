function [PairedIndices] = ...
    pairTimeStrings(TimeStrings, TimeStringOptions, Comparison)
%pairTimeStrings pairs time strings based on temporal proximity.
% This function selects indices from the input 'TimeStringOptions' that are
% closest in time to the timestrings provided in 'TimeStrings'.
%
% INPUTS:
%   TimeStrings: Cell array of timestrings formatted as would be output by
%                smi_helpers.genTimeString().
%   TimeStringOptions: Set of timestrings to be paired to 'TimeStrings'.
%   Comparison: Method of selection between the timestrings.
%               (Default = 'Before')
%               'Before': closest timestring that happened BEFORE
%                         'TimeStrings' is selected.
%               'After': closest timestring that happened AFTER
%                        'TimeStrings' is selected.
%               'Nearest': closest timestring is selected, whether it's
%                          before or after in time.
%
% OUTPUTS:
%   PairedIndices: Indices of 'TimeStringOptions' corresponding to pairings
%                  with 'TimeStrings'.  If any of 'TimeStrings' remains
%                  unpaired, it will be given an index of NaN.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('Comparison', 'var') || isempty(Comparison))
    Comparison = 'Before';
end

% Loop through the timestrings and compare.
TimeNums = smi_helpers.convertTimeStringToNum(TimeStrings);
TimeNumOptions = smi_helpers.convertTimeStringToNum(TimeStringOptions);
PairedIndices = NaN(size(TimeStrings));
NStrings = numel(TimeStrings);
NOptions = numel(TimeStringOptions);
for ii = 1:NStrings
    TimeDiff = TimeNums(ii) - TimeNumOptions;
    switch lower(Comparison)
        case 'before'
            ValidInd = find(TimeDiff >= 0);
        case 'after'
            ValidInd = find(TimeDiff <= 0);
        case 'nearest'
            ValidInd = 1:NOptions;
        otherwise
            error('Invalid input ''Comparison'' = %s', Comparison)
    end
    [~, BestIndex] = min(TimeDiff(ValidInd));
    if ~isempty(BestIndex)
        PairedIndices(ii) = ValidInd(BestIndex);
    end
end


end