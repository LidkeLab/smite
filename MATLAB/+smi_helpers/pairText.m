function [PairedList1, PairedList2] = pairText(List1, List2, ExcludeText)
%pairText creates paired lists of text from two sets of text lists.
% This method will remove the 'ExcludeText' from each char/string array in
% 'List1' and 'List2' and then attempt to match the two lists based on the
% remaining char/strings.
%
% EXAMPLE:
%   List1 = {'file1_ch1.mat', 'file3_ch1.mat', 'file2_ch1.mat'};
%   List2 = {'file2_ch2.mat', 'file1_ch2.mat', 'file3_ch2.mat'};
%   [PairedList1, PairedList2] = smi_helpers.pairText(List1, List2, ...
%                                {'ch1', 'ch2'});
%
% INPUTS:
%   List1: Cell array of text. (cell array of char)
%   List2: Cell array of text.  Entries in 'List2' will (ideally) be
%          matched one-to-one with files in 'List1' after erasing the
%          sub-strings contained in 'ExcludeText'. (cell array of char)
%   ExcludeText: Sub-strings which will be erased from 'List1' and 'List2'
%                before attempting to match those lists. 
%                (char/cell array of char)
%
% OUTPUTS:
%   PairedList1: Text entries from 'List1' which were paired to entries in
%                'List2'. (cell array of char)
%   PairedList2: Text entries from 'List2' which were paired to entries in
%                'List1'.  The indexing of 'PairedList1' and 'PairedList2'
%                should match, so that PairedList1{ii} is matched to the
%                text entry in PairedList2{ii}. (cell array of char)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Remove the 'ExcludeText' from each of the lists.
List1Erase = erase(List1, ExcludeText);
List2Erase = erase(List2, ExcludeText);

% Loop through the lists and search for matching strings in the opposite
% list.
NList1 = numel(List1);
NList2 = numel(List2);
PairedList1 = {};
PairedList2 = {};
for ii = 1:NList1
    for jj = 1:NList2
        if strcmp(List1Erase{ii}, List2Erase{jj})
            PairedList1 = [PairedList1; List1{ii}];
            PairedList2 = [PairedList2; List2{jj}];
        end
    end
end

% Ensure output arrays match the "shape" of the input arrays.
if isrow(List1)
    PairedList1 = PairedList1.';
end
if isrow(List2)
    PairedList2 = PairedList2.';
end


end