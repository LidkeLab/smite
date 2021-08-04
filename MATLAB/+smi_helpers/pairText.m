function [PairedList1, PairedList2] = pairText(List1, List2, ExcludeText)
%pairText creates paired lists of text from two sets of text lists.
% This method will remove the 'ExcludeText' from each char/string array in
% 'List1' and 'List2' and then attempt to match the two lists based on the
% remaining char/strings.


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


end