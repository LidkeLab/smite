function [lia] = my_ismemberBuiltinTypes(a,b)
% Extracted and simplified from MATLAB's ismember.m for the above usage.
% General handling.
% Use FIND method for very small sizes of the input vector to avoid SORT.
% Handle empty arrays and scalars.  
numelA = numel(a);
numelB = numel(b);
if numelA == 0 || numelB <= 1
    if numelA > 0 && numelB == 1
        lia = (a == b);
    else
        lia = false(size(a));
    end
    return
end

scalarcut = 5;
if numelA <= scalarcut
    lia = false(size(a));
        for i=1:numelA
            lia(i) = any(a(i)==b(:));   % ANY returns logical.
        end
else
    % Use method which sorts list, then performs binary search.
    % Convert to full to work in C helper.
    
%       % Find out whether list is presorted before sort
%       % If the list is short enough, SORT will be faster than ISSORTED
%       % If the list is longer, ISSORTED can potentially save time
%       checksortcut = 1000;
%       if numelB > checksortcut
%           sortedlist = issorted(b(:));
%       else
%           sortedlist = 0;
%       end
%       if ~sortedlist
%           b = sort(b(:));
%       end
    
    % Use builtin helper function ISMEMBERHELPER:
    % [LIA,LOCB] = ISMEMBERHELPER(A,B) Returns logical array LIA indicating
    % which elements of A occur in B and a double array LOCB with the
    % locations of the elements of A occuring in B. If multiple instances
    % occur, the first occurence is returned. B must be already sorted.
    
            lia = builtin('_ismemberhelper',a,b);
end
