function [String] = stringMUX(StringOptions, Select)
%stringMUX is a string multiplexer for selecting from a string array.
% This method is a 2-bit multiplexer for strings, meaning that one of two 
% 'StringOptions' will be selected and output as String based on the 
% selector signal 'Select'.
%
% EXAMPLE USAGE:
%   One usage of this would be to output a descriptive string based on a
%   UnitFlag, which is 0 for camera units and 1 for physical units. You
%   might use this code as folows:
%       UnitType = smi_helpers.stringMUX({'camera'; 'physical'}, UnitFlag);
%       If UnitFlag = 0, UnitType = 'camera'
%       If UnitFlag = 1, UnitType = 'physical'
%
% INPUTS:
%   StringOptions: A string array or cell array of char containing the
%                  string options. 
%                  (2x1 string array, or 2x1 cell array of char)
%   Select: Logical (boolean) indicating which of 'StringOptions' to output
%           in 'String'.  Select = 0 outputs StringOptions(1) and 
%           Select = 1 outputs StringOptions(2).
%
% OUTPUTS:
%   String: An element of StringOptions selected by Select.
%           (string or char)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Attempt to select the desired string.
if isstring(StringOptions)
    String = StringOptions(logical(Select) + 1);
elseif iscell(StringOptions)
    String = StringOptions{logical(Select) + 1};
else
    error(['Input ''StringOptions'' must be either a string array ', ...
        'or a cell array of char'])
end


end