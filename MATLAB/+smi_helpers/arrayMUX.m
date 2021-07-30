function [Output] = arrayMUX(OutputOptions, Select)
%arrayMUX is a numel(OutputOptions)-bit multiplexer.
% This function is a multiplexer intended to generalize for most arrays.
% This function works by selecting one of the input 'OutputOptions' based
% on the input 'Select', which defines the index of 'OutputOptions' to be
% returned (with the index starting at 0 following multiplexer convention)
%
% EXAMPLE USAGE:
%   One usage of this is to output a char. array defining options for some
%   other input function, e.g., to change the Display parameter of
%   optimset based on a verbosity level:
%       optimset('fmincon', 'Display', ...
%       	smi_helpers.arrayMUX({'none', 'iter', 'final'}, Verbose)) 
%
% INPUTS:
%   OutputOptions: Array which will be indexed based on 'Select'.
%   Select: Signal which defines which of 'OutputOptions' will be returned,
%           with indexing starting at 0 to maintain multiplexer
%           conventions.
%
% OUTPUTS:
%   Output: The element of 'OutputOptions' selected by 'Select'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Attempt to select the desired output. 
% NOTE: If a cell array is provided, I'll assume we want the contents of
%       the specified cell and not the cell itself.  This choice was made
%       based on my intended usage of this function at the time of writing.
if iscell(OutputOptions)
    Output = OutputOptions{Select + 1};
else
    Output = OutputOptions(Select + 1);
end


end