function [TimeString] = genTimeString(Delimiter)
%genTimeString creates a char array with the current time.
% This method creates a char array which contains the current time as
% determined by clock().
%
% INPUTS:
%   Delimiter: Delimiter between the time values in the output 'TimeString'
%              (e.g., the underscore in '2021_2_17_19_3_26')
%              (char array)(Default = '_')
%
% OUTPUTS:
%   TimeString: char array containing the current time. (char array)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('Delimiter', 'var') || isempty(Delimiter))
    Delimiter = '_';
end

% Create the time string.
CurrentTime = cellstr(num2str(round(clock().')));
TimeString = erase(strjoin(CurrentTime, Delimiter), ' ');
    

end