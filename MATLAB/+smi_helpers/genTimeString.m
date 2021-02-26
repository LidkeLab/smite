function [TimeString] = genTimeString(Delimiter, MinFieldWidth)
%genTimeString creates a char array with the current time.
% This method creates a char array which contains the current time as
% determined by clock().
%
% INPUTS:
%   Delimiter: Delimiter between the time values in the output 'TimeString'
%              (e.g., the underscore in '2021_2_17_19_3_26')
%              (char array)(Default = '_')
%   MinFieldWidth: Minimum field width.  For example, if MinFieldWidth = 2,
%                  the output would look like '2021_02_17_19_03_26',
%                  whereas MinFieldWidth = 1 would result in 
%                  '2021_2_17_19_3_26'. (Default = 2)
%
% OUTPUTS:
%   TimeString: char array containing the current time. (char array)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('Delimiter', 'var') || isempty(Delimiter))
    Delimiter = '_';
end
if (~exist('MinFieldWidth', 'var') || isempty(MinFieldWidth))
    MinFieldWidth = 2;
end

% Create the time string.
CurrentTime = round(clock());
TimeString = sprintf(['%0', num2str(MinFieldWidth), 'i', Delimiter], ...
    CurrentTime);
TimeString = TimeString(1:end-1);
    

end