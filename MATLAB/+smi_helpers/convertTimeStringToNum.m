function [TimeNum] = convertTimeStringToNum(TimeString, Delimiter)
%convertTimeStringToNum converts a time string to a number.
% This function will take a timestring (e.g., the output from
% smi_helpers.genTimeString()) and convert it to a number which can be used
% for comparison to other times of the same format.
%
% INPUTS:
%   TimeString: A time string given in the format output by genTimeString()
%               (char/string, or cell arraay of char/string)
%               NOTE: Changing the input 'MinFieldWidth' of genTimeString()
%                     might cause issues in this function!
%   Delimiter: A delimiter between the year, month, day, etc. in TimeString
%              (Default = {'_'; '-'; ','; '.'})
%
% OUTPUTS:
%   TimeNum: A number corresponding to TimeString.
%            CAUTION: I haven't been careful about this conversion, so use
%                     caution if trying to use this number for anything
%                     other than, e.g., file sorting, where the absolute
%                     time doesn't really matter.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults and revise inputs if needed.
if (~exist('Delimiter', 'var') || isempty(Delimiter))
    Delimiter = {'_'; '-'; ','; '.'; ';'; ':'};
end
if ~iscell(TimeString)
    TimeString = {TimeString};
end

% Loop through the input timestrings and convert.
NStrings = numel(TimeString);
TimeNum = NaN(size(TimeString));
for ii = 1:NStrings
    % Break up the time string and convert to a number.
    if iscell(TimeString{ii})
        % Some methods output a cell array of cells (e.g., strsplit()) so
        % we need to access the character array with {1}.
        SplitTimeStamp = strsplit(TimeString{ii}{1}, Delimiter);
    else
        SplitTimeStamp = strsplit(TimeString{ii}, Delimiter);
    end
    for jj = 1:numel(SplitTimeStamp)
        % Ensure we pad with zeros wherever needed, e.g., 2018-2-9-17-41
        % should be 2018-02-09-17-41
        if mod(numel(SplitTimeStamp{jj}), 2)
            SplitTimeStamp{jj} = ['0', SplitTimeStamp{jj}];
        end
    end
    
    % Convert the cell array to a double.
    TimeNum(ii) = str2double(cell2mat(SplitTimeStamp));
end


end