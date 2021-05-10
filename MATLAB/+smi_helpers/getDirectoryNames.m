function DirNames = getDirectoryNames(BaseDir, NameString)
%getDirectoryNames generates directory names matching a NameString.
% This function generates a cell array of directory names matching 
% NameString (which can contain a wild card '*') within a directory
% BaseDir.
%
% INPUTS:
%   BaseDir: Character array/string containing the parent directory which
%            contains the directories of interest.
%   NameString: Character array/string describing the names of 
%               sub-directories of interest. (Default = '*')
%
% OUTPUTS:
%   DirNames: Cell array of directory names matching the NameString.

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Set defaults if needed.
if (~exist('NameString', 'var') || isempty(NameString))
    NameString = '*';
end

% Create a list of the contentes within BaseDir.
DirContents = dir(fullfile(BaseDir, NameString));

% Grab the boolean array from DirContents which tells us which of the
% structure of contents are actually directories.
SubDirBoolean = [DirContents.isdir];

% Create a boolean describing the "valid" directories (i.e., not the . and
% .. current and parent directories).
ValidDirBoolean = ~(strcmp({DirContents.name}, '.') ...
    | strcmp({DirContents.name}, '..'));

% Create a cell array containing the names of the sub-directories.
DirNames = {DirContents(SubDirBoolean & ValidDirBoolean).name}.';


end