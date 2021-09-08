function [FileNames] = getFileNames(FileDir, NameString)
%getFileNames creates a list of filenames in FileDir.
% This function generates a cell array of filenames matching 'NameString'
% (which can contain a wild card '*', but is otherwise not a regexp [see
% usage of MATLAB built-in dir()]) within a directory 'FileDir'.
%
% INPUTS:
%   FileDir: Directory containing the files of interest.
%            (char array)(Default = pwd())
%   NameString: Pattern that files in 'FileDir' must match to be listed.
%               (char array)(Default = '*')
%
% OUTPUTS:
%   FileNames: List of the filenames in 'FileDir' matching 'NameString'.
%              (cell array)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% If 'FileDir' is empty, stop now (this is convenient for some methods that
% use this function).
if (~exist('FileDir', 'var') || isempty(FileDir) || ~isfolder(FileDir))
    FileNames = {};
    return
end

% Set defaults if needed.
if (~exist('NameString', 'var') || isempty(NameString))
    NameString = '*';
end

% Determine the contents in 'FileDir'.
DirContents = dir(fullfile(FileDir, NameString));

% Remove the directories listed in 'DirContents'.
SubDirBoolean = [DirContents.isdir];
FileNames = {DirContents(~SubDirBoolean).name}.';


end