function [FilePaths] = getFilePaths(FileDir, NameString)
%getFilePaths creates a list of filepaths to files in FileDir.
% This function generates a cell array of filepaths matching 'NameString'
% (which can contain a wild card '*', but is otherwise not a regexp [see
% usage of MATLAB built-in dir()]) within a directory 'FileDir'.
%
% INPUTS:
%   FileDir: Directory containing the files of interest. (char array)
%   NameString: Pattern that files in 'FileDir' must match to be listed.
%               (char array)(Default = '*')
%
% OUTPUTS:
%   FilePaths: List of the full paths to files in 'FileDir' matching
%              'NameString'. (cell array)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('NameString', 'var') || isempty(NameString))
    NameString = '*';
end

% Determine the contents in 'FileDir'.
DirContents = dir(fullfile(FileDir, NameString));

% Remove the directories listed in 'DirContents'.
SubDirBoolean = [DirContents.isdir];
FilePaths = fullfile(FileDir, {DirContents(~SubDirBoolean).name}).';


end