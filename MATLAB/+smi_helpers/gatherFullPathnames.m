function fullpaths = gatherFullPathnames(paths, files, pattern)
%gatherFullPathnames combines or searches for directory paths and filenames
% matching a pattern.
%
% Matches paths to files.  If neither paths nor files are provided, call up a
% GUI to let the user interactively choose.  If files are provided, but paths
% are not, take the current directory as the starting point.  If paths are
% provided, but files are not, look for filenames matching the pattern.
%
% INPUTS:
%    paths      cell array of directory paths
%    files      cell array of filenames
%    pattern    MATLAB file pattern (e.g., '*_Results.mat')
%
% OUTPUT:
%   fullpaths   cell array of full pathnames of combined directory paths with
%               filenames

% Created by
%    Michael J. Wester (Lidkelab 2022)

   if ~exist('pattern', 'var')
      pattern = '*_Results.mat';
   end
   if ~exist('files', 'var')
      files = {};
   end
   if ~exist('paths', 'var')
      paths = {};
   end

   if ~iscell(paths)
      paths = { paths };
   end
   if ~iscell(files)
      files = { files };
   end

   n_paths = numel(paths);
   n_files = numel(files);

   if n_paths == 0 && n_files == 0
      % If no paths or files are provided, interactively prompt for those
      % matching a pattern.  First, transform a MATLAB pattern into a regular
      % expression pattern.
      REpattern = regexprep(pattern, '*', '.*');
      fullpaths = uipickfiles('FilterSpec', '.', 'REFilter', REpattern, ...
                              'Prompt', [pattern, ' files']);
      fullpaths = fullpaths';
      return;
   end

   if n_files > 0
      n_fullpaths = max(n_paths, n_files);
      fullpaths = cell(n_fullpaths, 1);

      if n_paths == 0
         % Only files are provided.
         for i = 1 : n_fullpaths
            fullpaths{i} = files{i};
         end
      elseif n_paths == n_files
         % Paths and files are both provided in a matched set.
         for i = 1 : n_fullpaths
            fullpaths{i} = fullfile(paths{i}, files{i});
         end
      else
         % Numbers don't match---error!
         error(sprintf('n_paths (%d) != n_files (%d)', n_paths, n_files));
      end
   else   % n_files == 0
      % Look for the pattern in the filenames on each directory path.
      n_fullpaths = 0;
      for i = 1 : n_paths
          D = dir(fullfile(paths{i}, pattern));
          for j = 1 : numel(D)
             n_fullpaths = n_fullpaths + 1;
             fullpaths{n_fullpaths} = fullfile(paths{i}, D(j).name);
          end
      end
      fullpaths = fullpaths';
   end

end
