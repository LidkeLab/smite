function [pathname, files] = selectFiles(startdatadir, txt, pattern)
% Select a list of files.
%
% INPUTS:
%    startdatadir   the directory to start looking from
%    txt            [OPTIONAL] a message to display on the GUI
%    pattern        [OPTIONAL] default pattern for files chosen individually
%
% OUTPUTS:
%    pathname       directory path where files reside
%    files          filenames of selected files

% Created by
%    Samantha Schwartz (originally); modified by Michael Wester (2021)

   if ~exist('txt', 'var')
      txt = '*.mat files';
   end
   if ~exist('pattern', 'var') | isempty(pattern)
      pattern = '*.mat';
   end

   [files, pathname] = uigetfile(fullfile(startdatadir, pattern), txt, ...
                                 'Multiselect', 'on');

   % The pathname is 0 when the user aborts the selction.
   if ~exist('pathname', 'var') | pathname == 0
      error('Directory pathname is empty!');
   end

   % Just make sure files is a cell (in case only one file was selected).
   if ~iscell(files)
      files = {files};
   end

end
