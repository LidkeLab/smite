function [n_ROIs, RoI] = filterROIs(pathname, files, filter)
% ---------- Possibly, filter some ROIs out
% Filter out ROIs that meet certain criteria on their contents.
%
% INPUTS:
%    pathname          path to where the files are located
%    files             _ROIs.mat files to collect ROI info from
%    filter            structure that contains allowable limits
%       minLocROI      minimum number of localizations allowed in a ROI
%       maxLocROI      maximum number of localizations allowed in a ROI
%
% OUTPUTS:
%    _filtered_ROIs.mat files containing updated n_ROIs and RoI structure

   n_files = numel(files);
   for i = 1 : n_files
      data = load(fullfile(pathname, files{i}));
      toDelete = [];
      for j = 1 : data.n_ROIs
         n_locs = numel(data.RoI{j}.X{1});
         % Collect ROIs to delete.
         if (isfield(filter, 'minLocROI') && n_locs < filter.minLocROI) || ...
            (isfield(filter, 'maxLocROI') && n_locs > filter.maxLocROI)
            toDelete = [toDelete, j];
         end
      end
      if ~isempty(toDelete)
         fprintf('Deleting ROIs');
         for k = 1 : numel(toDelete)
            fprintf(' %d', toDelete(k));
         end
         fprintf(' from %s\n', files{i});
      end
      % Delete ROIs that exceeded the limits in filter.
      n_ROIs = data.n_ROIs - numel(toDelete);
      RoI = data.RoI;
      RoI(toDelete) = [];
      % Generate a filtered _ROIs.mat file.
      filename = files{i};
      filename = regexprep(filename, '_ROI', '_filtered_ROI');
      ResultsFile = data.ResultsFile;
      Pixel2nm    = data.Pixel2nm;
      XYsize      = data.XYsize;
      save(fullfile(pathname, filename), ...
           'ResultsFile', 'Pixel2nm', 'XYsize', 'n_ROIs', 'RoI');
   end

   fprintf('Done filtering ROIs.\n');

end
