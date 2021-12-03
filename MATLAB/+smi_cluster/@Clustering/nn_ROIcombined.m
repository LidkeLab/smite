function nn_dists = nn_ROIcombined(obj, base_name, n_ROIs, RoI)
%nn_ROIcombined plots the mean particle nearest neighbor distances for a series
% of ROIs.
%
% INPUTS:
%    obj         various properties used by the algorithms
%       Font_props      [{'FontSize', 15, 'FontWeight', 'bold'}]
%       Fig_ext         ['png']  figure extension
%       ResultsDir      ['.']    directory to store results
%       Xlim            []       x-axis limits if defined
%       Ylim            []       y-axis limits if defined
%    base_name   name to identify saved plots, which can have various
%                descriptive suffixes attached 
%    n_ROIs      number of ROIs to combine
%    RoI         n_ROIs cell array containing the following fields:
%       X,Y{,Z}  x, y {, z} coordinates
%
% OUTPUT:
%    nn_dists    nearest nighbor distances over all particles for each ROI

% Created by
%    Michael J. Wester (2020)

   nn_dists = cell(n_ROIs, 1);
   if isfield(RoI{1}, 'Z') && ~isempty(RoI{1}.Z)
      for i = 1 : n_ROIs
         nn_dists{i} = ...
            Clustering.nn_distances([RoI{i}.X, RoI{i}.Y, RoI{i}.Z]);
      end
   else
      for i = 1 : n_ROIs
         nn_dists{i} = Clustering.nn_distances([RoI{i}.X, RoI{i}.Y]);
      end
   end

   if isempty(base_name)
      return;
   end

   % Histogram
   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on
   histogram(arrayfun(@(i) mean(nn_dists{i}), 1 : n_ROIs), 25);
   if ~isempty(obj.Xlim)
       xlim(obj.Xlim);
   end
   title(base_name);
   xlabel('nearest neighbor distance (nm)');
   ylabel('frequency');
   hold off
   name = fullfile(obj.ResultsDir, [base_name, '_nn_RC']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

   % PDF
   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on
   histogram(arrayfun(@(i) mean(nn_dists{i}), 1 : n_ROIs), 25, ...
             'Normalization', 'probability');
   ylim([0, 1]);
   if ~isempty(obj.Xlim)
       xlim(obj.Xlim);
   end
   if ~isempty(obj.Ylim)
       ylim(obj.Ylim);
   end
   title([base_name, ' [all ROIs]']);
   xlabel('nearest neighbor distance (nm)');
   ylabel('probability');
   hold off
   name = fullfile(obj.ResultsDir, [base_name, '_nn_PDF_RC']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

end
