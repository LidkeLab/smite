function results_ls = ...
   doLocSep2(n_ROIs, RoI, desc, particles, results_dir, plotting)
% Find the nearest neighbor of each label 2 localization to each label 1
% localization.
%
% INPUTS:
%    n_ROIs        number of ROIs
%    RoI           data structure containing x and y coordinates per ROI
%    desc          string identifying the analysis
%    particles     string array describing the two particles
%    results_dir   output directory
%    plotting      true if producing plots
%
% OUTPUTS:
%    results_ls      localization separation between labels
%                       results_ls{1:n_ROIs}s
%                          .indx   index of nearest neighbor to each loc.
%                          .dist   NN localization distances
%       Also, figures *_ls_ROI*   histogram per ROI of loc. L1+L2 separations
%             *_ls_RC_PDF/CDF     ROI combined PDF and CDF over each cell

% Created by
%    Michael J. Wester (2022)

   xtxt = sprintf('%s-%s localization separation distances (nm)', ...
                  particles{1}, particles{2});
   results_ls = cell(n_ROIs, 1);
   bin_width = 5;
   dist_all = [];
   for i = 1 : n_ROIs
      label1 = [RoI{i}.X{1}, RoI{i}.Y{1}];
      label2 = [RoI{i}.X{2}, RoI{i}.Y{2}];
      [indx, dist] = knnsearch(label1, label2);
      results_ls{i}.indx = indx;
      results_ls{i}.dist = dist;
      dist_all = [dist_all; dist];

      if plotting
         figure('Visible', 'off');
         hold on
         histogram(results_ls{i}.dist, 'BinWidth', bin_width);
         xlabel(xtxt);
         ylabel('frequency');
         txt = sprintf('%s_ls_ROI%d', desc, i);
         title(regexprep(txt, '_', '\\_'));
         hold off
         saveas(gcf, fullfile(results_dir, txt));
         print(fullfile(results_dir, txt), '-dpng');
      end
   end

   % Combined plot over all ROIs.
   SC = smi_cluster.StatisticsClustering();
   SC.BaseName = desc;
   SC.ResultsDir = results_dir;
   SC.PlotDo = 'pC';
   P = SC.plotCombined({dist_all}, bin_width, xtxt, {}, '_ls_RC');

end
