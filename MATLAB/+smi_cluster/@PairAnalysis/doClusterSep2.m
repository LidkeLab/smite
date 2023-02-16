function results_cs = ...
   doClusterSep2(n_ROIs, results_c, desc, particles, results_dir, plotting)
% Find the nearest neighbor of each label 2 cluster to each label 1 cluster
% using center-to-center distances.
%
% INPUTS:
%    n_ROIs        number of ROIs
%    results_c     clustering results data structure produced by doClustering
%                  invoking smi_cluster.Clustering.clusterStats
%    desc          string identifying the analysis
%    particles     string array describing the two particles
%    results_dir   output directory
%    plotting      true if producing plots
%
% OUTPUTS:
%    results_cs      cluster separation between labels
%                       results_cs{1:n_ROIs}s
%                          .indx   index of nearest neighbor to each cluster
%                          .dist   NN cluster c2c distances
%       Also, figures *_cs_ROI*   histogram per ROI of cluster L1+L2 seps.
%             *_cs_RC_PDF/CDF     ROI combined PDF and CDF over each cell

% Created by
%    Michael J. Wester (2022)

   xtxt = sprintf('%s-%s cluster separation distances (nm)', ...
                  particles{1}, particles{2});
   results_cs = cell(n_ROIs, 1);
   dist_all = [];
   for i = 1 : n_ROIs
      centers1 = results_c{i}{1}.centers';
      centers2 = results_c{i}{2}.centers';
      [indx, dist] = knnsearch(centers1, centers2);
      results_cs{i}.indx = indx;
      results_cs{i}.dist = dist;
      dist_all = [dist_all; dist];

      if plotting
         figure('Visible', 'off');
         hold on
         histogram(results_cs{i}.dist, 'BinWidth', 50);
         xlabel(xtxt);
         ylabel('frequency');
         txt = sprintf('%s_cs_ROI%d', desc, i);
         title(regexprep(txt, '_', '\\_'));
         hold off
         saveas(gcf, fullfile(results_dir, txt));
         print(fullfile(results_dir, txt), '-dpng');
      end
   end

%  figure('Visible', 'off');
%  hold on
%  histogram(dist_all, 'BinWidth', 50);
%  xlabel(xtxt);
%  ylabel('frequency');
%  txt = sprintf('%s_cs_RC', desc);
%  title(regexprep(txt, '_', '\\_'));
%  hold off
%  saveas(gcf, fullfile(results_dir, txt));
%  print(fullfile(results_dir, txt), '-dpng');

   % Combined plot over all ROIs.
   SC = smi_cluster.StatisticsClustering();
   SC.BaseName = desc;
   SC.ResultsDir = results_dir;
   SC.PlotDo = 'pC';
   P = SC.plotCombined({dist_all}, 50, xtxt, {}, '_cs_RC');

end
