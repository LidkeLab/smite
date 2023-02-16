function results_c = doClustering(n_ROIs, RoI, desc, results_dir, plotting, ...
                                  PixelSize, E, minPts)
% Clustering for each label in each ROI given epsilon (E) and minPts.
%
% INPUTS:
%    n_ROIs        number of ROIs
%    RoI           data structure containing x and y coordinates per ROI
%    desc          string identifying the analysis
%    results_dir   output directory
%    plotting      true if producing plots
%    PixelSize     conversion factor conversion factor (nm)
%    E             2-element array defining epsilon for clustering each label
%    minPts        2-element array defining min cluster size for each label
%
% OUTPUTS:
%    results_c       clustering of the 2 labels
%                       results_c{1:n_ROIs}{1:2}
%                          output of smi_cluster.Clustering.clusterStats
%                          .C         cluster indices
%                          .centers   coordinates of cluster centers
%                          .ptsI      indices of isolated points
%       Also, figures *_ROI_L1/2_algorithm

% Created by
%    Michael J. Wester (2022)

   %E = 50;      % epsilon: max distance between 2 adjacent points in a cluster
   %minPts = 3;  % minimum number of points in a cluster

   c = smi_cluster.Clustering();
   c.PixelSize = PixelSize;
   c.Plotting = true;
   c.Alpha = 2;
   c.Valgorithm = 2;
   c.ShrinkFactor = 0.5;   % used to make cluster boundaries convex or concave
   algorithm = 'DBSCAN';   % clustering algorithm
   options = 'O';

   XY = cell(2, 1);
   results_c = cell(n_ROIs, 1);
   for i = 1 : n_ROIs
      XY{1} = [ RoI{i}.X{1}, RoI{i}.Y{1} ];
      XY{2} = [ RoI{i}.X{2}, RoI{i}.Y{2} ];

      for j = 1 : 2
         [nC, C, centers, ptsI] = c.cluster(algorithm, XY{j}, E(j), minPts(j));
         fprintf('%s number of clusters ROI %d label %d = %d\n', ...
                 algorithm, i, j, nC);

         results_c{i}{j} = c.clusterStats(XY{j}, C, centers);
         results_c{i}{j}.C = C;
         results_c{i}{j}.centers = centers;
         results_c{i}{j}.ptsI = ptsI;

         if plotting
            txt = sprintf('%s ROI%d L%d', desc, i, j);
            clusterFig = c.plotClusters(XY{j}, C, centers, ptsI, txt, options);
            %showm(clusterFig);
            txt = sprintf('%s_ROI%d_L%d_%s', desc, i, j, algorithm);
            saveas(clusterFig, fullfile(results_dir, sprintf('%s.png', txt)));
         end
      end
   end

end
