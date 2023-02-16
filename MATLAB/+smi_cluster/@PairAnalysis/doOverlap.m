function results_o = doOverlap(n_ROIs, RoI, results_c, l12, desc,      ...
                               particles, results_dir, PlotNonOverlap, ...
                               Color, plotting)
% Find overlaps between label 1/2 clusters and label 2/1 localizations.
%
% INPUTS:
%    n_ROIs        number of ROIs
%    RoI           data structure containing x and y coordinates per ROI; see
%                  defineROIs2
%    results_c     results from ckustering (see doClustering)
%    l12           true  (1) for label 1 clusters/label 2 localizations
%                  false (0) for label 2 clusters/label 1 localizations
%    desc          string identifying the analysis
%    particles     string array describing the two particles
%    results_dir   output directory
%    PlotNonOverlap if true, plot non-overlapping locs as black dots, otherwise
%                  omit them from the plot
%    Color         label 1 and label 2 colors on display
%    plotting      true if producing plots
%
% OUTPUTS:
%    results_o       overlap between label l1 cluster/label l2 localizations
%                       results_o{1:n_ROIs}
%                          .CluLabel   clusters label
%                          .LocLabel   localizations label
%                          .Ltotal     total Loclabel localizations in ROI
%                          .Lin        overlap locs. wrt CluLabel clusters
%       Also, figures *_cLl1_lLl2 Color(l1) Ll1 clusters, Color(l2) Ll2 locs.
%       where l1 identifies clusters and l2 identifies localizations (see
%       l12 in the inputs).

% Created by
%    Michael J. Wester (2022)

   results_o = cell(n_ROIs, 1);
   if l12
      l1 = 1;   % clusters label
      l2 = 2;   % localizations label
   else
      l1 = 2;   % clusters label
      l2 = 1;   % localizations label
   end
   for i = 1 : n_ROIs
      ROI = RoI{i}.ROI;

      % Collect label l1 cluster information for ROI i.
      nC    = results_c{i}{l1}.nC;   % number of clusters found
      C     = results_c{i}{l1}.C;    % indices of cluster points per cluster
      xyROI = results_c{i}{l1}.XY;   % (x, y) coordinates of locs in ROI i
      ihull = results_c{i}{l1}.indices_hull;   % boundary indices per cluster
      % Note that the last boundary index is the same as the first for drawing
      % purposes, so the number of actual boundary vertices is one less per
      % cluster.
      nBoundaryPtsPerC = cellfun(@numel, ihull) - 1;   % perC = per Cluster
      % Construct xv, yv from all the cluster boundary points.  This will be
      % arranged like:
      %    [b11, b12, b13, NaN, b21, b22, b23, b24];
      % where bij is the j boundary point of the ith cluster.  NaN's separate
      % clusters.
      xyv = zeros(sum(nBoundaryPtsPerC) - 1, 2);
      j = 1;
      k = 0;
      nB = nBoundaryPtsPerC(j);
      xyv(k + 1 : k + nB, :) = xyROI(C{j}(ihull{j}(1 : nB)), :);
      k = k + nBoundaryPtsPerC(j);
      for j = 2 : nC
         k = k + 1;
         xyv(k, :) = NaN;
         nB = nBoundaryPtsPerC(j);
         xyv(k + 1 : k + nB, :) = xyROI(C{j}(ihull{j}(1 : nB)), :);
         k = k + nBoundaryPtsPerC(j);
      end
      xv = xyv(:, 1);
      yv = xyv(:, 2);

      % Collect label l2 localizations inside ROI i.
      xq = RoI{i}.X{l2};
      yq = RoI{i}.Y{l2};

      % Here's the magic that does the hard work.
      in = inpolygon(xq, yq, xv, yv);
      n_total = numel(in);
      n_in = sum(in);

      results_o{i}.CluLabel = l1;      % clusters label
      results_o{i}.LocLabel = l2;      % localizations label
      results_o{i}.Ltotal = n_total;   % total l2 localizations in ROI i
      results_o{i}.Lin = n_in;         % overlapping l2 locs wrt l1 clusters

      if plotting
         figure;
         axes('FontSize', 15, 'FontWeight', 'bold');
         hold on
         % Plot the ROI boundaries.
         plot([ROI(1), ROI(2), ROI(2), ROI(1), ROI(1)], ...
              [ROI(3), ROI(3), ROI(4), ROI(4), ROI(3)], 'r-', 'LineWidth', 2);
         % Plot the l1 clusters.
         for j = 1 : nC
            xyc = xyROI(C{j}(ihull{j}), :);
            plot(xyc(:, 1), xyc(:, 2), [Color(l1), '-']);
         end
         % Plot l2 localizations outside clusters.
         if PlotNonOverlap
            plot(xq(~in), yq(~in), 'k.');
         end
         % Plot l2 localizations inside clusters (or on their boundary).
         plot(xq(in),  yq(in),  [Color(l2), 'x']);
         xlabel('x (nm)');
         ylabel('y (nm)');
         title(sprintf(['ROI %d L%d clusters/L%d locs, ', ...
                        'overlap %d/%d = %.2f%%'],        ...
                       i, l1, l2, n_in, n_total, 100 * n_in/n_total));
         hold off

         txt = sprintf('%s_ROI%d_cL%d_lL%d', desc, i, l1, l2);
         print(gcf, fullfile(results_dir, sprintf('%s.png', txt)), ...
                    '-dpng', '-r1200');
         savefig(gcf, fullfile(results_dir, sprintf('%s.fig', txt)), ...
                      'compact');
         close
      end
   end

end
