function [results_pcc, resultsRC_pcc, results_ss, results_c, results_cs, ...
          results_ls, results_o1, results_o2] =                          ...
   doAnalysis(n_ROIs, RoI, ROI_sizes, desc, particles, results_dir,      ...
              options, PixelSize, HistBinSize, RmaxAxisLimit, E, minPts, ...
              PlotNonOverlap, Color)
% Dispatches analyses to various helper functions depending options provided.
%
% INPUTS:
%    n_ROIs        number of ROIs
%    RoI           data structure containing x and y coordinates per ROI
%    ROI_sizes     ROI sizes (nm)
%    desc          string identifying the analysis
%    particles     string array describing the two particles
%    results_dir   output directory
%    options       string array of strings specifying analyses to be performed:
%                      'combined'    the analysis is combined over all ROIs
%                      'plotting'    plots are to be produced.
%                      'SimpleStats' Pearson's correlation and Manders' split
%                                    coefficients per ROI
%                      'BiStats'     pairwise mutual distances and bivariate
%                                    Ripley's statistics for each ROI
%                      'Clustering'  clusters for each label per ROI
%                      'Clustering2' C2C nearest neighbor distances between
%                                    label 1/label 2 clusters per ROI and
%                                    combined over all ROIs
%                      'LocSep2'     NN localization distances between L1 & L2
%                      'Overlap1'    overlap between L1 clusters/L2 locs.
%                      'Overlap2'    overlap between L2 clusters/L1 locs.
%                      'PairCorr'    pair correlation per ROI and combined over
%                                    all ROIs
%                      'Plot2'       2D plot per ROI
%    PixelSize     conversion factor conversion factor (nm)
%    HistBinSize   number of pixels per bin to collect correlation statistics;
%                  a good guess is the PixelSize
%    RmaxAxisLimit sets r axis limit for pair correlation plots if > 0 (nm)
%    E             2-element array defining epsilon for clustering each label
%    minPts        2-element array defining min cluster size for each label
%    PlotNonOverlap if true, plot non-overlapping locs as black dots, otherwise
%                  omit them from the plot
%    Color         label 1 and label 2 colors on display
%
% OUTPUTS:
%    desc_results.mat containing the various results_ cell arrays for one cell
%    Various results containers from the called helper functions in case the
%    user wants to have more details:
%       results_bi      bivariate statistics
%          Also, figures *_ROI*_L1/2_pairwiseCDF/PDF compared to a random dist.
%                *_ROI*_L1,L2_pairwisePDF2/CDF2      2-label PDF/CDF
%                *_ROI*_bivripley                    bivariate Ripley
%       results_pcc     pair cross-correlation
%                          results_pcc{1:n_ROIs}
%                             .G       2D pair-correlation function values
%                             .r       radius values (nm)
%                             .g       angularly averaged pair corr. function
%                             .dg      errors on angularly averaged g
%                             .model   model calculated at estimated value
%                             ....     various model results
%          Also, figures *_ROI*_crosscorr (ROIwise pairwise cross-correlations)
%       resultsRC_pcc   ROIs combined pair cross-correlation
%                          results_pcc{1:n_ROIs}
%                             see results_pcc (ROI combined pairwise_crosscorr)
%          Also, figures *_RC_crosscorrR
%       results_c       clustering of the 2 labels
%                          results_c{1:n_ROIs}{1:2}
%                             output of smi_cluster.Clustering.clusterStats
%                             .C         cluster indices
%                             .centers   coordinates of cluster centers
%                             .ptsI      indices of isolated points
%          Also, figures *_ROI*_L1/2_algorithm
%       results_cs      cluster separation between labels
%                          results_cs{1:n_ROIs}s
%                             .indx   index of nearest neighbor to each cluster
%                             .dist   NN cluster c2c distances
%          Also, figures *_cs_ROI*   histogram per ROI of cluster L1+L2 seps.
%                *_cs_RC_PDF/CDF     ROI combined PDF and CDF over each cell
%       results_ls      nearest neighbor localization separation between labels
%                          results_ls{1:n_ROIs}s
%                             .indx   index of nearest neighbor to each cluster
%                             .dist   NN localization distances
%          Also, figures *_ls_ROI*   histogram per ROI of loc. L1+L2 seps.
%                *_ls_RC_PDF/CDF     ROI combined PDF and CDF over each cell
%       results_o1      overlap between label 1 cluster/label 2 localizations
%                          results_o1{1:n_ROIs}
%                             .CluLabel   clusters label
%                             .LocLabel   localizations label
%                             .Ltotal     total Loclabel localizations in ROI
%                             .Lin        overlap locs. wrt CluLabel clusters
%          Also, figures *_cL1_lL2 Color(1) L1 clusters, Color(2) L2 locs.
%       results_o2      overlap between label 2 cluster/label 1 localizations
%                          results_o2{1:n_ROIs}
%                              as resultso1
%          Also, figures *_cL2_lL1 Color(1) L2 clusters, Color(2) L1 locs.
%                        plotting
%          Also, figures *_ROI*_L1+L2   Color(1) L1 and Color(2) L2 locs.

% Created by
%    Michael J. Wester (2022)

   PA = smi_cluster.PairAnalysis();

   % combined   indicates the analysis is combined over all ROIs.
   % plotting   indicates if plots are to be produced.
   combined = false;
   plotting = false;

   if any(contains(options, "combined"))
      combined = true;
   end
   if any(contains(options, "plotting"))
      plotting = true;
   end

   results_ss = [];
   results_bi = [];
   results_pcc = [];
   resultsRC_pcc = [];
   results_c = [];
   results_cs = [];
   results_ls = [];
   results_o1 = [];
   results_o2 = [];

   % Pearson's correlation and Manders' split coefficients per ROI.
   if any(contains(options, "SimpleStats"))
      results_ss = PA.doSimpleStats(n_ROIs, RoI, PixelSize, desc, particles, ...
                                    results_dir);
      fprintf("Done SimpleStats\n");
   end

   % Pairwise mutual distances and bivariate Ripley's per ROI and the latter
   % also combined over all ROIs.
   if any(contains(options, "BiStats"))
      results_bi = PA.doBiStats(n_ROIs, RoI, desc, particles, results_dir, ...
                                combined);
      fprintf("Done BiStats\n");
   end

   % Clusters for each label per ROI.
   if any(contains(options, "Clustering"))
      results_c = PA.doClustering(n_ROIs, RoI, desc, results_dir, plotting,...
                                  PixelSize, E, minPts);
      fprintf("Done Clustering\n");
   end

   % C2C nearest neighbor distances between label 1/label 2 clusters per ROI
   % and combined over all ROIs.
   if any(contains(options, "Clustering2"))
      results_cs = PA.doClusterSep2(n_ROIs, results_c, desc, particles, ...
                                    results_dir, plotting);
      fprintf("Done Clustering2\n");
   end

   % Nearest neighbor distances between label 1/label 2 localizations per ROI
   % and combined over all ROIs.
   if any(contains(options, "LocSep2"))
      results_ls = PA.doLocSep2(n_ROIs, RoI, desc, particles, results_dir, ...
                                plotting);
      fprintf("Done LocSep2\n");
   end

   % Overlaps between label 1 clusters and label 2 localizations.
   if any(contains(options, "Overlap1"))
      l12 = true;
      results_o1 = PA.doOverlap(n_ROIs, RoI, results_c, l12, desc,      ...
                                particles, results_dir, PlotNonOverlap, ...
                                Color, plotting);
      fprintf("Done Overlap1\n");
   end

   % Overlaps between label 2 clusters and label 1 localizations.
   if any(contains(options, "Overlap2"))
      l12 = false;
      results_o2 = PA.doOverlap(n_ROIs, RoI, results_c, l12, desc,      ...
                                particles, results_dir, PlotNonOverlap, ...
                                Color, plotting);
      fprintf("Done Overlap2\n");
   end

   % Pair correlation per ROI and combined over all ROIs.
   if any(contains(options, "PairCorr"))
      [results_pcc, resultsRC_pcc] =                                        ...
         PA.doPairCorr(n_ROIs, RoI, ROI_sizes, desc, results_dir, combined, ...
                       plotting, HistBinSize, RmaxAxisLimit);
      fprintf("Done PairCorr\n");
   end

   % 2D plot per ROI.
   if any(contains(options, "Plot2"))
      PA.doPlot2(n_ROIs, RoI, desc, particles, results_dir, Color, plotting);
      fprintf("Done Plot2\n");
   end

   % Save results.
   save(fullfile(results_dir, sprintf('%s_results.mat', desc)), 'n_ROIs', ...
        'RoI', 'results_ss', 'results_c', 'results_cs',  'results_ls',    ...
        'results_o1', 'results_o2', 'results_pcc', 'resultsRC_pcc');
   fprintf("Done saving results\n");

end
