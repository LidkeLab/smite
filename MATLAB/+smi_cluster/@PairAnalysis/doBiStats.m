function results_bi = ...
            doBiStats(n_ROIs, RoI, desc, particles, results_dir, combined)
% Pairwise mutual distances and bivariate Ripley's statistics for each ROI.
%
% INPUTS:
%    n_ROIs        number of ROIs
%    RoI           data structure containing x and y coordinates per ROI
%    desc          string identifying the analysis
%    particles     string array describing the two particles
%    results_dir   output directory
%    combined      produce combined bivariate Ripley over all ROIs
%
% OUTPUTS:
%    Figures *_ROI*_L1/2_pairwiseCDF/PDF    compared to a random distribution
%            *_ROI*_L1,L2_pairwisePDF2/CDF2 2-label PDF/CDF
%            *_ROI*_bivripley               bivariate Ripley

% Created by
%    Michael J. Wester (2022)

   results_bi = [];

   SC = smi_cluster.StatisticsClustering();
   SC.ResultsDir = results_dir;

   if ~combined
      %results_birip = cell(n_ROIs, 1);
      for i = 1 : n_ROIs
         SC.BaseName = sprintf('%s_ROI%d', desc, i);
         % Compare the actual distribution of nearest neighbor distances to a
         % random distribution per ROI.
         SC.pairwiseDist(particles, RoI);
         % Produce a plot of the distribution of nearest neighbor distances
         % between the localizations in the two labels.
         SC.pairwiseMutualDist(particles, RoI);
         % Bivariate Ripley plot of the localizations per ROI.
         SC.bivariateRipley(particles, RoI);
      end
   else
      % Combined bivariate Ripley over all ROIs.
      SC.BaseName = desc;
      SC.bivariateRipley_ROIcombined(particles, n_ROIs, RoI);
   end

end
