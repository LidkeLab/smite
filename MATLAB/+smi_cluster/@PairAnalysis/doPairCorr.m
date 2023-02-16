function [results_pcc, resultsRC_pcc] = ...
   doPairCorr(n_ROIs, RoI, ROI_sizes, desc, results_dir, combined, ...
              plotting, HistBinSize, RmaxAxisLimit)
% Pair cross-correlation for each ROI and combined ROIs.
%
% INPUTS:
%    n_ROIs        number of ROIs
%    RoI           data structure containing x and y coordinates per ROI
%    ROI_sizes     ROI sizes (nm)
%    desc          string identifying the analysis
%    results_dir   output directory
%    combined      analysis is combined over all ROIs
%    plotting      true if producing plots
%    HistBinSize   number of pixels per bin to collect correlation statistics
%    RmaxAxisLimit sets r axis limit for pair correlation plots if > 0 (nm)
%
% OUTPUTS:
%    results_pcc     pair cross-correlation
%                       results_pcc{1:n_ROIs}
%                          .G       2D pair-correlation function values
%                          .r       radius values (nm)
%                          .g       angularly averaged pair corr. function
%                          .dg      errors on angularly averaged g
%                          .model   model calculated at estimated value
%                          ....     various model results
%       Also, figures *_ROI*_crosscorr (ROIwise pairwise cross-correlations)
%    resultsRC_pcc   ROIs combined pair cross-correlation
%                       results_pcc{1:n_ROIs}
%                          see results_pcc (ROI combined pairwise_crosscorr)
%       Also, figures *_RC_crosscorrR

% Created by
%    Michael J. Wester (2022)

   pc = smi_cluster.PairCorrelation();

   pc.ResultsDir = results_dir;    % results directory
   pc.Fig_ext = 'png';             % figure extension
   pc.Rmax_axis = RmaxAxisLimit;   % sets plotting limit if > 0
   % Histogram bin size for pairwise correlation---this is the number of pixels
   % per bin over which correlation statistics are collected.
   pc.HistBinSize = HistBinSize;

   ROIs_combined = cell(n_ROIs, 1);
   results_pcc   = cell(n_ROIs, 1);
   resultsRC_pcc = [];
   for i = 1 : n_ROIs
      if plotting
         fprintf('ROI %d\n\n', i);
      end

      txt = sprintf('%s_ROI%d', desc, i);
      pc.BaseName = txt;
      pc.ROI = RoI{i}.ROI;
      XY1 = [ RoI{i}.X{1}, RoI{i}.Y{1} ];
      XY2 = [ RoI{i}.X{2}, RoI{i}.Y{2} ];

      ROIs_combined{i} = RoI{i};

      if plotting || combined
         results_pcc{i} = pc.pair_correlation(XY1, XY2);
      end
   end

   % Combined plot over all the ROIs.
   if combined && n_ROIs > 0
      fprintf('ROI combined\n\n');
      txt = sprintf('%s_RC', desc);
      pc.BaseName = txt;
      resultsRC_pcc = ...
         pc.pair_correlation_ROIcombined(2, n_ROIs, ROIs_combined, 1)
   end

   if combined
      gmax = zeros(n_ROIs, 1);
      den1 = zeros(n_ROIs, 1);
      den2 = zeros(n_ROIs, 1);
      ROI_area = prod(ROI_sizes);
      for i = 1 : n_ROIs
         gmax(i) = max(results_pcc{i}.g);
         den1(i) = numel(ROIs_combined{i}.X{1}) ./ ROI_area;
         den2(i) = numel(ROIs_combined{i}.X{2}) ./ ROI_area;
      end

      h = figure;
      axes('FontSize', 15, 'FontWeight', 'bold');
      hold on
      plot(den1, gmax, 'k.', 'MarkerSize', 10);
      xlabel('density Label 1 (#/nm^2)');
      ylabel('g_{max}');
      txt = sprintf('%s_gmax_den1', desc);
      title(regexprep(txt, '_', '\\_'));
      hold off
      saveas(h, fullfile(results_dir, sprintf('%s.png', txt)));
      close

      h = figure;
      axes('FontSize', 15, 'FontWeight', 'bold');
      hold on
      plot(den2, gmax, 'k.', 'MarkerSize', 10);
      xlabel('density Label 2 (#/nm^2)');
      ylabel('g_{max}');
      txt = sprintf('%s_gmax_den2', desc);
      title(regexprep(txt, '_', '\\_'));
      hold off
      saveas(h, fullfile(results_dir, sprintf('%s.png', txt)));
      close
   end

end
