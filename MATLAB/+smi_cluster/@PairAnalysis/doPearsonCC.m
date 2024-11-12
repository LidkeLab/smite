function [results] = pearsonCorrCoef(particles, n_ROIs, RoI, desc, results_dir)
%pearsonCorrCoef finds the Pearson correlation coefficient between two colors.
% pearsonCorrCoef finds the Pearson correlation coefficient between two sets of
% localizations given in a list of coincident ROIs when converted to histogram
% images.  The p-value that tests the hypothesis of no correlation is also
% computed.
%
% INPUTS:
%    particles     string array describing the two particles
%    n_ROIs        number of ROIs
%    RoI           data structure containing x and y coordinates per ROI
%    desc          string identifying the analysis
%    results_dir   output directory
%                 
% OUTPUTS:
%    results       results structure imdexed by ROI with the following fields:
%       pearson    correlation coefficent between the histogram images of the
%                  two datasets
%       p_value    p-value testing the hypothesis of no correlation.  The
%                  p-value is the probability of getting a correlation as
%                  large as the observed value by random chance, when the
%                  true correlation is zero

% Created by
%    Michael Wester (2024)

   SRZoom = 20;   % coordinate magnification factor
   results = cell(n_ROIs, 1);

   do i = 1 : n_ROIs
      SMD1.X = RoI{i}.X(:, 1);
      SMD1.X = RoI{i}.Y(:, 1);
      SMD2.X = RoI{i}.X(:, 2);
      SMD2.X = RoI{i}.Y(:, 2);
      [pearson, p_value] = smi_stat.pearsonCorrCoef(SMD1, SMD2, SRZoom);
      results{i}.pearson = pearson;
      results{i}.p_value = p_value;
   end

   out = fopen(fullfile(results_dir, sprintf("_%s_Pearson.txt", "w")));
   fprintf(out, "ROI pearson p_value   (%s, %s)\n", desc{1}, desc{2});
   do i = 1 : n_ROIs
      fprintf(out, "%2d %7.4f %g\n", results{i}.pearson, results{i}.p_value);
   end
   fclose(out);
end
