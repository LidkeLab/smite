function results_ss = doSimpleStats(n_ROIs, RoI, PixelSize, desc, particles,...
                                   results_dir)
%doSimpleStats finds the Pearson and Manders' coefficients between two colors.
% pearsonCorrCoef finds the Pearson correlation coefficient between two sets of
% localizations given in a list of coincident ROIs when converted to histogram
% images.  The p-value that tests the hypothesis of no correlation is also
% computed.
% Manders' coefficients M1 and M2 measure colocalization by computing the
% intensity fraction of pixels from one label that overlap with pixels from
% another.
%
% INPUTS:
%    n_ROIs        number of ROIs
%    RoI           data structure containing x and y coordinates per ROI
%    PixelSize     nm per pixel conversion factor
%    desc          string identifying the analysis
%    particles     string array describing the two particles
%    results_dir   output directory
%                 
% OUTPUTS:
%    results_ss    results structure imdexed by ROI with the following fields:
%       pearson    correlation coefficent between the histogram images of the
%                  two datasets
%       p_value    p-value testing the hypothesis of no correlation.  The
%                  p-value is the probability of getting a correlation as
%                  large as the observed value by random chance, when the
%                  true correlation is zero
%       M1         intensity fraction of pixels in label 1 that overlap with
%                  pixels in label 2
%       M2         intensity fraction of pixels in label 2 that overlap with
%                  pixels in label 1

% Created by
%    Michael Wester (2024)

   SRZoom = 20;   % coordinate magnification factor
   results_ss = cell(n_ROIs, 1);

   out = fopen(fullfile(results_dir, sprintf('%s_stats.txt', desc)), 'w');
   fprintf(out, 'PixelSize = %f nm/pixel, SRZoom = %g\n\n', PixelSize, SRZoom);

   fprintf(out, 'ROI pearson p_value  M1      M2     (%s, %s)\n', ...
                particles{1}, particles{2});
   for i = 1 : n_ROIs
      SMD1 = RoI{i}.SMD{1};
      SMD2 = RoI{i}.SMD{2};
      ROI  = RoI{i}.ROI / PixelSize;

      % Pearson correlation coefficient.
      [pearson, p_value] = smi_stat.pearsonCorrCoef(SMD1, SMD2, SRZoom, ROI);
      results_ss{i}.pearson = pearson;
      results_ss{i}.p_value = p_value;

      % Manders' split coefficients.
      [M1, M2] = smi_stat.mandersSplitCoefs(SMD1, SMD2, SRZoom);
      results_ss{i}.M1 = M1;
      results_ss{i}.M2 = M2;

      fprintf(out, '%2d  %7.4f %7.2g %7.4f %7.4f\n', ...
                   i, pearson, p_value, M1, M2);
   end
   fclose(out);

end
