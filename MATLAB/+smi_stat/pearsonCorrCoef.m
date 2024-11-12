function [pearson, p_value] = pearsonCorrCoef(SMD1, SMD2, SRZoom)
%pearsonCorrCoef finds the Pearson correlation coefficient between two SMDs.
% pearsonCorrCoef finds the Pearson correlation coefficient between two sets of
% localizations given in SMD1 and SMD2 when converted to Gaussian blob images.
% p-value that tests the hypothesis of no correlation is also returned.
%
% INPUTS:
%    SMD1, SMD2   single molecule data structures containing 2D localization
%                 coordinates in fields X, Y as well as the image sizes in
%                 XSize and YSize (the two datasets are assumed to have come
%                 from the same size images).  If XSize and YSize are not
%                 provided, they will be estimated from the localizations
%    SRZoom       magnification factor for the SMD coordinates in order to get
%                 a better estimate of the Pearson Correlation Coefficient
%                 (default: 29)
%
% OUTPUTS:
%    pearson      correlation coefficent between the Gaussian blob images of
%                 the two datasets
%    p_value      p-value testing the hypothesis of no correlation.  The
%                 p-value is the probability of getting a correlation as large
%                 as the observed value by random chance, when the true
%                 correlation is zero

% Created by
%    Michael Wester (2024)

   if ~exist('SRZoom', 'var')
      SRZoom = 20;
   end

   % Make grayscale Gaussian blob images from the coordinates in SMD1 and SMD2.
   G1 = smi_vis.GenerateImages.grayscaleImage(SMD1, SRZoom, 0);
   G2 = smi_vis.GenerateImages.grayscaleImage(SMD2, SRZoom, 0);

   [r, p] = corrcoef(G1, G2);
   pearson = r(1, 2);
   p_value = p(1, 2);
   %fprintf('pearson  = %f\n', pearson);

   %H1 = H1(:);
   %H2 = H2(:);
   %pearsonA = sum((H1 - mean(H1)) .* (H2 - mean(H2))) ...
   %           ./ ((numel(H1) - 1) * std(H1) * std(H2));
   %fprintf('pearsonA = %f\n', pearsonA);

end
