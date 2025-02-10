function [pearson, p_value] = pearsonCorrCoef(SMD1, SMD2, SRZoom, ROI)
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
%                 (default: 20)
%    ROI          for a ROI, corner coordinates [xmin, xmax, ymin, ymax] pixels
%                 (default: [0, 256, 0, 256] pixels)
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
   if ~exist('ROI', 'var')
      ROI = [0, 256, 0, 256];
   end

   SMD1a = SMD1;
   SMD2a = SMD2;
   if ROI ~= [0, 256, 0, 256]
      % If the ROI is not the whole image ... move the ROI to (x, y) = (0, 0)
      % corner.  This is done to perform the correlation only within the ROI,
      % excluding the zero pixels outside the ROI,
      SMD1a.X = SMD1a.X - ROI(1);
      SMD2a.X = SMD2a.X - ROI(1);
      SMD1a.Y = SMD1a.Y - ROI(3);
      SMD2a.Y = SMD2a.Y - ROI(3);
      % delta_x/y should be integer numbers of pixels.
      delta_x = floor(ROI(2) - ROI(1));
      delta_y = floor(ROI(4) - ROI(3));
      SMD1a.XSize = delta_x;
      SMD2a.XSize = delta_x;
      SMD1a.YSize = delta_y;
      SMD2a.YSize = delta_y;
   end

   % Make grayscale Gaussian blob images from the coordinates in SMD1 and SMD2.
   % 0s below are to omit the scalebar in the generated images.
   G1 = smi_vis.GenerateImages.grayscaleImage(SMD1a, SRZoom, 0);
   G2 = smi_vis.GenerateImages.grayscaleImage(SMD2a, SRZoom, 0);

   [r, p] = corrcoef(G1, G2);
   pearson = r(1, 2);
   p_value = p(1, 2);

   verbose = false;
   if verbose
      fprintf('pearson  = %f\n', pearson);

      % Code to calculate pearson by hand.  The results should be the same as
      % using the MATLAB function corrcoef.
      H1 = G1(:);
      H2 = G2(:);
      pearsonA = sum((H1 - mean(H1)) .* (H2 - mean(H2))) ...
                 ./ ((numel(H1) - 1) * std(H1) * std(H2));
      fprintf('pearsonA = %f\n', pearsonA);
   end

end
