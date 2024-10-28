function [pearson, p_value] = pearsonCorrCoef(SMD1, SMD2)
%pearsonCorrCoef finds the Pearson correlation coefficient between two SMDs.
% pearsonCorrCoef finds the Pearson correlation coefficient between two sets of
% localizations given in SMD1 and SMD2 when converted to histogram images.  Thw
% p-value that tests the hypothesis of no correlation is also returned.
%
% INPUTS:
%    SMD1, SMD2   single molecule data structures containing 2D localization
%                 coordinates in fields X, Y as well as the image sizes in
%                 XSize and YSize (the two datasets are assumed to have come
%                 from teh same size images)
%                 
% OUTPUTS:
%    pearson      correlation coefficent between the histogram images of the
%                 two datasets
%    p_value      p-value testing the hypothesis of no correlation.  The
%                 p-value is the probability of getting a correlation as large
%                 as the observed value by random chance, when the true
%                 correlation is zero

% Created by
%    Michael Wester (2024)

   % Make histogram images from the coordinates in SMD1 and SMD2.
   nx = SMD1.XSize;   % assumed the same as SMD2.XSize
   ny = SMD1.YSize;   % assumed the same as SMD2.YSize
   H1 = zeros(nx, ny, 'single');
   H2 = zeros(nx, ny, 'single');
   n1 = numel(SMD1.X);   % number of localizations
   n2 = numel(SMD2.X);   % number of localizations
   for l = 1 : n1
      ii = ceil(SMD1.X(l));
      jj = ceil(SMD1.Y(l));
      if ii >= 1 && ii <= nx && jj >= 1 && jj <= ny
         H1(ii, jj) = H1(ii, jj) + 1;
      end
   end
   for l = 1 : n2
      ii = ceil(SMD2.X(l));
      jj = ceil(SMD2.Y(l));
      if ii >= 1 && ii <= nx && jj >= 1 && jj <= ny
         H2(ii, jj) = H2(ii, jj) + 1;
      end
   end

   [r, p] = corrcoef(H1, H2);
   pearson = r(1, 2);
   p_value = p(1, 2);

end
