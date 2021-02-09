function [dist1, rmse1, dist2, rmse2, nnfig] = ...
   calcDCRMSE(SMD, X_True, Y_True, Z_True,     ...
              DriftX_True, DriftY_True, DriftZ_True)
% calcDCRMSE calculates the RMSE of SMD relative to true coordinates/curves.
% (1) considers differences in localization coordinates wrt the true values.
% (2) considers differences in drift curves wrt to the true values.
% (1) and (2) are typically close in value.
%
% INPUTS:
%   SMD              A structure with fields:
%      X                  x-coordinates (Nx1) where N is total number of points
%      Y                  y-coordinates (Nx1) [pixels, also for X]
%      Z                  z-coordinates (Nx1) [um]
%      DriftX             x-drifts found (Nframes x Ndatasets) [pixels]
%      DriftY             y-drifts found (Nframes x Ndatasets) [pixels]
%      DriftZ             z-drifts found (Nframes x Ndatasets) [um]
%      PixelSizeZUnit     um per pixel
%   X_True           true x-coordinates (Nx1) [pixels]
%   Y_True           true y-coordinates (Nx1) [pixels]
%   Z_True           true z-coordinates (Nx1) [um]
%   DriftX_True      true x-drift curve (Nframes x Ndatasets) [pixels]
%   DriftY_True      true y-drift curve (Nframes x Ndatasets) [pixels]
%   DriftZ_True      true z-drift curve (Nframes x Ndatasets) [um]
% OUTPUTS:
%   dist1   average distance between true and drift corrected coordinates (nm)
%   dist2   average distance between true and drift corrected curves (nm)
%   rmse1   root-mean-square-error between true and drift corrected coords (nm)
%   rmse2   root-mean-square-error between true and drift corrected curves (nm)
%   nnfig   histogram of nearest neighbor distances between drift corrected and
%           true coordinates
%
% Created by:
%   Michael J. Wester (Lidke Lab 2020)

   SRImageZoom = 10;

   % Remove any NaNs.
   nans = find(isnan(SMD.X) | isnan(SMD.Y) | isnan(X_True) | isnan(Y_True));
   n_nans = numel(nans);
   if n_nans > 0
      fprintf('calcDCRMSE: %d NaNs removed!\n', n_nans);
      SMD.X(nans) = [];
      SMD.Y(nans) = [];
      X_True(nans) = [];
      Y_True(nans) = [];
   end
   N = numel(SMD.X);
   Nframes = SMD.Nframes;

   % --------------------------------------------------------------------------

%  % Drift corrected image.
%  DCimage = SMA_DriftCorrect.histogramImage(SMD, SRImageZoom);

%  % True image.
%  SMD_True = SMD;
%  SMD_True.X = X_True;
%  SMD_True.Y = Y_True;
%  TrueImage = SMA_DriftCorrect.histogramImage(SMD_True, SRImageZoom);

%  % Reconstructed original image.
%  X_unDC = zeros(N, 1, 'single');
%  Y_unDC = zeros(N, 1, 'single');
%  for k = 1:N
%     i = SMD.FrameNum(k);
%     j = SMD.DatasetNum(k);
%     X_unDC(k) = SMD.X(k) + SMD.DriftX(i, j);
%     Y_unDC(k) = SMD.Y(k) + SMD.DriftY(i, j);
%  end
%  SMD_unDC = SMD;
%  SMD_unDC.X = X_unDC;
%  SMD_unDC.Y = Y_unDC;
%  unDCimage = SMA_DriftCorrect.histogramImage(SMD_unDC, SRImageZoom);

%  % Residual image.
%  residual = abs(DCimage - TrueImage);

   % --------------------------------------------------------------------------

   P2nm = SMD.PixelSizeZUnit * 1000;   % nm per pixel

   % Drift corrected and true coordinates (nm).
   XX = SMD.X * P2nm;
   YY = SMD.Y * P2nm;
   XX_True = X_True * P2nm;
   YY_True = Y_True * P2nm;
   if ~isempty(Z_True)
      ZZ      = SMD.Z  * 1000;
      ZZ_True = Z_True * 1000;
   end

   % Average distance between the true and drift corrected coordinates.
   if ~isempty(Z_True)
      dists1 = sqrt((XX_True - XX).^2 + (YY_True - YY).^2 + (ZZ_True - ZZ).^2);
   else
      dists1 = sqrt((XX_True - XX).^2 + (YY_True - YY).^2);
   end
   dist1 = sum(dists1) / N;
   % RMSE1 between the true and drift corrected coordinates (nm).
   rmse1 = sqrt(sum(dists1.^2) / N);

   % Histogram of nearest neighbor distances between the true and the drift
   % corrected coordinates.
   if ~isempty(Z_True)
      [~, D] = knnsearch([XX_True, YY_True, ZZ_True], [XX, YY, ZZ]);
   else
      [~, D] = knnsearch([XX_True, YY_True], [XX, YY]);
   end
   nnfig = figure('Visible', 'off');
   hold on
   histogram(D, 100);
   title('true vs. drift corrected image');
   xlabel('nearest neighbor distances (nm)');
   ylabel('frequency');
   hold off

   % --------------------------------------------------------------------------

   % Compare found versus true drift.
   % Convert drift from pixels per frame to nm per frame.
   x_drift_true_nm = DriftX_True(:) * P2nm;
   y_drift_true_nm = DriftY_True(:) * P2nm;
   x_drift_found_nm = SMD.DriftX(:) * P2nm;
   y_drift_found_nm = SMD.DriftY(:) * P2nm;
   if ~isempty(DriftZ_True)
      % Convert drift from um per frame to nm per frame for z.
      z_drift_true_nm = DriftZ_True(:) * 1000; 
      z_drift_found_nm = SMD.DriftZ(:) * 1000;
   end

   % Average distance between the true and the found drift corrected curves.
   if ~isempty(DriftZ_True)
      dists2 = ...
         sqrt((x_drift_true_nm - x_drift_found_nm).^2 + ...
              (y_drift_true_nm - y_drift_found_nm).^2 + ...
              (z_drift_true_nm - z_drift_found_nm).^2);
   else
      dists2 = ...
         sqrt((x_drift_true_nm - x_drift_found_nm).^2 + ...
              (y_drift_true_nm - y_drift_found_nm).^2);
   end
   dist2 = sum(dists2) / numel(x_drift_found_nm);
   % RMSE for true versus found drift curves.
   rmse2 = sqrt(sum(dists2.^2) / numel(x_drift_found_nm));

end
