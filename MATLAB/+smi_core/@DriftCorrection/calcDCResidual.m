function [residual, dist, rmse, nnfig] = ...
   calcDCResidual(SMD, X_True, Y_True, Z_True)
%calcDCResidual calculates the residual of 2D SMD relative to true coordinates
% also returning averaged statistics comparing true vs. drift corrected
% results for 2D or 3D.
%
% INPUTS:
%   SMD              A structure with fields:
%      X                  x-coordinates (Nx1) where N is total number of points
%      Y                  y-coordinates (Nx1) (pixels, also for X)
%      Z                  z-coordinates (Nx1) (um) [OPTIONAL]
%      PixelSizeZUnit     um per pixel
%   X_True           true x-coordinates (pixels)
%   Y_True           true y-coordinates (pixels)
%   Z_True           true z-coordinates (um) [OPTIONAL]
%
% OUTPUTS:
%   residual     residual image
%   dist         average distance between true and drift corrected coordinates
%                (nm)
%   rmse         root-mean-square-error between true and drift corrected coords
%                (nm)
%   nnfig        histogram of nearest neighbor distances between drift
%                corrected and true coordinates

% Created by
%   Michael J. Wester (Lidke Lab 2018)

   SRImageZoom = 10;

   % Remove any NaNs.
   nans = find(isnan(SMD.X) | isnan(SMD.Y) | isnan(X_True) | isnan(Y_True));
   n_nans = numel(nans);
   if n_nans > 0
      fprintf('calcDCResidual: %d NaNs removed!\n', n_nans);
      SMD.X(nans) = [];
      SMD.Y(nans) = [];
      X_True(nans) = [];
      Y_True(nans) = [];
   end
   N = numel(SMD.X);

   % Drift corrected image.
   DCimage = SMA_Vis.histogramImage(SMD, SRImageZoom);

   % True image.
   SMD_True = SMD;
   SMD_True.X = X_True;
   SMD_True.Y = Y_True;
   TrueImage = SMA_Vis.histogramImage(SMD_True, SRImageZoom);

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
%  unDCimage = SMA_Vis.histogramImage(SMD_unDC, SRImageZoom);

   P2nm = SMD.PixelSizeZUnit * 1000;

   % Drift corrected and true coordinates.
   X = SMD.X * P2nm;
   Y = SMD.Y * P2nm;
   X_True = X_True * P2nm;
   Y_True = Y_True * P2nm;
   if exist('Z_True', 'var')
      Z      = SMD.Z  * 1000;
      Z_True = Z_True * 1000;
   end

   % Average distance between the true and drift corrected coordinates.
   if exist('Z_True', 'var')
      dists = sqrt((X_True - X).^2 + (Y_True - Y).^2 + (Z_True - Z).^2);
   else
      dists = sqrt((X_True - X).^2 + (Y_True - Y).^2);
   end
   dist = sum(dists) / N;
   % RMSE between the true and drift corrected coordinates.
   rmse = sqrt(sum(dists.^2) / N);
   clear dists

   % Histogram of nearest neighbor distances between the true and the drift
   % corrected coordinates.
   if exist('Z_True', 'var')
      [~, D] = knnsearch([X_True, Y_True, Z_True], [X, Y, Z]);
   else
      [~, D] = knnsearch([X_True, Y_True], [X, Y]);
   end
   nnfig = figure('Visible', 'off');
   hold on
   histogram(D, 100);
   title('true vs. drift corrected image');
   xlabel('nearest neighbor distances (nm)');
   ylabel('frequency');
   hold off

   residual = abs(DCimage - TrueImage);

end
