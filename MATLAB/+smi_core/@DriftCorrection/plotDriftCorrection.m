function DC_fig = plotDriftCorrection(obj, SMD, option)
%plotDriftCorrection plots the computed drift correction stored in SMD
% stucture for 2D or 3D data.  The plot is color coded so as to indicate the
% drift correction as a function of time.
%
% INPUTS:
%    SMD:         A structure with fields:
%       X             x-coordinates (Nx1) where N is total number of points
%       Y             y-coordinates (Nx1)
%       Z             z-coordinates (Nx1) [OPTIONAL]
%                  Note: X and Y are typically in pixels, Z in um
%       NFrames       number of frames in each dataset
%       NDatasets     number of datasets
%       DriftX        found x-drift (NFrames x NDatasets)
%       DriftY        found y-drift (NFrames x NDatasets)
%       DriftZ        found z drift (NFrames x NDatasets) [OPTIONAL]
%    option       [OPTIONAL] type of drift correction plot to make:
%                     'A' - absolute, 'R' - relative, '1' - initial absolute
%                     (Default = 'A')
%    obj:         [class properties]
%       PixelSizeZUnit pixel size in um (Default = 0.1)
%       PDegree       degree of the intra-dataset drift rate fitting polynomial
%       Init_inter    inter-dataset initialization with respect to the previous
%                     dataset (see driftCorrectKNN)
%
% OUTPUTS:
%    DC_fig       figure handle

% Created by
%    Michael J. Wester (Lidke Lab 2018)

   if ~isfield(SMD, 'DriftX') | isempty(SMD.DriftX)
      error('SMD.DriftX either missing or empty, so DriftCorrection plot unavailable!');
   end

   if ~exist('option', 'var')
      option = 'A';
   end

   PixelSizeZUnit = obj.PixelSizeZUnit;
%  DriftParams = obj.DriftParams;
%  if exist('DriftParams', 'var') & isfield('DriftParams', 'PixelSizeZUnit')
%     PixelSizeZUnit = DriftParams.PixelSizeZUnit;
%  else
%     PixelSizeZUnit = 0.1;
%  end
   P2nm = PixelSizeZUnit * 1000;

   SMZ.DriftX = SMD.DriftX * P2nm;
   SMZ.DriftY = SMD.DriftY * P2nm;
   if isfield(SMD, 'Z') && numel(SMD.Z) == numel(SMD.X)
      Ndims = 3;
      SMZ.DriftZ = SMD.DriftZ * 1000;
   else
      Ndims = 2;
   end

   N = SMD.NDatasets * SMD.NFrames;

   % Drift correction plot.
   DC_fig = figure('Visible', 'off');
   hold on
   % Color code the plot from blue to red to indicate the passage of time.
   cm = colormap(jet);
   n_cm = size(cm, 1);

   % Within each dataset, taper the segment width from large to small to
   % indicate the direction of time.
   M = 500;
   sz = arrayfun(@(i) max(1, round(M * i / SMD.NFrames)), ...
                 SMD.NFrames : -1 : 1);
   for i = 1 : SMD.NDatasets
      k = (i - 1)*SMD.NFrames;
      indx = max(1, ceil(n_cm * (k + [1 : SMD.NFrames]) / N));
      if Ndims == 2
         if option == 'R' & i > 1
            SMZ.DriftX(:, i) = SMZ.DriftX(:, i) - ...
               (SMZ.DriftX(1, i) - SMZ.DriftX(SMD.NFrames, i - 1));
            SMZ.DriftY(:, i) = SMZ.DriftY(:, i) - ...
               (SMZ.DriftY(1, i) - SMZ.DriftY(SMD.NFrames, i - 1));
         end
         if option == '1'
            scatter(SMZ.DriftX(1, i), SMZ.DriftY(1, i), sz(1), ...
                    cm(indx(1), :), '.');
         else
            scatter(SMZ.DriftX(:, i), SMZ.DriftY(:, i), sz, cm(indx, :), '.');
         end
      else
         if option == 'R' & i > 1
            SMZ.DriftX(:, i) = SMZ.DriftX(:, i) - ...
               (SMZ.DriftX(1, i) - SMZ.DriftX(SMD.NFrames, i - 1));
            SMZ.DriftY(:, i) = SMZ.DriftY(:, i) - ...
               (SMZ.DriftY(1, i) - SMZ.DriftY(SMD.NFrames, i - 1));
            SMZ.DriftZ(:, i) = SMZ.DriftZ(:, i) - ...
               (SMZ.DriftZ(1, i) - SMZ.DriftZ(SMD.NFrames, i - 1));
         end
         if option == '1'
            scatter3(SMZ.DriftX(1, i), SMZ.DriftY(1, i), SMZ.DriftZ(1, i), ...
                     sz(1), cm(indx(1), :), '.');
         else
            scatter3(SMZ.DriftX(:, i), SMZ.DriftY(:, i), SMZ.DriftZ(:, i), ...
                     sz, cm(indx, :), '.');
         end
      end
   end
   cb = colorbar;
   cb.Label.String = sprintf('fraction of total frames (= %d x %d)', ...
                             SMD.NDatasets, SMD.NFrames);

   % Root mean square error computed from the drift in the first frame of each
   % dataset which should be zero (in theory) for registered data.
   RMSE = sum(SMZ.DriftX(1, :).^2)/SMD.NDatasets + ...
          sum(SMZ.DriftY(1, :).^2)/SMD.NDatasets;
   if Ndims == 3
      RMSE = RMSE + sum(SMZ.DriftZ(1, :).^2)/SMD.NDatasets;
   end
   RMSE = sqrt(RMSE);

   if option == 'R'
      prefix = 'relative ';
   else
      prefix = '';
   end
%  if exist('DriftParams', 'var')
%     title(sprintf([prefix, 'drift correction: RMSE_{reg} = %f nm\n', ...
%                    '(PDegree = %d, Init\\_inter = %d)'],          ...
%                   RMSE, DriftParams.PDegree, DriftParams.Init_inter));
%  else
      title(sprintf([prefix, 'drift correction: RMSE_{reg} = %f nm\n'], RMSE));
%  end
   xlabel('x (nm)');
   ylabel('y (nm)');
   if Ndims == 3
      zlabel('z (nm)');
   end
   
   hold off

end
