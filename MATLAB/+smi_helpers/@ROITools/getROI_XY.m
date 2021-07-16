function [n_ROIs, RoI] = getROI_XY(obj, XY, XY_SE, x_size, y_size, ...
                                   txt, XYsize, SMR)
% Let the user select ROIs from (x, y)-coordinates.
%
% INPUTS:
%    XY               cell array of (x, y) coordinates of points labeled in
%                     image (nm); XY{i} is the coordinates for label i
%    XY_SE            cell array of (x, y) standard deviations of the
%                     coordinate values (nm)
%    x_size, y_size   box diameters used when clicking the left mouse button
%                     (nm)
%    txt              text to label the ROI figure
%    XYsize           [OPTIONAL] (x, y) image size (nm) [1 x 2]; needed for
%                     displaying coordinates like DIPimage does (XYvsDIP)
%    SMR              [OPTIONAL] SMR structure from an SMA_SR data structure
%                     needed for gaussianImage in get_ROI (GaussIm)
%
% OUTPUTS:
%    n_ROIs           number of ROIs created
%    RoI              cell array for each ROI of
%                        ROI          [xmin, xmax, ymin, ymax] of ROI
%                        X, Y         (x, y) coordinates of points inside
%                        X_SE, Y_SE   (x, y) localization errors for the above
%                     NOTE: X{i}{j} is the x-coordinates of label j in ROI i,
%                     etc.

% Created by
%    Michael J. Wester (2021)

   if ~exist('SMR', 'var')
      SMR = [];
   end

   if obj.XYvsYX
      ix = 1;
      iy = 2;
   else
      ix = 2;
      iy = 1;
   end

   %if obj.XYvsDIP
      ix = 1;
      iy = 2;
   %else   % Below needed for SRtest, but not for SMR.
   %   ix = 2;
   %   iy = 1;
   %end

   n_labels = numel(XY);

   for j = 1 : n_labels
      n = size(XY{j}, 1);
      X_SE{j} = NaN(n, 1);
      Y_SE{j} = NaN(n, 1);

      if obj.XYvsDIP
         X{j} = XY{j}(:, ix);
         Y{j} = XY{j}(:, iy);
      else
         X{j} = XY{j}(:, ix);
         Y{j} = XYsize(2) - XY{j}(:, iy);
      end
      if ~isempty(XY_SE{j})
         X_SE{j} = XY_SE{j}(:, ix);
         Y_SE{j} = XY_SE{j}(:, iy);
      end
   end

   if obj.GaussIm && exist('SMR', 'var') && ~isempty(SMR)
      [n_ROIs, ROI, index_ROI] = ...
         obj.get_ROI_GaussIm(X, Y, x_size, y_size, txt, SMR);
   else
      [n_ROIs, ROI, index_ROI] = obj.get_ROI(X, Y, x_size, y_size, txt);
   end
   RoI = cell(1, n_ROIs);
   for i = 1 : n_ROIs
      fprintf('ROI %d = %.3f %.3f %.3f %.3f\n', i, ROI{i});
      RoI{i}.ROI   = ROI{i};
      for j = 1 : n_labels
         k = index_ROI{i}{j};
         RoI{i}.X{j}    = X{j}(k);
         RoI{i}.Y{j}    = Y{j}(k);
         RoI{i}.X_SE{j} = X_SE{j}(k);
         RoI{i}.Y_SE{j} = Y_SE{j}(k);
      end
   end

end
