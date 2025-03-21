function [n_ROIs, RoI, XYsize] = getROI(obj, src, txt)
% Let the user select ROIs from
%    a BGL, SMA_SR or SMD data structure,
%    point (x, y) coordinates [provided as N x 2 matrices, or if X_SE and
%       Y_SE are appended, then N x 4 matrices],
%    super-resolution image files containing a BGL, SMA_SR or SMD data
%       structure,
%    or directly from a data structure with fields X and Y (and X_SE, Y_SE or
%       X_STD, Y_STD).
%
% INPUTS:
%    obj   various properties used by the algorithms
%       Pixel2nm    conversion factor from pixels to nm 
%       ROI_sizes   box diameters used when clicking the left mouse button (nm)
%       Color       label colors
%       Order       plotting order
%       Msize       marker size
%       XYvsYX      coordinate order
%       OriginLLvsUL coordinate origin (llower left vs. upper left)
%       GaussIm     use Gaussian image for display
%       EM          if true, data is in EM format
%    src   cell array of (x, y) coordinate sources, although it is permissible
%          to omit the cell array for a single source
%    txt   [OPTIONAL] text to label the ROI figure
%
% OUTPUTS (units in nm):
%    n_ROIs      number of ROIs created
%    RoI         cell array for each ROI of
%                   ROI            [xmin, xmax, ymin, ymax] of ROI
%                   X, Y           (x, y) coordinates of points inside
%                   X_SE, Y_SE     (x, y) localization errors for the above
%                   SMD            SMD structure defining the ROI
%    XYsize      (x, y) image size (nm) [1 x 2]; needed for
%                displaying coordinates where the origin in the UL corner
%                (OriginLLvsUL)

% Created by:
%    Michael J. Wester (2021)

   if ~exist('txt', 'var')
      txt = '';
   end

   x_size = obj.ROI_sizes(1);
   y_size = obj.ROI_sizes(2);

   n_ROIs = 0;
   RoI = [];

   fmt = '';
   if obj.EM
      fmt = 'EM';
   end

   if ~iscell(src)
      n_labels = 1;
      [XY{1}, XY_SE{1}, XYsize, SMD{1}] = ...
         obj.import_XY(src, obj.Pixel2nm, fmt);
   else
      n_labels = numel(src);
      for i = 1 : n_labels
         [XY{i}, XY_SE{i}, XYsize, SMD{i}] = ...
            obj.import_XY(src{i}, obj.Pixel2nm, fmt);
      end
   end
   [n_ROIs, RoI] = ...
      obj.getROI_XY(XY, XY_SE, x_size, y_size, txt, XYsize, SMD);

   % Gather the appropriate SMD ROI to stash into RoI.
   % NOTE:
   %    RoI.ROI [xmin, xmax, ymin, ymax] of ROI in nm
   %    isolateSubROI wants [YStart, XStart, YEnd, XEnd] in pixels
   if ~isempty(SMD)
      for i = 1 : n_labels
         for j = 1 : n_ROIs
            ROIpx = RoI{j}.ROI ./ obj.Pixel2nm;
            ROIpx = ROIpx([3, 1, 4, 2]);
            RoI{j}.SMD{i} = ...
               smi_core.SingleMoleculeData.isolateSubROI(SMD{i}, ROIpx);
         end
      end
   end

end
