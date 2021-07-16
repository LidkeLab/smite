function [XY, XY_SE, XYsize, SMDimport] = import_XY(src, pixel2nm, fmt)
% Import N x 2 (x, y) coordinates and standard deviations from ...
%    a BGL, SMA_SR or SMD data structure,
%    point (x, y) coordinates [provided as N x 2 matrices, or if x_SE and
%       y_SE are appended, then N x 4 matrices],
%    super-resolution image files containing a BGL, SMA_SR or SMD data
%       structure,
%    or directly from a data structure with fields X and Y (and X_SE, Y_SE or
%       X_STD, Y_STD).
%
% INPUTS:
%    src              (x, y) coordinate source
%    pixel2nm         conversion factor from pixels to nm.
%                        NOTE: overridden by existing value in SMD
%    fmt              special data format (e.g., 'EM');
%
% OUTPUTS (units in nm):
%    XY               array of (x, y) point coordinates [N x 2]
%    XY_SE            array of (x, y) standard deviations of the coordinate
%                     values [N x 2]
%    XYsize           (x, y) sizes of the src image [1 x 2]
%    SMDimport        SMD structure either returned or contructed in order for
%                     get_ROI_GaussIm to properly produce a Gaussian image

% Created by:
%    Michael J. Wester (2021)

   XY_SE = [];
   SMDimport = [];

   if strcmp(class(src), 'SMA_SR')
      % src is an SMA_SR data structure
      XY = [ double(src.SMR.X) .* pixel2nm, ...
             double(src.SMR.Y) .* pixel2nm ];
      XY_SE = [ double(src.SMR.X_SE) .* pixel2nm, ...
                double(src.SMR.Y_SE) .* pixel2nm ];
      XYsize = [src.XSize, src.YSize] .* pixel2nm;
      SMDimport = src.SMR;
      SMDimport.PixelSize = pixel2nm / 1000;

   elseif strcmp(class(src), 'BaGoL') || strcmp(class(src), 'smi.BaGoL')
      % src is a BGL.SMD data structure (from BaGoL)
      XY = [ double(src.MAPN.X), double(src.MAPN.Y) ];
      XY_SE = [ double(src.MAPN.X_SE), double(src.MAPN.Y_SE) ];
      %XYsize = [floor(min(min(XY(:, 1:2)))), ceil(max(max(XY(:, 1:2))))];
      XYsize = [256, 256] .* pixel2nm;
      SMDimport = src.MAPN;
      SMDimport.X = SMDimport.X ./ pixel2nm;
      SMDimport.Y = SMDimport.Y ./ pixel2nm;
      SMDimport.X_SE = SMDimport.X_SE ./ pixel2nm;
      SMDimport.Y_SE = SMDimport.Y_SE ./ pixel2nm;
      SMDimport.XSize = XYsize(1) ./ pixel2nm;
      SMDimport.YSize = XYsize(2) ./ pixel2nm;
      SMDimport.PixelSize = pixel2nm / 1000;
      n_locs = numel(src.MAPN.X);
      SMDimport.Bg = zeros(n_locs, 1);
      SMDimport.Photons = 1000 * ones(n_locs, 1);
      SMDimport.FrameNum = ones(n_locs, 1);

   elseif ismatrix(src) && ~ischar(src) && ~isstruct(src)
      % src is a matrix
      n_cols = size(src, 2);
      if n_cols == 2 || n_cols == 4
         XY = src(:, 1:2) .* pixel2nm;
      else
         error('src matrix has %d rather than 2 or 4 columns!', n_cols);
      end
      if n_cols == 4
         XY_SE = src(:, 3:4) .* pixel2nm;
      end
      XYsize = [floor(min(src(:, 1:2))), ceil(max(src(:, 1:2)))] .* pixel2nm;

   elseif ischar(src)
      % src is a filename
      if isempty(fmt)
         load(src);
         if exist('SMD', 'var')
            [XY, XY_SE, XYsize, SMDimport] = ...
               smi_helpers.ROITools.import_XY(SMD, pixel2nm);
         elseif exist('SMASR', 'var')
            [XY, XY_SE, XYsize, SMDimport] = ...
               smi_helpers.ROITools.import_XY(SMASR, pixel2nm);
         elseif exist('SMR', 'var') & ~isempty(SMR)
            [XY, XY_SE, XYsize] = ...
               smi_helpers.ROITools.import_XY(SMR, pixel2nm);
         elseif exist('BGL', 'var')
            [XY, XY_SE, XYsize, SMDimport] = ...
               smi_helpers.ROITools.import_XY(BGL, pixel2nm);
         else
            error('No BGL, SMASR, SMD or SMR object found in %s!', src);
         end
      elseif strcmp(fmt, 'EM')
         % (x, y) data is in columns (2, 3).
         [x, y] = textread(src, '%*u %u %u %*u', 'headerlines', 1);
         XY = [x, y] .* pixel2nm;
         XYsize = [256, 256] .* pixel2nm;
      else
         error('Unknown fmt (%s)!', fmt);
      end

   elseif isstruct(src)
      % src should be a data structure with fields X, Y and optionally
      % X_SE, Y_SE or X_STD, Y_STD as well as XSize, YSize
      if isfield(src, 'X') & isfield(src, 'Y')
         % If PixelSize is available, use it!
         if ~isempty(pixel2nm) && isfield(src, 'PixelSize') ...
            && pixel2nm ~= 1000 * src.PixelSize
            warning(['Incompatible pixel2nm/PixelSize specification:\n', ...
                     '(PixelSize from input source = %f overrides pixel2nm = %f)\n'], ...
                    pixel2nm, 1000 * src.PixelSize);
            pixel2nm = 1000 * src.PixelSize;
         end
         XY = [ double(src.X) .* pixel2nm, double(src.Y) .* pixel2nm ];
         if isfield(src, 'X_SE') & isfield(src, 'Y_SE')
            XY_SE = [ double(src.X_SE) .* pixel2nm, ...
                      double(src.Y_SE) .* pixel2nm ];
         elseif isfield(src, 'X_STD') & isfield(src, 'Y_STD')
            XY_SE = [ double(src.X_STD) .* pixel2nm, ...
                      double(src.Y_STD) .* pixel2nm ];
         end
         if isfield(src, 'XSize') & isfield(src, 'YSize')
            XYsize = [src.XSize, src.YSize] .* pixel2nm;
         else
            %XYsize = [ceil(max(src.X)), ceil(max(src.Y))] .* pixel2nm;
            XYsize = [256, 256] .* pixel2nm;
         end
         SMDimport = src;
         if ~isfield(src, 'PixelSize')
            SMDimport.PixelSize = pixel2nm / 1000;
         end
      else
         error('Fields X, Y not found in src!');
      end

   else
      error('Incorrect argument type for src!');
   end

end
