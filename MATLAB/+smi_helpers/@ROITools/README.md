### +smi_helpers/@ROITools

These are tools to select ROIs from an image and save them into a data
structure.  The main interface function is getROI.  All other methods are
called by this routine.  getROI can take many types of sources (src):

- a BGL, SMA_SR or SMD data structure,
- point (x, y) coordinates [provided as N x 2 matrices, or if X_SE and
     Y_SE are appended, then N x 4 matrices],
- super-resolution image files containing a BGL, SMA_SR or SMD data
     structure,
- or directly from a data structure with fields X and Y (and X_SE, Y_SE or
     X_STD, Y_STD).

---

A simple invocation of getROI for a single label image is as follows:
```
   RT = ROITools;
   [n_ROIs, RoI, XYsize] = RT.getROI(Label1, 'Image Name');
```
This produces a cell array of dimension n_ROIs called RoI, defined as follows
(all units are in nm):
```
   RoI      cell array for each ROI of
               ROI          [xmin, xmax, ymin, ymax] of ROI
               X, Y         (x, y) coordinates of points inside
               X_SE, Y_SE   (x, y) localization errors for the above
```
Note that X, Y, X_SE, Y_SE are actually cell arrays.  For example, RoI{1}
will contain:
```
   X: {[n1 x 1 double]}, etc.
```
If src is a cell array, then each element will be taken to be a separate
label.  For example,
```
   RT = ROITools;
   [n_ROIs, RoI, XYsize] = RT.getROI({Label1, Label2}, 'Image Name');
```
will combine the images for Label1 and Label2, and allow selection of ROIs
that are compatible with both labels.  In this case, RoI{1} will contain:
```
   X: {[n1 x 1 double], [n2 x 1 double]}, etc.,
```
where n1, n2 are the number of localizations from Label1, Label2,
respectively.

By default, point images will be displayed from which to choose ROIs.  To
use Gaussian images (currently only available for one label input), set
```
   RT.GaussIm = true;   RT.OriginLLvsUL = false;
```
before invoking RT.getROI.  As a final comment, XYsize is the (x, y) image
size (nm) [1 x 2] needed for displaying coordinates where the origin is in
the UL corner.  This only is needed for BaGoL files when GaussIm is true---
see code fragment:
```
  Xnm = BGL.MAPN.X;
  if RT.GaussIm   % if GaussIm, the origin is in the upper left corner
     Ynm = XYsize(2) - BGL.MAPN.Y;
  else
     Ynm = BGL.MAPN.Y;
  end
```

---

```
properties:
   % Box sizes used when clicking the left mouse button (nm).
   ROI_sizes = [3000, 3000];   % nm
   Color = ['g', 'm', 'k', 'c', 'r', 'y'];   % Label colors.
   Order = 1 : 6;     % Plotting order.
   Msize = 7;         % Marker size.
   XYvsYX  = true;    % Coordinate order.
   % GaussIm = true indicates that gaussianImage will be used for the ROI
   %    selection display.
   % OriginLLvsUL = true says to use lower left origin coordinates rather than
   %    upper left origin coordinates in the ROI selection display.
   % SRzoom is the zoom factor for displaying Gaussian images.
   %
   % If GaussIm is true, make OriginLLvsUL false for consistency.  The default
   % is: GaussIm = false and OriginLLvsUL = true.  Currently, GaussIm only
   % works for single labeled molecules.
   GaussIm = false;       % Use gaussianImage for ROI selection display.
   OriginLLvsUL = true;   % Use LL origin coordinates rather than UL.
   SRzoom = 1;            % Zoom factor for gaussianImage.
   EM = false;            % EM data format for import_XY.
   Transform = {};        % Coordinate transform from label i -> 1.
   Mask = {};             % Boolean image mask to be used with an SMD source
```
