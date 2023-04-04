### +smi_core/@ChannelRegistration

This class contains methods for performing channel registration and
methods used to interpret/visualize the results.

REQUIRES:
- MATLAB 2019b or later (some newer method inputs are used, e.g.,
      size(Image, [1, 2]) wasn't allowed pre-2019b).
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

---
    
```
properties:
   % Single molecule fitting structure (see SingleMoleculeFitting)
   % This SMF structure is used to find localizations in the fiducial
   % files specified by FiducialFilePath.  If
   % TransformationBasis = 'images' this is not needed.
   SMF
   
   % Coords used to compute transforms (cell array of numeric array)
   % These coordinates are ordered as follows: the coordinates used to
   % find a transform from fiducial m to fiducial n will be stored in
   % Coordinates{m}.  The coordinates from fiducial m will be in
   % Coordinates{m}(:, :, m), and those from fiducial n will be in
   % Coordinates{m}(:, :, n), with each of these being organized as
   % two-column arrays [X, Y].
   Coordinates cell
   
   % Indicates fiducial images set manually by user. (Default = false)
   % This can be useful if the user already has the set of fiducial
   % images loaded into their workspace (e.g., if they were simulated
   % and not saved to a file).
   % NOTE: If using this setting, the fiducials must be
   %       "pre-processed" before setting to obj.FiducialImages (e.g.,
   %       truncated to proper ROI, averaged over time, ...).
   %       Furthermore, obj.FiducialROI should be set to the
   %       appropriate value so that the correct ROI is saved along
   %       with the transform. obj.SplitFormat cannot be used with
   %       this setting.
   ManualSetFiducials(1, 1) logical = false;
   
   % Fiducial images (numeric array, MxP(xNFiducials))
   FiducialImages {mustBeNumeric(FiducialImages)}
   
   % Format guiding the fiducial ROIs to be used (Default = [1])
   % (see obj.convertSplitFormatToROIs() for a more complete
   % description).
   % NOTE: If you wish to manually set obj.FiducialROI, you must set
   %       obj.SplitFormat = [].
   SplitFormat {mustBeInteger(SplitFormat)} = 1;
   
   % Fiducial ROIs ([YStart, XStart, YEnd, XEnd, ZStart, ZPeriod])
   % NOTE: In its current usuage, FiducialROI is set automatically in
   %       findTransform(), or manually by the user.
   % NOTE: If you wish to manually define this array, you must set
   %       obj.SplitFormat = [].  Otherwise, the ROI splitting scheme
   %       defined by obj.SplitFormat takes precedence.
   % OPTIONS:
   %   If size(FiducialROI, 1) == 1, each image in FiducialImages will
   %       be truncated to the ROI specified by FiducialROI before
   %       computing the transform.
   %   If size(FiducialROI, 1) > 1, the image found in in
   %       SMF.Data.FileName{1} will be split up into
   %       the ROIs defined by each row of FiducialROI.  Each row of
   %       FiducialROI must specify an equal size ROI.
   %       FiducialROI(1, :) will be treated as the "reference" (or
   %       "fixed") fiducial, meaning all other ROIs will be
   %       transformed w.r.t. FiducialROI(1, :).  Furthermore, the
   %       properties 'Coordinates' and 'RegistrationTransform' will
   %       follow the same ordering as FiducialROI.
   %   If FiducialROI is not set by the user, it will be given a
   %       default depending on how many files are specified by
   %       obj.SMF.Data.FileName.  If there is only one file,
   %       FiducialROI will be set by default s.t. the image in the
   %       one file will be split in two along its columns.  If there
   %       are multiple files,
   %       FiducialROI = [1, 1, size(FiducialImages(:, :, 1))]
   %       where FiducialImages will contain the image stored in the
   %       file obj.SMF.Data.FileName.
   FiducialROI(:, 6) {mustBeInteger(FiducialROI)}
   
   % Data used to compute transform (char)(Default = 'coordinates')
   % OPTIONS:
   %   'coordinates': localizations (defined by (x, y) coordinates)
   %                  are used to find the transform.
   %   'images': images are used directly to find the transform.
   TransformationBasis char {mustBeMember(TransformationBasis, ...
       {'coordinates', 'images'})} = 'coordinates';
   
   % Type of transform to be computed (char array)(Default = 'lwm')
   % OPTIONS:
   %   If TransformationBasis = 'coords', this can be set to any of
   %       the transformationType options defined in doc fitgeotrans
   %   If TransformationBasis = 'images', this can be set to any of
   %       the transformType options defined in doc imregtform.
   TransformationType char = 'lwm';
   
   % Threshold for pairing localizations (Pixels)(Default = inf)
   % This only matters when TransformationBasis = 'coords'.
   SeparationThreshold(1, 1) = inf;
   
   % # of neighbor points used to compute transform (Default = 12)
   % This is only used when TransformationType = 'lwm'
   NNeighborPoints(1, 1) {mustBeInteger, ...
       mustBeGreaterThan(NNeighborPoints, 5)} = 12;
   
   % Degree of polynomial for 'polynomial' tform (Default = 2)
   % This is only used when TransformationType = 'polynomial'.
   PolynomialDegree(1, 1) ...
       {mustBeMember(PolynomialDegree, [2, 3, 4])} = 2;
   
   % Auto-scale fiducial images (boolean)(Default = true)
   % This flag lets this class do a somewhat arbitrary scaling of the
   % fiducial images in an attempt to simplify the code usage.  This
   % allows us to avoid gain/offset correcting the data, which might
   % be annoying in some cases (as in, it's nice to just use the
   % default SMF instead of having to tweak parameters just for this
   % code).
   AutoscaleFiducials(1, 1) logical = true;
   
   % Manually cull localization pairs (boolean)(Default = true)
   % This flag lets the user manually cull the paired localizations
   % used to produce the transform (this is only applicable for
   % TransformationBasis = 'coords').
   ManualCull(1, 1) logical = true;
   
   % Verbosity level for standard workflow. (Default = 1)
   %   0: Command Window updates will be supressed where possible and
   %      reasonable.
   %   1: Some updates may appear in Command Window
   %   2: More detailed updates in Command Window
   %   3: Lot's of info. may be passed to Command Window. This mode
   %      may be useful for debugging large workflows encompassing
   %      this class.
   Verbose = 1;
```

---

methods:
- **[convertSplitFormatToROIs](convertSplitFormatToROIs.m)**:
  converts a split format to an array of ROIs
- **[estimateRegErrorLOO](estimateRegErrorLOO.m)**:
  estimates registration error by leave-one-out analysis
- **[estimateRegistrationError](estimateRegistrationError.m)**:
  estimates the registration error
- **[exportTransform](exportTransform.m)**:
  exports transform information into a .mat file
- **[findTransform](findTransform.m)**:
  finds a channel registration transform
- **[gui](gui.m)**:
  is the GUI method for the ChannelRegistration class
- **[loadFiducials](loadFiducials.m)**:
  loads fiducial files and sets associated class properties
- **[pairCoordinates](pairCoordinates.m)**:
  pairs sets of coordinates with each other
- **[performManualCull](performManualCull.m)**:
  performs the interactive (graphical) culling process
- **[plotCoordsOnData](plotCoordsOnData.m)**:
  plots coordinates in Coordinates on top of ScaledData
- **[rescaleFiducials](rescaleFiducials.m)**:
  rescales the images in Fiducials as needed
- **[simFiducials](simFiducials.m)**:
  simulates two channel fiducials to test channel registration
- **[transformCoords](transformCoords.m)**:
  transforms a set of coordinates with the given transform
- **[transformCoordsDirect](transformCoordsDirect.m)**:
  transforms a set of coordinates directly
- **[transformImages](transformImages.m)**:
  transforms a set of images with the given transform
- **[transformSMD](transformSMD.m)**:
  transforms SMD structures using the specified transform
- **[transformTR](transformTR.m)**:
  transforms a TR structure using the specified transform
- **[unitTest](unitTest.m)**:
  tests vital functionality of smi_core.ChannelRegistration
- **[visualizeCoordTransform](visualizeCoordTransform.m)**:
  creates visuals for a coordinate transform
- **[visualizeImageTransform](visualizeImageTransform.m)**:
  creates visuals for an image transform
- **[visualizeRegistrationError](visualizeRegistrationError.m)**:
  visualizes the error in RegistrationTransform
- **[visualizeRegistrationResults](visualizeRegistrationResults.m)**:
  shows registration results on the fiducials
