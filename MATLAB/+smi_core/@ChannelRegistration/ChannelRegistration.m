classdef ChannelRegistration < handle
    %ChannelRegistration contains methods for channel registration.
    % This class contains methods for performing channel registration and
    % methods used to interpret/visualize the results.
    %
    % REQUIRES:
    %   MATLAB 2019b or later (some newer method inputs are used, e.g.,
    %       size(Image, [1, 2]) wasn't allowed pre-2019b).
    %   MATLAB Image Processing Toolbox 2013b or later
    
    % Created by:
    %   David J. Schodt (Lidke Lab, 2020)
    
    
    properties
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
    end
    
    properties (SetAccess = protected)
        % Computed transformation(s) (cell array of tform objects)
        % RegistrationTransform will be organized differently depending on
        % the usage of FiducialROI:
        % If size(FiducialROI, 1) == 1, then each fiducial was provided as
        %   a separate image, in which case each element corresponds to a
        %   transform w.r.t. the fiducial in obj.SMF.Data.FileName{1}. For
        %   example, RegistrationTransform{3} is a registration between
        %   obj.SMF.Data.FileName{3} and obj.SMF.Data.FileName{1}.
        % If size(FiducialROI, 1) > 1, only one fiducial image was
        %   provided, but it will be split into different ROIs.  In this
        %   case, RegistrationTransform{n} will be the transform between
        %   the ROI defined by FiducialROI(n, :) and FiducialROI(1, :).
        RegistrationTransform cell
    end
    
    properties (Hidden, SetAccess = protected)
        % These properties are used for convenience internally (e.g., for
        % use in the GUI) and shouldn't be modified.
        
        TransformationBasisOptions cell = {'coordinates', 'images'};
        CoordTransformOptions cell = {'nonreflectivesimilarity', ...
            'similarity', 'affine', 'projective', 'polynomial', 'pwl', ...
            'lwm'};
        ImageTransformOptions cell = {'translation', 'rigid', ...
            'similarity', 'affine'};
        SplitFormatOptions cell = {[1, 2], [1; 2], [1, 3; 2, 4]};
        SplitFormatOptionsChar cell = ...
            {'[1, 2]', '[1; 2]', '[1, 3; 2, 4]'};
        
    end
    
    methods
        function [obj] = ChannelRegistration(...
                FiducialFileDir, FiducialFileNames, SMF, Verbose)
            %ChannelRegistration class constructor.
            % The inputs can be used to set class properties if desired.
            
            % Set inputs to class properties if needed.
            if (exist('Verbose', 'var') && ~isempty(Verbose))
                obj.Verbose = Verbose;
            end
            if (exist('SMF', 'var') && ~isempty(SMF))
                obj.SMF = SMF;
                if (obj.Verbose > 2)
                    fprintf(['\tChannelRegistration constructor: ', ...
                        'Input SMF structure stored as a class ', ...
                        'property\n'])
                end
            else
                % Set a (mostly) default SMF, with a few tweaks that tend
                % to help out for several types of fiducial images.
                if (obj.Verbose > 2)
                    fprintf(['\tChannelRegistration constructor: ', ...
                        'Using default SMF structure.\n'])
                end
                obj.SMF = smi_core.SingleMoleculeFitting;
                obj.SMF.BoxFinding.MinPhotons = 100;
                obj.SMF.Fitting.FitType = 'XYNBS';
            end
            if (exist('FiducialFileDir', 'var') ...
                    && ~isempty(FiducialFileDir))
                obj.SMF.Data.FileDir = FiducialFileDir;
                if (obj.Verbose > 2)
                    fprintf(['\tChannelRegistration constructor: ', ...
                        'Input \n\tFiducialFileDir = ''%s'' stored ', ...
                        'as a class property.\n'], FiducialFileDir)
                end
            end
            if (exist('FiducialFileNames', 'var') ...
                    && ~isempty(FiducialFileNames))
                obj.SMF.Data.FileName = FiducialFileNames;
                if (obj.Verbose > 2)
                    fprintf(['\tChannelRegistration constructor: ', ...
                        'Input FiducialFileNames stored as a ', ...
                        'class property.\n'])
                end
            end
            if (obj.Verbose > 1)
                fprintf(['\tChannelRegistration: constructor ran ', ...
                    'succesfully.\n'])
            end
            
        end
        
        function set.SMF(obj, SMFInput)
            % This is a set method for the SMF to ensure a valid SMF is
            % provided.
            obj.SMF = smi_core.SingleMoleculeFitting.reloadSMF(SMFInput);
        end
        
        [RegistrationTransform] = findTransform(obj);
        loadFiducials(obj)
        [FilePath] = exportTransform(obj, FileDir);
        gui(obj, GUIParent)
        
    end
    
    methods (Static)
        [PlotAxes, LineHandles] = ...
            plotCoordsOnData(PlotAxes, RawData, Coordinates);
        [SMDMoving] = transformSMD(RegistrationTransform, SMDMoving);
        [TRMoving] = transformTR(RegistrationTransform, TRMoving);
        [MovingCoordsTransformed] = ...
            transformCoords(RegistrationTransform, MovingCoords);
        [TransformedImages] = transformImages(...
            RegistrationTransform, Images)
        [PlotFigure] = visualizeCoordTransform(PlotFigure, ...
            RegistrationTransform, FrameSize, GridSpacing);
        [PlotAxes] = visualizeImageTransform(PlotAxes, ...
            RegistrationTransform, FrameSize, GridSpacing);
        [SquaredError] = estimateRegistrationError(...
            RegistrationTransform, MovingCoordinates, FixedCoordinates);
        [SquaredError] = estimateRegErrorLOO(...
            TransformationType, TransformationParams, ...
            MovingCoordinates, FixedCoordinates);
        [PlotAxes] = visualizeRegistrationError(PlotAxes, ...
            RegistrationTransform, MovingCoordinates, FixedCoordinates, ...
            FOV, GridSpacing)
        [PlotAxes] = visualizeRegistrationResults(PlotAxes, ...
            RegistrationTransform, MovingCoordinates, FixedCoordinates, ...
            MovingImage, FixedImage);
        unitTest()
    end
    
    methods (Static, Hidden)
        % These methods aren't expected to be used directly by the user, so
        % it's nice to hide them from view (so they don't distract the
        % user).
        [PairMap12, PairMap21] = pairCoordinates(Coords1, Coords2, ...
            SeparationThreshold);
        [CulledCoordinates] = performManualCull(RawData, Coordinates);
        [TransformedCoordinates] = transformCoordsDirect(...
            RegistrationTransform, Coordinates);
        [SplitROIs] = convertSplitFormatToROIs(FullROI, SplitFormat)
        [ScaledFiducials] = rescaleFiducials(Fiducials, SMF, AutoScale);
    end
    
    
end