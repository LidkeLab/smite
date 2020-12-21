classdef ChannelRegistration < handle
    %ChannelRegistration contains methods for channel registration.
    % This class contains methods for performing channel registration and
    % methods used to interpret/visualize the results.
    %
    % REQUIRES:
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
        Coordinates cell
        
        % Fiducial images (numeric array, MxPxNFiducials)
        FiducialImages {mustBeFloat}
        
        % Type of data used to compute transform (char)(Default = 'coords')
        % OPTIONS: 
        %   'coords': localizations (defined by (x, y) coordinates) are
        %             used to find the transform.
        %   'images': images are used directly to find the transform.
        TransformationBasis char = 'coords';
        
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
        NNeighborPoints(1, 1) {mustBeInteger} = 12;
        
        % Degree of polynomial for 'polynomial' tform (Default = 2)
        % This is only used when TransformationType = 'polynomial'.
        PolynomialDegree(1, 1) ...
            {mustBeMember(PolynomialDegree, [2, 3, 4])} = 2;
        
        % Auto-scale fiducial images (boolean)(Default = 1)
        % This flag lets this class do a somewhat arbitrary scaling of the
        % fiducial images in an attempt to simplify the code usage.  This
        % allows us to avoid gain/offset correcting the data, which might
        % be annoying in some cases (as in, it's nice to just use the
        % default SMF instead of having to tweak parameters just for this
        % code).
        AutoscaleFiducials(1, 1) logical = 1;
        
        % Manually cull localization pairs (boolean)(Default = 1)
        % This flag lets the user manually cull the paired localizations
        % used to produce the transform (this is only applicable for
        % TransformationBasis = 'coords').
        ManualCull(1, 1) logical = 1;
    end
    
    properties (SetAccess = protected)
        % Computed transformation(s) (cell array of tform objects)
        % Each element corresponds to a transform w.r.t. the fiducial in
        % obj.SMF.Data.FileName{1}, as in, RegistrationTransform{3} is a
        % registration between obj.SMF.Data.FileName{3} and 
        % obj.SMF.Data.FileName{1}.
        % RegistrationTransform{1} will either be empty or meaningless.
        RegistrationTransform cell
    end
    
    methods
        function [obj] = ChannelRegistration(...
                FiducialFileDir, FiducialFileNames, SMF)
            %ChannelRegistration class constructor.
            % The inputs can be used to set class properties if desired.
            %
            % INPUTS:
            %   FiducialFileDir: Name of directory containing the fiducial
            %                    files. (char array/string)
            %   FiducialFileNames: Filename(s) of the files containing the
            %                      fiducial images. 
            %                      (cell array of char/string).
            
            % Set inputs to class properties if needed.
            if (exist('SMF', 'var') && ~isempty(SMF))
                obj.SMF = SMF;
            else
                % Set a (mostly) default SMF, with a few tweaks that tend
                % to help out for several types of fiducial images.
                obj.SMF = smi_core.SingleMoleculeFitting;
                obj.SMF.Fitting.FitType = 'XYNBS';
            end
            if (exist('FiducialFileDir', 'var') ...
                    && ~isempty(FiducialFileDir))
                obj.SMF.Data.FileDir = FiducialFileDir;
            end
            if (exist('FiducialFileNames', 'var') ...
                    && ~isempty(FiducialFileNames))
                obj.SMF.Data.FileName = FiducialFileNames;
            end

            
        end
        
        [RegistrationTransform] = findTransform(obj);
        exportTransform(obj)
        gui(obj)
    end
    
    methods (Static)
        [PlotAxes, LineHandles] = ...
            plotCoordsOnData(PlotAxes, RawData, Coordinates);
        [SMDMoving, SMDFixed] = transformSMD(...
            RegistrationTransform, SMDMoving, SMDFixed);
        [MovingCoordinates, FixedCoordinates] = transformCoords(...
            RegistrationTransform, MovingCoordinates, FixedCoordinates);
        [TransformedImages] = transformImages(...
            RegistrationTransform, Images)
        [PlotFigure] = visualizeCoordTransform(PlotFigure, ...
            RegistrationTransform, FrameSize, GridSpacing);
        [PlotAxes] = visualizeImageTransform(PlotAxes, ...
            RegistrationTransform, FrameSize, GridSpacing);
        [SquaredError] = estimateRegistrationError(...
            RegistrationTransform, Coords1, Coords2);
        [PlotAxes] = visualizeRegistrationError(PlotAxes, ...
            RegistrationTransform, MovingCoordinates, FixedCoordinates, ...
            FOV, GridSpacing)
        [PlotAxes] = visualizeRegistrationResults(PlotAxes, ...
            RegistrationTransform, MovingCoordinates, FixedCoordinates, ...
            MovingImage, FixedImage);
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
    end
    
    
end