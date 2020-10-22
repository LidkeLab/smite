classdef ChannelRegistration < handle
    %ChannelRegistration contains methods for channel registration.
    % This class contains methods for performing channel registration and
    % methods used to interpret/visualize the results.
    %
    % REQUIRES:
    %   MATLAB Image Processing Toolbox
    
    % Created by:
    %   David J. Schodt (Lidke Lab, 2020)
    
    
    properties
        % Variable in FiducialFileNames containing the raw data (char)
        DataVariable = 'sequence';
        
        % Directory containing the fiducial file(s) (char)
        FiducialFileDir
        
        % Filenames of the fiducial data files (cell array of char)
        % FiducialFileNames{1} is always assumed to be the "fixed" or
        % "reference" file, meaning that FiducialFileNames{n>1} will all be
        % registered with respepct to FiducialFileNames{1}.
        FiducialFileNames
        
        % Single molecule fitting structure (see SingleMoleculeFitting)
        % This SMF structure is used to find localizations in the fiducial
        % files specified by FiducialFilePath.  If 
        % TransformationBasis = 'images' this is not needed.
        SMF struct
        
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
        %       the transformType options defined in doc imregtform
        TransformationType char = 'lwm';
    end
    
    properties (SetAccess = protected)
        % Transform computed from the fiducial files (tform object)
        RegistrationTransform
    end
    
    methods
        function [obj] = ChannelRegistration(...
                FiducialFileDir, FiducialFileNames, SMF)
            %ChannelRegistration class constructor.
            % The inputs can be used to set class properties if desired.
            
            % Set inputs to class properties if needed.
            if (exist('FiducialFileDir', 'var') ...
                    && ~isempty(FiducialFileDir))
                obj.FiducialFileDir = FiducialFileDir;
            end
            if (exist('FiducialFileNames', 'var') ...
                    && ~isempty(FiducialFileNames))
                obj.FiducialFileNames = FiducialFileNames;
            end
            if (exist('SMF', 'var') && ~isempty(SMF))
                obj.SMF = SMF;
            end
            
        end
        
        [RegistrationTransform] = findTransform(obj);
        
    end
    
    methods (Static)
        transformImages()
        transformCoords()
        visualizeTransform()
    end
    
    methods (Static, Hidden)
        % These methods aren't expected to be used directly by the user, so
        % it's nice to hide them from view (so they don't distract the
        % user).
        
        findImageTransform()
        findCoordTransform()
        visualizeImageTransform()
        visualizeCoordTransform()
        
    end
    
    
end