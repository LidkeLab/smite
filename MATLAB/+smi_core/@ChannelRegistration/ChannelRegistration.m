classdef ChannelRegistration < handle
    %ChannelRegistration contains methods for channel registration.
    % This class contains methods for performing channel registration and
    % methods used to interpret/visualize the results.
    
    % Created by:
    %   David J. Schodt (Lidke Lab, 2020)
    
    properties
        % Directory containing the fiducial file(s) (cell array of char)
        FiducialFileDir
        
        % Filename(s) of the fiducial data file(s) (cell array of char)
        FiducialFilePath
        
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
        function [obj] = ChannelRegistration()
            %ChannelRegistration class constructor.
        end
        
        findTransform(obj)
        
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