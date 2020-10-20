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
        applyImageTransform()
        applyCoordTransform()
        visualizeImageTransform()
        visualizeCoordTransform()
    end
    
    methods (Static, Hidden)
        % These methods aren't expected to be used directly by the user, so
        % it's nice to hide them from view (so they don't distract the
        % user).
        
        findImageTransform()
        findCoordTransform()
    end
    
    
end