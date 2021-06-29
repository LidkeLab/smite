classdef GenerateMovies
    %GenerateMovies contains methods for generating movies.
    % This class consists of methods useful for generating movies,
    % particularly for single-particle tracking data.
    
    % Created by:
    %   David J. Schodt (Lidke Lab, 2021)
    
    
    properties
%         % Clip data to range explored by trajectories (Default = false)
%         AutoClip = false;
%         FrameRate = 10;
%         MaxTrajLength = inf;
%         Resolution = 0;
%         SavePath = '';
%         UnitFlag = false;
    end
    
    properties (SetAccess = 'protected')
        VideoObject
    end
    
    methods
        function obj = GenerateMovies()
            %GenerateMovies is the class constructor.
        end
        
        makeMovie(obj)
        
    end
    
    methods (Static)
        playMovie(PlotAxes, TR, RawData, SMD, DisplayParams)
    end
    
    methods (Static, Hidden)
        % These methods are Hidden because I don't expect the user to
        % access these directly, however there's no harm in leaving them
        % unrestricted.
        [LineHandles] = plotTrajectories(PlotAxes, ...
            TR, FrameRange, MaxTrajLength, Color, varargin);
    end
    
    
end