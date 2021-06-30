classdef GenerateMovies < handle
    %GenerateMovies contains methods for generating movies.
    % This class consists of methods useful for generating movies,
    % particularly for single-particle tracking data.
    
    % Created by:
    %   David J. Schodt (Lidke Lab, 2021)
    
    
    properties
        % Structure of parameters enforced in the generated movie.
        % See prepDefaults() for a description of the parameter options and
        % playMovie() for usage.
        Params
    end
    
    properties (SetAccess = 'protected')
        VideoObject
    end
    
    methods
        function obj = GenerateMovies(Params)
            %GenerateMovies is the class constructor.
            if (exist('Params', 'var') && ~isempty(Params))
                obj.Params = Params;
            end
        end
        
        function set.Params(obj, ParamsInput)
            %set method for the class property 'Params'.
            % This method ensures that the class property 'Params' is
            % complete, i.e., that it has all necessary fields set.
            DefaultParams = obj.prepDefaults();
            obj.Params = smi_helpers.padParams(ParamsInput, DefaultParams);
        end
        
        generateMovie(obj)
        
    end
    
    methods (Static)
        playMovie(PlotAxes, TR, RawData, Params, SMD, VideoObject)
        Params = prepDefaults();
    end
    
    methods (Static, Hidden)
        % These methods are Hidden because I don't expect the user to
        % access these directly, however there's no harm in leaving them
        % unrestricted.
        [LineHandles] = plotTrajectories(PlotAxes, ...
            TR, FrameRange, MaxTrajLength, Color, varargin);
    end
    
    
end