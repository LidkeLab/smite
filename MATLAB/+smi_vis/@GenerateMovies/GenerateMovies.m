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
        
        % Raw data to be displayed under trajectories.
        RawData

        % Single Molecule Fitting structure with pixel size and framerate.
        SMF = smi_core.SingleMoleculeFitting;
        
        % Optional SMD containing points to mark in the movie.
        SMD = smi_core.SingleMoleculeData.createSMD();
        
        % Tracking Results structure for the trajectories.
        TR = smi_core.TrackingResults.createTR();
    end
    
    properties (SetAccess = 'protected')
        % Graphics handles to trajectory lines displayed in movie frames.
        % NOTE: This is index the same way as obj.TR, i.e., obj.TR(n)
        %       corresponds to obj.LineHandles(n).
        LineHandles
        
        % Axes in which the movie is prepared if using generateMovie().
        MovieAxes
                
        % Rescaled/cropped version of property RawData (see rescaleData())
        % NOTE: I've made this a protected property so that we can ensure
        %       some of obj.Params are updated when this is written to
        %       (which I only do in obj.rescaleData())
        ScaledData
        
        % MATLAB VideoWriter object used to write a movie to a file.
        VideoObject
    end
    
    properties (Hidden)
        LengthUnitOptions = {'pixels', '$\mu m$'};
        TimeDimensionOptions = {'Frame', 'Time'};
        TimeUnitOptions = {'frames', 's'};
    end
    
    properties (Hidden, Dependent)
        LengthUnitString
        TimeDimensionString
        TimeUnitString
    end
    
    methods
        function obj = GenerateMovies(Params)
            %GenerateMovies is the class constructor.
            if (exist('Params', 'var') && ~isempty(Params))
                obj.Params = Params;
            else
                obj.Params = obj.prepDefaults();
            end
        end
        
        function set.Params(obj, ParamsInput)
            % Ensure that the class property 'Params' is complete (i.e., 
            % that it has all necessary fields set) and that they are typed
            % correctly (e.g., logical, char, ...).
            DefaultParams = obj.prepDefaults();
            obj.Params = smi_helpers.padParams(ParamsInput, DefaultParams);
        end
        
        function [LengthUnitString] = get.LengthUnitString(obj)
            LengthUnitString = smi_helpers.stringMUX(...
                obj.LengthUnitOptions, obj.Params.UnitFlag);
        end
        
        function [TimeDimensionString] = get.TimeDimensionString(obj)
            TimeDimensionString = smi_helpers.stringMUX(...
                obj.TimeDimensionOptions, obj.Params.UnitFlag);
        end
        
        function [TimeUnitString] = get.TimeUnitString(obj)
            TimeUnitString = smi_helpers.stringMUX(...
                obj.TimeUnitOptions, obj.Params.UnitFlag);
        end
        
        rescaleData(obj)
        saveMovie(obj, SavePath)
        generateMovie(obj)
        gui(obj)
        
    end
    
    methods (Hidden)
        % These methods are Hidden because I don't expect the user to
        % access these directly, however there's no harm in leaving them
        % unrestricted.
        setVitalParams(obj)
        prepAxes(obj)
    end
    
    methods (Static)
        playMovie(PlotAxes, TR, ScaledData, Params, SMD, VideoObject)
        [Params] = prepDefaults();
    end
    
    methods (Static, Hidden)
        % These methods are Hidden because I don't expect the user to
        % access these directly, however there's no harm in leaving them
        % unrestricted.
        [LineHandles] = plotTrajectories(PlotAxes, ...
            TR, FrameRange, MaxTrajLength, Color, varargin);
        [LineHandles] = makeFrame(PlotAxes, TR, ScaledData, ...
            Params, SMD, Frame);
    end
    
    
end