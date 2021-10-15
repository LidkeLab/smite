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
        Params struct = struct([])
        
        % Raw data displayed under trajectories. (YSizexXSizex3xNFrames)
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
        
        % Axes in which the movie is prepared if using obj.gui().
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
        DispPlotsOptions = {'DatasetNum', 'FrameNum', 'X', 'Y', 'Z', ...
            'X_SE', 'Y_SE', 'Z_SE', ...
            'Photons', 'Photons_SE', 'Bg', 'Bg_SE', ...
            'PSFSigma', 'PSFSigmaX', 'PSFSigmaY', ...
            'PSFSigma_SE', 'PSFSigmaX_SE', 'PSFSigmaY_SE', ...
            'PValue', 'LogLikelihood', 'ThreshFlag', 'ConnectID'}
        LengthUnitOptions = {'pixels', '$\mu m$'};
        TimeDimensionOptions = {'Frame', 'Time'};
        TimeUnitOptions = {'frames', 's'};
    end
    
    properties (Hidden, Dependent)
        LengthUnitString
        TimeDimensionString
        TimeUnitString
    end
    
    properties (Hidden, GetAccess = 'protected')
        % Internally modified version of the user set field 'TR'.
        TRInternal
    end
    
    methods
        function obj = GenerateMovies()
            %GenerateMovies is the class constructor.
            obj.Params = obj.prepDefaults();
        end
        
        function set.Params(obj, ParamsInput)
            % Ensure that the class property 'Params' is complete (i.e.,
            % that it has all necessary fields set) and that they are typed
            % correctly (e.g., logical, char, ...).
            DefaultParams = obj.prepDefaults();
            obj.Params = smi_helpers.padStruct(ParamsInput, DefaultParams);
        end
        
        function set.TR(obj, TRInput)
            % Update TRMod whenever the user resets the TR.
            obj.TR = TRInput;
            obj.TRInternal = TRInput;
        end
        
        function [LengthUnitString] = get.LengthUnitString(obj)
            LengthUnitString = smi_helpers.arrayMUX(...
                obj.LengthUnitOptions, obj.Params.UnitFlag);
        end
        
        function [TimeDimensionString] = get.TimeDimensionString(obj)
            TimeDimensionString = smi_helpers.arrayMUX(...
                obj.TimeDimensionOptions, obj.Params.UnitFlag);
        end
        
        function [TimeUnitString] = get.TimeUnitString(obj)
            TimeUnitString = smi_helpers.arrayMUX(...
                obj.TimeUnitOptions, obj.Params.UnitFlag);
        end
        
        prepRawData(obj)
        saveMovie(obj, SavePath)
        generateMovie(obj)
        gui(obj)
        
    end
    
    methods (Hidden)
        % These methods are Hidden because I don't expect the user to
        % access these directly, however there's no harm in leaving them
        % unrestricted.
        setVitalParams(obj)
        prepAxes(obj, PlotAxes)
    end
    
    methods (Static)
        playMovie(PlotAxes, TR, ScaledData, Params, SMF, SMD, VideoObject)
        [Params] = prepDefaults();
    end
    
    methods (Static, Hidden)
        % These methods are Hidden because I don't expect the user to
        % access these directly, however there's no harm in leaving them
        % unrestricted.
        [LineHandles] = plotTrajectories(PlotAxes, ...
            TR, FrameRange, MaxTrajLength, Color, varargin);
        [LineHandles] = makeFrame(PlotAxes, TR, ScaledData, ...
            Params, SMF, SMD, Frame);
        [Params] = defineCropROI(TR, Params);
        [RawData] = cropRawData(RawData, Params);
        addTimeStamp(PlotAxes, Frame, FrameRate, Params)
    end
    
    
end