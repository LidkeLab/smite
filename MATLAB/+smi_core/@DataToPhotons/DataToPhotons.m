classdef DataToPhotons < handle
    %DataToPhotons contains static methods to convert raw data to photons.
    %
    % This class contains static methods associated with the gain and
    % offset correction needed to convert raw data from the camera (arrays
    % given in Analog to Digital Units (ADU)) to units of photons.  The
    % main usage of this class is shown in the EXAMPLE USAGE section below.
    %
    % EXAMPLE USAGE:
    %   Given RawData in units of ADU, and an SMF structure with the fields
    %   SMF.Data.CameraGain, SMF.Data.CameraOffset, and
    %   SMF.Data.CameraReadNoise (see SingleMoleculeFitting class for
    %   details), you can convert RawData and CameraReadNoise to units of
    %   photons and photons^2, respectively, as follows:
    %       [~, RawDataConverted, CameraReadNoiseConverted] = ...
    %           smi_core.DataToPhotons(SMF, ...
    %           RawData, RawDataROI, CalibrationROI, true);
    %   Alternatively, you can prepare the class for usage (and set class
    %   parameters immediately) as follows:
    %       DTP = smi_core.DataToPhotons(SMF, ...
    %           RawData, RawDataROI, CalibrationROI);
    %       [Data, ReadNoise] = DTP.convertData();
    %
    % REQUIRES:
    %   Image Processing Toolbox
    %   Statistics and Machine Learning Toolbox
    
    
    properties
        % obj.RawData converted to units of photons (float array)
        CorrectedData {mustBeNumeric(CorrectedData)}
        
        % obj.CameraReadNoise converted to units of photons (float array)
        CorrectedReadNoise {mustBeNumeric(CorrectedReadNoise)}
        
        % Data that is to be gain/offset corrected (float array)
        RawData {mustBeNumeric(RawData)}
        
        % Read noise that is to be gain/offset corrected (float array)
        ReadNoise {mustBeNumeric(ReadNoise)}
        
        % Region of interest of the raw data (float array)
        % (see obj.convertToPhotons() for details/usage)
        RawDataROI {mustBeNumeric(RawDataROI)}
        
        % Gain of the camera used to collect RawData (float array)(ADU/e-)
        CameraGain {mustBeNumeric(CameraGain)}
        
        % Offset of the camera used to collect RawData (float array)(ADU)
        CameraOffset {mustBeNumeric(CameraOffset)}
        
        % Read noise of the camera used to collect Raw Data (float array)(ADU^2)
        CameraReadNoise {mustBeNumeric(CameraReadNoise)}
        
        % Region of interest of the gain/offset arrays (float array)
        % (see obj.convertToPhotons() for details/usage)
        CalibrationROI {mustBeNumeric(CalibrationROI)}
        
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
    
    methods
        function [obj, Data, ReadNoise] = DataToPhotons(SMF, ...
                RawData, RawDataROI, CalibrationROI, Verbose, AutoRun)
            %DataToPhotons is the class constructor.
            % This constructor has several optional inputs which are set to
            % class properties.  If you provide all inputs (e.g., SMF,
            % RawData, ...) and set AutoRun = 1, this class will
            % automatically run obj.convertData() on your data.  The
            % results will be stored in the appropriate class properties,
            % but you can also retrieve them by specifying 2nd and 3rd
            % outputs from this constructor. Note that the input
            % CalibrationROI is unused when SMF.Data.CalibrationFilePath is
            % defined. Instead, we'll attempt to set that field with the
            % value stored in that file.
            
            % Set defaults if needed.
            if (~exist('AutoRun', 'var') || isempty(AutoRun))
                AutoRun = 0;
            end
            
            % Set class properties based on the inputs.
            AllFieldsSet = true;
            if (exist('Verbose', 'var') && ~isempty(Verbose))
                obj.Verbose = Verbose;
            end
            if (exist('CalibrationROI', 'var') && ~isempty(CalibrationROI))
                % NOTE: This may be overwritten below if a calibration file
                %       is present and the camera type is specified as
                %       'SCMOS'. If that happens, it was done on purpose!
                %       It's best to use the value in the calibration file
                %       rather than what the user has entered.
                obj.CalibrationROI = CalibrationROI;
                if (obj.Verbose > 2)
                    fprintf(['\tDataToPhotons constructor: ', ...
                        'Input CalibrationROI structure stored as a ', ...
                        'class property.\n'])
                end
            end
            if (exist('SMF', 'var') && ~isempty(SMF))
                % Reload the SMF to ensure it has all required properties.
                SMF = smi_core.SingleMoleculeFitting.reloadSMF(SMF);
                
                % Set the desired SMF fields to class properties.
                if (isempty(SMF.Data.CalibrationFilePath) ...
                        || strcmpi(SMF.Data.CameraType, 'EMCCD'))
                    % If no filepath is given, or if the camera is
                    % specified as an EMCCD, we don't need to load a
                    % calibration file.
                    obj.CameraGain = SMF.Data.CameraGain;
                    obj.CameraOffset = SMF.Data.CameraOffset;
                    obj.CameraReadNoise = SMF.Data.CameraReadNoise;
                    if (obj.Verbose > 2)
                        fprintf(['\tDataToPhotons constructor: ', ...
                            'Camera gain, offset, and read-noise ', ...
                            'provided in input SMF have been stored ', ...
                            'as class properties.\n'])
                    end
                else
                    % Attempt to load the calibration data.
                    [obj.CameraGain, obj.CameraOffset, ...
                        obj.CameraReadNoise, obj.CalibrationROI] = ...
                        smi_core.LoadData.loadDataCalibration(SMF);
                    if (obj.Verbose > 2)
                        fprintf(['\tDataToPhotons constructor: ', ...
                            'Camera gain, offset, and read-noise ', ...
                            'loaded from file specified by ', ...
                            'SMF.Data.CalibrationFilePath.\n'])
                    end
                end
            else
                AllFieldsSet = false;
            end
            if (exist('RawData', 'var') && ~isempty(RawData))
                obj.RawData = RawData;
                if (obj.Verbose > 2)
                    fprintf(['\tDataToPhotons constructor: ', ...
                        'Input RawData stored as a class property.\n'])
                end
            else
                AllFieldsSet = false;
            end
            if (exist('RawDataROI', 'var') && ~isempty(RawDataROI))
                obj.RawDataROI = RawDataROI;
                if (obj.Verbose > 2)
                    fprintf(['\tDataToPhotons constructor: ', ...
                        'Input RawDataROI stored as a class property.\n'])
                end
            else
                AllFieldsSet = false;
            end
            
            % Run obj.convertToPhotons() if requested.
            Data = [];
            ReadNoise = [];
            if (AutoRun && AllFieldsSet)
                if (obj.Verbose > 2)
                    fprintf(['\tDataToPhotons constructor: ', ...
                        'Auto-running obj.convertData()...\n'])
                end
                obj.convertData();
            end
        end
        
        [CorrectedData, CorrectedReadNoise] = convertData(obj);
        
    end
    
    methods (Static)
        [CorrectedData, CorrectedReadNoise] = convertToPhotons(RawData, ...
            CameraGain, CameraOffset, CameraReadNoise, ...
            RawDataROI, CalibrationROI);
        [Success] = unitTest();
    end
    
    
end