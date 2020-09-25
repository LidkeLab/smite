classdef SingleMoleculeFittingTest < handle
% SingleMoleculeFitting A class defining the Single Molecule Fitting structure
%
% The SMF structure is a structure of structures that collectively contain
% all parameters required to go from raw data to an SMD results structure.
% The SMF structure is an input of many smi methods. It
% intended to be extensible to enable new analysis tools and methods.
% The SMF class implements tools for working with SMF structures,
% but the data structure itself is not an object of the class.
%
% Parameters of sub-structures are explained in more detail in
% the classes and methods that use them.  An incomplete list of classes
% that use each sub-structure is listed in {}.
%
% The SMF structure has the following sub-structures and fields:
%
% SMF:  Fields that are structures: 
%
% Data:             {LoadData}
%   FileName:       File name (char array or cell array of char array)
%   FileDir:        File directory (char array or cell array of char array)
%   ResultsDir:     Results directory (char array)(Default='FileDir/Results')
%   CameraType:     'EMCCD','SCMOS' (Default='EMCCD')
%   CameraGain:     Camera Gain, scalar or image (Default=1)
%   CameraOffset:   Camera Offset, scalar or image (Default=0)
%   CameraNoise:    Camera readnoise, scalar or image (Default=0)
%   FrameRate:      Data Collection Frame Rate (1/s) 
%   PixelSize:      Camera back-projected pixel size (micron)   
%
% BoxFinding:       {FindROI}
%   BoxSize:        Linear box size for fitting (Pixels)(Default=7)
%   BoxOverlap:     Overlap of boxes allowed (Pixels)(Default=2)
%   MinPhotons:     Minimum number of photons from emitter (Default=200)
%
% Fitting           {GaussMLE}
%   PSFSigma:   Initial or fixed Sigma of 2D Gaussian PSF Model (Pixels)(Default=1)
%   FitType:    See fit class for options  (Default='XYNB')
%   Iterations: Newton Raphson iterations (Default=20)
%   ZFitStruct: Structure for astigmatic fitting:
%       Ax:         Asigmatism fit parameter (see GaussMLE)         
%       Ay:         Asigmatism fit parameter (see GaussMLE)
%       Bx:         Asigmatism fit parameter (see GaussMLE)
%       By:         Asigmatism fit parameter (see GaussMLE)
%       Gamma:      Asigmatism fit parameter (see GaussMLE)
%       D:          Asigmatism fit parameter (see GaussMLE)
%
% Thresholding      {ThresholdFits,SRA}
%   On              Perform thresholding? (Default=true)
%   MaxSE_XY:       Maximum allowed precision in x,y (Pixels)(Default=.2)
%   MaxZ_SE:        Maximum allowed precision in z (Microns)(Default=.5)
%   LoS:            Minimum accepted p-value from fit (Default=.01)
%   MinPSFSigma:    Minimum PSF Sigma from fit (Pixels)(Default=.5);
%   MaxPSFSigma:    Maximum PSF Sigma from fit (Pixels)(Default=2);
%   MinPhotons:     Minimum accepted photons from fit (Default=100)
%   MaxBg:          Maximum background accepted from fit (Default=Inf)
%
% FrameConnection:  {FrameConnect,SRA}
%   On              Perform frame connection? (Default=true)
%   MaxSeparation:  Maximum separation for connection (Pixels)(Default=1)
%   MaxFrameGap:    Maximum frame gap for connection (Frames)(Default=4)
%   LoS:            Minimum accepted p-value for connection (Default=.01)
%
% DriftCorrection   {DriftCorrection,SRA}
%  On               Perform drift correction? (Default=true)
%  L_intra          Intra-dataset threshold (Pixel)(Default=1)
%  L_inter          Inter-dataset threshold (Pixel)(Default=2)
%  PixelSizeZUnit   X/Y pixel size (3D drift correction) (um)(Default=0.1)
%  PDegree          Degree intra-dataset fitting poly for drift rate (Default=1)
%  Init_inter       Init wrt previous dataset for inter-dataset DC (Default=0)
% 
% Tracking          {SPT}
%   Method:         Type of method used for tracking (Default='smi_spt')
%   D:              Diffusion Constant (Pixels^2/Frame) (Default=1)
%   K_on:           Off to On Rate (Frame^-1) (Default=.1)
%   K_off:          On to Off Rate (Frame^-1) (Default=.1)
%   MaxFrameGap:    Maximum frame gap for Gap Closing (Pixels) (Default=10)
%   MaxDist:        Maximum distance gap for Gap Closing (Pixels) (Default=10)
%   MinTrackLength  Minimum track length of trajectory (Frames) (Default=3)
%

% created by:
% Keith Lidke, Hanieh Mazloom-Farsibaf, David Schodt. Lidke Lab 2018
    
    properties
        Data struct {mustBeNonempty} = struct();
        BoxFinding struct {mustBeNonempty} = struct();
        Fitting struct {mustBeNonempty} = struct();
        Thresholding struct {mustBeNonempty} = struct();
        FrameConnection struct {mustBeNonempty} = struct();
        DriftCorrection struct {mustBeNonempty} = struct();
        Tracking struct {mustBeNonempty} = struct();
    end
    
    methods
        
        function obj = SingleMoleculeFittingTest()
            % Class constructor used to set default property values.
            
            %Data
            obj.Data.FileName='';
            obj.Data.FileDir='';
            obj.Data.ResultsDir='';
            obj.Data.CameraType='EMCCD';
            obj.Data.CameraGain=1;
            obj.Data.CameraOffset=0;
            obj.Data.CameraReadNoise=0;
            obj.Data.FrameRate=1;
            obj.Data.PixelSize=0.1;
            
            %BoxFinding
            obj.BoxFinding.BoxSize=7;
            obj.BoxFinding.BoxOverlap=2;
            obj.BoxFinding.MinPhotons=200;
            
            %Fitting
            obj.Fitting.PSFSigma=1;
            obj.Fitting.FitType='XYNB';
            obj.Fitting.Iterations=20;
            obj.Fitting.ZFitStruct.Ax=[];
            obj.Fitting.ZFitStruct.Ay=[];
            obj.Fitting.ZFitStruct.Bx=[];
            obj.Fitting.ZFitStruct.By=[];
            obj.Fitting.ZFitStruct.Gamma=[];
            obj.Fitting.ZFitStruct.D=[];

            %Thresholding
            obj.Thresholding.On=true;
            obj.Thresholding.MaxXY_SE=.2;
            obj.Thresholding.MaxZ_SE=.05;
            obj.Thresholding.MinPValue=.01;
            obj.Thresholding.MinPSFSigma=0.5;
            obj.Thresholding.MaxPSFSigma=2;
            obj.Thresholding.MinPhotons=100;
            obj.Thresholding.MaxBg=Inf;

            %FrameConnection
            obj.FrameConnection.On=true;
            obj.FrameConnection.MaxSeparation=1; % pixels 
            obj.FrameConnection.MaxFrameGap=4; % frames
            obj.FrameConnection.LoS=.01;

            %DriftCorrection
            obj.DriftCorrection.On = true;
            obj.DriftCorrection.L_intra = 1; % pixel
            obj.DriftCorrection.L_inter = 2; % pixel
            obj.DriftCorrection.PixelSizeZUnit = 0.1; % um
            obj.DriftCorrection.PDegree = 1;
            obj.DriftCorrection.Init_inter = 0;
            
            %Tracking
            obj.Tracking.TrackMethods='SMA_SPT';
            obj.Tracking.D=1;
            obj.Tracking.K_on=.1;
            obj.Tracking.K_off=.1;
            obj.Tracking.MaxFrameGap=10;
            obj.Tracking.MaxDist=10;
            obj.Tracking.MinTrackLength=3;
            
        end
        
        function set.Data(obj, DataInput)
            % This is a set method for the class property Data.
            
            % Validate the inputs given in DataInput.
            if isfield(DataInput, 'FileName')
                if ~(ischar(DataInput.FileName) ...
                        || isstring(DataInput.FileName) ...
                        || iscell(DataInput.FileName))
                    error(['''SMF.Data.FileName'' must be of type ', ...
                        'char, string, or cell.'])
                end
            end
            if isfield(DataInput, 'FileDir')
                if ~(ischar(DataInput.FileDir) ...
                        || isstring(DataInput.FileDir) ...
                        || iscell(DataInput.FileDir))
                    error(['''SMF.Data.FileDir'' must be of type ', ...
                        'char, string, or cell.'])
                elseif (isfield(obj.Data, 'ResultsDir') ...
                            && isempty(obj.Data.ResultsDir))
                        obj.Data.ResultsDir = fullfile(...
                            DataInput.FileDir, 'Results');
                end
            end
            if isfield(DataInput, 'ResultsDir')
                if ~(ischar(DataInput.ResultsDir) ...
                        || isstring(DataInput.ResultsDir))
                    error(['''SMF.Data.ResultsDir'' must be of type ', ...
                        'char or string.'])
                end
            end
            if isfield(DataInput, 'CameraType')
                if ~ismember(DataInput.CameraType, {'EMCCD', 'SCMOS'})
                    error(['''SMF.Data.CameraType'' must be either ', ...
                        '''EMCCD'' or ''SCMOS'''])
                end
            end
            if isfield(DataInput, 'CameraGain')
                if ~isnumeric(DataInput.CameraGain)
                    error('''SMF.Data.CameraGain'' must be numeric.')
                end
            end
            if isfield(DataInput, 'CameraOffset')
                if ~isnumeric(DataInput.CameraOffset)
                    error('''SMF.Data.CameraOffset'' must be numeric.')
                end
            end
            if isfield(DataInput, 'CameraReadNoise')
                if ~isnumeric(DataInput.CameraReadNoise)
                    error('''SMF.Data.CameraReadNoise'' must be numeric.')
                end
            end
            if isfield(DataInput, 'FrameRate')
                if ~isnumeric(DataInput.FrameRate)
                    error('''SMF.Data.FrameRate'' must be numeric.')
                end
            end
            if isfield(DataInput, 'PixelSize')
                if ~isnumeric(DataInput.PixelSize)
                    error('''SMF.Data.PixelSize'' must be numeric.')
                end
            end
            
            % Set the input fields as class properties.
            InputFields = fieldnames(DataInput);
            for ff = 1:numel(InputFields)
                obj.Data.(InputFields{ff}) = DataInput.(InputFields{ff});
            end
        end
        
        function set.BoxFinding(obj, BoxFindingInput)
            % This is a set method for the class property BoxFinding.
            
            % Validate the inputs given in BoxFindingInput.
            if isfield(BoxFindingInput, 'BoxSize')
                if mod(BoxFindingInput.BoxSize, 1)
                    error(['''SMF.BoxFinding.BoxSize'' ', ...
                        'must be an integer.'])
                end
            end
            if isfield(BoxFindingInput, 'BoxOverlap')
                if mod(BoxFindingInput.BoxOverlap, 1)
                    error(['''SMF.BoxFinding.BoxOverlap'' ', ...
                        'must be an integer.'])
                end
            end
            if isfield(BoxFindingInput, 'MinPhotons')
                if ~isnumeric(BoxFindingInput.MinPhotons)
                    error(['''SMF.BoxFinding.MinPhotons'' ', ...
                        'must be numeric.'])
                end
            end
            
            % Set the input fields as class properties.
            InputFields = fieldnames(BoxFindingInput);
            for ff = 1:numel(InputFields)
                obj.BoxFinding.(InputFields{ff}) = ...
                    BoxFindingInput.(InputFields{ff});
            end
        end
        
        function set.Fitting(obj, FittingInput)
            % This is a set method for the class property Fitting.
            
            % Validate the inputs given in FittingInput.
            if isfield(FittingInput, 'PSFSigma')
                if ~isnumeric(FittingInput.PSFSigma)
                    error('''SMF.Fitting.PSFSigma'' must be numeric.')
                end
            end
            if (isfield(FittingInput, 'FitType') ...
                    && ~ismember(FittingInput.FitType, ...
                    {'XYNB', 'XYNBS', 'XYNBSXSY', 'XYZNB'}))
                error(['SMF.Fitting.FitType must be one of ', ...
                    'XYNB, XYNBS, XYNBSXSY, or XYZNB'])
            end
            if isfield(FittingInput, 'Iterations')
                if mod(FittingInput.Iterations, 1)
                    error('''SMF.Fitting.Iterations'' must be an integer.')
                end
            end
            if isfield(FittingInput, 'ZFitStruct')
                if ~isstruct(FittingInput.ZFitStruct)
                    error(['''SMF.Fitting.ZFitStruct'' must be of ', ...
                        'type struct.'])
                end
            end
            
            % Set the input fields as class properties.
            InputFields = fieldnames(FittingInput);
            for ff = 1:numel(InputFields)
                obj.Fitting.(InputFields{ff}) = ...
                    FittingInput.(InputFields{ff});
            end
        end
        
        function set.Thresholding(obj, ThresholdingInput)
            % This is a set method for the class property Thresholding.
            
            % Validate the inputs given in ThresholdingInput.
            if isfield(ThresholdingInput, 'On')
                if ~(islogical(ThresholdingInput.On) ...
                        || isnumeric(ThresholdingInput.On))
                    error(['''SMF.Thresholding.On'' must be logical ', ...
                        'or interpretable as logical (numeric).'])
                end
            end
            if isfield(ThresholdingInput, 'MaxXY_SE')
                if ~isnumeric(ThresholdingInput.MaxXY_SE)
                    error('''SMF.Thresholding.MaxXY_SE'' must be numeric.')
                end
            end
            if isfield(ThresholdingInput, 'MaxZ_SE')
                if ~isnumeric(ThresholdingInput.MaxZ_SE)
                    error('''SMF.Thresholding.MaxZ_SE'' must be numeric.')
                end
            end
            if isfield(ThresholdingInput, 'MinPValue')
                if ~isnumeric(ThresholdingInput.MinPValue)
                    error(['''SMF.Thresholding.MinPValue'' ', ...
                        'must be numeric.'])
                end
            end
            if isfield(ThresholdingInput, 'MinPSFSigma')
                if ~isnumeric(ThresholdingInput.MinPSFSigma)
                    error(['''SMF.Thresholding.MinPSFSigma'' ', ...
                        'must be numeric.'])
                end
            end
            if isfield(ThresholdingInput, 'MaxPSFSigma')
                if ~isnumeric(ThresholdingInput.MaxPSFSigma)
                    error(['''SMF.Thresholding.MaxPSFSigma'' ', ...
                        'must be numeric.'])
                end
            end
            if isfield(ThresholdingInput, 'MinPhotons')
                if ~isnumeric(ThresholdingInput.MinPhotons)
                    error(['''SMF.Thresholding.MinPhotons'' ', ...
                        'must be numeric.'])
                end
            end
            if isfield(ThresholdingInput, 'MaxBg')
                if ~isnumeric(ThresholdingInput.MaxBg)
                    error('''SMF.Thresholding.MaxBg'' must be numeric.')
                end
            end
            
            % Set the input fields as class properties.
            InputFields = fieldnames(ThresholdingInput);
            for ff = 1:numel(InputFields)
                obj.Thresholding.(InputFields{ff}) = ...
                    ThresholdingInput.(InputFields{ff});
            end
        end
        
        function set.FrameConnection(obj, FCInput)
            % This is a set method for the class property FrameConnection.
            
            % Validate the inputs given in FCInput.
            if isfield(FCInput, 'On')
                if ~(islogical(FCInput.On) || isnumeric(FCInput.On))
                    error(['''SMF.FrameConnection.On'' must be ', ...
                        'logical or interpretable as logical (numeric).'])
                end
            end
            if isfield(FCInput, 'MaxSeparation')
                if ~isnumeric(FCInput.MaxSeparation)
                    error(['''SMF.FrameConnection.MaxSeparation'' ', ...
                        'must be numeric.'])
                end
            end
            if isfield(FCInput, 'MaxFrameGap')
                if mod(FCInput.MaxFrameGap, 1)
                    error(['''SMF.FrameConnection.MaxFrameGap'' ', ...
                        'must be an integer.'])
                end
            end
            if isfield(FCInput, 'LoS')
                if ~isnumeric(FCInput.LoS)
                    error('''SMF.FrameConnection.LoS'' must be numeric.')
                end
            end
            
            % Set the input fields as class properties.
            InputFields = fieldnames(FCInput);
            for ff = 1:numel(InputFields)
                obj.FrameConnection.(InputFields{ff}) = ...
                    FCInput.(InputFields{ff});
            end
        end
        
        function set.DriftCorrection(obj, DCInput)
            % This is a set method for the class property DriftCorrection.
            
            % Validate the inputs given in DCInput.
            if isfield(DCInput, 'On')
                if ~(islogical(DCInput.On) || isnumeric(DCInput.On))
                    error(['''SMF.DriftCorrection.On'' must be ', ...
                        'logical or interpretable as logical (numeric).'])
                end
            end
            if isfield(DCInput, 'L_intra')
                if ~isnumeric(DCInput.L_intra)
                    error(['''SMF.DriftCorrection.L_intra'' ', ...
                        'must be numeric.'])
                end
            end
            if isfield(DCInput, 'L_inter')
                if ~isnumeric(DCInput.L_inter)
                    error(['''SMF.DriftCorrection.L_inter'' ', ...
                        'must be numeric.'])
                end
            end
            if isfield(DCInput, 'PixelSizeZUnit')
                if ~isnumeric(DCInput.PixelSizeZUnit)
                    error(['''SMF.DriftCorrection.PixelSizeZUnit'' ', ...
                        'must be numeric.'])
                end
            end
            if isfield(DCInput, 'PDegree')
                if mod(DCInput.PDegree, 1)
                    error(['''SMF.DriftCorrection.PDegree'' ', ...
                        'must be an integer.'])
                end
            end
            if isfield(DCInput, 'Init_inter')
                if ~isnumeric(DCInput.Init_inter)
                    error(['''SMF.DriftCorrection.Init_inter'' ', ...
                        'must be numeric.'])
                end
            end
            
            % Set the input fields as class properties.
            InputFields = fieldnames(DCInput);
            for ff = 1:numel(InputFields)
                obj.DriftCorrection.(InputFields{ff}) = ...
                    DCInput.(InputFields{ff});
            end
        end
        
        function set.Tracking(obj, TrackingInput)
            % This is a set method for the class property Tracking.
            
            % Validate the inputs given in TrackingInput.
            if (isfield(TrackingInput, 'Method') ...
                    && ~ismember(FittingInput.Method, {'SMA_SPT'}))
                error('SMF.Tracking.Method must be one of SMA_SPT,')
            end
            if isfield(TrackingInput, 'D')
                if ~isnumeric(TrackingInput.D)
                    error('''SMF.Tracking.D'' must be numeric.')
                end
            end
            if isfield(TrackingInput, 'K_on')
                if ~isnumeric(TrackingInput.K_on)
                    error('''SMF.Tracking.K_on'' must be numeric.')
                end
            end
            if isfield(TrackingInput, 'K_off')
                if ~isnumeric(TrackingInput.K_off)
                    error('''SMF.Tracking.K_off'' must be numeric.')
                end
            end
            if isfield(TrackingInput, 'MaxFrameGap')
                if ~isnumeric(TrackingInput.MaxFrameGap)
                    error(['''SMF.DriftCorrection.MaxFrameGap'' ', ...
                        'must be logical or interpretable as logical ', ...
                        '(numeric).'])
                end
            end
            if isfield(TrackingInput, 'MaxDist')
                if ~isnumeric(TrackingInput.MaxDist)
                    error('''SMF.Tracking.MaxDist'' must be numeric.')
                end
            end
            if isfield(TrackingInput, 'MinTrackLength')
                if ~isnumeric(TrackingInput.MinTrackLength)
                    error(['''SMF.DriftCorrection.MinTrackLength'' ', ...
                        'must be logical or interpretable as logical ', ...
                        '(numeric).'])
                end
            end
            
            % Set the input fields as class properties.
            InputFields = fieldnames(TrackingInput);
            for ff = 1:numel(InputFields)
                obj.Tracking.(InputFields{ff}) = ...
                    TrackingInput.(InputFields{ff});
            end
        end
        
    end
end

