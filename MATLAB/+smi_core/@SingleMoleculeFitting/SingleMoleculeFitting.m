classdef SingleMoleculeFitting < matlab.mixin.Copyable
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
    %   FileName:       File name (cell array of char array)
    %   FileDir:        File directory (char array)
    %   ResultsDir:     Results directory (char array)(Default='FileDir/Results')
    %   AnalysisID:     ID tagged onto saved results (char array)(Default='')
    %   FileType:       Type of data specified by FileName. If using a custom
    %                   extension, you must set this field manually to the true
    %                   underlying file type (e.g., if using a .mat file saved
    %                   as exFile.spt, set obj.Data.FileType = 'mat')
    %                   (char array)(Default set to extension of FileName{1})
    %   DataVariable:   Name of variable saved in FileName which contains the
    %                   raw data. (char array)(Default='sequence')
    %   DatasetList:    List of datasets of the raw data to be analyzed.
    %                   (array of int32)(Default=int32([]))
    %   DatasetMods: Cell array containing datasets to be used/excluded from
    %                analysis (Mods <-> modifiers). This is meant to be the
    %                user-facing lists which define DatasetList, meaning that
    %                this is what would be set in the GUI. DatasetMods{1} will
    %                contain an array of the "inclusion" dataset numbers and
    %                DatasetMods{2} will contain an array of the "exclusion"
    %                datasets. DatasetList will be set elsewhere (e.g.,
    %                smi_core.LoadData) to include the set
    %                   intersect(intersect(1:NDatasets, DatasetMods{1}), ...
    %                   setdiff(1:NDatasets, DatasetMods{2}))
    %                unless DatasetMods{1} is empty, in which case the first
    %                parantheses term is dropped. For example, if
    %                NDatasets = 20, and you only want to analyze datasets 1:5,
    %                you can set DatasetMods{1} = 1:5. If you further decide to
    %                exclude datsaets 2 and 4, you could set
    %                DatasetMods{2} = [2, 4].
    %                (cell array of int32 arrays)(Default={[]; []})
    %   CameraType:     'EMCCD','SCMOS' (Default='EMCCD')
    %                   NOTE: The next 3 quantities (CameraGain, CameraOffset,
    %                   CameraReadNoise) should be scalars if the CameraType is
    %                   'EMCCD' while if 'SCMOS', square arrays taken from the
    %                   CalibrationFile located at the CalibrationFilePath.
    %   CameraGain:     Camera Gain, scalar or image (Default=1)
    %   CameraOffset:   Camera Offset, scalar or image (Default=0)
    %   CameraNoise:    Camera readnoise, scalar or image (Default=0)
    %   CalibrationFilePath: Path to the camera calibration file (Default='')
    %   RegistrationFilePath: Path to channel registration file (Default='')
    %   DataROI:        Region of interest of data file to be used (Default=[])
    %   FrameRate:      Data Collection Frame Rate (1/s)
    %   PixelSize:      Camera back-projected pixel size (micrometers)
    %   SEAdjust:       Standard error inflation per localization (Pixels)(Default=0)
    %
    % BoxFinding:       {FindROI}
    %   BoxSize:        Linear box size for fitting (Pixels)(Default=7)
    %   BoxOverlap:     Overlap of boxes allowed (Pixels)(Default=2)
    %   MinPhotons:     Minimum number of photons from emitter (Default=200)
    %
    % Fitting           {GaussMLE}
    %   PSFSigma:   Initial or fixed Sigma of 2D Gaussian PSF Model (Pixels)(Default=1)
    %   FitType:    See fit class for options  (Default='XYNB')
    %   NParams:    Number of fitting parameters (auto-set based on FitType)
    %   Iterations: Newton Raphson iterations (Default=20)
    %   ZFitStruct: Structure for astigmatic fitting:
    %       Ax:         Astigmatism fit parameter (see GaussMLE)
    %       Ay:         Astigmatism fit parameter (see GaussMLE)
    %       Bx:         Astigmatism fit parameter (see GaussMLE)
    %       By:         Astigmatism fit parameter (see GaussMLE)
    %       Gamma:      Astigmatism fit parameter (see GaussMLE)
    %       D:          Astigmatism fit parameter (see GaussMLE)
    %
    % Thresholding      {ThresholdFits,SRA}
    %   On              Perform thresholding? (Default=true)
    %   MaxXY_SE:       Maximum allowed precision in x,y (Pixels)(Default=.2)
    %   MaxZ_SE:        Maximum allowed precision in z (Microns)(Default=.5)
    %   MinPValue:      Minimum accepted p-value from fit (Default=.01)
    %   AutoThreshLogL: Automatically threshold on LogL and ignore MinPValue (Default = false)
    %   AutoThreshPrctile: Extrema percentile thrown out when computing LogL auto-threshold (Default = 1e-4)
    %   MinPSFSigma:    Minimum PSF Sigma from fit (Pixels)(Default=.5);
    %   MaxPSFSigma:    Maximum PSF Sigma from fit (Pixels)(Default=2);
    %   MinPhotons:     Minimum accepted photons from fit (Default=100)
    %   MaxBg:          Maximum background accepted from fit (Default=Inf)
    %   InMeanMultiplier:   Determines maximum intensity accepted (Default=Inf)
    %   NNMedianMultiplier: Nearest neighbor acceptance region (Default=3)
    %   MinNumNeighbors:    Minimum number of neighbors in above (Default=0)
    %
    % FrameConnection:  {FrameConnect,SRA}
    %   On              Perform frame connection? (Default=true)
    %   Method:         Frame connection method being used (Default='LAP-FC')
    %   MaxSeparation:  Maximum separation for connection (Pixels)(Default=1)
    %   LoS:            Minimum accepted p-value for connection (Default=.01)
    %   MaxFrameGap:    Maximum frame gap for connection (Frames)(Default=5)
    %   NSigmaDev:      SE multiplier for pre-cluster distance threshold (Default=5)
    %   NNearestClusters: Number of clusters used in density estimates (Default=2)
    %   NIterations:    Number of iterative FC attempts when Method=lap-fc (Default=1)
    %   MinNFrameConns  Minimum accepted number of frame connections (Default=1)
    %
    % DriftCorrection   {DriftCorrection,SRA}
    %  On               Perform drift correction? (Default=true)
    %  Method:          Drift correction method being used (Default='DC-KNN')
    %  BFRegistration   Was brightfield registration performed? (Default=true)
    %  L_intra          Intra-dataset threshold (Pixel)(Default=1)
    %  L_inter          Inter-dataset threshold (Pixel)(Default=2)
    %  PixelSizeZUnit   X/Y pixel size (3D drift correction) (um)(Default=0.1)
    %  PDegree          Degree intra-dataset fitting poly for drift rate (Default=1)
    %
    % Tracking          {SPT}
    %   Method:         Type of method used for tracking (Default='CostMatrix')
    %   D:              Diffusion Constant (Pixels^2/Frame) (Default=0.01)
    %   TrajwiseD:      Use traj.-wise value for D (logical)(Default=true)
    %   K_on:           Off to On Rate (Frame^-1) (Default=.9)
    %   K_off:          On to Off Rate (Frame^-1) (Default=.1)
    %   MaxDistFF:      Maximum distance gap for frame-to-frame connection (Pixels)(Default=5)
    %   MaxDistGC:      Maximum distance gap for Gap Closing (Pixels) (Default=10)
    %   MaxFrameGap:    Maximum frame gap for Gap Closing (Pixels) (Default=10)
    %   MinTrackLength: Minimum track length of trajectory (Frames) (Default=3)
    %   NIterMax: Max. number of iterative tracking attempts (Integer)(Default = 5)
    %   NIterMaxBatch: Max. number of batch tracking iterations (Integer)(Default = 5)
    %   MaxRelativeChange: Max. relative param. change to end iterations (Default = 1e-5)
    %   MaxZScoreDist:  Max. abs(z-score) x/y jump size (Default=inf)
    %   MaxZScorePhotons: Max. abs(z-score) for photon diffs. (Default=inf)
    %   TryLowPValueLocs: Try to incorporate low p-val. locs. (Default=false)
    
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
    
    properties (Access = protected)
        % Current version of the SingleMoleculeFitting class/structs.
        % NOTE: Should stay 1.0 until we decide as a group to change it!
        SMFVersion = 1.0;
        
        % This is a cell array of class property names to be present in GUI
        SMFPropertyNames = {'Data', 'BoxFinding', 'Fitting', ...
            'Thresholding', 'FrameConnection', 'DriftCorrection', ...
            'Tracking'};
        
        % This is a structure similar to SMF but field entries define units
        % This is basically an SMF structure but the sub-fields are all
        % char arrays defining the units/datatype/other notes. For example,
        % SMFFieldNotes.Data.CameraGain.Units = 'ADU'. This can also be
        % used to specify a datatype where the unit isn't relevant, e.g.,
        % SMFFieldNotes.Data.FileName.Units = 'cell array'.
        SMFFieldNotes struct = struct();
    end
    
    methods
        
        function obj = SingleMoleculeFitting()
            % Class constructor used to set default property values.
            
            %Data
            obj.Data.FileName={''};
            obj.Data.FileDir='';
            obj.Data.ResultsDir='';
            obj.Data.AnalysisID='';
            obj.Data.FileType='';
            obj.Data.DataVariable='sequence';
            obj.Data.DatasetList=[];
            obj.Data.DatasetMods={[]; []};
            obj.Data.CameraType='EMCCD';
            obj.Data.CameraGain=1;
            obj.Data.CameraOffset=0;
            obj.Data.CameraReadNoise=0;
            obj.Data.CalibrationFilePath='';
            obj.Data.RegistrationFilePath='';
            obj.Data.DataROI=[];
            obj.Data.FrameRate=1;
            obj.Data.PixelSize=0.1;
            obj.Data.SEAdjust=0;
            
            %BoxFinding
            obj.BoxFinding.BoxSize=7;
            obj.BoxFinding.BoxOverlap=2;
            obj.BoxFinding.MinPhotons=200;
            
            %Fitting
            obj.Fitting.PSFSigma=1;
            obj.Fitting.FitType='XYNB';
            obj.Fitting.NParams=4;
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
            obj.Thresholding.AutoThreshLogL=false;
            obj.Thresholding.AutoThreshPrctile=1e-4;
            obj.Thresholding.MinPSFSigma=0.5;
            obj.Thresholding.MaxPSFSigma=2;
            obj.Thresholding.MinPhotons=100;
            obj.Thresholding.MaxBg=Inf;
            obj.Thresholding.InMeanMultiplier=Inf;
            obj.Thresholding.NNMedianMultiplier=3;
            obj.Thresholding.MinNumNeighbors=0;
            
            %FrameConnection
            obj.FrameConnection.On=true;
            obj.FrameConnection.Method='LAP-FC';
            obj.FrameConnection.MaxSeparation=1; % pixels
            obj.FrameConnection.LoS=.01;
            obj.FrameConnection.MaxFrameGap=5; % frames
            obj.FrameConnection.NSigmaDev=5;
            obj.FrameConnection.NNearestClusters=2;
            obj.FrameConnection.NIterations=1;
            obj.FrameConnection.MinNFrameConns=1;
            
            %DriftCorrection
            obj.DriftCorrection.On = true;
            obj.DriftCorrection.Method = 'DC-KNN';
            obj.DriftCorrection.BFRegistration = true;
            obj.DriftCorrection.L_intra = 1; % pixel
            obj.DriftCorrection.L_inter = 2; % pixel
            obj.DriftCorrection.PixelSizeZUnit = 0.1; % um
            obj.DriftCorrection.PDegree = 1;
            
            %Tracking
            obj.Tracking.Method='CostMatrix';
            obj.Tracking.D=0.01;
            obj.Tracking.TrajwiseD=true;
            obj.Tracking.K_on=.9;
            obj.Tracking.K_off=.1;
            obj.Tracking.MaxDistFF=5;
            obj.Tracking.MaxDistGC=10;
            obj.Tracking.MaxFrameGap=10;
            obj.Tracking.MinTrackLength=3;
            obj.Tracking.MaxZScoreDist=inf;
            obj.Tracking.MaxZScorePhotons=inf;
            obj.Tracking.NIterMax=5;
            obj.Tracking.NIterMaxBatch=5;
            obj.Tracking.ParamsHistory={};
            obj.Tracking.MaxRelativeChange=1e-5;
            obj.Tracking.TryLowPValueLocs=false;
            
            % Store a note about the unit/type of various sub-fields.
            % NOTE: Sub-structs (e.g., SMF.Fitting.ZFitStruct, should only
            %       have one note pertaining to the overall structure!).
            obj.SMFFieldNotes.Data.FileName.Units = 'cell array';
            obj.SMFFieldNotes.Data.FileDir.Units = 'char array';
            obj.SMFFieldNotes.Data.ResultsDir.Units = 'char array';
            obj.SMFFieldNotes.Data.AnalysisID.Units = 'char array';
            obj.SMFFieldNotes.Data.FileType.Units = 'char array';
            obj.SMFFieldNotes.Data.DataVariable.Units = 'char array';
            obj.SMFFieldNotes.Data.DatasetList.Units = 'integer array';
            obj.SMFFieldNotes.Data.DatasetMods.Units = 'integer array';
            obj.SMFFieldNotes.Data.CameraType.Units = 'EMCCD, SCMOS';
            obj.SMFFieldNotes.Data.CameraGain.Units = 'ADU / e-';
            obj.SMFFieldNotes.Data.CameraOffset.Units = 'ADU';
            obj.SMFFieldNotes.Data.CameraReadNoise.Units = 'ADU^2';
            obj.SMFFieldNotes.Data.CalibrationFilePath.Units = ...
                'char array';
            obj.SMFFieldNotes.Data.RegistrationFilePath.Units = ...
                'char array';
            obj.SMFFieldNotes.Data.DataROI.Units = ...
                'pixels, [YStart, XStart, YEnd, XEnd, ZStart, ZPeriod]';
            obj.SMFFieldNotes.Data.FrameRate.Units = ...
                'frames / second';
            obj.SMFFieldNotes.Data.PixelSize.Units = 'micrometers';
            obj.SMFFieldNotes.Data.SEAdjust.Units = 'pixels';
            obj.SMFFieldNotes.BoxFinding.BoxSize.Units = 'pixels';
            obj.SMFFieldNotes.BoxFinding.BoxOverlap.Units = 'pixels';
            obj.SMFFieldNotes.BoxFinding.MinPhotons.Units = 'photons';
            obj.SMFFieldNotes.Fitting.PSFSigma.Units = 'pixels';
            obj.SMFFieldNotes.Fitting.FitType.Units = ...
                'XYNB, XYNBS, XYNBSXSY, XYZNB';
            obj.SMFFieldNotes.Fitting.NParams.Units = '(not set by user)';
            obj.SMFFieldNotes.Fitting.Iterations.Units = '';
            obj.SMFFieldNotes.Fitting.ZFitStruct.Units = 'see GaussMLE';
            obj.SMFFieldNotes.Thresholding.On.Units = 'logical';
            obj.SMFFieldNotes.Thresholding.MaxXY_SE.Units = 'pixels';
            obj.SMFFieldNotes.Thresholding.MaxZ_SE.Units = 'pixels';
            obj.SMFFieldNotes.Thresholding.MinPValue.Units = ...
                'number between 0 and 1';
            obj.SMFFieldNotes.Thresholding.AutoThreshLogL.Units = 'logical';
            obj.SMFFieldNotes.Thresholding.AutoThreshPrctile.Units = ...
                'number between 0 and 100';
            obj.SMFFieldNotes.Thresholding.MinPSFSigma.Units = 'pixels';
            obj.SMFFieldNotes.Thresholding.MaxPSFSigma.Units = 'pixels';
            obj.SMFFieldNotes.Thresholding.MinPhotons.Units = 'photons';
            obj.SMFFieldNotes.Thresholding.MaxBg.Units = 'photons';
            obj.SMFFieldNotes.Thresholding.InMeanMultiplier.Units = ...
                'positive number';
            obj.SMFFieldNotes.Thresholding.NNMedianMultiplier.Units = ...
                'positive number';
            obj.SMFFieldNotes.Thresholding.MinNumNeighbors.Units = ...
                'non-negative integer';
            obj.SMFFieldNotes.FrameConnection.On.Units = 'logical';
            obj.SMFFieldNotes.FrameConnection.Method.Units = '';
            obj.SMFFieldNotes.FrameConnection.MaxSeparation.Units = ...
                'pixels';
            obj.SMFFieldNotes.FrameConnection.LoS.Units = ...
                'number between 0 and 1';
            obj.SMFFieldNotes.FrameConnection.MaxFrameGap.Units = 'frames';
            obj.SMFFieldNotes.FrameConnection.NSigmaDev.Units = ...
                'non-negative number';
            obj.SMFFieldNotes.FrameConnection.NNearestClusters.Units = ...
                'positive integer';
            obj.SMFFieldNotes.FrameConnection.NIterations.Units = ...
                'positive integer';
            obj.SMFFieldNotes.FrameConnection.MinNFrameConns.Units = ...
                'positive integer';
            obj.SMFFieldNotes.DriftCorrection.On.Units = 'logical';
            obj.SMFFieldNotes.DriftCorrection.Method.Units = '';
            obj.SMFFieldNotes.DriftCorrection.BFRegistration.Units = ...
                'logical';
            obj.SMFFieldNotes.DriftCorrection.L_intra.Units = 'pixels';
            obj.SMFFieldNotes.DriftCorrection.L_inter.Units = 'pixels';
            obj.SMFFieldNotes.DriftCorrection.PixelSizeZUnit.Units = ...
                'micrometers';
            obj.SMFFieldNotes.DriftCorrection.PDegree.Units = '';
            obj.SMFFieldNotes.Tracking.Method.Units = ...
                'must be CostMatrix for now';
            obj.SMFFieldNotes.Tracking.D.Units = 'pixel^2 / frame';
            obj.SMFFieldNotes.Tracking.TrajwiseD.Units = 'logical';
            obj.SMFFieldNotes.Tracking.K_on.Units = '1 / frame';
            obj.SMFFieldNotes.Tracking.K_off.Units = '1 / frame';
            obj.SMFFieldNotes.Tracking.MaxDistFF.Units = 'pixels';
            obj.SMFFieldNotes.Tracking.MaxDistGC.Units = 'pixels';
            obj.SMFFieldNotes.Tracking.MaxZScoreDist.Units = '';
            obj.SMFFieldNotes.Tracking.MaxZScorePhotons.Units = '';
            obj.SMFFieldNotes.Tracking.MaxFrameGap.Units = 'frames';
            obj.SMFFieldNotes.Tracking.MinTrackLength.Units = ...
                'observations';
            obj.SMFFieldNotes.Tracking.NIterMax.Units = '';
            obj.SMFFieldNotes.Tracking.NIterMaxBatch.Units = '';
            obj.SMFFieldNotes.Tracking.ParamsHistory.Units = '';
            obj.SMFFieldNotes.Tracking.MaxRelativeChange.Units = '';
            obj.SMFFieldNotes.Tracking.TryLowPValueLocs.Units = 'logical';
            
            % Store a 'tip' for certain sub-fields (intended to be
            % displayed in the GUI as a tooltip when hovering over
            % associated GUI elements).
            % NOTE: Sub-structs (e.g., SMF.Fitting.ZFitStruct, should only
            %       have one tip pertaining to the overall structure!).
            % NOTE: If you want multi-line tooltips displayed in the GUI,
            %       use sprintf() here with the \n control.
            obj.SMFFieldNotes.Data.FileName.Tip = ...
                sprintf(['Name of the raw data file(s) you wish to\n', ...
                'analyze. Note that selecting a specific file in\n', ...
                'the popup-menu has no effect.']);
            obj.SMFFieldNotes.Data.FileDir.Tip = ...
                'Name of the directory containing the raw data file(s)';
            obj.SMFFieldNotes.Data.ResultsDir.Tip = ...
                'Name of the directory in which results will be saved';
            obj.SMFFieldNotes.Data.AnalysisID.Tip = ...
                sprintf(['Optional identifier to be tagged onto\n', ...
                'the filenames of saved results.']);
            obj.SMFFieldNotes.Data.FileType.Tip = ...
                sprintf(['Underlying type of the raw data file,\n' ...
                'e.g., ''mat'', ''h5'', ...']);
            obj.SMFFieldNotes.Data.DataVariable.Tip = ...
                'Name of variable in raw data file(s) containing the data';
            obj.SMFFieldNotes.Data.DatasetList.Tip = ...
                sprintf(['Array specifying the dataset number(s) to be\n', ...
                'analyzed---if empty, defaults to all FileNames chosen']);
            obj.SMFFieldNotes.Data.DatasetMods.Tip = ...
                sprintf(['DatasetMods{1} (Include Datasets) is a set\n',...
                'of integers specifying which datasets should be\n', ...
                'included in the analysis. This can be set to empty\n', ...
                'to include all datasets.\n', ...
                'DatasetMods{2} (Exclude Datasets) is a set\n', ...
                'of integers specifying which datasets should be\n', ...
                'strictly excluded from the analysis. This set takes\n',...
                'precedence over the inclusion set']);
            obj.SMFFieldNotes.Data.CameraType.Tip = ...
                sprintf(['Type of camera used to collect the raw data.\n', ...
                'CameraGain, CameraOffset, CameraReadNoise) should be\n', ...
                'scalars if the CameraType is EMCCD while if SCMOS,\n', ...
                'square arrays taken from the file located at the ', ...
                'CalibrationFilePath.']);
            obj.SMFFieldNotes.Data.CameraGain.Tip = ...
                'Gain of the camera used to collect the raw data';
            obj.SMFFieldNotes.Data.CameraOffset.Tip = ...
                'Offset of the camera used to collect the raw data';
            obj.SMFFieldNotes.Data.CameraReadNoise.Tip = ...
                sprintf(['Variance of the read-noise of the camera\n', ...
                'used to collect the raw data']);
            obj.SMFFieldNotes.Data.CalibrationFilePath.Tip = ...
                'Path to the camera calibration file to be used.';
            obj.SMFFieldNotes.Data.RegistrationFilePath.Tip = ...
                sprintf(['Path to a channel registration file ', ...
                'containing a\n transform that will be applied to the data']);
            obj.SMFFieldNotes.Data.DataROI.Tip = ...
                sprintf(['ROI of data to be analyzed (currently only ', ...
                'used in smi.SPT).\nThis should be organized as ', ...
                '[YStart, XStart, YEnd, XEnd, ZStart, ZPeriod]. \nThe ', ...
                'elements ZStart and ZPeriod are used to unmix data ', ...
                'shuffled along their\nthird dimension (e.g., for ', ...
                'stacked two-channel data, ZStart = 1 and\n', ...
                'ZPeriod = 2 would select images 1:2:end).']);
            obj.SMFFieldNotes.Data.FrameRate.Tip = ...
                sprintf(['Acquisition frame rate of the camera used\n', ...
                'to collect the raw data']);
            obj.SMFFieldNotes.Data.PixelSize.Tip = ...
                sprintf(['Pixel size of the camera used to collect\n', ...
                'the raw data, back-projected to the objective\n', ...
                'focal plane.']);
            obj.SMFFieldNotes.Data.SEAdjust.Tip = ...
                'Standard error inflation to be applied to each localization';
            obj.SMFFieldNotes.BoxFinding.BoxSize.Tip = ...
                'See FindROI for usage';
            obj.SMFFieldNotes.BoxFinding.BoxOverlap.Tip = ...
                'See FindROI for usage';
            obj.SMFFieldNotes.BoxFinding.MinPhotons.Tip = ...
                'See FindROI for usage';
            obj.SMFFieldNotes.Fitting.PSFSigma.Tip = ...
                'See GaussMLE for usage';
            obj.SMFFieldNotes.Fitting.FitType.Tip = ...
                'See GaussMLE for usage';
            obj.SMFFieldNotes.Fitting.NParams.Tip = ...
                sprintf(['Number of fitting parameters.  This is\n', ...
                'set automatically based on the FitType.']);
            obj.SMFFieldNotes.Fitting.Iterations.Tip = ...
                'See GaussMLE for usage';
            obj.SMFFieldNotes.Fitting.ZFitStruct.Tip = ...
                'See GaussMLE for usage';
            obj.SMFFieldNotes.Thresholding.On.Tip = ...
                'Indicates whether or not thresholding will be performed';
            obj.SMFFieldNotes.Thresholding.MaxXY_SE.Tip = ...
                sprintf(['Maximum value of X or Y standard error of\n', ...
                'localizations retained after thresholding']);
            obj.SMFFieldNotes.Thresholding.MaxZ_SE.Tip = ...
                sprintf(['Maximum value of Z standard error of\n', ...
                'localizations retained after thresholding']);
            obj.SMFFieldNotes.Thresholding.MinPValue.Tip = ...
                sprintf(['Minimum value of the p-value of\n', ...
                'localizations retained after thresholding. In this\n', ...
                'context, the p-value is bigger for good\n', ...
                'localizations and smaller for bad localizations']);
            obj.SMFFieldNotes.Thresholding.AutoThreshLogL.Tip = ...
                sprintf(['Automatically select a threshold for the ', ...
                'log-likelihood \n using the triangle method and use ', ...
                'it in place of the \nMinPValue threshold.']);
            obj.SMFFieldNotes.Thresholding.AutoThreshPrctile.Tip = ...
                sprintf(['Percentile of extrema of log-likelihood ', ...
                'thrown out\nbefore computing the auto-threshold when ', ...
                'using\nAutoThreshLogL. ' ...
                '(see smi_helpers.triangle_threshold())']);
            obj.SMFFieldNotes.Thresholding.MinPSFSigma.Tip = ...
                sprintf(['Minimum value of PSF sigma of\n', ...
                'localizations retained after thresholding']);
            obj.SMFFieldNotes.Thresholding.MaxPSFSigma.Tip = ...
                sprintf(['Maximum value of PSF sigma of\n', ...
                'localizations retained after thresholding']);
            obj.SMFFieldNotes.Thresholding.MinPhotons.Tip = ...
                sprintf(['Minimum number of photons in\n', ...
                'localizations retained after thresholding']);
            obj.SMFFieldNotes.Thresholding.MaxBg.Tip = ...
                sprintf(['Maximum number of background photons of\n', ...
                'localizations retained after thresholding']);
            obj.SMFFieldNotes.Thresholding.InMeanMultiplier.Tip = ...
                sprintf(['Intensity mean multiplier defining maximum\n', ...
                'photons allowed to retain a localization']);
            obj.SMFFieldNotes.Thresholding.NNMedianMultiplier.Tip = ...
                sprintf(['Standard error mean multiplier defining the\n'...
                'acceptance region for counting numbers of neighbors\n', ...
                '(Do not use on dSTORM data)']);
            obj.SMFFieldNotes.Thresholding.MinNumNeighbors.Tip = ...
                sprintf(['In conjunction with ''NNMedianMultiplier'',\n', ...
                'the minumum number of neighbors that must be in the\n', ...
                'acceptance region to retain a localization\n', ...
                '(Do not use on dSTORM data)']);
            obj.SMFFieldNotes.FrameConnection.On.Tip = ...
                sprintf(['Indicates whether or not frame connection\n', ...
                'will be performed']);
            obj.SMFFieldNotes.FrameConnection.Method.Tip = ...
                sprintf(['Method used to connect repeated\n', ...
                'localizations of the same emitter on event']);
            obj.SMFFieldNotes.FrameConnection.MaxSeparation.Tip = ...
                sprintf(['Maximum separation between two\n', ...
                'localizations such that they can still be\n', ...
                'considered candidates for frame connection']);
            obj.SMFFieldNotes.FrameConnection.LoS.Tip = ...
                sprintf(['Level of Significance compared to\n', ...
                'p-values computed in the frame-connection\n', ...
                'procedure. In this context, a lower value of LoS\n', ...
                'corresponds to a more liberal connection of\n', ...
                'localizations.']);
            obj.SMFFieldNotes.FrameConnection.MaxFrameGap.Tip = ...
                sprintf(['Maximum number of frames separating two\n', ...
                'localizations in time such that they can still be\n', ...
                'considered candidates for frame connection']);
            obj.SMFFieldNotes.FrameConnection.NSigmaDev.Tip = ...
                sprintf(['Localization error multiplier used to set\n', ...
                'pre-clustering separation threshold when\n', ...
                'Method = ''LAP-FC''.']);
            obj.SMFFieldNotes.FrameConnection.NNearestClusters.Tip = ...
                sprintf(['Number of nearest clusters used in local\n', ...
                'emitter density estimation when Method=''LAP-FC''']);
            obj.SMFFieldNotes.FrameConnection.NIterations.Tip = ...
                sprintf(['Number of iterative calls to the LAP-FC\n', ...
                'algorithm when Method=''LAP-FC''.  With each\n', ...
                'iteration, internal parameters are re-estimated\n', ...
                'from the previous iteration''s results.']);
            obj.SMFFieldNotes.FrameConnection.MinNFrameConns.Tip = ...
                sprintf(['Do not retain localizations representing\n', ...
                'this or fewer frame connection sequences\n', ...
                '(Do not use on dSTORM data)']);
            obj.SMFFieldNotes.DriftCorrection.On.Tip = ...
                sprintf(['Indicates whether or not drift correction\n', ...
                'will be performed']);
            obj.SMFFieldNotes.DriftCorrection.Method.Tip = ...
                sprintf(['Method used to correct for emitter drift']);
            obj.SMFFieldNotes.DriftCorrection.BFRegistration.Tip = ...
                sprintf(['Indicates whether or not brightfield\n', ...
                'registration was performed during data collection']);
            obj.SMFFieldNotes.DriftCorrection.L_intra.Tip = ...
                'See DriftCorrection for usage';
            obj.SMFFieldNotes.DriftCorrection.L_inter.Tip = ...
                'See DriftCorrection for usage';
            obj.SMFFieldNotes.DriftCorrection.PixelSizeZUnit.Tip = ...
                'See DriftCorrection for usage';
            obj.SMFFieldNotes.DriftCorrection.PDegree.Tip = ...
                'See DriftCorrection for usage';
            obj.SMFFieldNotes.Tracking.Method.Tip = ...
                sprintf(['Method used for tracking. Note that\n', ...
                'currently there is only one option; this is just\n', ...
                'a placeholder!']);
            obj.SMFFieldNotes.Tracking.D.Tip = ...
                sprintf(['Known/anticipated diffusion constant of\n', ...
                'emitters present in the raw data.']);
            obj.SMFFieldNotes.Tracking.TrajwiseD.Tip = ...
                sprintf(['Use trajectory-wise value for D when', ...
                'iteratively tracking.\n''NIterMax'' or ''NIterMaxBatch''', ...
                'must beat least 2 to use this property']);
            obj.SMFFieldNotes.Tracking.K_on.Tip = ...
                sprintf(['Known/anticipated rate at which dark\n', ...
                'emitters become fluorescent/return from a dark state.']);
            obj.SMFFieldNotes.Tracking.K_off.Tip = ...
                sprintf(['Known/anticipated rate at which emitters\n', ...
                'transition to a dark state.']);
            obj.SMFFieldNotes.Tracking.MaxDistFF.Tip = ...
                sprintf(['Maximum separation between localizations\n', ...
                'such that they can still be considered candidates\n', ...
                'for the frame-to-frame connection procedure']);
            obj.SMFFieldNotes.Tracking.MaxDistGC.Tip = ...
                sprintf(['Maximum separation between localizations\n', ...
                'such that they can still be considered candidates\n', ...
                'for the gap closing procedure']);
            obj.SMFFieldNotes.Tracking.MaxZScoreDist.Tip = ...
                sprintf(['Maximum z-score for jump sizes allowed\n', ...
                'for trajectory connections.']);
            obj.SMFFieldNotes.Tracking.MaxZScorePhotons.Tip = ...
                sprintf(['Maximum z-score for photon differences\n', ...
                'allowed for trajectory connections.']);
            obj.SMFFieldNotes.Tracking.MaxFrameGap.Tip = ...
                sprintf(['Maximum number of frames elapsed between\n', ...
                'localizations such that they can still be\n', ...
                'considered candidates for the gap closing procedure']);
            obj.SMFFieldNotes.Tracking.MinTrackLength.Tip = ...
                sprintf(['Minimum number of observations\n', ...
                '(localizations) a trajectory must have to not be\n', ...
                'culled after tracking']);
            obj.SMFFieldNotes.Tracking.NIterMax.Tip = ...
                'Max. # of iterations permitted if iteratively tracking.';
            obj.SMFFieldNotes.Tracking.NIterMaxBatch.Tip = ...
                ['Max. # of iterations permitted if iteratively ', ...
                'tracking over batches.'];
            obj.SMFFieldNotes.Tracking.ParamsHistory.Tip = ...
                sprintf(['History of parameters used if iteratively ', ...
                'tracking.\nNot to be edited by the user.  If using ', ...
                'smi.SPT.batchTrack(),\nthis will show the history over ', ...
                'the ''NIterMaxBatch''\niterations.  Otherwise, it''s ', ...
                'the history over ''NIterMax''\niterations.']);
            obj.SMFFieldNotes.Tracking.MaxRelativeChange.Tip = ...
                sprintf(['Max. relative change in parameters allowed\n', ...
                'before ending iterative tracking']);
            obj.SMFFieldNotes.Tracking.TryLowPValueLocs.Tip = ...
                sprintf(['Attempt to track with localizations that\n', ...
                'were thresholded based on their p-value.  If those\n', ...
                'localizations are incorporated into trajectories,\n', ...
                'they are kept, otherwise they are discarded.']);
            
            % Check that the protected property 'SMFFieldNames' makes
            % sense, i.e., it doesn't have entries that don't exist as
            % properties.
            PropertyNames = fieldnames(obj);
            if ~all(ismember(obj.SMFPropertyNames, PropertyNames))
                warning(['smi_core.SingleMoleculeFitting: Not all ', ...
                    'entries in the protected property ', ...
                    '''SMFPropertyNames'' exist as class properties. ', ...
                    'Please revise in SingleMoleculeFitting.m'])
            end
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
                elseif (ischar(DataInput.FileName) ...
                        || isstring(DataInput.FileName))
                    DataInput.FileName = {DataInput.FileName};
                elseif (size(DataInput.FileName, 1) ...
                        < size(DataInput.FileName, 2))
                    % I prefer column arrays so I'll transpose if needed.
                    DataInput.FileName = DataInput.FileName.';
                end
                
                % Attempt to set obj.Data.FileType, if needed.
                if (isfield(obj.Data, 'FileType') ...
                        && isempty(obj.Data.FileType) ...
                        && isempty(DataInput.FileType))
                    [~, ~, FileType] = fileparts(DataInput.FileName{1});
                    DataInput.FileType = FileType(2:end);
                end
            end
            if isfield(DataInput, 'FileDir')
                if ~(ischar(DataInput.FileDir) ...
                        || isstring(DataInput.FileDir) ...
                        || iscell(DataInput.FileDir))
                    error(['''SMF.Data.FileDir'' must be of type ', ...
                        'char, string, or cell.'])
                elseif ~isempty(DataInput.FileDir)
                    % Note that if the FileDir isn't being changed, we
                    % don't want to update ResultsDir since that might
                    % overwrite a change to ResultsDir intended by the
                    % user. Similarly, if the current ResultsDir isn't the
                    % expected default value, we should assume the user has
                    % specified some custom directory that shouldn't be
                    % changed.
                    if (isfield(obj.Data, 'FileDir') ...
                            && ~strcmp(DataInput.FileDir, ...
                            obj.Data.FileDir))
                        DefaultResultsDir = ...
                            fullfile(DataInput.FileDir, 'Results');
                        if isfield(obj.Data, 'ResultsDir')
                            OldDefaultResultsDir = ...
                                fullfile(obj.Data.FileDir, 'Results');
                            if (strcmp(obj.Data.ResultsDir, ...
                                    OldDefaultResultsDir) ...
                                    || isempty(obj.Data.ResultsDir))
                                % If this condition isn't met, I'm assuming
                                % the ResultsDir is set to a custom
                                % directory and shouldn't be changed.
                                DataInput.ResultsDir = DefaultResultsDir;
                            end
                        else
                            DataInput.ResultsDir = DefaultResultsDir;
                        end
                    end
                end
            end
            if isfield(DataInput, 'ResultsDir')
                if ~(ischar(DataInput.ResultsDir) ...
                        || isstring(DataInput.ResultsDir))
                    error(['''SMF.Data.ResultsDir'' must be of type ', ...
                        'char or string.'])
                end
            end
            if isfield(DataInput, 'AnalysisID')
                if ~(ischar(DataInput.AnalysisID) ...
                        || isstring(DataInput.AnalysisID))
                    error(['''SMF.Data.AnalysisID'' must be of type ', ...
                        'char or string.'])
                end
            end
            if isfield(DataInput, 'FileType')
                if ~(ischar(DataInput.FileType) ...
                        || isstring(DataInput.FileType))
                    error(['''SMF.Data.FileType'' must be of type ', ...
                        'char or string.'])
                end
            end
            if isfield(DataInput, 'DataVariable')
                if ~(ischar(DataInput.DataVariable) ...
                        || isstring(DataInput.DataVariable))
                    error(['''SMF.Data.DataVariable'' must be of type ',...
                        'char or string.'])
                end
            end
            if isfield(DataInput, 'DatasetList')
                if (~isempty(DataInput.DatasetList) ...
                        && any(mod(DataInput.DatasetList, 1)))
                    % NOTE: This condition allows obj.Data.DatasetList to
                    %       be set to empty, which may cause issues
                    %       depending on how it's used.
                    error(['''SMF.Data.DatasetList must be an array ', ...
                        'of integers.'])
                end
            end
            if isfield(DataInput, 'DatasetMods')
                if (~isempty(DataInput.DatasetMods{1}) ...
                        && any(mod(DataInput.DatasetMods{1}, 1)))
                    error(['''SMF.Data.DatasetMods{1}'' ', ...
                        'must be an array of integers']);
                elseif (size(DataInput.DatasetMods{1}, 1) ...
                        > size(DataInput.DatasetMods{1}, 2))
                    % Row arrays work nicely for use in for loops.
                    DataInput.DatasetMods{1} = DataInput.DatasetMods{1}.';
                end
                if (~isempty(DataInput.DatasetMods{2}) ...
                        && any(mod(DataInput.DatasetMods{2}, 1)))
                    error(['''SMF.Data.DatasetMods{2}'' ', ...
                        'must be an array of integers']);
                elseif (size(DataInput.DatasetMods{2}, 1) ...
                        > size(DataInput.DatasetMods{2}, 2))
                    % Row arrays work nicely for use in for loops.
                    DataInput.DatasetMods{2} = DataInput.DatasetMods{2}.';
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
            if isfield(DataInput, 'CalibrationFilePath')
                if ~(ischar(DataInput.CalibrationFilePath) ...
                        || isstring(DataInput.CalibrationFilePath))
                    error(['''SMF.Data.CalibrationFilePath'' must ', ...
                        'be of type char or string'])
                end
            end
            if isfield(DataInput, 'RegistrationFilePath')
                if ~(ischar(DataInput.RegistrationFilePath) ...
                        || isstring(DataInput.RegistrationFilePath))
                    error(['''SMF.Data.RegistrationFilePath'' must ', ...
                        'be of type char or string'])
                end
            end
            if isfield(DataInput, 'DataROI')
                if (~isempty(DataInput.DataROI) ...
                        && any(mod(DataInput.DataROI, 1)))
                    % NOTE: This condition allows obj.Data.DataROI to
                    %       be set to empty, which may cause issues
                    %       depending on how it's used.
                    error(['''SMF.Data.DataROI must be an array ', ...
                        'of integers.'])
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
            if isfield(DataInput, 'SEAdjust')
                if ~isnumeric(DataInput.SEAdjust)
                    error('''SMF.Data.SEAdjust'' must be numeric.')
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
            if isfield(FittingInput, 'FitType')
                if ismember(FittingInput.FitType, ...
                        {'XYNB', 'XYNBS', 'XYNBSXSY', 'XYZNB'})
                    % Set the NParams field automatically based on this
                    % selection.
                    switch FittingInput.FitType
                        case 'XYNB'
                            FittingInput.NParams = 4;
                        case {'XYNBS', 'XYZNB'}
                            FittingInput.NParams = 5;
                        case 'XYNBSXSY'
                            FittingInput.NParams = 6;
                    end
                else
                    error(['SMF.Fitting.FitType must be one of ', ...
                        'XYNB, XYNBS, XYNBSXSY, or XYZNB'])
                end
            end
            if isfield(FittingInput, 'NParams')
                % Set the NParams field automatically based on the FitType,
                % regardless of what the user tried to enter.
                switch FittingInput.FitType
                    case 'XYNB'
                        FittingInput.NParams = 4;
                    case {'XYNBS', 'XYZNB'}
                        FittingInput.NParams = 5;
                    case 'XYNBSXSY'
                        FittingInput.NParams = 6;
                end
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
                elseif isnumeric(ThresholdingInput.On)
                    ThresholdingInput.On = logical(ThresholdingInput.On);
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
                if ~(isnumeric(ThresholdingInput.MinPValue) ...
                        && (ThresholdingInput.MinPValue<=1) ...
                        && (ThresholdingInput.MinPValue>=0))
                    error(['''SMF.Thresholding.MinPValue'' ', ...
                        'must be a number in the interval [0, 1].'])
                end
            end
            if isfield(ThresholdingInput, 'AutoThreshLogL')
                if ~(islogical(ThresholdingInput.AutoThreshLogL) ...
                        || isnumeric(ThresholdingInput.AutoThreshLogL))
                    error(['''SMF.Thresholding.AutoThreshLogL'' must be logical ', ...
                        'or interpretable as logical (numeric).'])
                elseif isnumeric(ThresholdingInput.AutoThreshLogL)
                    ThresholdingInput.AutoThreshLogL = ...
                        logical(ThresholdingInput.AutoThreshLogL);
                end
            end
            if isfield(ThresholdingInput, 'AutoThreshPrctile')
                if ~(isnumeric(ThresholdingInput.AutoThreshPrctile) ...
                        && (ThresholdingInput.AutoThreshPrctile<=100) ...
                        && (ThresholdingInput.AutoThreshPrctile>=0))
                    error(['''SMF.Thresholding.AutoThreshPrctile'' ', ...
                        'must be a number in the interval [0, 100].'])
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
            if isfield(ThresholdingInput, 'InMeanMultiplier')
                if ~isnumeric(ThresholdingInput.InMeanMultiplier)
                    error('''SMF.Thresholding.InMeanMultiplier'' must be numeric.')
                end
            end
            if isfield(ThresholdingInput, 'NNMedianMultiplier')
                if ~isnumeric(ThresholdingInput.NNMedianMultiplier)
                    error('''SMF.Thresholding.NNMedianMultiplier'' must be numeric.')
                end
            end
            if isfield(ThresholdingInput, 'MinNumNeighbors')
                if (mod(ThresholdingInput.MinNumNeighbors, 1) || (ThresholdingInput.MinNumNeighbors<0))
                    error(['''SMF.Thresholding.MinNumNeighbors'' ', ...
                        'must be a non-negative integer.'])
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
                elseif isnumeric(FCInput.On)
                    FCInput.On = logical(FCInput.On);
                end
            end
            if isfield(FCInput, 'Method')
                if ~ismember(lower(FCInput.Method), ...
                        {'hypothesis test', 'lap-fc', ...
                        'classical', 'revised classical'})
                    error(['''SMF.FrameConnection.Method'' must be ', ...
                        '''Hypothesis test'', ''LAP-FC'', ', ...
                        '''Classical'', or ''Revised classical''.'])
                end
            end
            if isfield(FCInput, 'MaxSeparation')
                if ~isnumeric(FCInput.MaxSeparation)
                    error(['''SMF.FrameConnection.MaxSeparation'' ', ...
                        'must be numeric.'])
                end
            end
            if isfield(FCInput, 'LoS')
                if ~(isnumeric(FCInput.LoS) ...
                        && (FCInput.LoS<=1) ...
                        && (FCInput.LoS>=0))
                    error(['''SMF.FrameConnection.LoS'' ', ...
                        'must be a number in the interval [0, 1].'])
                end
            end
            if isfield(FCInput, 'MaxFrameGap')
                if mod(FCInput.MaxFrameGap, 1)
                    error(['''SMF.FrameConnection.MaxFrameGap'' ', ...
                        'must be an integer.'])
                end
            end
            if isfield(FCInput, 'NNearestClusters')
                if (mod(FCInput.NNearestClusters, 1) ...
                        || (FCInput.NNearestClusters<1))
                    error(['''SMF.FrameConnection.NNearestClusters'' ', ...
                        'must be a positive integer.'])
                end
            end
            if isfield(FCInput, 'NIterations')
                if (mod(FCInput.NIterations, 1) || (FCInput.NIterations<1))
                    error(['''SMF.FrameConnection.NIterations'' ', ...
                        'must be a positive integer.'])
                end
            end
            if isfield(FCInput, 'MinNFrameConns')
                if (mod(FCInput.MinNFrameConns, 1) || (FCInput.MinNFrameConns<1))
                    error(['''SMF.FrameConnection.MinNFrameConns'' ', ...
                        'must be a positive integer.'])
                end
            end
            if isfield(FCInput, 'NSigmaDev')
                if ~(isnumeric(FCInput.NSigmaDev) ...
                        && (FCInput.NSigmaDev>=0))
                    error(['''SMF.FrameConnection.NSigmaDev'' ', ...
                        'must be a non-negative number.'])
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
                elseif isnumeric(DCInput.On)
                    DCInput.On = logical(DCInput.On);
                end
            end
            if isfield(DCInput, 'Method')
                if ~ismember(lower(DCInput.Method), ...
                        {'dc-knn', 'dc-bf'})
                    error(['''SMF.DriftCorrection.Method'' must be ', ...
                        '''DC-KNN'' or ''DC-BF''.'])
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
                    && ~ismember(TrackingInput.Method, {'CostMatrix'}))
                error('SMF.Tracking.Method must be one of CostMatrix,')
            end
            if isfield(TrackingInput, 'D')
                if ~isnumeric(TrackingInput.D)
                    error('''SMF.Tracking.D'' must be numeric.')
                end
            end
            if isfield(TrackingInput, 'TrajwiseD')
                if ~(islogical(TrackingInput.TrajwiseD) ...
                        || isnumeric(TrackingInput.TrajwiseD))
                    error(['''SMF.Tracking.TrajwiseD'' must be ', ...
                        'logical or interpretable as logical (numeric).'])
                elseif isnumeric(TrackingInput.TrajwiseD)
                    TrackingInput.TrajwiseD = ...
                        logical(TrackingInput.TrajwiseD);
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
            if isfield(TrackingInput, 'MaxDistFF')
                if ~isnumeric(TrackingInput.MaxDistFF)
                    error('''SMF.Tracking.MaxDistFF'' must be numeric.')
                end
            end
            if isfield(TrackingInput, 'MaxDistGC')
                if ~isnumeric(TrackingInput.MaxDistGC)
                    error('''SMF.Tracking.MaxDistGC'' must be numeric.')
                end
            end
            if isfield(TrackingInput, 'MaxZScoreDist')
                if ~isnumeric(TrackingInput.MaxZScoreDist)
                    error('''SMF.Tracking.MaxZScoreDist'' must be numeric.')
                end
            end
            if isfield(TrackingInput, 'MaxZScorePhotons')
                if ~isnumeric(TrackingInput.MaxZScorePhotons)
                    error(['''SMF.Tracking.MaxZScorePhotons'' ', ...
                        'must be numeric.'])
                end
            end
            if isfield(TrackingInput, 'MaxFrameGap')
                if mod(TrackingInput.MaxFrameGap, 1)
                    error(['''SMF.Tracking.MaxFrameGap'' ', ...
                        'must be an integer.'])
                end
            end
            if isfield(TrackingInput, 'MinTrackLength')
                if mod(TrackingInput.MinTrackLength, 1)
                    error(['''SMF.Tracking.MinTrackLength'' ', ...
                        'must be an integer.'])
                end
            end
            if isfield(TrackingInput, 'NIterMax')
                if mod(TrackingInput.NIterMax, 1)
                    error(['''SMF.Tracking.NIterMax'' ', ...
                        'must be an integer.'])
                end
            end
            if isfield(TrackingInput, 'NIterMaxBatch')
                if mod(TrackingInput.NIterMaxBatch, 1)
                    error(['''SMF.Tracking.NIterMaxBatch'' ', ...
                        'must be an integer.'])
                end
            end
            if isfield(TrackingInput, 'ParamsHistory')
                if ~iscell(TrackingInput.ParamsHistory)
                    error(['''SMF.Tracking.ParamsHistory'' ', ...
                        'must be a cell array.'])
                end
            end
            if isfield(TrackingInput, 'MaxRelativeChange')
                if ~isnumeric(TrackingInput.MaxRelativeChange)
                    error(['''SMF.Tracking.MaxRelativeChange'' ', ...
                        'must be numeric'])
                end
            end
            if isfield(TrackingInput, 'TryLowPValueLocs')
                if ~(islogical(TrackingInput.TryLowPValueLocs) ...
                        || isnumeric(TrackingInput.TryLowPValueLocs))
                    error(['''SMF.Tracking.TryLowPValueLocs'' must be ', ...
                        'logical or interpretable as logical (numeric).'])
                elseif isnumeric(TrackingInput.TryLowPValueLocs)
                    TrackingInput.TryLowPValueLocs = ...
                        logical(TrackingInput.TryLowPValueLocs);
                end
            end
            
            % Set the input fields as class properties.
            InputFields = fieldnames(TrackingInput);
            for ff = 1:numel(InputFields)
                obj.Tracking.(InputFields{ff}) = ...
                    TrackingInput.(InputFields{ff});
            end
        end
        
        
        function [SMFStruct] = packageSMF(obj)
            %packageSMF converts class instance obj to a structure array.
            % This method will convert the class instance into a structure
            % array by setting all class properties as fields in the
            % structure array.
            
            % Generate the output SMFStruct.
            ClassFields = fieldnames(obj);
            for ff = 1:numel(ClassFields)
                SMFStruct.(ClassFields{ff}) = obj.(ClassFields{ff});
            end
            SMFStruct.SMFVersion = obj.SMFVersion;
        end
        
        function importSMF(obj, SMFStruct)
            %importSMF imports fields from SMF into the class instance obj.
            % This method will take the fields given in the structure array
            % SMFStruct and set them to the corresponding class properties
            % in obj.
            
            % Update class properties based on fields in SMFStruct.
            InputFields = fieldnames(SMFStruct);
            ClassFields = fieldnames(obj);
            ValidFields = InputFields(ismember(InputFields, ClassFields));
            for ff = 1:numel(ValidFields)
                obj.(ValidFields{ff}) = SMFStruct.(ValidFields{ff});
            end
        end
        
        function resetSMF(obj)
            %resetSMF resets all SMF fields to their default values.
            % This method will overwrite all properties of the class
            % instance obj with their default values defined by the class
            % constructor. This is done by creating a new instance of the
            % SMF class and then importing those properties to the current
            % instance obj using importSMF().
            SMFDefault = smi_core.SingleMoleculeFitting;
            obj.importSMF(SMFDefault);
        end
        
        [GUIParent] = gui(obj, GUIParent);
        
    end
    
    methods (Static)
        function [SMF] = reloadSMF(SMFStruct)
            %reloadSMF loads a struct type SMF into a class instance SMF
            % This method will take the fields given in the structure array
            % SMFStruct and set them to the corresponding class properties,
            % i.e., it creates an SMF class instance from an input
            % SMFStruct.
            % NOTE: This is basically a static wrapper for the non-static
            %       method importSMF(), which needed to be non-static due
            %       to its intended usage in the GUI.
            
            % Create an instance of the SingleMoleculeFitting class.
            SMF = smi_core.SingleMoleculeFitting;
            
            % Update class properties based on fields in SMFStruct.
            SMF.importSMF(SMFStruct)
        end
        
        [SMFPadded, PaddedFields] = padSMF(SMF, SMFPadding, ...
            DisplayMessages);
    end
    
end
