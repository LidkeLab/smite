classdef DiffusionEstimator < handle
    %DiffusionEstimator contains methods useful for diffusion estimation.
    % This class contains several methods used to estimate diffusion
    % constants from single-particle tracking data.
    %
    % REQUIRES:
    %   Optimization Toolbox (to use fmincon() for certain models)
    
    % Created by:
    %   David J. Schodt (Lidke lab, 2021)
    
    
    properties
        % ID of the diffusion model to be considered (char array/string)
        % OPTIONS:
        %   'Brownian'
        DiffusionModel{mustBeMember(DiffusionModel, {'Brownian'})} = ...
            'Brownian';
        
        % Number of diffusing components. (Default = 1)
        % This can only be changed for FitTarget = CDFOfJumps or
        % LikelihoodOfJumps.
        NComponents = 1;
        
        % Fit method for fitting data (char array/string)
        FitMethod{mustBeMember(FitMethod, {'WeightedLS', 'LS'})} = ...
            'WeightedLS';
        
        % Target data that will be fit (char array/string)
        FitTarget{mustBeMember(FitTarget, ...
            {'MSD', 'CDFOfJumps', 'LikelihoodOfJumps'})} = 'MSD';

        % Range of frame lags used to estimate D (Default = [1, 5])
        FrameLagRange = [1, 5];
        
        % Number of MSD points to be fit (scalar, integer)(Default = 5)
        NFitPoints = 5;
        
        % Flag to estimate standard errors (Default = true)
        % NOTE: For some of the fit methods/fit targets (e.g., FitTarget =
        %       'CDFOfJumps') we have to use a bootstrap, in which case
        %       estimating SEs can be unreasonably slow.
        EstimateSEs = true;
        
        % Flag to fit individual trajectory MSDs/CDFs (Default = true)
        FitIndividualTrajectories = true;
        
        % Number of spatial dimensions (scalar, integer)(Default = 2)
        NDimensions = 2;
        
        % Directory in which results will be saved by saveResults().
        SaveDir = pwd();
        
        % Base name of saved results. Default defined in obj.saveDir().
        BaseSaveName
        
        % Tracking results structure.
        TR
        
        % Single Molecule Fitting structure, for pixel size and framerate.
        SMF = smi_core.SingleMoleculeFitting;
        
        % Boolean flag to indicate units of outputs (boolean)(Default = 0)
        % 1 (true) will make the outputs of estimateDiffusionConstant()
        %   micrometers and seconds.
        % 0 (false) will make the outputs of estimateDiffusionConstant()
        %   pixels and frames.
        % NOTE: Most methods of this class will use pixels and frames
        %       regardless of obj.UnitFlag. This property will only affect
        %       the "user-facing" wrapper methods, such as
        %       estimateDiffusionConstant() and saveResults()
        UnitFlag = false;
        
        % Verbosity level for estimateDiffusionConstant() (Default = 0)
        %   Verbose 0: no Command Window outputs
        %   Verbose 1: General progress printed to Command Window
        %   Verbose 2: Display some intermediate results
        %   Verbose 3: Debugging mode, extensive outputs
        Verbose = 0;
    end
    
    properties (SetAccess = protected)
        % Structure array containing diffusion estimates.
        DiffusionStruct
        
        % Structure array containing trajectory-wise MSDs.
        MSDSingleTraj
        
        % Structure array containing the ensemble MSD.
        MSDEnsemble
    end
    
    methods
        function [obj, DiffusionStruct] = ...
                DiffusionEstimator(TR, SMF, Verbose, AutoRun)
            %DiffusionEstimator is the class constructor.
            % Several optional inputs can be provided to directly set class
            % properties.  'AutoRun' is a boolean flag which specifies
            % whether or not this constructor should call
            % obj.estimateDiffusionConstant() if all requisite class
            % properties were provided.
            
            % Set defaults if needed.
            if (~exist('AutoRun', 'var') || isempty(AutoRun))
                AutoRun = false;
            end
            
            % Set class properties based on the inputs.
            AllFieldsSet = true;
            if (exist('TR', 'var') && ~isempty(TR))
                obj.TR = TR;
            else
                AllFieldsSet = false;
            end
            if (exist('SMF', 'var') && ~isempty(SMF))
                obj.SMF = SMF;
            end
            if (exist('Verbose', 'var') && ~isempty(Verbose))
                obj.Verbose = Verbose;
            end
            
            % Run obj.estimateDiffusionConstant() if requested.
            if (AutoRun && AllFieldsSet)
                [DiffusionStruct] = obj.estimateDiffusionConstant();
            end
        end
        
        [DiffusionStruct] = estimateDiffusionConstant(obj);
        saveResults(obj, SaveParams)
        
    end
    
    methods (Static)
        [MSDSingleTraj, MSDEnsemble] = ...
            computeMSD(TR, FrameLagRange, Verbose);
        [MSDStruct] = computeCDFOfJumps(MSDStruct, FrameLagRange);
        [FitParams, FitParamsSE] = ...
            fitMSD(MSDStruct, FitMethod, NFitPoints, ...
            DiffusionModel, Verbose);
        [FitParams, FitParamsSE] = fitCDFOfJumps(MSDStruct, ...
            FitMethod, NComponents, DiffusionModel, Verbose);
        [MLEParams, MLEParamsSE] = ...
            mleOfJumps(MSDStruct, NComponents, DiffusionModel, Verbose);
        [PlotAxes] = plotEnsembleMSD(PlotAxes, ...
            MSDEnsemble, DiffusionStruct, DiffusionModel, UnitFlag);
        [PlotAxes] = plotEnsembleCDFOfJumps(PlotAxes, ...
            MSDEnsemble, DiffusionStruct, UnitFlag);
    end
    
    methods (Static, Hidden)
        % These methods are 'Hidden' because they are primarily intended
        % for use within other more user-centric methods (i.e., we don't
        % want to distract the user with these options, but if they need
        % them they are still accessible).
        
        [MSDSingleTraj] = computeSingleTrajMSD(TR, FrameLagRange, Verbose);
        [FitParams, FitParamsSE] = ...
            fitMSDBrownian(FrameLags, MSD, NPoints, FitMethod);
        [FitParams, FitParamsSE] = ...
            fitCDFOfJumpsBrownian(SortedJumps, CDFOfJumps, ...
            FrameLags, NPoints, LocVarianceSum, NComponents, ...
            Weights, FitMethod, FitOptions);
        [CDFOfJumps] = brownianJumpCDF(MotionParams, ...
            SortedSquaredDisp, FrameLags, NPoints, LocVarianceSum);
        [MLEParams, MLEParamsSE] = mleOfJumpsBrownian(...
            SquaredDisplacement, FrameLagsAll, ...
            LocVarianceSum, NComponents, FitOptions);
        [LogLikelihood] = brownianJumpLikelihood(MotionParams, ...
            SquaredDisplacement, FrameLagsAll, LocVarianceSum);
    end
    
    
end