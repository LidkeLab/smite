function [DiffusionStruct] = estimateDiffusionConstant(obj)
%estimateDiffusionConstant estimates the diffusion constant from an MSD.
% This method will fit the mean squared displacement data in 'MSD' to make
% an estimate of the diffusion constant.
%
% OUTPUTS:
%   DiffusionStruct: Structure array containing the fit diffusion
%                    constants and their standard errors. The units will be
%                    specified by the property obj.UnitFlag, with 1
%                    specifying physical units (micrometers, seconds) and 0
%                    specifying camera units (pixels, frames).
%
% REQUIRES:
%   Statistics and Machine Learning Toolbox (for ecdf())

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Compute the MSDs.
if (obj.Verbose > 0)
    fprintf('estimateDiffusionConstant(): computing MSDs...\n');
end
[obj.MSDSingleTraj, obj.MSDEnsemble] = ...
    obj.computeMSD(obj.TR, obj.MaxFrameLag, obj.Verbose);

% Fit the results and convert units if necessary.
if (obj.Verbose > 1)
    fprintf(['estimateDiffusionConstant(): fitting trajectory-wise ', ...
        'data...\n']);
elseif (obj.Verbose > 0)
    fprintf(['estimateDiffusionConstant(): estimating diffusion ', ...
        'constants...\n']);
end
DUnitConversion = ~obj.UnitFlag ...
    + obj.UnitFlag*(obj.TR(1).PixelSize^2)*obj.TR(1).FrameRate;
FitParamsSingleTraj = NaN;
FitParamsSingleTrajSE = NaN;
DiffusionConstantSingleTraj = NaN;
DiffusionConstantSingleTrajSE = NaN;
switch obj.FitTarget
    case 'MSD'
        % Fit the trajectory-wise MSDs (if requested).
        if obj.FitIndividualTrajectories
            [FitParamsSingleTraj, FitParamsSingleTrajSE] = ...
                obj.fitMSD(obj.MSDSingleTraj, ...
                obj.FitMethod, obj.DiffusionModel, obj.Verbose);
            
            % Compute the diffusion constants.
            DiffusionConstantSingleTraj = DUnitConversion ...
                * FitParamsSingleTraj(:, 2) / (2*obj.NDimensions);
            DiffusionConstantSingleTrajSE = DUnitConversion ...
                * FitParamsSingleTrajSE(:, 2) / (2*obj.NDimensions);
        end
        
        % Fit the ensemble MSD.
        if (obj.Verbose > 1)
            fprintf(['estimateDiffusionConstant(): fitting ensemble ', ...
                'MSD...\n']);
        end
        [FitParamsEnsemble, FitParamsEnsembleSE] = ...
            obj.fitMSD(obj.MSDEnsemble, ...
            obj.FitMethod, obj.DiffusionModel, obj.Verbose);
        
        % Compute the ensemble diffusion constant.
        DiffusionConstantEnsemble = DUnitConversion ...
            * FitParamsEnsemble(:, 2) / (2*obj.NDimensions);
        DiffusionConstantEnsembleSE = DUnitConversion ...
            * FitParamsEnsembleSE(:, 2) / (2*obj.NDimensions);
    case 'CDFOfJumps'
        % Compute the CDF (cumulative distribution function, a.k.a. 
        % cumulative probability distribution, CPD) of the trajectory-wise
        % displacements.
        obj.MSDSingleTraj = obj.computeCDFOfJumps(obj.MSDSingleTraj);
        obj.MSDEnsemble = obj.computeCDFOfJumps(obj.MSDEnsemble);
        
        % Fit the trajectory-wise CDFs (if requested).
        if obj.FitIndividualTrajectories
            [FitParamsSingleTraj, FitParamsSingleTrajSE] = ...
                obj.fitCDFOfJumps(obj.MSDSingleTraj, ...
                obj.FitMethod, obj.DiffusionModel, obj.Verbose);
            
            % Compute the diffusion constants.
            DiffusionConstantSingleTraj = DUnitConversion ...
                * FitParamsSingleTraj / (2*obj.NDimensions);
            DiffusionConstantSingleTrajSE = DUnitConversion ...
                * FitParamsSingleTrajSE / (2*obj.NDimensions);
        end
        
        % Fit the ensemble CDF of jumps.
        if (obj.Verbose > 1)
            fprintf(['estimateDiffusionConstant(): fitting ensemble ', ...
                'CDF of jumps...\n']);
        end
        [FitParamsEnsemble, FitParamsEnsembleSE] = ...
            obj.fitCDFOfJumps(obj.MSDEnsemble, ...
            obj.FitMethod, obj.DiffusionModel, obj.Verbose);
        
        % Compute the ensemble diffusion constant.
        DiffusionConstantEnsemble = DUnitConversion ...
            * FitParamsEnsemble / (2*obj.NDimensions);
        DiffusionConstantEnsembleSE = DUnitConversion ...
            * FitParamsEnsembleSE / (2*obj.NDimensions);
end

% Store the results in the DiffusionStruct.
DiffusionStruct(1).Name = 'trajectory';
DiffusionStruct(1).Units = ...
    smi_helpers.stringMUX({'pixels, frames', 'micrometers, seconds'}, ...
    obj.UnitFlag);
DiffusionStruct(1).FitParams = FitParamsSingleTraj;
DiffusionStruct(1).FitParamsSE = FitParamsSingleTrajSE;
DiffusionStruct(1).PixelSize = obj.TR(1).PixelSize;
DiffusionStruct(1).FrameRate = obj.TR(1).FrameRate;
DiffusionStruct(1).DiffusionConstant = DiffusionConstantSingleTraj;
DiffusionStruct(1).DiffusionConstantSE = DiffusionConstantSingleTrajSE;
DiffusionStruct(2).Name = 'ensemble';
DiffusionStruct(2).Units = ...
    smi_helpers.stringMUX({'pixels, frames', 'micrometers, seconds'}, ...
    obj.UnitFlag);
DiffusionStruct(2).FitParams = FitParamsEnsemble;
DiffusionStruct(2).FitParamsSE = FitParamsEnsembleSE;
DiffusionStruct(2).PixelSize = obj.TR(1).PixelSize;
DiffusionStruct(2).FrameRate = obj.TR(1).FrameRate;
DiffusionStruct(2).DiffusionConstant = DiffusionConstantEnsemble;
DiffusionStruct(2).DiffusionConstantSE = DiffusionConstantEnsembleSE;
obj.DiffusionStruct = DiffusionStruct;


end