function [DiffusionStruct] = estimateDiffusionConstant(obj)
%estimateDiffusionConstant estimates the diffusion constant.
% This method will fit the MSD/CDF/etc. data to estimate the diffusion
% constant(s). This is meant to be the main run method of the
% DiffusionEstimator class.
%
% OUTPUTS:
%   DiffusionStruct: Structure array containing the fit diffusion
%                    constants and their standard errors. The units will be
%                    specified by the property obj.UnitFlag, with 1
%                    specifying physical units (micrometers, seconds) and 0
%                    specifying camera units (pixels, frames).

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Compute the MSDs.
if (obj.Verbose > 0)
    fprintf('estimateDiffusionConstant(): computing MSDs...\n');
end
[obj.MSDSingleTraj, obj.MSDEnsemble] = ...
    obj.computeMSD(obj.TR, obj.MaxFrameLag, obj.Verbose);

% Determine how many diffusion components are used in the desired model.
switch lower(obj.DiffusionModel)
    case 'brownian1c'
        NComponents = 1;
    case 'brownian2c'
        NComponents = 2;
    otherwise
        error('Unknown ''DiffusionModel'' = %s', DiffusionModel)
end

% Fit the results and convert units if necessary.
if (obj.Verbose > 1)
    fprintf(['estimateDiffusionConstant(): fitting trajectory-wise ', ...
        'data...\n']);
elseif (obj.Verbose > 0)
    fprintf(['estimateDiffusionConstant(): estimating diffusion ', ...
        'constants...\n']);
end
JumpUnitConversion = ~obj.UnitFlag + obj.UnitFlag*obj.TR(1).PixelSize;
FrameUnitConversion = ~obj.UnitFlag + obj.UnitFlag*obj.TR(1).FrameRate;
FitParamsSingleTraj = NaN;
FitParamsSingleTrajSE = NaN;
DiffusionConstantSingleTraj = NaN;
DiffusionConstantSingleTrajSE = NaN;
PopulationRatiosSingleTraj = NaN;
PopulationRatiosSingleTrajSE = NaN;
switch obj.FitTarget
    case 'MSD'
        if obj.FitIndividualTrajectories
            % Fit the trajectory-wise MSDs.
            % NOTE: fitMSD() is written such that the SEs are only computed
            %       when the second output is requested.
            if obj.EstimateSEs
                [FitParamsSingleTraj, FitParamsSingleTrajSE] = ...
                    obj.fitMSD(obj.MSDSingleTraj, ...
                    obj.FitMethod, obj.NFitPoints, ...
                    obj.DiffusionModel, obj.Verbose);
            else
                FitParamsSingleTraj = obj.fitMSD(obj.MSDSingleTraj, ...
                    obj.FitMethod, obj.NFitPoints, ...
                    obj.DiffusionModel, obj.Verbose);
                FitParamsSingleTrajSE = NaN(size(FitParamsSingleTraj));
            end
            FitParamsSingleTraj = FitParamsSingleTraj ...
                .* (JumpUnitConversion.^2) .* [1, FrameUnitConversion];
            FitParamsSingleTrajSE = FitParamsSingleTrajSE ...
                .* (JumpUnitConversion.^2) .* [1, FrameUnitConversion];
            
            % Compute the diffusion constants.
            DiffusionConstantSingleTraj = ...
                FitParamsSingleTraj(:, 2) / (2*obj.NDimensions);
            DiffusionConstantSingleTrajSE = ...
                FitParamsSingleTrajSE(:, 2) / (2*obj.NDimensions);
            PopulationRatiosSingleTraj = 1;
            PopulationRatiosSingleTrajSE = 0;
        end
        
        % Fit the ensemble MSD.
        if (obj.Verbose > 1)
            fprintf(['estimateDiffusionConstant(): fitting ensemble ', ...
                'MSD...\n']);
        end
        if obj.EstimateSEs
            [FitParamsEnsemble, FitParamsEnsembleSE] = ...
                obj.fitMSD(obj.MSDEnsemble, ...
                obj.FitMethod, obj.NFitPoints, ...
                obj.DiffusionModel, obj.Verbose);
        else
            FitParamsEnsemble = obj.fitMSD(obj.MSDEnsemble, ...
                obj.FitMethod, obj.NFitPoints, ...
                obj.DiffusionModel, obj.Verbose);
            FitParamsEnsembleSE = NaN(size(FitParamsEnsemble));
        end
        FitParamsEnsemble = FitParamsEnsemble ...
            .* (JumpUnitConversion.^2) .* [1, FrameUnitConversion];
        FitParamsEnsembleSE = FitParamsEnsembleSE ...
            .* (JumpUnitConversion.^2) .* [1, FrameUnitConversion];
        
        % Compute the ensemble diffusion constant.
        DiffusionConstantEnsemble = ...
            FitParamsEnsemble(:, 2) / (2*obj.NDimensions);
        DiffusionConstantEnsembleSE = ...
            FitParamsEnsembleSE(:, 2) / (2*obj.NDimensions);
        PopulationRatiosEnsemble = 1;
        PopulationRatiosEnsembleSE = 0;
    case 'CDFOfJumps'
        % Compute the CDF (cumulative distribution function, a.k.a.
        % cumulative probability distribution, CPD) of the trajectory-wise
        % displacements.
        obj.MSDSingleTraj = obj.computeCDFOfJumps(obj.MSDSingleTraj);
        obj.MSDEnsemble = obj.computeCDFOfJumps(obj.MSDEnsemble);
        
        % Fit the trajectory-wise CDFs.
        if obj.FitIndividualTrajectories
            if obj.EstimateSEs
                [FitParamsSingleTraj, FitParamsSingleTrajSE] = ...
                    obj.fitCDFOfJumps(obj.MSDSingleTraj, ...
                    obj.FitMethod, obj.NFitPoints, ...
                    obj.DiffusionModel, obj.Verbose);
            else
                FitParamsSingleTraj = ...
                    obj.fitCDFOfJumps(obj.MSDSingleTraj, ...
                    obj.FitMethod, obj.NFitPoints, ...
                    obj.DiffusionModel, obj.Verbose);
                FitParamsSingleTrajSE = NaN(size(FitParamsSingleTraj));
            end
            FitParamsSingleTraj(1:NComponents) = ...
                FitParamsSingleTraj(1:NComponents) ...
                .* (JumpUnitConversion.^2) .* FrameUnitConversion;
            FitParamsSingleTrajSE(1:NComponents) = ...
                FitParamsSingleTrajSE(1:NComponents) ...
                .* (JumpUnitConversion.^2) .* FrameUnitConversion;
            
            % Compute the diffusion constants.
            DiffusionConstantSingleTraj = ...
                FitParamsSingleTraj(1:NComponents) / 2;
            DiffusionConstantSingleTrajSE(1:NComponents) = ...
                FitParamsSingleTrajSE / 2;
            
            % Extract the population ratios.
            PopulationRatiosSingleTraj = ...
                [FitParamsSingleTraj((NComponents+1):end); ...
                1 - sum(FitParamsSingleTraj((NComponents+1):end))];
            PopulationRatiosSingleTrajSE = ...
                [FitParamsSingleTrajSE((NComponents+1):end); ...
                NaN];
        end
        
        % Fit the ensemble CDF of jumps.
        if (obj.Verbose > 1)
            fprintf(['estimateDiffusionConstant(): fitting ensemble ', ...
                'CDF of jumps...\n']);
        end
        if obj.EstimateSEs
            [FitParamsEnsemble, FitParamsEnsembleSE] = ...
                obj.fitCDFOfJumps(obj.MSDEnsemble, ...
                obj.FitMethod, obj.DiffusionModel, obj.Verbose);
        else
            FitParamsEnsemble = ...
                obj.fitCDFOfJumps(obj.MSDEnsemble, ...
                obj.FitMethod, obj.DiffusionModel, obj.Verbose);
            FitParamsEnsembleSE = NaN(size(FitParamsEnsemble));
        end
        FitParamsEnsemble(1:NComponents) = ...
            FitParamsEnsemble(1:NComponents) ...
            .* (JumpUnitConversion.^2) .* FrameUnitConversion;
        FitParamsEnsembleSE(1:NComponents) = ...
            FitParamsEnsembleSE(1:NComponents) ...
            .* (JumpUnitConversion.^2) .* FrameUnitConversion;
        
        % Compute the ensemble diffusion constant.
        DiffusionConstantEnsemble = ...
            FitParamsEnsemble(1:NComponents) / 2;
        DiffusionConstantEnsembleSE = ...
            FitParamsEnsembleSE(1:NComponents) / 2;
        
        % Extract the population ratios.
        PopulationRatiosEnsemble = ...
            [FitParamsEnsemble((NComponents+1):end); ...
            1 - sum(FitParamsEnsemble((NComponents+1):end))];
        PopulationRatiosEnsembleSE = ...
            [FitParamsEnsembleSE((NComponents+1):end); ...
            NaN];
end

% Store the results in the DiffusionStruct.
DiffusionStruct(1).Name = 'trajectory';
JumpUnit = smi_helpers.stringMUX({'pixels', 'micrometers'}, obj.UnitFlag);
TimeUnit = smi_helpers.stringMUX({'frames', 'seconds'}, obj.UnitFlag);
DiffusionStruct(1).Units = {JumpUnit; TimeUnit};
DiffusionStruct(1).FitParams = FitParamsSingleTraj;
DiffusionStruct(1).FitParamsSE = FitParamsSingleTrajSE;
DiffusionStruct(1).PixelSize = obj.TR(1).PixelSize;
DiffusionStruct(1).FrameRate = obj.TR(1).FrameRate;
DiffusionStruct(1).DiffusionConstant = DiffusionConstantSingleTraj;
DiffusionStruct(1).DiffusionConstantSE = DiffusionConstantSingleTrajSE;
DiffusionStruct(1).PopulationRatios = PopulationRatiosSingleTraj;
DiffusionStruct(1).PopulationRatiosSE = PopulationRatiosSingleTrajSE;
DiffusionStruct(1).NFitPoints = obj.NFitPoints;
DiffusionStruct(2).Name = 'ensemble';
DiffusionStruct(2).Units = {JumpUnit; TimeUnit};
DiffusionStruct(2).FitParams = FitParamsEnsemble;
DiffusionStruct(2).FitParamsSE = FitParamsEnsembleSE;
DiffusionStruct(2).PixelSize = obj.TR(1).PixelSize;
DiffusionStruct(2).FrameRate = obj.TR(1).FrameRate;
DiffusionStruct(2).DiffusionConstant = DiffusionConstantEnsemble;
DiffusionStruct(2).DiffusionConstantSE = DiffusionConstantEnsembleSE;
DiffusionStruct(2).PopulationRatios = PopulationRatiosEnsemble;
DiffusionStruct(2).PopulationRatiosSE = PopulationRatiosEnsembleSE;
DiffusionStruct(2).NFitPoints = obj.NFitPoints;
obj.DiffusionStruct = DiffusionStruct;


end