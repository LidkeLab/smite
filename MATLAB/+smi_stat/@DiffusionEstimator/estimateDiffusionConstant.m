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
    obj.computeMSD(obj.TR, obj.FrameLagRange, obj.Verbose);

% Fit the results and convert units where necessary.
if (obj.Verbose > 1)
    fprintf(['estimateDiffusionConstant(): fitting trajectory-wise ', ...
        'data...\n']);
elseif (obj.Verbose > 0)
    fprintf(['estimateDiffusionConstant(): estimating diffusion ', ...
        'constants...\n']);
end
JumpUnitConversion = ~obj.UnitFlag + obj.UnitFlag*obj.SMF.Data.PixelSize;
FrameUnitConversion = ~obj.UnitFlag + obj.UnitFlag*obj.SMF.Data.FrameRate;
ParamsSingleTraj = NaN;
ParamsSingleTrajSE = NaN;
DiffusionConstantSingleTraj = NaN;
DiffusionConstantSingleTrajSE = NaN;
PopulationRatiosSingleTraj = NaN;
PopulationRatiosSingleTrajSE = NaN;
switch obj.FitTarget
    case 'MSD'
        % Update NComponents to make sure it's 1 (we can't fit multiple
        % components with MSD).
        obj.NComponents = 1;
        
        % Fit the trajectory-wise MSDs.
        % NOTE: fitMSD() is written such that the SEs are only computed
        %       when the second output is requested.
        if obj.FitIndividualTrajectories
            if obj.EstimateSEs
                [ParamsSingleTraj, ParamsSingleTrajSE] = ...
                    obj.fitMSD(obj.MSDSingleTraj, ...
                    obj.FitMethod, obj.NFitPoints, ...
                    obj.DiffusionModel, obj.Verbose);
            else
                ParamsSingleTraj = obj.fitMSD(obj.MSDSingleTraj, ...
                    obj.FitMethod, obj.NFitPoints, ...
                    obj.DiffusionModel, obj.Verbose);
                ParamsSingleTrajSE = NaN(size(ParamsSingleTraj));
            end
            ParamsSingleTraj = ParamsSingleTraj ...
                .* (JumpUnitConversion.^2) .* [1, FrameUnitConversion];
            ParamsSingleTrajSE = ParamsSingleTrajSE ...
                .* (JumpUnitConversion.^2) .* [1, FrameUnitConversion];
            
            % Compute the diffusion constants.
            DiffusionConstantSingleTraj = ...
                ParamsSingleTraj(:, 2) / (2*obj.NDimensions);
            DiffusionConstantSingleTrajSE = ...
                ParamsSingleTrajSE(:, 2) / (2*obj.NDimensions);
            PopulationRatiosSingleTraj = 1;
            PopulationRatiosSingleTrajSE = 0;
        end
        
        % Fit the ensemble MSD.
        if (obj.Verbose > 1)
            fprintf(['estimateDiffusionConstant(): fitting ensemble ', ...
                'MSD...\n']);
        end
        if obj.EstimateSEs
            [ParamsEnsemble, ParamsEnsembleSE] = ...
                obj.fitMSD(obj.MSDEnsemble, ...
                obj.FitMethod, obj.NFitPoints, ...
                obj.DiffusionModel, obj.Verbose);
        else
            ParamsEnsemble = obj.fitMSD(obj.MSDEnsemble, ...
                obj.FitMethod, obj.NFitPoints, ...
                obj.DiffusionModel, obj.Verbose);
            ParamsEnsembleSE = NaN(size(ParamsEnsemble));
        end
        ParamsEnsemble = ParamsEnsemble ...
            .* (JumpUnitConversion.^2) .* [1, FrameUnitConversion];
        ParamsEnsembleSE = ParamsEnsembleSE ...
            .* (JumpUnitConversion.^2) .* [1, FrameUnitConversion];
        
        % Compute the ensemble diffusion constant.
        DiffusionConstantEnsemble = ...
            ParamsEnsemble(:, 2) / (2*obj.NDimensions);
        DiffusionConstantEnsembleSE = ...
            ParamsEnsembleSE(:, 2) / (2*obj.NDimensions);
        PopulationRatiosEnsemble = 1;
        PopulationRatiosEnsembleSE = 0;
    case 'CDFOfJumps'
        % Modify obj.NFitPoints to reflect the fact that we're fitting the
        % entire CDF, no matter what the user has requested!
        obj.NFitPoints = numel(obj.MSDEnsemble.FrameLagsAll);
        
        % Compute the CDF (cumulative distribution function, a.k.a.
        % cumulative probability distribution, CPD) of the trajectory-wise
        % displacements.
        obj.MSDSingleTraj = ...
            obj.computeCDFOfJumps(obj.MSDSingleTraj, obj.FrameLagRange);
        obj.MSDEnsemble = ...
            obj.computeCDFOfJumps(obj.MSDEnsemble, obj.FrameLagRange);
        
        % Fit the trajectory-wise CDFs.
        if obj.FitIndividualTrajectories
            if obj.EstimateSEs
                [ParamsSingleTraj, ParamsSingleTrajSE] = ...
                    obj.fitCDFOfJumps(obj.MSDSingleTraj, ...
                    obj.FitMethod, obj.NComponents, ...
                    obj.DiffusionModel, obj.Verbose);
            else
                ParamsSingleTraj = ...
                    obj.fitCDFOfJumps(obj.MSDSingleTraj, ...
                    obj.FitMethod, obj.NComponents, ...
                    obj.DiffusionModel, obj.Verbose);
                ParamsSingleTrajSE = NaN(size(ParamsSingleTraj));
            end
            ParamsSingleTraj(1:obj.NComponents) = ...
                ParamsSingleTraj(1:obj.NComponents) ...
                .* (JumpUnitConversion.^2) .* FrameUnitConversion;
            ParamsSingleTrajSE(1:obj.NComponents) = ...
                ParamsSingleTrajSE(1:obj.NComponents) ...
                .* (JumpUnitConversion.^2) .* FrameUnitConversion;
            
            % Extract the diffusion constants and population ratios.
            DiffusionConstantSingleTraj = ...
                ParamsSingleTraj(:, 1:obj.NComponents);
            DiffusionConstantSingleTrajSE = ...
                ParamsSingleTrajSE(:, 1:obj.NComponents);
            PopulationRatiosSingleTraj = ...
                ParamsSingleTraj(:, (obj.NComponents+1):end);
            PopulationRatiosSingleTrajSE = ...
                ParamsSingleTrajSE(:, (obj.NComponents+1):end);
        end
        
        % Fit the ensemble CDF of jumps.
        if (obj.Verbose > 1)
            fprintf(['estimateDiffusionConstant(): fitting ensemble ', ...
                'CDF of jumps...\n']);
        end
        if obj.EstimateSEs
            [ParamsEnsemble, ParamsEnsembleSE] = ...
                obj.fitCDFOfJumps(obj.MSDEnsemble, ...
                obj.FitMethod, obj.NComponents, ...
                obj.DiffusionModel, obj.Verbose);
        else
            ParamsEnsemble = ...
                obj.fitCDFOfJumps(obj.MSDEnsemble, ...
                obj.FitMethod, obj.NComponents, ...
                obj.DiffusionModel, obj.Verbose);
            ParamsEnsembleSE = NaN(size(ParamsEnsemble));
        end
        ParamsEnsemble(1:obj.NComponents) = ...
            ParamsEnsemble(1:obj.NComponents).' ...
            .* (JumpUnitConversion.^2) .* FrameUnitConversion;
        ParamsEnsembleSE(1:obj.NComponents) = ...
            ParamsEnsembleSE(1:obj.NComponents).' ...
            .* (JumpUnitConversion.^2) .* FrameUnitConversion;
        
        % Extract the diffusion constants and population ratios.
        DiffusionConstantEnsemble = ...
            ParamsEnsemble(1:obj.NComponents);
        DiffusionConstantEnsembleSE = ...
            ParamsEnsembleSE(1:obj.NComponents);
        PopulationRatiosEnsemble = ...
            ParamsEnsemble((obj.NComponents+1):end);
        PopulationRatiosEnsembleSE = ...
            ParamsEnsembleSE((obj.NComponents+1):end);
    case 'LikelihoodOfJumps'
        % Modify obj.NFitPoints to reflect the fact that we're using all of
        % the data in the MLE.
        obj.NFitPoints = numel(obj.MSDEnsemble.FrameLagsAll);
        
        % Compute the trajectory-wise MLEs.
        % NOTE: This is the same model as the 'CDFOfJumps', but making an
        %       MLE from the PDF instead of fitting the CDF.
        if obj.FitIndividualTrajectories
            if obj.EstimateSEs
                [ParamsSingleTraj, ParamsSingleTrajSE] = ...
                    obj.mleOfJumps(obj.MSDSingleTraj, ...
                    obj.NComponents, obj.DiffusionModel, obj.Verbose);
            else
                ParamsSingleTraj = ...
                    obj.mleOfJumps(obj.MSDSingleTraj, ...
                    obj.NComponents, obj.DiffusionModel, obj.Verbose);
                ParamsSingleTrajSE = NaN(size(ParamsSingleTraj));
            end
            ParamsSingleTraj(1:obj.NComponents) = ...
                ParamsSingleTraj(1:obj.NComponents) ...
                .* (JumpUnitConversion.^2) .* FrameUnitConversion;
            ParamsSingleTrajSE(1:obj.NComponents) = ...
                ParamsSingleTrajSE(1:obj.NComponents) ...
                .* (JumpUnitConversion.^2) .* FrameUnitConversion;
            
            % Extract the diffusion constants and population ratios.
            DiffusionConstantSingleTraj = ...
                ParamsSingleTraj(:, 1:obj.NComponents);
            DiffusionConstantSingleTrajSE = ...
                ParamsSingleTrajSE(:, 1:obj.NComponents);
            PopulationRatiosSingleTraj = ...
                ParamsSingleTraj(:, (obj.NComponents+1):end);
            PopulationRatiosSingleTrajSE = ...
                ParamsSingleTrajSE(:, (obj.NComponents+1):end);
        end
        
        % Compute the ensemble MLEs.
        if (obj.Verbose > 1)
            fprintf(['estimateDiffusionConstant(): ', ...
                'computing MLE of jumps...\n']);
        end
        if obj.EstimateSEs
            [ParamsEnsemble, ParamsEnsembleSE] = obj.mleOfJumps(...
                obj.MSDEnsemble, obj.NComponents, ...
                obj.DiffusionModel, obj.Verbose);
        else
            ParamsEnsemble = obj.mleOfJumps(...
                obj.MSDEnsemble, obj.NComponents, ...
                obj.DiffusionModel, obj.Verbose);
            ParamsEnsembleSE = NaN(size(ParamsEnsemble));
        end
        ParamsEnsemble(1:obj.NComponents) = ...
            ParamsEnsemble(1:obj.NComponents) ...
            .* (JumpUnitConversion.^2) .* FrameUnitConversion;
        ParamsEnsembleSE(1:obj.NComponents) = ...
            ParamsEnsembleSE(1:obj.NComponents) ...
            .* (JumpUnitConversion.^2) .* FrameUnitConversion;
        
        % Extract the diffusion constants and population ratios.
        DiffusionConstantEnsemble = ...
            ParamsEnsemble(1:obj.NComponents);
        DiffusionConstantEnsembleSE = ...
            ParamsEnsembleSE(1:obj.NComponents);
        PopulationRatiosEnsemble = ...
            ParamsEnsemble((obj.NComponents+1):end);
        PopulationRatiosEnsembleSE = ...
            ParamsEnsembleSE((obj.NComponents+1):end);
end

% Store the results in the DiffusionStruct.
DiffusionStruct(1).Name = 'trajectory';
JumpUnit = smi_helpers.arrayMUX({'pixels', 'micrometers'}, obj.UnitFlag);
TimeUnit = smi_helpers.arrayMUX({'frames', 'seconds'}, obj.UnitFlag);
DiffusionStruct(1).Units = {JumpUnit; TimeUnit};
DiffusionStruct(1).FitParams = ParamsSingleTraj;
DiffusionStruct(1).FitParamsSE = ParamsSingleTrajSE;
DiffusionStruct(1).PixelSize = obj.SMF.Data.PixelSize;
DiffusionStruct(1).FrameRate = obj.SMF.Data.FrameRate;
DiffusionStruct(1).DiffusionConstant = DiffusionConstantSingleTraj;
DiffusionStruct(1).DiffusionConstantSE = DiffusionConstantSingleTrajSE;
DiffusionStruct(1).PopulationRatios = PopulationRatiosSingleTraj;
DiffusionStruct(1).PopulationRatiosSE = PopulationRatiosSingleTrajSE;
DiffusionStruct(1).NFitPoints = obj.NFitPoints;
DiffusionStruct(2).Name = 'ensemble';
DiffusionStruct(2).Units = {JumpUnit; TimeUnit};
DiffusionStruct(2).FitParams = ParamsEnsemble;
DiffusionStruct(2).FitParamsSE = ParamsEnsembleSE;
DiffusionStruct(2).PixelSize = obj.SMF.Data.PixelSize;
DiffusionStruct(2).FrameRate = obj.SMF.Data.FrameRate;
DiffusionStruct(2).DiffusionConstant = DiffusionConstantEnsemble;
DiffusionStruct(2).DiffusionConstantSE = DiffusionConstantEnsembleSE;
DiffusionStruct(2).PopulationRatios = PopulationRatiosEnsemble;
DiffusionStruct(2).PopulationRatiosSE = PopulationRatiosEnsembleSE;
DiffusionStruct(2).NFitPoints = obj.NFitPoints;
obj.DiffusionStruct = DiffusionStruct;


end