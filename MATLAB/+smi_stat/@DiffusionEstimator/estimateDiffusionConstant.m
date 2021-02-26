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

% Created by:
%   David J. Schodt (Lidke lab, 2021) 


% Compute the MSDs.
if (obj.Verbose > 0)
    fprintf('estimateDiffusionConstant(): computing MSDs...\n');
end
[obj.MSDSingleTraj, obj.MSDEnsemble] = ...
    obj.computeMSD(obj.TR, obj.MaxFrameLag, obj.Verbose);

% Fit the MSDs and convert units if necessary.
if (obj.Verbose > 1)
    fprintf(['estimateDiffusionConstant(): fitting trajectory-wise ', ...
        'MSDs...\n']);
elseif (obj.Verbose > 0)
    fprintf('estimateDiffusionConstant(): fitting MSDs...\n');
end
[FitParams, FitParamsSE] = ...
    obj.fitMSD(obj.MSDSingleTraj, obj.FitMethod, obj.Verbose);
DConversionFactor = ~obj.UnitFlag ...
    + obj.UnitFlag*(obj.TR(1).PixelSize^2)*obj.TR(1).FrameRate;
DiffusionStruct(1).Name = 'trajectory';
DiffusionStruct(1).Units = ...
    smi_helpers.stringMUX({'pixels, frames', 'micrometers, seconds'}, ...
    obj.UnitFlag);
DiffusionStruct(1).FitParams = FitParams;
DiffusionStruct(1).FitParamsSE = FitParamsSE;
DiffusionStruct(1).PixelSize = obj.TR(1).PixelSize;
DiffusionStruct(1).FrameRate = obj.TR(1).FrameRate;
DiffusionStruct(1).DiffusionConstant = DConversionFactor ...
    * FitParams(:, 2) / (2*2);
DiffusionStruct(1).DiffusionConstantSE = DConversionFactor ...
    * FitParamsSE(:, 2) / (2*2);
if (obj.Verbose > 1)
    fprintf(['estimateDiffusionConstant(): fitting ensemble ', ...
        'MSD...\n']);
end
[FitParams, FitParamsSE] = ...
    obj.fitMSD(obj.MSDEnsemble, obj.FitMethod, obj.Verbose);
DiffusionStruct(2).Name = 'ensemble';
DiffusionStruct(2).Units = ...
    smi_helpers.stringMUX({'pixels, frames', 'micrometers, seconds'}, ...
    obj.UnitFlag);
DiffusionStruct(2).FitParams = FitParams;
DiffusionStruct(2).FitParamsSE = FitParamsSE;
DiffusionStruct(2).PixelSize = obj.TR(1).PixelSize;
DiffusionStruct(2).FrameRate = obj.TR(1).FrameRate;
DiffusionStruct(2).DiffusionConstant = DConversionFactor ...
    * FitParams(2) / (2*2);
DiffusionStruct(2).DiffusionConstantSE = DConversionFactor ...
    * FitParamsSE(2) / (2*2);
obj.DiffusionStruct = DiffusionStruct;


end