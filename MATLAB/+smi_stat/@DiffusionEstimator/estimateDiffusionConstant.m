function [DiffusionStruct] = estimateDiffusionConstant(obj, SaveFlag)
%estimateDiffusionConstant estimates the diffusion constant from an MSD.
% This method will fit the mean squared displacement data in 'MSD' to make
% an estimate of the diffusion constant.
%
% INPUTS:
%   SaveFlag: Boolean indicating whether or not obj.saveResults() gets
%             called within this method.  (boolean)(Default = false)
%
% OUTPUTS:
%   DiffusionStruct: Structure array containing the fit diffusion
%                    constants and their standard errors. The units will be
%                    specified by the property obj.UnitFlag, with 1
%                    specifying physical units (micrometers, seconds) and 0
%                    specifying camera units (pixels, frames).

% Created by:
%   David J. Schodt (Lidke lab, 2021) 


% Set defaults if needed.
if (~exist('SaveFlag', 'var') || isempty(SaveFlag))
    SaveFlag = false;
end

% Compute the MSDs.
[obj.MSDSingleTraj, obj.MSDEnsemble] = ...
    obj.computeMSD(obj.TR, obj.MaxFrameLag);

% Fit the MSDs and convert units if necessary.
[FitParams, FitParamsSE] = obj.fitMSD(obj.MSDSingleTraj, obj.FitMethod);
DConversionFactor = ~obj.UnitFlag ...
    + obj.UnitFlag*(obj.TR(1).PixelSize^2)*obj.TR(1).FrameRate;
DiffusionStruct(1).Name = 'trajectory';
DiffusionStruct(1).Units = ...
    smi_helpers.stringMUX({'pixels, frames', 'micrometers, seconds'}, ...
    obj.UnitFlag);
DiffusionStruct(1).DiffusionConstant = DConversionFactor ...
    * FitParams(:, 1) / (2*2);
DiffusionStruct(1).DiffusionConstantSE = DConversionFactor ...
    * FitParamsSE(:, 1) / (2*2);
[FitParams, FitParamsSE] = obj.fitMSD(obj.MSDEnsemble, obj.FitMethod);
DiffusionStruct(2).Name = 'ensemble';
DiffusionStruct(2).Units = ...
    smi_helpers.stringMUX({'pixels, frames', 'micrometers, seconds'}, ...
    obj.UnitFlag);
DiffusionStruct(2).DiffusionConstant = DConversionFactor ...
    * FitParams(1) / (2*2);
DiffusionStruct(2).DiffusionConstantSE = DConversionFactor ...
    * FitParamsSE(1) / (2*2);
obj.DiffusionStruct = DiffusionStruct;

% Save the results if requested.
if SaveFlag
    obj.saveResults();
end


end