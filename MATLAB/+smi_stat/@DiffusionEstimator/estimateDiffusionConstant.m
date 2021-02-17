function [DiffusionStruct] = estimateDiffusionConstant(obj)
%estimateDiffusionConstant estimates the diffusion constant from an MSD.
% This method will fit the mean squared displacement data in 'MSD' to make
% an estimate of the diffusion constant.
%
% OUTPUTS:
%   DiffusionStruct: Structure array containing the fit diffusion
%                    constants and their standard errors.  All units are
%                    camera units (pixels, frames).

% Created by:
%   David J. Schodt (Lidke lab, 2021) 


% Compute the MSDs.
[obj.MSDSingleTraj, obj.MSDEnsemble] = ...
    obj.computeMSD(obj.TR, obj.MaxFrameLag);

% Fit the MSDs.
[FitParams, FitParamsSE] = obj.fitMSD(obj.MSDSingleTraj, obj.FitMethod);
DiffusionStruct(1).Name = 'trajectory';
DiffusionStruct(1).DiffusionConstant = FitParams(:, 1) / (2 * 2);
DiffusionStruct(1).DiffusionConstantSE = FitParamsSE(:, 1) / (2 * 2);
[FitParams, FitParamsSE] = obj.fitMSD(obj.MSDEnsemble, obj.FitMethod);
DiffusionStruct(2).Name = 'ensemble';
DiffusionStruct(2).DiffusionConstant = FitParams(1) / (2 * 2);
DiffusionStruct(2).DiffusionConstantSE = FitParamsSE(1) / (2 * 2);


end