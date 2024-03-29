function [DiffusionStruct] = ...
    estimateDiffCoeffs(TR, DiffusionEstimator, DReset)
%estimateDiffCoeffs estimates diffusion coefficients from TR.
% This method estimates diffusion coefficients for the trajectories in 'TR'
% using the provided instance of the smi_stat.DiffusionEstimator class.
%
% INPUTS:
%   TR: Tracking Results structure. (see smi_core.TrackingResults)
%   DiffusionEstimator: Instance of the smi_stat.DiffusionEstimator class.
%                       (Default = smi_stat.DiffusionEstimator)
%   DReset: If 2 elements: Min. and max. diffusion coefficients.
%           If 1 element: Value used to reset all "bad" values.
%           If empty: A bad ensemble diffusion coefficient is reset to a
%                     value of 0.01.  Individual trajectory values are
%                     reset to the ensemble value. (Default = [])
%
% OUTPUTS:
%   DiffusionStruct: DiffusionEstimator.DiffusionStruct with DReset
%                    applied.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults.
if (~exist('DiffusionEstimator', 'var') || isempty(DiffusionEstimator))
    DiffusionEstimator = smi_stat.DiffusionEstimator;
end
if ~exist('DReset', 'var')
    DReset = [];
end

% Estimate the diffusion coefficients.
DiffusionEstimator.UnitFlag = false;
DiffusionEstimator.FitIndividualTrajectories = true;
DiffusionEstimator.TR = TR;
DiffusionEstimator.estimateDiffusionConstant();
DiffusionStruct = DiffusionEstimator.DiffusionStruct;

% Filter the diffusion coefficients.
if isempty(DReset)
    % An empty DReset means we'll set a bad ensemble value to 1e-5 and bad
    % individual values to the ensemble value.
    ResetBool = ((DiffusionStruct(2).DiffusionConstant<=0) ...
        || isnan(DiffusionStruct(2).DiffusionConstant) ...
        || isinf(DiffusionStruct(2).DiffusionConstant));
    DiffusionStruct(2).DiffusionConstant(ResetBool) = 0.01;
    DiffusionStruct(2).DiffusionConstantSE(ResetBool) = inf;
    ResetBool = ((DiffusionStruct(1).DiffusionConstant<=0) ...
        | isnan(DiffusionStruct(1).DiffusionConstant) ...
        | isinf(DiffusionStruct(1).DiffusionConstant));
    DiffusionStruct(1).DiffusionConstant(ResetBool) = ...
        DiffusionStruct(2).DiffusionConstant;
    DiffusionStruct(1).DiffusionConstantSE(ResetBool) = inf;
elseif (numel(DReset) > 1)
    % DReset represents a minimum and maximum allowed value.
    SetToMin = ((DiffusionStruct(1).DiffusionConstant<DReset(1)) ...
        | isnan(DiffusionStruct(1).DiffusionConstant));
    SetToMax = (DiffusionStruct(1).DiffusionConstant > DReset(2));
    DiffusionStruct(1).DiffusionConstant(SetToMin) = DReset(1);
    DiffusionStruct(1).DiffusionConstant(SetToMax) = DReset(2);
    DiffusionStruct(1).DiffusionConstant(SetToMin | SetToMax) = inf;
    SetToMin = ((DiffusionStruct(2).DiffusionConstant<DReset(1)) ...
        || isnan(DiffusionStruct(2).DiffusionConstant));
    SetToMax = (DiffusionStruct(2).DiffusionConstant > DReset(2));
    DiffusionStruct(2).DiffusionConstant(SetToMin) = DReset(1);
    DiffusionStruct(2).DiffusionConstant(SetToMax) = DReset(2);
    DiffusionStruct(2).DiffusionConstant(SetToMin || SetToMax) = inf;
else
    % DReset is a scalar (assumed) that we'll use to set all bad values.
    ResetBool = ((DiffusionStruct(1).DiffusionConstant<=0) ...
        | isnan(DiffusionStruct(1).DiffusionConstant) ...
        | isinf(DiffusionStruct(1).DiffusionConstant));
    DiffusionStruct(1).DiffusionConstant(ResetBool) = DReset;
    DiffusionStruct(1).DiffusionConstantSE(ResetBool) = inf;
    ResetBool = ((DiffusionStruct(2).DiffusionConstant<=0) ...
        || isnan(DiffusionStruct(2).DiffusionConstant) ...
        || isinf(DiffusionStruct(2).DiffusionConstant));
    DiffusionStruct(2).DiffusionConstant(ResetBool) = DReset;
    DiffusionStruct(2).DiffusionConstantSE(ResetBool) = inf;
end


end