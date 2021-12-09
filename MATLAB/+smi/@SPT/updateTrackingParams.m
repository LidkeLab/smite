function [] = updateTrackingParams(obj, SMD)
%updateTrackingParams updates tracking parameters based on SMD and obj.TR.
% This method will attempt to update as many SMF.Tracking parameters as
% possible from the provided SMD and obj.TR, as well as computing
% DiffusionCoefficients for the trajectories in SMD.  The intention is that
% this can simplify iterative tracking by making all parameter updates in
% one methods.
%
% INPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%        containing localizations associated by SMD.ConnectID.  Note that
%        this is provided as an optional input in case some larger set of
%        trajectory data is being provided than what is stored in obj.SMD
%        (e.g., if batch-tracking over multiple files, obj.SMD might be for
%        only one file, but we want to update rate parameter estimates
%        based on the batch results provided by the input SMD).
%        (Default = obj.SMD)
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults.
if (~exist('SMD', 'var') || isempty(SMD))
    SMD = obj.SMD;
end

% If SMD has the field ConnectID populated, estimate diffusion constants
% and rate parameters.
if isempty(SMD.ConnectID)
    % Make sure SMF.Tracking.D matches the size of SMD.
    obj.SMF.Tracking.D = padarray(obj.SMF.Tracking.D, ...
        [max(0, numel(SMD.FrameNum)-numel(obj.SMF.Tracking.D)), 0], ...
        median(obj.SMF.Tracking.D), 'post');
else
    % Update the diffusion coefficients, storing results so that indices
    % match indices of SMD.
    obj.DiffusionEstimator.FitIndividualTrajectories = ...
        obj.SMF.Tracking.TrajwiseD;
    DiffusionStruct = obj.estimateDiffCoeffs(obj.TR, ...
        obj.DiffusionEstimator, median(obj.SMF.Tracking.D));
    obj.SMF.Tracking.D = ...
        DiffusionStruct(2).DiffusionConstant * ones(size(SMD.FrameNum));
    if obj.SMF.Tracking.TrajwiseD
        for ii = 1:numel(obj.TR)
            obj.SMF.Tracking.D(obj.TR(ii).IndSMD, 1) = ...
                DiffusionStruct(1).DiffusionConstant(ii);
        end
    end
    
    % Update our rate parameter estimates.
    [obj.SMF.Tracking.K_on, obj.SMF.Tracking.K_off] = ...
        obj.estimateRateParameters(SMD);
end

% Estimate the density of off emitters.
obj.RhoOff = obj.estimateDensities(SMD, obj.SMF);


end