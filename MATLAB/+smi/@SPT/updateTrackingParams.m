function [] = updateTrackingParams(obj, IsBatch)
%updateTrackingParams updates tracking parameters based on previous results
% This method will attempt to update as many SMF.Tracking parameters as
% possible from obj.SMD, obj.SMDBatch, obj.TR, ..., as well as computing
% DiffusionCoefficients for the trajectories in obj.SMD.  The intention is 
% that this can simplify iterative tracking by making all parameter updates
% in one methods.
%
% INPUTS:
%   IsBatch: Flag indicating parameters should be estimated for the entire
%            batch if possible.  Note that only the rate parameters will be
%            affected by this flag, as the density and diffusion
%            coefficients should always be estimated for the current file
%            being tracked.  (Default = true)
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults.
if (~exist('IsBatch', 'var') || isempty(IsBatch))
    IsBatch = true;
end

% If SMD has the field ConnectID populated, estimate diffusion constants
% and rate parameters.
if isempty(obj.SMD.ConnectID)
    % Make sure SMF.Tracking.D matches the size of SMD.
    obj.SMF.Tracking.D = padarray(obj.SMF.Tracking.D, ...
        [max(0, numel(obj.SMD.FrameNum)-numel(obj.SMF.Tracking.D)), 0], ...
        median(obj.SMF.Tracking.D(:)), 'post');
else
    % Update the diffusion coefficients, storing results so that indices
    % match indices of obj.SMD.
    obj.DiffusionEstimator.FitIndividualTrajectories = ...
        obj.SMF.Tracking.TrajwiseD;
    DiffusionStruct = obj.estimateDiffCoeffs(obj.TR, ...
        obj.DiffusionEstimator, median(obj.SMF.Tracking.D));
    obj.SMF.Tracking.D = ...
        DiffusionStruct(2).DiffusionConstant * ones(size(obj.SMD.FrameNum));
    if obj.SMF.Tracking.TrajwiseD
        for ii = 1:numel(obj.TR)
            obj.SMF.Tracking.D(obj.TR(ii).IndSMD, 1) = ...
                DiffusionStruct(1).DiffusionConstant(ii);
        end
    end
    
    % Update our rate parameter estimates.  If batch tracking, we should
    % try to update 
    if (~isempty(obj.SMDBatch) && IsBatch)
        [obj.SMF.Tracking.K_on, obj.SMF.Tracking.K_off] = ...
            obj.estimateRateParameters(obj.SMDBatch);
    else
        [obj.SMF.Tracking.K_on, obj.SMF.Tracking.K_off] = ...
            obj.estimateRateParameters(obj.SMD);
    end
end

% Estimate the density of off emitters for the current file.
obj.RhoOff = obj.estimateDensities(obj.SMD, obj.SMF);


end