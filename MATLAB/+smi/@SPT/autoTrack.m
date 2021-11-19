function [] = autoTrack(obj)
%autoTrack performs iterative tracking by updating tracking parameters.
% This method will iteratively track the data pointed to by obj.SMF making
% new estimates of some parameters (e.g., diffusion constant) after each
% iteration.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Estimate the density of off emitters (if requested).
if obj.SMF.Tracking.EstimateRho
    RhoOnMean = ...
        mean(smi_core.SingleMoleculeData.computeDensity(obj.SMD));
    obj.SMF.Tracking.Rho_off = RhoOnMean ...
        * (obj.SMF.Tracking.K_off/obj.SMF.Tracking.K_on);
end
    
% Generate trajectories from the localizations in obj.SMD.
% NOTE: On the first pass, the diffusion constant in obj.SMF.Tracking.D
%       will be used for all cost matrices.
obj.DiffusionConstant = [];
obj.SMF.Tracking.ParamsHistory = {obj.SMF.Tracking};
obj.generateTrajectories()

% If needed, perform the iterative tracking.
OnesArray = ones(numel(obj.SMD.FrameNum), 1);
obj.DiffusionEstimator.FitIndividualTrajectories = ...
    obj.SMF.Tracking.TrajwiseD;
for rr = 1:(obj.SMF.Tracking.NIterationsMax-1)
    % Send an update to the command window.
    if (obj.Verbose > 1)
        fprintf(['\tsmi.spt.performFullAnalysis(): ', ...
            'Tracking iteration %i...\n'], rr)
    end
    
    % Estimate the diffusion constants from the previous tracking results.
    obj.DiffusionEstimator.TR = obj.TR;
    DiffusionStruct = obj.DiffusionEstimator.estimateDiffusionConstant();
    obj.SMF.Tracking.D = DiffusionStruct(2).DiffusionConstant;
    if obj.SMF.Tracking.TrajwiseD
        % Store the trajectory-wise diffusion constants in the TR
        % structures.
        DiffusionConstantCurrent = ...
            [DiffusionStruct(1).DiffusionConstant, ...
            DiffusionStruct(1).DiffusionConstantSE];
        obj.DiffusionConstant = OnesArray ...
            * [DiffusionStruct(2).DiffusionConstant, ...
            DiffusionStruct(2).DiffusionConstantSE];
        for ii = 1:numel(obj.TR)
            obj.DiffusionConstant(obj.TR(ii).IndSMD, 1) = ...
                DiffusionConstantCurrent(ii, 1);
            obj.DiffusionConstant(obj.TR(ii).IndSMD, 2) = ...
                DiffusionConstantCurrent(ii, 2);
        end
        
        % Filter the diffusion constants, setting those which are invalid
        % to the previously estimated ensemble value. For the corresponding
        % SEs, we'll set the bad values to an SE of inf (since we don't
        % actually know the value, and having a low SE might prevent those
        % localizations from being connected to others).
        BadValueBool = ((obj.DiffusionConstant(:, 1)<0) ...
            | isnan(obj.DiffusionConstant(:, 1)) ...
            | isinf(obj.DiffusionConstant(:, 1)));
        if isnan(DiffusionStruct(2).DiffusionConstant)
            % This can happen when no good trajectories were formed, but we
            % may still wish to iterate.
            obj.DiffusionConstant(BadValueBool, 1) = obj.SMF.Tracking.D;
        else
            obj.DiffusionConstant(BadValueBool, 1) = ...
                DiffusionStruct(2).DiffusionConstant;
        end
        obj.DiffusionConstant(BadValueBool, 2) = inf;
    end
    
    % Estimate the blinking rates from the previous tracking results.
    [KOn, KOff] = obj.estimateRateParameters(obj.SMD);
    obj.SMF.Tracking.K_on = KOn;
    obj.SMF.Tracking.K_off = KOff;
    
    % Estimate the density of off emitters (if requested).
    if obj.SMF.Tracking.EstimateRho
        RhoOnMean = ...
            mean(smi_core.SingleMoleculeData.computeDensity(obj.SMD));
        obj.SMF.Tracking.Rho_off = RhoOnMean ...
            * (obj.SMF.Tracking.K_off/obj.SMF.Tracking.K_on);
    end
    
    % Check if the parameters have changed below the requested relative
    % change.  If so, we can stop early.
    PreviousParams = obj.SMF.Tracking.ParamsHistory{rr};
    DChange = abs((PreviousParams.D-obj.SMF.Tracking.D) ...
        / PreviousParams.D);
    KOnChange = abs((PreviousParams.K_on-obj.SMF.Tracking.K_on) ...
        / PreviousParams.K_on);
    KOffChange = abs((PreviousParams.K_off-obj.SMF.Tracking.K_off) ...
        / PreviousParams.K_off);
    RhoOffChange = abs((PreviousParams.Rho_off-obj.SMF.Tracking.Rho_off) ...
        / PreviousParams.Rho_off);
    if  all([DChange, KOnChange, KOffChange, RhoOffChange] ...
            <= obj.SMF.Tracking.MaxRelativeChange)
        return
    end
    % Store the current set of tracking parameters.
    obj.SMF.Tracking.ParamsHistory{rr+1, 1} = obj.SMF.Tracking;
    
    % Re-track the data.
    obj.SMD.ConnectID = [];
    obj.generateTrajectories();
end


end