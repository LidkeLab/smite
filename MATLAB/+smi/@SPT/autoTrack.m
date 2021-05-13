function [] = autoTrack(obj)
%autoTrack performs recursive tracking by updating tracking parameters.
% This method will recursively track the data pointed to by obj.SMF making
% new estimates of some parameters (e.g., diffusion constant) after each
% recursion.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Generate trajectories from the localizations in obj.SMD.
% NOTE: On the first pass, the diffusion constant in obj.SMF.Tracking.D
%       will be used for all cost matrices.
obj.DiffusionConstant = [];
obj.ParamsHistory = cell(0);
obj.generateTrajectories()

% If needed, perform the recursive tracking.
OnesArray = ones(numel(obj.SMD.FrameNum), 1);
obj.DiffusionEstimator.FitIndividualTrajectories = obj.UseTrackByTrackD;
for rr = 1:(obj.NRecursions-1)
    % Send an update to the command window.
    if (obj.Verbose > 1)
        fprintf(['\tsmi.spt.performFullAnalysis(): ', ...
            'tracking recursion iteration %i...\n'], rr)
    end
    
    % Estimate the diffusion constants from the previous tracking results.
    obj.DiffusionEstimator.TR = obj.TR;
    DiffusionStruct = obj.DiffusionEstimator.estimateDiffusionConstant();
    obj.SMF.Tracking.D = DiffusionStruct(2).DiffusionConstant;
    if obj.UseTrackByTrackD
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
    
    % Re-track the data.
    obj.SMD.ConnectID = [];
    obj.generateTrajectories();
end


end