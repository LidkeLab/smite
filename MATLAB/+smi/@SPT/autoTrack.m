function [] = autoTrack(obj)
%autoTrack performs iterative tracking by updating tracking parameters.
% This method will iteratively track the data pointed to by obj.SMF making
% new estimates of some parameters (e.g., diffusion constant) after each
% iteration.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Perform the iterative tracking.
obj.DiffusionCoefficients = [];
NoBatchIter = (obj.SMF.Tracking.NIterMaxBatch == 1);
obj.updateTrackingParams();
ParamsHistory = {obj.SMF.Tracking};
for ii = 1:obj.SMF.Tracking.NIterMax
    % Send an update to the command window.
    if (obj.Verbose > 1)
        fprintf(['\tsmi.spt.performFullAnalysis(): ', ...
            'Tracking iteration %i...\n'], ii)
    end
    
    % Track the data.
    obj.SMD.ConnectID = [];
    obj.generateTrajectories();
    
    % Update the tracking parameters.
    obj.updateTrackingParams(obj.SMD)
    
    % Store the current set of tracking parameters.
    ParamsHistory{ii+1, 1} = obj.SMF.Tracking;
    
    % Check if the parameters have changed below the requested relative
    % change.  If so, we can stop early.
    PreviousParams = ParamsHistory{ii, 1};
    DChange = abs((PreviousParams.D-obj.SMF.Tracking.D) ...
        / PreviousParams.D);
    KOnChange = abs((PreviousParams.K_on-obj.SMF.Tracking.K_on) ...
        / PreviousParams.K_on);
    KOffChange = abs((PreviousParams.K_off-obj.SMF.Tracking.K_off) ...
        / PreviousParams.K_off);
    if (all([DChange, KOnChange, KOffChange] ...
            <= obj.SMF.Tracking.MaxRelativeChange) ...
            || (ii==obj.SMF.Tracking.NIterMax))
        % Store the parameter history (if needed) and stop iterating.
        if NoBatchIter
            obj.SMF.Tracking.ParamsHistory = ParamsHistory(1:ii);
        end
        return
    end
end


end