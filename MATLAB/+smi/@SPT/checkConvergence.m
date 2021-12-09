function [Converged] = checkConvergence(NewParams, PreviousParams, ...
    MaxRelativeChange)
%checkConvergence checks if tracking parameters have converged.
% This method checks if the set of parameters in NewParams differ from
% PreviousParams within the provided tolerance to claim convergence.
%
% INPUTS:
%   NewParams: New set of SMF.Tracking parameters.
%   PreviousParams: Previous set of SMF.Tracking parameters.
%   MaxRelativeChange: Maximum relative change between checked parameters
%                      allowed for convergence.
%
% OUTPUTS:
%   Converged: Flag indicating convergence (true) or failure to converge
%              (false).

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Check relative differences in parameters and decide if they converged.
DChange = (median(PreviousParams.D)-median(NewParams.D)) ...
    / median(PreviousParams.D);
KOnChange = abs((PreviousParams.K_on-NewParams.K_on) ...
    / PreviousParams.K_on);
KOffChange = abs((PreviousParams.K_off-NewParams.K_off) ...
    / PreviousParams.K_off);
Converged = all([DChange, KOnChange, KOffChange] <= MaxRelativeChange);


end