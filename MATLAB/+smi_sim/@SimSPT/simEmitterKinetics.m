function [TrajectoryStruct] = simEmitterKinetics(TrajectoryStruct, ...
    SimParams)
%simEmitterKinetics simulates blinking and photobleaching.
% This method modifies the trajectories in 'TrajectoryStruct' to simulate
% blinking (transitions from visible to dark, dark to visible) and
% photobleaching.
%
% INPUTS:
%   TrajectoryStruct: Structure containing trajectory data (see
%                     smi_sim.SimSPT.simTrajectories())
%   SimParams: Structure of simulation parameters (see
%              smi_sim.SimSPT.defineDefaultParams())
%
% OUTPUTS:
%   TrajectoryStruct: Input 'TrajectoryStruct' updated based on the
%                     simulated emitter kinetics.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Convert some parameters to units of subframes.
KOnToOffSub = SimParams.KOnToOff / SimParams.SubframeDensity;
KOffToOnSub = SimParams.KOffToOn / SimParams.SubframeDensity;
KOnToBleachSub = SimParams.KOnToBleach / SimParams.SubframeDensity;
IntensitySub = SimParams.Intensity / SimParams.SubframeDensity;

% Loop through all of the trajectories and simulate blinking.
NTraj = size(TrajectoryStruct.Trajectories, 1);
NSubframes = SimParams.NFrames * SimParams.SubframeDensity;
IsOn = logical(KOnToOffSub*zeros(NTraj, NSubframes) ...
    + ~KOnToOffSub*ones(NTraj, NSubframes));
IsBleached = zeros(NTraj, NSubframes, 'logical');
if ((KOnToOffSub>0) || (KOnToBleachSub>0))
    % Initialize the (existing) particles to an equilibrium state.
    PhotonsSub = zeros(NTraj, NSubframes, 'single');
    IsOn(:, 1) = (rand(NTraj, 1) ...
        < (KOffToOnSub/(KOnToOffSub+KOffToOnSub)));
    
    % Loop through frames and simulate the blinking on and blinking off
    % processes using an SR latch.
    for ff = 2:NSubframes
        % Define the set and reset signals.  The set signal is based on the
        % the probability of turning on sometime in frame ff (given that
        % it's not already on or bleached).  The reset signal is based on
        % the probability of blinking off or bleaching in frame ff (given
        % that it's not already off or bleached).
        RandomArray = rand(NTraj, 1);
        SetSignal = (~(IsOn(:, ff-1)|IsBleached(:, ff-1)) ...
            & (RandomArray<(1-exp(-KOffToOnSub))));
        ResetSignal = ((IsOn(:, ff-1)|IsBleached(:, ff-1)) ...
            & (RandomArray<(1-exp(-(KOnToOffSub+KOnToBleachSub)))));
        
        % Update the IsOn and IsBleached arrays for the ff-th frame. IsOn
        % is given by the standard SR latch equation in terms of the set
        % and reset signals.  IsBleached is either propagating a previous
        % bleaching signal, or is determining if the "Reset" above was a
        % blinking off or a bleaching event.  Note that we can't use the
        % same random array as above for the bleaching calculation,
        % because we're now performing a new "coin flip".
        IsOn(:, ff) = (SetSignal | (IsOn(:, ff-1)&(~ResetSignal)));
        IsBleached(:, ff) = ((IsOn(:, ff-1) ...
            & (rand(NTraj, 1)<(KOnToBleachSub/(KOnToOffSub+KOnToBleachSub)))) ...
            | IsBleached(:, ff-1));
        
        % For emitters that turned on or off, determine how many photons
        % were emitted during the frame they came on/turned off. 
        % NOTE: An emitter turning on or off can be considered a Poisson
        %       arrival event, and as such the event time within the
        %       subframe is just a uniform random number between 0 and 1 (0
        %       is the start of the frame and 1 is the end of the frame).
        TOn = rand(sum(SetSignal), 1);
        TurnedOff = (ResetSignal & IsOn(:, ff));
        TOff = rand(sum(TurnedOff), 1);
        PhotonsSub(:, ff) = IntensitySub * IsOn(:, ff);
        PhotonsSub(SetSignal, ff) = IntensitySub * (1-TOn);
        PhotonsSub(TurnedOff, ff) = IntensitySub * (1-TOff);
    end
else
    PhotonsSub = IntensitySub * ones(NTraj, NSubframes, 'single');
end

% Remove trajectories that were never visible.
% NOTE: TrajectoryStruct.IsOn must be checked so that we aren't allowing
%       trajectories to blink on before their birth from a periodic 
%       boundary.
IsOn = (IsOn & TrajectoryStruct.IsOn);
NotAlwaysOff = ~all(~IsOn, 2);
TrajectoryStruct.IsOn = IsOn(NotAlwaysOff, :);
TrajectoryStruct.DSub = TrajectoryStruct.DSub(NotAlwaysOff);
TrajectoryStruct.ConnectionMapT = ...
    TrajectoryStruct.ConnectionMapT(NotAlwaysOff, :);
TrajectoryStruct.Trajectories = ...
    TrajectoryStruct.Trajectories(NotAlwaysOff, :, :);
TrajectoryStruct.PhotonsSub = PhotonsSub(NotAlwaysOff, :);


end