function [TrajStruct] = simEmitterKinetics(TrajStruct, SimParams)
%simEmitterKinetics simulates blinking and photobleaching.
% This method modifies the trajectories in 'TrajStruct' to simulate
% blinking (transitions from visible to dark, dark to visible) and
% photobleaching.
%
% INPUTS:
%   TrajStruct: Structure containing trajectory data (see
%               smi_sim.SimSPT.simTrajectories())
%   SimParams: Structure of simulation parameters (see
%              smi_sim.SimSPT.defineDefaultParams())
%
% OUTPUTS:
%   TrajStruct: Input 'TrajStruct' updated based on the simulated emitter
%               kinetics.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Convert some parameters to units of subframes.
KOnToOffSub = SimParams.KOnToOff / SimParams.SubframeDensity;
KOffToOnSub = SimParams.KOffToOn / SimParams.SubframeDensity;
KOnToBleachSub = SimParams.KOnToBleach / SimParams.SubframeDensity;
IntensitySub = SimParams.Intensity / SimParams.SubframeDensity;

% Loop through all of the trajectories and simulate blinking.
NTraj = size(TrajStruct.Trajectories, 1);
NSubframes = SimParams.NFrames * SimParams.SubframeDensity;
IsOn = KOnToOffSub*zeros(NTraj, NSubframes, 'logical') ...
    | ~KOnToOffSub*ones(NTraj, NSubframes, 'logical');
IsBleached = zeros(NTraj, NSubframes, 'logical');
BlinkOn = zeros(NTraj, NSubframes, 'logical');
BlinkOff = zeros(NTraj, NSubframes, 'logical');
if ((KOnToOffSub>0) || (KOnToBleachSub>0))
    % Initialize the (existing) particles to an equilibrium state.
    IsOn(:, 1) = (rand(NTraj, 1) ...
        < (KOffToOnSub/(KOnToOffSub+KOffToOnSub)));
    PhotonsSub = single(IntensitySub * IsOn);
    
    % Loop through frames and simulate the emitter kinetics using SR
    % latches.
    KBRatio = KOnToBleachSub / (KOnToOffSub+KOnToBleachSub);
    for ff = 2:NSubframes
        % Define the set and reset signals.  The set signal is based on the
        % the probability of turning on sometime in frame ff (given that
        % it's not already on or bleached).  The reset signal is based on
        % the probability of blinking off or bleaching in frame ff (given
        % that it's not already off or bleached).
        RandomArray = rand(NTraj, 1);
        SetOn = (~(IsOn(:, ff-1)|IsBleached(:, ff-1)) ...
            & (RandomArray<(1-exp(-KOffToOnSub))));
        ResetOff = (IsOn(:, ff-1) ...
            & (RandomArray<(1-exp(-(KOnToOffSub+KOnToBleachSub)))));
        BlinkOn(:, ff) = SetOn;
        BlinkOff(:, ff) = ResetOff;
        
        % Define the set signal for the bleaching latch. The reset signal
        % is always false since bleaching isn't reversible.
        % NOTE: We can't use the previously sampled RandomArray here since
        %       this is a separate set of "coin flips" from those above.
        SetBleach = (IsOn(:, ff-1) ...
            & ResetOff ...
            & (rand(NTraj, 1)<KBRatio));
        
        % Update the IsOn and IsBleached arrays for the ff-th frame. IsOn
        % is given by the standard SR latch equation in terms of the set
        % and reset signals.  IsBleached is either propagating a previous
        % bleaching signal, or is determining if the "Reset" above was a
        % blinking off or a bleaching event.
        IsOn(:, ff) = (SetOn | (IsOn(:, ff-1)&(~ResetOff)));
        IsBleached(:, ff) = (SetBleach | IsBleached(:, ff-1));
        
        % For emitters that turned on or off, determine how many photons
        % were emitted during the frame they came on/turned off. 
        % NOTE: An emitter turning on or off can be considered a Poisson
        %       arrival event, and as such the event time within the
        %       subframe is just a uniform random number between 0 and 1 (0
        %       is the start of the frame and 1 is the end of the frame).
        TOn = rand(sum(SetOn), 1);
        TurnedOff = (ResetOff & IsOn(:, ff));
        TOff = rand(sum(TurnedOff), 1);
        PhotonsSub(:, ff) = IntensitySub * IsOn(:, ff);
        PhotonsSub(SetOn, ff) = IntensitySub * (1-TOn);
        PhotonsSub(TurnedOff, ff) = IntensitySub * (1-TOff);
    end
else
    PhotonsSub = IntensitySub * ones(NTraj, NSubframes, 'single');
end

% Remove trajectories that were never visible.
% NOTE: TrajectoryStruct.IsOnSub must be checked so that we aren't allowing
%       trajectories to blink on before their birth from a periodic 
%       boundary.
IsOn = (IsOn & TrajStruct.IsOn);
NotAlwaysOff = ~all(~IsOn, 2);
TrajStruct.IsOn = IsOn(NotAlwaysOff, :);
TrajStruct.IsBleached = IsBleached(NotAlwaysOff, :);
TrajStruct.BlinkOn = BlinkOn(NotAlwaysOff, :);
TrajStruct.BlinkOff = BlinkOff(NotAlwaysOff, :);
TrajStruct.D = TrajStruct.D(NotAlwaysOff);
TrajStruct.Photons = PhotonsSub(NotAlwaysOff, :);
TrajStruct.Photons_SE = TrajStruct.Photons_SE(NotAlwaysOff, :);
TrajStruct.Bg = TrajStruct.Bg(NotAlwaysOff, :);
TrajStruct.Bg_SE = TrajStruct.Bg_SE(NotAlwaysOff, :);
TrajStruct.ConnectionMapT = ...
    TrajStruct.ConnectionMapT(NotAlwaysOff, :);
TrajStruct.Trajectories = ...
    TrajStruct.Trajectories(NotAlwaysOff, :, :);
TrajStruct.Trajectories_SE = ...
    TrajStruct.Trajectories_SE(NotAlwaysOff, :, :);


end