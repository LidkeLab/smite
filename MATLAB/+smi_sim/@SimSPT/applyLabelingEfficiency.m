function [TrajectoryStruct, KeepInd] = applyLabelingEfficiency(...
    TrajectoryStruct, LabelingEfficiency)
%applyLabelingEfficiency randomly discards trajectories to mimic labeling.
% This method will randomly discard trajectories from TrajectoryStruct
% based on the labeling efficiency.  Note that not all fields of 
% TrajectoryStruct will be updated (e.g., the PeriodicityMapT, which will
% not be updated due to complications in doing so).
%
% INPUTS:
%   TrajectoryStruct: Structure with information about the trajectories
%                     (see, e.g., usage in
%                     smi_sim.SimSPT.simTrajectories())
%   LabelingEfficiency: Labeling efficiency (i.e., the probability that a
%                       trajectory will be labeled). 
%                       (scalar between 0 and 1)(Default = 1)
%
% OUTPUTS:
%   TrajectoryStruct: Input 'TrajectoryStruct' after discarding unlabeled
%                     trajectories.
%   KeepInd: Indices from the original input 'TrajectoryStruct'
%            trajectories that were kept after labeling.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults.
if (~exist('LabelingEfficiency', 'var') || isempty(LabelingEfficiency))
    LabelingEfficiency = 1;
end

% Simulate labeling.
NTrajectories = size(TrajectoryStruct.Trajectories, 1);
KeepInd = sort(randperm(NTrajectories, round(LabelingEfficiency*NTrajectories)));
TrajectoryStruct.Trajectories = TrajectoryStruct.Trajectories(KeepInd, :, :);
TrajectoryStruct.IsOn = TrajectoryStruct.IsOn(KeepInd, :);
TrajectoryStruct.DSub = TrajectoryStruct.DSub(KeepInd);
TrajectoryStruct.ConnectionMapT = TrajectoryStruct.ConnectionMapT(KeepInd, :);


end