function [SMD, SMDModel, SMDLabeled, SMDTrue, OligomerStruct] = ...
    createSimulation(obj)
%createSimulation creates simulated trajectories.
% This method is the primary user-focused method in the smi_sim.SimSPT
% class, meaning that most users will only need to use this method.  The
% intention is that this method can generate all of the vital results
% needed when simulating trajectories.  Several class properties are copied
% as outputs accessible to the user for convenience (although they can
% still be accessed in the class instance obj)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Simulate diffusing targets and, if needed, oligomerization between those
% diffusing targets.
obj.TrajectoryStruct = obj.simTrajectories(obj.SimParams);

% Simulate the effect of labeling efficiency.
obj.TrajectoryStruct = obj.applyLabelingEfficiency(obj.TrajectoryStruct, ...
    obj.SimParams.LabelingEfficiency);

% Simulate the emitter kinetics (e.g., blinking and bleaching).
obj.TrajectoryStruct = obj.simEmitterKinetics(obj.TrajectoryStruct, ...
    obj.SimParams);

% Simulate measurement effects (e.g., motion blur, camera noise, etc.)
obj.SMD = obj.applyMeasurementModel(obj.SMDModel, obj.SimParams);
SMD = obj.SMD;

% If no outputs were requested, clear them from the workspace (this is nice
% to do so the user doesn't clutter their Command Window with unrequested
% outputs).
if ~nargout
    clearvars()
end


end