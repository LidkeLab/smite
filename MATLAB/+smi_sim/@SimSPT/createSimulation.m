function [] = createSimulation(obj)
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
obj.TrajStructSubTrue = obj.simTrajectories(obj.SimParams);

% Simulate the effect of labeling efficiency.
obj.TrajStructSubLabeled = obj.applyLabelingEfficiency(...
    obj.TrajStructSubTrue, obj.SimParams.LabelingEfficiency);

% Simulate the emitter kinetics (e.g., blinking and bleaching).
obj.TrajStructSubModel = obj.simEmitterKinetics(...
    obj.TrajStructSubLabeled, obj.SimParams);

% Simulate measurement effects (e.g., motion blur, camera noise, etc.)
[obj.TrajStruct, obj.TrajStructModel] = obj.applyMeasurementModel(...
    obj.TrajStructSubModel, obj.SimParams);


end