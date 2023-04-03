### +smi_sim/@SimSPT

SimSPT is a collection of methods that can be used to simulate
single-particle tracking (SPT) data, including realistic raw data and
ground-truth trajectories.
Class structure is modeled on smi_sim.SimSMLM .

properties:
- SimParams = struct(); % Structure of parameters used in the simulation.

methods:
- **applyBoundaryCondition**: applies a boundary condition to Trajectories
- **applyCoordMask**: applies a (discrete) mask to the provided coordinates
- **applyLabelingEfficiency**: randomly discards trajectories to mimic labeling
- **applyMeasurementModel**: simulates measurement effects in trajectories
- **convertTrajToSMD**: converts a TrajStruct to an SMD
- **createSimulation**: creates simulated trajectories
- **defineDefaultParams**: creates a ParamStruct with all default values set
  (see the actual for further details)
- **enforcePeriodicBoundary**: breaks trajectories that hit a periodic boundary
- **makeExampleSim**: prepares a basic two channel SPT simulation
- **simEmitterKinetics**: simulates blinking and photobleaching
- **simOligomers**: simulates oligomerization between trajectories
- **simTrajBrownian**: simulates Brownian trajectories
- **simTrajectories**: simulates trajectories with oligomerization

Examples using this class can be found in MATLAB/examples:
- Example_DiffusionEstimator.m
- Example_HMM.m
- Example_SPT.m
- Example_SPTBatch.m
For example,
```
   % Simulate and save some SPT data in the format expected for real data.
   % Simulate some diffusing blobs.
   rng(12)
   SPTSim = smi_sim.SimSPT;
   SPTSim.SimParams.FrameSize = [128, 128];
   SPTSim.SimParams.ParticleDensity = 0.002; % particles / px^2
   SPTSim.SimParams.D = 0.1; % px^2 / s
   SPTSim.SimParams.KOffToOn = 0.9;
   SPTSim.SimParams.KOnToOff = 0.05;
   SPTSim.SimParams.Intensity = 1000;
   SPTSim.createSimulation()
```
See defineDefaultParams for definitions of the above parameters.
