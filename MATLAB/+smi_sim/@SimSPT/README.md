### +smi_sim/@SimSPT

SimSPT is a collection of methods that can be used to simulate
single-particle tracking (SPT) data, including realistic raw data and
ground-truth trajectories.
Class structure is modeled on smi_sim.SimSMLM .

properties:
- SimParams = struct(); % Structure of parameters used in the simulation.

---

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
See [defineDefaultParams](defineDefaultParams.m)
for definitions of the above parameters.

----

methods:
- **[applyBoundaryCondition](applyBoundaryCondition.m)**:
  applies a boundary condition to Trajectories
- **[applyCoordMask](applyCoordMask.m)**:
  applies a (discrete) mask to the provided coordinates
- **[applyLabelingEfficiency](applyLabelingEfficiency.m)**:
  randomly discards trajectories to mimic labeling
- **[applyMeasurementModel](applyMeasurementModel.m)**:
  simulates measurement effects in trajectories
- **[convertTrajToSMD](convertTrajToSMD.m)**:
  converts a TrajStruct to an SMD
- **[createSimulation](createSimulation.m)**:
  creates simulated trajectories
- **[defineDefaultParams](defineDefaultParams.m)**:
  creates a ParamStruct with all default values set
  (see the actual method for further details)
- **[enforcePeriodicBoundary](enforcePeriodicBoundary.m)**:
  breaks trajectories that hit a periodic boundary
- **[makeExampleSim](makeExampleSim.m)**:
  prepares a basic two channel SPT simulation
- **[simEmitterKinetics](simEmitterKinetics.m)**:
  simulates blinking and photobleaching
- **[simOligomers](simOligomers.m)**:
  simulates oligomerization between trajectories
- **[simTrajBrownian](simTrajBrownian.m)**:
  simulates Brownian trajectories
- **[simTrajectories](simTrajectories.m)**:
  simulates trajectories with oligomerization
