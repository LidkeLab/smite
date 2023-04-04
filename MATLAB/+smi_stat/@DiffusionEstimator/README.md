### +smi_stat/@DiffusionEstimator

DiffusionEstimator contains several methods used to estimate diffusion
constants from single-particle tracking data.

REQUIRES:
- Optimization Toolbox (to use fmincon() for certain models)

For a detailed example, see
[Example_DiffusionEstimator.m](../../examples/Example_DiffusionEstimator.m).

---

```
properties:
   % ID of the diffusion model to be considered (char array/string)
   % OPTIONS:
   %   'Brownian'
   DiffusionModel{mustBeMember(DiffusionModel, {'Brownian'})} = ...
       'Brownian';
   
   % Number of diffusing components. (Default = 1)
   % This can only be changed for FitTarget = CDFOfJumps or
   % LikelihoodOfJumps.
   NComponents = 1;
   
   % Fit method for fitting data (char array/string)
   FitMethod{mustBeMember(FitMethod, {'WeightedLS', 'LS'})} = ...
       'WeightedLS';
   
   % Target data that will be fit (char array/string)
   FitTarget{mustBeMember(FitTarget, ...
       {'MSD', 'CDFOfJumps', 'LikelihoodOfJumps'})} = 'MSD';

   % Range of frame lags used to estimate D (Default = [1, 5])
   FrameLagRange = [1, 5];
   
   % Number of MSD points to be fit (scalar, integer)(Default = 5)
   NFitPoints = 5;
   
   % Flag to estimate standard errors (Default = true)
   % NOTE: For some of the fit methods/fit targets (e.g., FitTarget =
   %       'CDFOfJumps') we have to use a bootstrap, in which case
   %       estimating SEs can be unreasonably slow.
   EstimateSEs = true;
   
   % Flag to fit individual trajectory MSDs/CDFs (Default = true)
   FitIndividualTrajectories = true;
   
   % Number of spatial dimensions (scalar, integer)(Default = 2)
   NDimensions = 2;
   
   % Directory in which results will be saved by saveResults().
   SaveDir = pwd();
   
   % Base name of saved results. Default defined in obj.saveDir().
   BaseSaveName
   
   % Tracking results structure.
   TR
   
   % Single Molecule Fitting structure, for pixel size and framerate.
   SMF = smi_core.SingleMoleculeFitting;
   
   % Boolean flag to indicate units of outputs (boolean)(Default = 0)
   % 1 (true) will make the outputs of estimateDiffusionConstant()
   %   micrometers and seconds.
   % 0 (false) will make the outputs of estimateDiffusionConstant()
   %   pixels and frames.
   % NOTE: Most methods of this class will use pixels and frames
   %       regardless of obj.UnitFlag. This property will only affect
   %       the "user-facing" wrapper methods, such as
   %       estimateDiffusionConstant() and saveResults()
   UnitFlag = false;
        
   % Verbosity level for estimateDiffusionConstant() (Default = 0)
   %   Verbose 0: no Command Window outputs
   %   Verbose 1: General progress printed to Command Window
   %   Verbose 2: Display some intermediate results
   %   Verbose 3: Debugging mode, extensive outputs
   Verbose = 0;
```

---

methods:
- **[brownianJumpCDF](brownianJumpCDF.m)**:
  generates a model of the CDF of Brownian jumps
- **[brownianJumpLikelihood](brownianJumpLikelihood.m)**:
  computes the likelihood of given Brownian jumps
- **[computeCDFOfJumps](computeCDFOfJumps.m)**:
  computes the CDF (CPD) of the jumps in 'MSDStruct'.
- **[computeMSD](computeMSD.m)**:
  computes the mean squared displacement from TR
- **[computeSingleTrajMSD](computeSingleTrajMSD.m)**:
  computes the mean squared displacement from TR
- **[estimateDiffusionConstant](estimateDiffusionConstant.m)**:
  estimates the diffusion constant
- **[fitCDFOfJumps](fitCDFOfJumps.m)**:
  fits the CDF of displacement data (jumps)
- **[fitCDFOfJumpsBrownian](fitCDFOfJumpsBrownian.m)**:
  fits the CDF of jumps to a Brownian motion model
- **[fitMSD](fitMSD.m)**:
  fits mean squared displacement data
- **[fitMSDBrownian](fitMSDBrownian.m)**:
  fits an MSD to a Brownian motion model
- **[mleOfJumps](mleOfJumps.m)**:
  finds the MLE from the likelihood of observed jumps
- **[mleOfJumpsBrownian](mleOfJumpsBrownian.m)**:
  finds the MLE for Brownian motion jumps
- **[plotEnsembleCDFOfJumps](plotEnsembleCDFOfJumps.m)**:
  plots the CDF of jumps and an associated fit
- **[plotEnsembleMSD](plotEnsembleMSD.m)**:
  plots an ensemble MSD and an associated fit
- **[saveResults](saveResults.m)**:
  saves useful results of diffusion estimation analysis
