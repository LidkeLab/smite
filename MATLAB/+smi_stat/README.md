### +smi_stat

+smi_stat is the namespace for various statistical functions and classes of
***smite***:
- [@ChangeDetection](@ChangeDetection/README.md):
  change detection analysis methods
- [@DiffusionEstimator](@DiffusionEstimator/README.md):
  diffusion estimation methods
- [@HMM](@HMM/README.md):
  hidden Markov model methods

---

Functions:
- **[bootstrapFit](bootstrapFit.m)**:
  minimizes CostFunction with constraints and performs a basic bootstrap
- **[bootstrapFitCon](bootstrapFitCon.m)**:
  minimizes CostFunction and performs a basic bootstrap
- **[computeHessian](computeHessian.m)**:
  computes the Hessian of FunctionHandle around ParamsHat
- **[findCoordAffine](findCoordAffine.m)**:
  finds an affine transform to transform Coords2 to Coords1
- **[findOffset](findOffset.m)**:
  estimates a sub-pixel offset between two stacks of images
- **[findOffsetIter](findOffsetIter.m)**:
  iteratively estimates the sub-pixel shift between images
- **[findZOffset](findZOffset.m)**:
  finds the offset between Image and Stack along Z
- **[frequencyMask](frequencyMask.m)**:
  prepares a boolean mask defining a frequency cutoff
- **[leastSquaresFit](leastSquaresFit.m)**:
  performs a least squares fit on the provided data
- **[shiftImage](shiftImage.m)**:
  shifts an image by the provided shift
