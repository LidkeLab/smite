### +smi_stat

+smi_stat is the namespace for various statistical functions and classes of
***smite***:
- [@ChangeDetection](@ChangeDetection/README.md):
- [@DiffusionEstimator](@DiffusionEstimator/README.md):
- [@HMM](@HMM/README.md):

Functions:
- bootstrapFit:
  minimizes CostFunction with constraints and performs a basic bootstrap
- bootstrapFitCon:
  minimizes CostFunction and performs a basic bootstrap
- computeHessian:
  computes the Hessian of FunctionHandle around ParamsHat
- findCoordAffine:
  finds an affine transform to transform Coords2 to Coords1
- findOffset:
  estimates a sub-pixel offset between two stacks of images
- findOffsetIter:
  iteratively estimates the sub-pixel shift between images
- findZOffset:
  finds the offset between Image and Stack along Z
- frequencyMask:
  prepares a boolean mask defining a frequency cutoff
- leastSquaresFit:
  performs a least squares fit on the provided data
- shiftImage:
  shifts an image by the provided shift
