### smi_core.DriftCorrection

driftCorrectKNN performs drift correction on 2D or 3D data
provided in an SMD structure using K nearest neighbor (KNN) searching,
returning an updated structure with drift corrected coordinates.  Plots of
the drift estimates can be produced with plotDriftCorrection and some
additional measures with calcDCResidual.

EXAMPLE USAGE (see also unitTest):
   DC = smi_core.DriftCorrection(SMF, SMDin);
   SMDIntra = [];
   for i = 1 : NDatasets
      [SMDIntra_i, StatisticsIntra] = DC.driftCorrectKNNIntra(SMDin_i, i, i);
      SMDIntra = smi_core.SingleMoleculeData.catSMD(SMDIntra, SMDIntra_i, false);
   end
   [SMDInter, StatisticsInter] = DC.driftCorrectKNNInter(SMDIntra);
   SMDout = SMDInter;

properties:
```
   % Intra-dataset threshold (pixel)
   L_intra        = 1;
   % Inter-dataset threshold (pixel)
   L_inter        = 2;
   % X/Y pixel size in um (only needed for 3D drift correction)
   PixelSizeZUnit = 0.1;
   % Degree of the intra-dataset fitting polynomial for drift rate
   PDegree        = 1;
   % Termination tolerance on the intra-dataset function value
   TolFun_intra   = 1e-2;
   % Termination tolerance on the intra-dataset fitting polynomial
   TolX_intra     = 1e-4;
   % Termination tolerance on the inter-dataset function value
   TolFun_inter   = 1e-2;
   % Termination tolerance on the inter-dataset fitting polynomial
   TolX_inter     = 1e-4;
   % Initialization wrt the previous dataset for inter-dataset drift correction
   % The value should be either 0 (no initial drift), 1 (initial drift of the
   % previous dataset) or SMD.NFrames (final drift); zero or initial drift
   % should work well with brightfield registration, while final drift works
   % well generally (but the optimization process may not converge quite as
   % quickly).
   Init_inter     = 0;
   % DriftCorrectKNNInterPair, instead of using Init_inter, uses the actual
   % values of P0_inter to initialize the minimizer search.  The default of []
   % indicates that the initial P0 value will be all zeros, otherwise the
   % values in P0_inter will be used directly.
   P0_inter       = [];
   % Semi-redundant variable, needed because SMD may not exist when the
   % constructor is invoked, so Init_inter has to be set in
   % driftCorrectKNNInter (only needed when breaking intra-dataset and
   % inter-dataset calculations up).
   BFRegistration = true;
   % If non-empty, override the collected value of number of datasets
   NDatasets      = [];
   % If non-empty, override the collected value of number of frames per dataset
   NFrames        = [];
   % Verbosity level
   Verbose        = 1;
   SMF            = [];
```
methods:
- **calcDCRMSE**:
  calculates the RMSE of SMD relative to true coordinates/curves
- **changeInterRef**:
  shifts SMD coordinates to new dataset reference
- **driftCorrectBF**:
  performs drift correction from brightfield images stored in an h5 file
- **driftCorrectBFInit**:
  initializes drift correction from brightfield images
- **driftCorrectBFInter**:
  computes inter-DS drift correction from brightfield images
- **driftCorrectBFIntra**:
  computes intra-DS drift correction from brightfield images
- **driftCorrectKNN**:
  calculates the drift directly from X,Y{,Z} coordinates
  by fitting a polynomial depending on time (i.e., frame number) to the frames
  with each dataset (intra-dataset), and fitting constant shifts between
  datasets (inter-dataset)
- **driftCorrectKNNInter**:
  Inter-dataset portion of **driftCorrectKNN**
- **driftCorrectKNNInterPair**:
  calculates inter-dataset drift directly from X,Y{,Z}
  coordinates (i.e., constant shifts between datasets)
- **driftCorrectKNNIntra**:
  Intra-dataset portion of **driftCorrectKNN**
- **minD**:
  Sum of nearest neighbor distances for intra/inter-dataset drift correction
- **plotCumDrift**:
  creates cumulative plots for Drift in any direction
- **plotDriftAndReg**:
  plots drift correction and brightfield corrections
- **plotDriftCorrection**:
  plots the computed drift correction stored in SMD
  structure for 2D or 3D data.  The plot is color coded so as to indicate the
  drift correction as a function of time
- **plotXYDriftParametric**:
  makes a parametric plot of the x,y drift model
- **regViaDC**:
  (registration via drift correction) takes two differently labeled
  data collections of the same biological phenomenon and attempts to align them
  using inter-dataset drift correction
- **unitTest**:
  tests smi_core.DriftCorrection.driftCorrectKNN
