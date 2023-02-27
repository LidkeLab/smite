Self-contained examples of using SMITE code:
   Example_ChangeDetection	change point detection
G  Example_ChannelRegstration	2-color channel registration
   Example_Clustering		various ways to invoke clustering routines
   Example_DiffusionEstimator	diffusion estimator
G  Example_GaussBlobs 		generate image stack of randomly located blobs
G  Example_HMM			hidden Markov model for dimer detection
G  Example_LocalizeData		find localizations in an image stack
   Example_PairCorrelation	compute auto- and cross-correlations of data
G  Example_SPT			single particle tracking
G  Example_SPTBatch		single particle tracking using channel reg.
   Example_StatisticsClustering	various statistical measures of clustering
   Example_simSMLM		generate synthetic data for a Siemen's star

G indicates a GPU is used.

SMITE code templates requiring user-supplied data:
   Example_Publish		generate results for a microscope experiment
   Example_Publish_generic	as above for multiple directories of data
   Example_SMLM_Basic		demonstrate basic SMLM functionality
   Example_SMLM_script		example of SMLM analysis
   hierBaGoL_wrapper		wrapper used to call BaGoL routines
   spt_resolft_track_demo	SPT-RESOLFT example

SMITE unit tests:
   smi.SMLM.unitTest
   smi.SPT.unitTestFFGC (frame-to-frame and gap closing processes)
   smi_core.ChannelRegistration.unitTest
   smi_core.DataToPhotons.unitTest
   smi_core.DriftCorrection.unitTest
   smi_core.FRC.unitTest
   smi_core.FrameConnection.unitTest
   smi_core.LocalizeData.unitTest
   smi_core.Threshold.unitTest
   smi_psf.PointSpreadFunction.unitTest (does the following unit tests)
      smi_psf.PointSpreadFunction.crlbPSFPupil_unitTest
      smi_psf.PointSpreadFunction.optimPSFZernike_unitTest
      smi_psf.PointSpreadFunction.oversamplePSFPupil_unitTest
      smi_psf.PointSpreadFunction.phaseRetrieval_unitTest
      smi_psf.PointSpreadFunction.psfROIStack_unitTest
      smi_psf.PointSpreadFunction.scalarPSFPrasadZone_unitTest
      smi_psf.PointSpreadFunction.zernikeImage_unitTest
   smi_psf.Zernike.unitTest
   smi_sim.SimSMLM.unitTest
   smi_vis.GenerateImages.blobColorOverlay_unitTest
   smi_vis.GenerateImages.circleImage_unitTest
   smi_vis.GenerateImages.colorImage_unitTest
   smi_vis.GenerateImages.driftImage_unitTest
   smi_vis.GenerateImages.gaussianImage_unitTest
   smi_vis.GenerateImages.histogramImage_unitTest
