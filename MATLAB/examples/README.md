### ***smite*** examples and unit tests

Self-contained examples of using SMITE code:

|| example code ||
---|---|---
&nbsp;| [Example_ChangeDetection](Example_ChangeDetection.m)                    | change point detection
G     | [Example_ChannelRegistration](Example_ChannelRegistration.m)              | 2-color channel registration
&nbsp;| [Example_Clustering](Example_Clustering.m)                              | various ways to invoke clustering routines
&nbsp;| [Example_DiffusionEstimator](Example_DiffusionEstimator.m)              | diffusion estimator
G     | [Example_GaussBlobs](Example_GaussBlobs.m)                              | generate image stack of randomly located blobs
G     | [Example_HMM](Example_HMM.m)                                            | hidden Markov model for dimer detection
G     | [Example_LocalizeData](Example_LocalizeData.m)                          | find localizations in an image stack
&nbsp;  | [Example_PairCorrelation](Example_PairCorrelation.m)                  | compute auto- and cross-correlations of data
G     | [Example_SPT](Example_SPT.m)                                            | single particle tracking
G     | [Example_SPTBatch](Example_SPTBatch.m)                                  | single particle tracking via batch processing
&nbsp;| [Example_StatisticsClustering](Example_StatisticsClustering.m)          | various statistical measures of clustering
&nbsp;| [Example_simSMLM](Example_simSMLM.m)                                    | generate synthetic data for a Siemen's star

G indicates a GPU is used.

SMITE code templates requiring user-supplied data:

| example code ||
---|---
[Example_Publish](Example_Publish.m)                                            | generate results for a microscope experiment
[Example_Publish_generic](Example_Publish_generic.m)                            | as above for multiple directories of data
[Example_SMLM_Basic](Example_SMLM_Basic.m)                                      | demonstrate basic SMLM functionality
[Example_SMLM_script](Example_SMLM_script.m)                                    | example of SMLM analysis
[hierBaGoL_wrapper](hierBaGoL_wrapper.m)                                        | wrapper used to call BaGoL routines
[plotROIDriver](plotROIDriver.m)                                                | plot dot, Gaussian or circle images of ROIs
[simplePairCorr](simplePairCorr.m)                                              | step-by-step script to choose 2-label ROIs & do various analyses for ROIs separate or combined
[simpleROIcluster](simpleROIcluster.m)                                          | step-by-step script to choose ROIs, cluster and analyze/compare conditions for 1-label data
[singleConditionDriver](singleConditionDriver.m)                                | batch cluster analysis for comparison of experimental conditions for 1-label data
[spt_resolft_track_demo](spt_resolft_track_demo.m)                              | SPT-RESOLFT example

SMITE unit tests (see also [ExpectedResults](../ExpectedResults/README.md)):
- [smi.SMLM.unitTest](../+smi/@SMLM/unitTest.m)
- [smi.SPT.unitTestFFGC](../+smi/@SPT/unitTestFFGC.m) (frame-to-frame and gap closing processes)
- [smi_cluster.Clustering.unitTest](../+smi_cluster/@Clustering/unitTest.m)
- [smi_cluster.PairCorrelation.unitTest](../+smi_cluster/@PairCorrelation/unitTest.m)
- [smi_cluster.StatisticsClustering.unitTest](../+smi_cluster/@StatisticsClustering/unitTest.m)
- [smi_core.ChannelRegistration.unitTest](../+smi_core/@ChannelRegistration/unitTest.m)
- [smi_core.DataToPhotons.unitTest](../+smi_core/@DataToPhotons/unitTest.m)
- [smi_core.DriftCorrection.unitTest](../+smi_core/@DriftCorrection/unitTest.m)
- [smi_core.FrameConnection.unitTest](../+smi_core/@FrameConnection/unitTest.m)
- [smi_core.LocalizeData.unitTest](../+smi_core/@LocalizeData/unitTest.m)
- [smi_core.Threshold.unitTest](../+smi_core/@Threshold/unitTest.m)
- [smi_psf.PointSpreadFunction.unitTest](../+smi_psf/@PointSpreadFunction/unitTest) (does the following unit tests)
  - [smi_psf.PointSpreadFunction.crlbPSFPupil_unitTest](../+smi_psf/@PointSpreadFunction/crlbPSFPupil_unitTest.m)
  - [smi_psf.PointSpreadFunction.optimPSFZernike_unitTest](../+smi_psf/@PointSpreadFunction/optimPSFZernike_unitTest.m)
  - [smi_psf.PointSpreadFunction.psfROIStack_unitTest](../+smi_psf/@PointSpreadFunction/psfROIStack_unitTest.m)
  - [smi_psf.PointSpreadFunction.scalarPSFPrasadZone_unitTest](../+smi_psf/@PointSpreadFunction/scalarPSFPrasadZone_unitTest.m)
  - [smi_psf.PointSpreadFunction.zernikeImage_unitTest](../+smi_psf/@PointSpreadFunction/zernikeImage_unitTest.m)
- [smi_psf.Zernike.unitTest](../+smi_psf/@Zernike/unitTest.m)
- [smi_sim.GaussBlobs.unitTest](../+smi_sim/@GaussBlobs.m)
- [smi_sim.SimSMLM.unitTest](../+smi_sim/@SimSMLM/unitTest.m)
- [smi_stat.ChangeDetection.unitTest](../+smi_stat/@ChangeDetection/unitTest.m)
- [smi_stat.DiffusionEstimator.unitTest](../+smi_stat/@DiffusionEstimator/unitTest.m)
- [smi_stat.HMM.unitTest](../+smi_stat/@HMM/unitTest.m)
- [smi_vis.GenerateImages.blobColorOverlay_unitTest](../+smi_vis/@GenerateImages/blobColorOverlay_unitTest.m)
- [smi_vis.GenerateImages.circleImage_unitTest](../+smi_vis/@GenerateImages/circleImage_unitTest.m)
- [smi_vis.GenerateImages.colorImage_unitTest](../+smi_vis/@GenerateImages/colorImage_unitTest.m)
- [smi_vis.GenerateImages.driftImage_unitTest](../+smi_vis/@GenerateImages/driftImage_unitTest.m)
- [smi_vis.GenerateImages.gaussianImage_unitTest](../+smi_vis/@GenerateImages/gaussianImage_unitTest.m)
- [smi_vis.GenerateImages.histogramImage_unitTest](smi_vis/@GenerateImages/histogramImage_unitTest.m)
