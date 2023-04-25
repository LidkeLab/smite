% Run various tests on the core functionality of SMITE.

% +smi

fprintf('smi.SMLM.unitTest');
smi.SMLM.unitTest()

fprintf('smi.SPT.unitTestFFGC (frame-to-frame and gap closing processes)');
smi.SPT.unitTestFFGC()

% +smi_cluster

fprintf('smi_cluster.Clustering.unitTest');
smi_cluster.Clustering.unitTest()

fprintf('smi_cluster.PairCorrelation.unitTest');
smi_cluster.PairCorrelation.unitTest()

fprintf('smi_cluster.StatisticsClustering.unitTest');
smi_cluster.StatisticsClustering.unitTest()

% +smi_core

fprintf('smi_core.ChannelRegistration.unitTest');
smi_core.ChannelRegistration.unitTest()

fprintf('smi_core.DataToPhotons.unitTest (convert raw data in ADUs to photons');
smi_core.DataToPhotons.unitTest()

fprintf('smi_core.DriftCorrection.unitTest');
smi_core.DriftCorrection.unitTest()

fprintf('smi_core.FRC.unitTest (Fourier Ring Correlation)');
smi_core.FRC.unitTest()

fprintf('smi_core.FrameConnection.unitTest');
smi_core.FrameConnection.unitTest()

fprintf('smi_core.LocalizeData.unitTest (find localizations in an image stack');
smi_core.LocalizeData.unitTest()

fprintf('smi_core.Threshold.unitTest');
smi_core.Threshold.unitTest()

% +smi_psf

% smi_psf.PointSpreadFunction.crlbPSFPupil_unitTest()
% smi_psf.PointSpreadFunction.optimPSFZernike_unitTest()
% smi_psf.PointSpreadFunction.oversamplePSFPupil_unitTest()
% smi_psf.PointSpreadFunction.phaseRetrieval_unitTest()
% smi_psf.PointSpreadFunction.psfROIStack_unitTest()
% smi_psf.PointSpreadFunction.scalarPSFPrasadZone_unitTest()
% smi_psf.PointSpreadFunction.zernikeImage_unitTest()
fprintf('smi_psf.PointSpreadFunction.unitTest (several individual tests)');
smi_psf.PointSpreadFunction.unitTest()

fprintf('smi_psf.Zernike.unitTest');
smi_psf.Zernike.unitTest()

% +smi_sim

fprintf('smi_sim.GaussBlobs.unitTest');
smi_sim.GaussBlobs.unitTest()

fprintf('smi_sim.SimSMLM.unitTest');
smi_sim.SimSMLM.unitTest()

% +smi_stat

fprintf('smi_stat.ChangeDetection.unitTest');
smi_stat.ChangeDetection.unitTest()

fprintf('smi_stat.DiffusionEstimator.unitTest');
smi_stat.DiffusionEstimator.unitTest()

fprintf('smi_stat.HMM.unitTest');
smi_stat.HMM.unitTest()

% +smi_vis

fprintf('smi_vis.GenerateImages.blobColorOverlay_unitTest');
smi_vis.GenerateImages.blobColorOverlay_unitTest()

fprintf('smi_vis.GenerateImages.circleImage_unitTest');
smi_vis.GenerateImages.circleImage_unitTest()

fprintf('smi_vis.GenerateImages.colorImage_unitTest');
smi_vis.GenerateImages.colorImage_unitTest()

fprintf('smi_vis.GenerateImages.driftImage_unitTest');
smi_vis.GenerateImages.driftImage_unitTest()

fprintf('smi_vis.GenerateImages.gaussianImage_unitTest');
smi_vis.GenerateImages.gaussianImage_unitTest()

fprintf('smi_vis.GenerateImages.histogramImage_unitTest');
smi_vis.GenerateImages.histogramImage_unitTest()
