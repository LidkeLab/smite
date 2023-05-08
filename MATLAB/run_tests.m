% Run various tests on the core functionality of SMITE.

% +smi

fprintf('smi.SMLM.unitTest');
try
   smi.SMLM.unitTest()
end

fprintf('smi.SPT.unitTestFFGC (frame-to-frame and gap closing processes)');
try
   smi.SPT.unitTestFFGC()
end

% +smi_cluster

fprintf('smi_cluster.Clustering.unitTest');
try
   smi_cluster.Clustering.unitTest()
end

fprintf('smi_cluster.PairCorrelation.unitTest');
try
   smi_cluster.PairCorrelation.unitTest()
end

fprintf('smi_cluster.StatisticsClustering.unitTest');
try
   smi_cluster.StatisticsClustering.unitTest()
end

% +smi_core

fprintf('smi_core.ChannelRegistration.unitTest');
try
   smi_core.ChannelRegistration.unitTest()
end

fprintf('smi_core.DataToPhotons.unitTest (convert raw data in ADUs to photons');
try
   smi_core.DataToPhotons.unitTest()
end

fprintf('smi_core.DriftCorrection.unitTest');
try
   smi_core.DriftCorrection.unitTest()
end

fprintf('smi_core.FRC.unitTest (Fourier Ring Correlation) [DIPimage needed]');
try
   smi_core.FRC.unitTest()
end

fprintf('smi_core.FrameConnection.unitTest');
try
   smi_core.FrameConnection.unitTest()
end

fprintf('smi_core.LocalizeData.unitTest (find localizations in an image stack');
try
   smi_core.LocalizeData.unitTest()
end

fprintf('smi_core.Threshold.unitTest');
try
   smi_core.Threshold.unitTest()
end

% +smi_psf

% smi_psf.PointSpreadFunction.crlbPSFPupil_unitTest()
% smi_psf.PointSpreadFunction.optimPSFZernike_unitTest()
% smi_psf.PointSpreadFunction.oversamplePSFPupil_unitTest() [DIPimage needed]
% smi_psf.PointSpreadFunction.phaseRetrieval_unitTest() [DIPimage needed]
% smi_psf.PointSpreadFunction.psfROIStack_unitTest()
% smi_psf.PointSpreadFunction.scalarPSFPrasadZone_unitTest()
% smi_psf.PointSpreadFunction.zernikeImage_unitTest()
fprintf('smi_psf.PointSpreadFunction.unitTest (several individual tests)');
try
   smi_psf.PointSpreadFunction.unitTest()
end

fprintf('smi_psf.Zernike.unitTest');
try
   smi_psf.Zernike.unitTest()
end

% +smi_sim

fprintf('smi_sim.GaussBlobs.unitTest');
try
   smi_sim.GaussBlobs.unitTest()
end

fprintf('smi_sim.SimSMLM.unitTest');
try
   smi_sim.SimSMLM.unitTest()
end

% +smi_stat

fprintf('smi_stat.ChangeDetection.unitTest');
try
   smi_stat.ChangeDetection.unitTest()
end

fprintf('smi_stat.DiffusionEstimator.unitTest');
try
   smi_stat.DiffusionEstimator.unitTest()
end

fprintf('smi_stat.HMM.unitTest');
try
   smi_stat.HMM.unitTest()
end

% +smi_vis

fprintf('smi_vis.GenerateImages.blobColorOverlay_unitTest');
try
   smi_vis.GenerateImages.blobColorOverlay_unitTest()
end

fprintf('smi_vis.GenerateImages.circleImage_unitTest');
try
   smi_vis.GenerateImages.circleImage_unitTest()
end

fprintf('smi_vis.GenerateImages.colorImage_unitTest');
try
   smi_vis.GenerateImages.colorImage_unitTest()
end

fprintf('smi_vis.GenerateImages.driftImage_unitTest');
try
   smi_vis.GenerateImages.driftImage_unitTest()
end

fprintf('smi_vis.GenerateImages.gaussianImage_unitTest');
try
   smi_vis.GenerateImages.gaussianImage_unitTest()
end

fprintf('smi_vis.GenerateImages.histogramImage_unitTest');
try
   smi_vis.GenerateImages.histogramImage_unitTest()
end
