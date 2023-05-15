% Run various tests on the core functionality of SMITE.  Much output will be
% saved in tempdir/smite/unitTest/name_of_test.  ExpectedResults are provided
% in this directory where very large files have been deleted so as to not bloat
% up the the SMITE distribution.

% +smi

fprintf('smi.SMLM.unitTest\n');
try
   smi.SMLM.unitTest()
end

fprintf('smi.SPT.unitTestFFGC (frame-to-frame and gap closing processes)\n');
try
   smi.SPT.unitTestFFGC()
end

% +smi_cluster

fprintf('smi_cluster.Clustering.unitTest\n');
try
   smi_cluster.Clustering.unitTest()
end

fprintf('smi_cluster.PairCorrelation.unitTest\n');
try
   smi_cluster.PairCorrelation.unitTest()
end

fprintf('smi_cluster.StatisticsClustering.unitTest\n');
try
   smi_cluster.StatisticsClustering.unitTest()
end

% +smi_core

fprintf('smi_core.ChannelRegistration.unitTest\n');
try
   smi_core.ChannelRegistration.unitTest()
end

fprintf( ...
   'smi_core.DataToPhotons.unitTest (convert raw data in ADUs to photons)\n');
try
   smi_core.DataToPhotons.unitTest()
end

fprintf('smi_core.DriftCorrection.unitTest\n');
try
   smi_core.DriftCorrection.unitTest()
end

%fprintf( ...
%   'smi_core.FRC.unitTest (Fourier Ring Correlation) [DIPimage needed]\n');
%try
%   smi_core.FRC.unitTest()
%end

fprintf('smi_core.FrameConnection.unitTest\n');
try
   smi_core.FrameConnection.unitTest()
end

fprintf( ...
   'smi_core.LocalizeData.unitTest (find localizations in an image stack\n');
try
   smi_core.LocalizeData.unitTest()
end

fprintf('smi_core.Threshold.unitTest\n');
try
   smi_core.Threshold.unitTest()
end

% +smi_psf

% smi_psf.PointSpreadFunction.crlbPSFPupil_unitTest()
% smi_psf.PointSpreadFunction.optimPSFZernike_unitTest()
% smi_psf.PointSpreadFunction.psfROIStack_unitTest()
% smi_psf.PointSpreadFunction.scalarPSFPrasadZone_unitTest()
% smi_psf.PointSpreadFunction.zernikeImage_unitTest()
fprintf('smi_psf.PointSpreadFunction.unitTest (several individual tests)\n');
try
   smi_psf.PointSpreadFunction.unitTest()
end

fprintf('smi_psf.Zernike.unitTest\n');
try
   smi_psf.Zernike.unitTest()
end

% +smi_sim

fprintf('smi_sim.GaussBlobs.unitTest\n');
try
   smi_sim.GaussBlobs.unitTest()
end

fprintf('smi_sim.SimSMLM.unitTest\n');
try
   smi_sim.SimSMLM.unitTest()
end

% +smi_stat

fprintf('smi_stat.ChangeDetection.unitTest\n');
try
   smi_stat.ChangeDetection.unitTest()
end

fprintf('smi_stat.DiffusionEstimator.unitTest\n');
try
   smi_stat.DiffusionEstimator.unitTest()
end

fprintf('smi_stat.HMM.unitTest\n');
try
   smi_stat.HMM.unitTest()
end

% +smi_vis

fprintf('smi_vis.GenerateImages.blobColorOverlay_unitTest\n');
try
   smi_vis.GenerateImages.blobColorOverlay_unitTest()
end

fprintf('smi_vis.GenerateImages.circleImage_unitTest\n');
try
   smi_vis.GenerateImages.circleImage_unitTest()
end

fprintf('smi_vis.GenerateImages.colorImage_unitTest\n');
try
   smi_vis.GenerateImages.colorImage_unitTest()
end

fprintf('smi_vis.GenerateImages.driftImage_unitTest\n');
try
   smi_vis.GenerateImages.driftImage_unitTest()
end

fprintf('smi_vis.GenerateImages.gaussianImage_unitTest\n');
try
   smi_vis.GenerateImages.gaussianImage_unitTest()
end

fprintf('smi_vis.GenerateImages.histogramImage_unitTest\n');
try
   smi_vis.GenerateImages.histogramImage_unitTest()
end
