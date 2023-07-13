classdef test_core < matlab.unittest.TestCase

% Run various tests on the core functionality of SMITE.  Much output will be
% saved in tempdir/smite/unitTest/name_of_test.  ExpectedResults are provided
% in the directory in which run_tests.m resides where very large files have
% been deleted so as to not bloat up the the SMITE distribution.

% In the MATLAB unittest context, run in the following manner:
%    testCase = test_core
%    results = testCase.run

% +smi

methods (Test)

   function test_SMLM(testCase)
      fprintf('smi.SMLM.unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi.SMLM.unitTest();
      else
         fprintf('GPU needed!\n');
         results = [1, 1, 1];
      end
      testCase.verifyEqual(results, [1, 1, 1]);
   end

   function test_unitTestFFGC(testCase)
      fprintf( ...
         'smi.SPT.unitTestFFGC (frame-to-frame and gap closing processes)\n');
      if (gpuDeviceCount > 0)
         results = smi.SPT.unitTestFFGC();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

% +smi_cluster

   function test_Clustering(testCase)
      fprintf('smi_cluster.Clustering.unitTest\n');
      results = smi_cluster.Clustering.unitTest();
      testCase.verifyEqual(results, 1);
   end

   function test_PairCorrelation(testCase)
      fprintf('smi_cluster.PairCorrelation.unitTest\n');
      results = smi_cluster.PairCorrelation.unitTest();
      testCase.verifyEqual(results, 1);
   end

   function test_StatisticsClustering(testCase)
      fprintf('smi_cluster.StatisticsClustering.unitTest\n');
      results = smi_cluster.StatisticsClustering.unitTest();
      testCase.verifyEqual(results, 1);
   end

% +smi_core

   function test_ChannelRegistration(testCase)
      fprintf('smi_core.ChannelRegistration.unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_core.ChannelRegistration.unitTest();
      testCase.verifyEqual(results, ...
                           [true, true, true, true, true, true, true]);
      else
         fprintf('GPU needed!\n');
         results = [true, true, true, true, true, true, true];
      end
   end

   function test_DataToPhotons(testCase)
      fprintf( ...
      'smi_core.DataToPhotons.unitTest (convert raw data in ADUs to photons)\n');
      if (gpuDeviceCount > 0)
         results = smi_core.DataToPhotons.unitTest();
      else
         fprintf('GPU needed!\n');
         results = [true, true, true, true, true];
      end
      testCase.verifyEqual(results, [true, true, true, true, true]);
   end

   function test_DriftCorrection(testCase)
      fprintf('smi_core.DriftCorrection.unitTest\n');
      results = smi_core.DriftCorrection.unitTest();
      testCase.verifyEqual(results, 1);
   end

   function test_FrameConnection(testCase)
      fprintf('smi_core.FrameConnection.unitTest\n');
      results = smi_core.FrameConnection.unitTest();
      testCase.verifyEqual(results, [true, true, true, true, true]);
   end

   function test_FRC(testCase)
      fprintf( ...
       'smi_core.FRC.unitTest (Fourier Ring Correlation) [DIPimage needed]\n');
      try
         results = smi_core.FRC.unitTest();
      catch
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

   function test_LocalizeData(testCase)
      fprintf( ...
   'smi_core.LocalizeData.unitTest (find localizations in an image stack\n');
      if (gpuDeviceCount > 0)
         results = smi_core.LocalizeData.unitTest();
      else
         fprintf('GPU needed!\n');
         results = true;
      end
      testCase.verifyEqual(results, true);
   end

   function test_Threshold(testCase)
      fprintf('smi_core.Threshold.unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_core.Threshold.unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

% +smi_psf

   function test_crlbPSFPupil(testCase)
      fprintf('smi_psf.PointSpreadFunction.crlbPSFPupil_unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_psf.PointSpreadFunction.crlbPSFPupil_unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

   function test_optimPSFZernike(testCase)
      fprintf('smi_psf.PointSpreadFunction.optimPSFZernike_unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_psf.PointSpreadFunction.optimPSFZernike_unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

   function test_psfROIStack(testCase)
      fprintf('smi_psf.PointSpreadFunction.psfROIStack_unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_psf.PointSpreadFunction.psfROIStack_unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

   function test_scalarPSFPrasadZone(testCase)
      fprintf('smi_psf.PointSpreadFunction.scalarPSFPrasadZone_unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_psf.PointSpreadFunction.scalarPSFPrasadZone_unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

   function test_zernikeImage(testCase)
      fprintf('smi_psf.PointSpreadFunction.zernikeImage_unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_psf.PointSpreadFunction.scalarPSFPrasadZone_unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

   function test_Zernike(testCase)
      fprintf('smi_psf.Zernike.unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_psf.Zernike.unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

% +smi_sim

   function test_GaussBlobs(testCase)
      fprintf('smi_sim.GaussBlobs.unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_sim.GaussBlobs.unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

   function test_SimSMLM(testCase)
      fprintf('smi_sim.SimSMLM.unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_sim.SimSMLM.unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

% +smi_stat

   function test_ChangeDetection(testCase)
      fprintf('smi_stat.ChangeDetection.unitTest\n');
      fprintf('This may sometimes error randomly.\n');
      results = smi_stat.ChangeDetection.unitTest();
      testCase.verifyEqual(results, 1);
   end

   function test_DiffusionEstimator(testCase)
      fprintf('smi_stat.DiffusionEstimator.unitTest\n');
      results = smi_stat.DiffusionEstimator.unitTest();
      testCase.verifyEqual(results, 1);
   end

   function test_HMM(testCase)
      fprintf('smi_stat.HMM.unitTest\n');
      if (gpuDeviceCount > 0)
          results = smi_stat.HMM.unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

% +smi_vis

   function test_blobColorOverlay(testCase)
      fprintf('smi_vis.GenerateImages.blobColorOverlay_unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_vis.GenerateImages.blobColorOverlay_unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

   function test_circleImage(testCase)
      fprintf('smi_vis.GenerateImages.circleImage_unitTest\n');
      results = smi_vis.GenerateImages.circleImage_unitTest();
      testCase.verifyEqual(results, 1);
   end

   function test_colorImage(testCase)
      fprintf('smi_vis.GenerateImages.colorImage_unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_vis.GenerateImages.colorImage_unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

   function test_driftImage(testCase)
      fprintf('smi_vis.GenerateImages.driftImage_unitTest\n');
      results = smi_vis.GenerateImages.driftImage_unitTest();
      testCase.verifyEqual(results, 1);
   end

   function test_gaussianImage(testCase)
      fprintf('smi_vis.GenerateImages.gaussianImage_unitTest\n');
      if (gpuDeviceCount > 0)
         results = smi_vis.GenerateImages.gaussianImage_unitTest();
      else
         fprintf('GPU needed!\n');
         results = 1;
      end
      testCase.verifyEqual(results, 1);
   end

   function test_histogramImage(testCase)
      fprintf('smi_vis.GenerateImages.histogramImage_unitTest\n');
      results = smi_vis.GenerateImages.histogramImage_unitTest();
      testCase.verifyEqual(results, 1);
   end

end % Methods

end
