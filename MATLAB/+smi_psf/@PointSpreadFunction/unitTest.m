function unitTest()
%unitTest tests various functionality.

   [CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil_unitTest();
   [Report]=smi_psf.PointSpreadFunction.optimPSFZernike_unitTest();
   [Report]=smi_psf.PointSpreadFunction.oversamplePSFPupil_unitTest();
   [PSFStruct]=smi_psf.PointSpreadFunction.phaseRetrieval_unitTest();
   [Report]=smi_psf.PointSpreadFunction.scalarPSFPrasadZone_unitTest();
   smi_psf.PointSpreadFunction.psfROIStack_unitTest();
   [Report]=smi_psf.PointSpreadFunction.zernikeImage_unitTest();

end
