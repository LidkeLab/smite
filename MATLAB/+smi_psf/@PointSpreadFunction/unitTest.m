function unitTest()
%unitTest tests various functionality.

   fprintf('crlbPSFPupil_unitTest ...\n');
   try
      [CRLB]=smi_psf.PointSpreadFunction.crlbPSFPupil_unitTest();
   end

   fprintf('optimPSFZernike_unitTest ...\n');
   try
      [Report]=smi_psf.PointSpreadFunction.optimPSFZernike_unitTest();
   end

   fprintf('oversamplePSFPupil_unitTest ...\n');
   try
      [Report]=smi_psf.PointSpreadFunction.oversamplePSFPupil_unitTest();
   end

   fprintf('phaseRetrieval_unitTest ...\n');
   try
      [PSFStruct]=smi_psf.PointSpreadFunction.phaseRetrieval_unitTest();
   end

   fprintf('scalarPSFPrasadZone_unitTest ...\n');
   try
      [Report]=smi_psf.PointSpreadFunction.scalarPSFPrasadZone_unitTest();
   end

   fprintf('psfROIStack_unitTest ...\n');
   try
      smi_psf.PointSpreadFunction.psfROIStack_unitTest();
   end

   fprintf('zernikeImage_unitTest ...\n');
   try
      [Report]=smi_psf.PointSpreadFunction.zernikeImage_unitTest();
   end

   fprintf('Done.\n');

end
