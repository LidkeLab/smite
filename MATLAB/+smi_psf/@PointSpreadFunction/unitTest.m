function unitTest()
%unitTest tests various functionality.

   [CRLB]=crlbPSFPupil_unitTest();
   [Report]=optimPSFZernike_unitTest();
   [Report]=oversamplePSFPupil_unitTest();
   [PSFStruct]=phaseRetrieval_unitTest();
   [Report]=scalarPSFPrasadZone_unitTest();
   psfROIStack_unitTest();
   [Report]=zernikeImage_unitTest();

end
