### +smi_core.DataToPhotons

This class contains static methods associated with the gain and
offset correction needed to convert raw data from the camera (arrays
given in Analog to Digital Units (ADU)) to units of photons.  The
main usage of this class is shown in the EXAMPLE USAGE section below.

EXAMPLE USAGE:
- Given RawData in units of ADU, and an SMF structure with the fields
  SMF.Data.CameraGain, SMF.Data.CameraOffset, and
  SMF.Data.CameraReadNoise (see SingleMoleculeFitting class for
  details), you can convert RawData and CameraReadNoise to units of
  photons and photons^2, respectively, as follows:
```
      [~, RawDataConverted, CameraReadNoiseConverted] = ...
          smi_core.DataToPhotons(SMF, ...
          RawData, RawDataROI, CalibrationROI, true);
```
  Alternatively, you can prepare the class for usage (and set class
  parameters immediately) as follows:
```
      DTP = smi_core.DataToPhotons(SMF, ...
          RawData, RawDataROI, CalibrationROI);
      [Data, ReadNoise] = DTP.convertData();
```

REQUIRES:
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
    
properties:
```
   % obj.RawData converted to units of photons (float array)
   CorrectedData {mustBeNumeric(CorrectedData)}
   
   % obj.CameraReadNoise converted to units of photons (float array)
   CorrectedReadNoise {mustBeNumeric(CorrectedReadNoise)}
   
   % Data that is to be gain/offset corrected (float array)
   RawData {mustBeNumeric(RawData)}

   % Region of interest of the raw data (float array)
   % (see obj.convertToPhotons() for details/usage)
   RawDataROI {mustBeNumeric(RawDataROI)}
   
   % Gain of the camera used to collect RawData (float array)(ADU/e-)
   CameraGain {mustBeNumeric(CameraGain)}
   
   % Offset of the camera used to collect RawData (float array)(ADU)
   CameraOffset {mustBeNumeric(CameraOffset)}
   
   % Read noise of the camera used to collect Raw Data (float array)(ADU^2)
   CameraReadNoise {mustBeNumeric(CameraReadNoise)}
   
   % Region of interest of the gain/offset arrays (float array)
   % (see obj.convertToPhotons() for details/usage)
   CalibrationROI {mustBeNumeric(CalibrationROI)}
   
   % Verbosity level for standard workflow. (Default = 1)
   %   0: Command Window updates will be supressed where possible and
   %      reasonable.
   %   1: Some updates may appear in Command Window
   %   2: More detailed updates in Command Window
   %   3: Lot's of info. may be passed to Command Window. This mode
   %      may be useful for debugging large workflows encompassing
   %      this class.
   Verbose = 1;
```
methods:
- **convertData**: performs gain/offset correction on data
- **convertToPhotons**: converts RawData to units of photons
- **unitTest**: checks vital functionality of the DataToPhotons class
