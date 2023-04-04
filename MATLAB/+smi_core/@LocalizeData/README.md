### +smi_core/@LocalizeData

This class contains method(s) to generate localizations from numeric
arrays of raw data in the form of images/stacks of images.

NOTE: All default class properties are set in the constructor.  These
      properties are extracted/constructed from either in constructor
      input SMF or from a default SMF created as
      ```SMF = smi_core.SingleMoleculeFitting.createSMF();```

EXAMPLE USAGE:
```
  LD = smi_core.LocalizeData(RawData, SMF);
  [SMD, SMDPreThresh] = LD.genLocalizations();
  or
  [~, SMD, SMDPreThresh] = smi_core.LocalizeData(RawData, SMF, 1);
```
REQUIRES:
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

---

```
properties:
  SMF             % see smi_core.SingleMoleculeFitting
  ScaledData      % (float array)(Photons) Gain/offset corrected data
  Verbose = 1     % verbosity level
  ResultsDir = [] % directory to save the color overlay if defined
```

---

methods:
- **[colorOverlay](colorOverlay.m)**:
  Displays an overlay of the (green) model with (red) data
- **[genLocalizations](genLocalizations.m)**:
  generates localizations from scaled data
- **[unitTest](unitTest.m)**:
  tests vital functionality of the LocalizeData class
