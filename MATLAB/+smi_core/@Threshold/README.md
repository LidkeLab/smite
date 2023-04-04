### +smi_core/Threshold

Threshold localizations based on various properties of the localizations.
This is done by creating a *ThreshFlag* field for the SMD structure using the
method **setThreshFlag** and the bounds in the SMF.MinMax (only those bounds
defined will be applied).  **applyThresh** then applies the *ThreshFlag* to the
localizations.  **rejectedLocalizations** will visually display what
localizations were eliminated and for what reasons.  **translateThreshFlag**
can also be used to produce human readable interpretations of the
**ThreshFlag**.

methods:
- **applyThresh**:
  applies ThreshFlag to perform thresholding on SMD,
  where SMD can be any appropriate smite coordinate containing structure.
- **rejectedLocalizations**:
  Produce plots of accepted and rejected localization fits,
  individually by reason rejected and combined by number of reasons rejected
  or by major reason rejected
- **setMinMax**:
  creates the MinMax structure used throughout the Threshold class
  from a provided SMF structure.  The MinMax fields contain a minimum and a
  maximum value on which to threshold
- **setThreshFlag**:
  Creates ThreshFlag field for SMD, the same size as SMD.X
- **translateThreshFlag**:
  is a static "wrapper" for translateThreshFlagNS
- **translateThreshFlagNS**:
  translates ThreshFlag to a human readable string
- **unitTest**:
  unitTest for setThreshFlag() in smi_core.Threshold class
