### +smi_core/@SingleMoleculeData

iSingleMoleculeData: A class defining the Single Molecule Data structure

This datatype is one of the primary results structures in the ***smite***
environment. The SMD structure is an input and output of many smi
methods. It intended to be extensible.
The SMD class implements tools for working with SMD structures,
but the data structure itself is not an object of the class.

The structure has these [properties](../../../doc/DataStructures/SMD.md).

---

methods:
- **[catSMD](catSMD.m)**:
  concatenates two SMD structures into one
- **[computeDensity](computeDensity.m)**:
  estimates the per frame density of observed emitters
- **[computeDensityImage](computeDensityImage.m)**:
  estimates the local density of observed emitters
- **[defineSMDMask](defineSMDMask.m)**:
  defines a boolean for masking SMD localizaions.
- **[extractDatasets](extractDatasets.m)**:
  extracts the 'Datasets' from 'SMD'
- **[isolateSubROI](isolateSubROI.m)**:
  isolates localizations within ROI
- **[isolateSubSMD](isolateSubSMD.m)**:
  isolates a subset of SMD defined by SubIndices
- **[maskSMD](maskSMD.m)**:
  masks the input SMD based on the image Mask
