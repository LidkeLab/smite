### +smi_core/@LoadData

LoadData has methods to load raw microsope dataset from .mat or
.h5 format files.

INPUTS:
- SMF - Single Molecule Fitting structure. The name and directory
  location of the data file is obtained from fields
  (SMF.Data.FileName and SMF.Data.FileDir) of SMF. If the fields are
  non-existent or are empty, data file selection will be required
  with a pop-up window.
- varargin - input arguments helping to locate the 'data' to be
  loaded from the data file. This varies for the different file
  extensions (.mat and .h5). For details please see
  documentation for the methods.

OUTPUTS:
- Data - data loaded from FileName, converted to type single
- SMF - SMF structure modified by adding fields FileName and FileDir,
  if not present.

---

methods:
- **[countNDatasets](countNDatasets.m)**:
  counts the number of datasets in SMF.Data.FileName .
- **[readH5File](readH5File.m)**:
  contents of an h5 file into H5Structure
- **[setSMFDatasetList](setSMFDatasetList.m)**:
  sets SMF.Data.DatasetList from DatasetMods and NDatasets
