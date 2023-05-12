### +smi_helpers

+smi_helpers is the namespace for helper functions and classes of ***smite***:
- [@Filters](@Filters/README.md):
  filters useful for BaGoL operating on SMDs
- [@ROITools](@ROITools/README.md):
  select ROIs from an image; save in a structure

---

functions:
- **[addBasicGUI](addBasicGUI.m)**:
  adds simple GUI controls based on the fields in 'ParamStruct'
- **[arrayMUX](arrayMUX.m)**:
  a multiplexer intended to generalize for most arrays
- **[compressToRange](compressToRange.m)**:
  compresses the 'IntegerArray' to a range of integers with no missing values
- **[convertTimeStringToNum](convertTimeStringToNum.m)**:
  converts a time string to a number
- **[filenameWindows2Unix](filenameWindows2Unix.m)**:
  converts the Windows filenameWindows to the Unix filename filenameUnix,
  optionally preceded by prefix
- **[findStartEndInds](findStartEndInds.m)**:
  finds the start and end indices of events in BoolArray
- **[gatherFullPathnames](gatherFullPathnames.m)**:
  combines or searches for directory paths and filenames matching a pattern
- **[genTimeString](genTimeString.m)**:
  creates a char array with the current time
- **[getDirectoryNames](getDirectoryNames.m)**:
  generates directory names matching a NameString
- **[getFileNames](getFileNames.m)**:
  creates a list of filenames in FileDir
- **[mkSMITETmpDir](mkSMITETmpDir.m)**:
  Create full pathnames and temporary directories for SMITE testing/examples
- **[nMODm](nMODm.m)**:
  modulus such that r is in [1, m] rather than [0, m - 1]
- **[padStruct](padStruct.m)**:
  pads the input 'Struct' with defaults in 'DefaultStruct'
- **[pairText](pairText.m)**:
  creates paired lists of text from two sets of text lists
- **[pairTimeStrings](pairTimeStrings.m)**:
  selects indices from the input 'TimeStringOptions' that are closest in time
  to the timestrings provided in 'TimeStrings'
- **[removeBorder](removeBorder.m)**:
  removes border pixels from the input image
- **[requiredToolboxes](requiredToolboxes.m)**:
  print out required toolboxes for each directory in the SMITE directory
  structure
- **[selectFiles](selectFiles.m)**:
  select a list of files
- **[subdivideImage](subdivideImage.m)**:
  divvys up an image into sub-ROIs of the image
- **[subdivideSMD](subdivideSMD.m)**:
  divvys up an SMD into sub-ROIs
- **[triangle_threshold](triangle_threshold.m)**:
  use the triangle method to find a threshold
- **[writeMPEG4](writeMPEG4.m)**:
  writes a .mp4 image constructed from the RGB image input on all OS
