### +smi_helpers

+smi_helpers is the namespace for helper functions and classes of ***smite***:
- (@Filters)[@Filters/README.md]:
  filters useful for BaGoL operating on SMDs
- (@ROITools)[@ROITools/README.md]:
  select ROIs from an image; save in a structure

Functions:
- addBasicGUI:
  adds simple GUI controls based on the fields in 'ParamStruct'
- arrayMUX:
  a multiplexer intended to generalize for most arrays
- compressToRange:
  compresses the 'IntegerArray' to a range of integers with no missing values
- convertTimeStringToNum:
  converts a time string to a number
- filenameWindows2Unix:
  converts the Windows filenameWindows to the Unix filename filenameUnix,
  optionally preceded by prefix
- findStartEndInds:
  finds the start and end indices of events in BoolArray
- gatherFullPathnames:
  combines or searches for directory paths and filenames matching a pattern
- genTimeString:
  creates a char array with the current time
- getDirectoryNames:
  generates directory names matching a NameString
- getFileNames:
  creates a list of filenames in FileDir
- nMODm:
  modulus such that r is in [1, m] rather than [0, m - 1]
- padStruct:
  pads the input 'Struct' with defaults in 'DefaultStruct'
- pairText:
  creates paired lists of text from two sets of text lists
- pairTimeStrings:
  selects indices from the input 'TimeStringOptions' that are closest in time
  to the timestrings provided in 'TimeStrings'
- removeBorder:
  removes border pixels from the input image
- requiredToolboxes:
  print out required toolboxes for each directory in the SMITE directory
  structure
- selectFiles:
  select a list of files
- subdivideImage:
  divvys up an image into sub-ROIs of the image
- subdivideSMD:
  divvys up an SMD into sub-ROIs
- triangle_threshold:
  use the triangle method to find a threshold
- writeMPEG4:
  writes a .mp4 image constructed from the RGB image input on all OS
