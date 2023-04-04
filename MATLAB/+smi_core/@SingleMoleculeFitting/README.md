### +smi_core/@SingleMoleculeFitting

SingleMoleculeFitting: A class defining the Single Molecule Fitting structure

The SMF structure is a structure of structures that collectively contain
all parameters required to go from raw data to an SMD results structure.
The SMF structure is an input of many smi methods. It
intended to be extensible to enable new analysis tools and methods.
The SMF class implements tools for working with SMF structures,
but the data structure itself is not an object of the class.

Parameters of sub-structures are explained in more detail in
the classes and methods that use them.  An incomplete list of classes
that use each sub-structure is listed in {}.

The SMF structure has these
[sub-structures and fields](../../../doc/DataStructures/SMF.md).

---

methods:
- **[gui](gui.m)**:
  the GUI method for the SingleMoleculeFitting class
- **[padSMF](padSMF.m)**:
  adds fields in SMFPadding to SMF that weren't already present
