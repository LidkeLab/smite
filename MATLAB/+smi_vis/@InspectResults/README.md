### +smi_vis/@InspectResults

InspectResults contains methods useful for inspecting SR data.
This class provides an interface for inspecting super-resolution
results.  The primary intention is that this class can be used to
associate localizations from a reconstruction image (e.g., a Gaussian
SR image) and the entries of a Single Molecule Data structure
(see smi_core.SingleMoleculeData).

The suggested usage of this class is through the GUI, which will be
opened by default when creating an instance of this class:
```
   Inspect = smi_vis.InspectResults();``
```
From the GUI, you should load a super-resolution image reconstruction
(saved as, e.g., a .png) as well as a *_Results.mat file containing
the SMD structure (see usage of smi.SMLM, which generates a 
*_Results.mat file).  Once the image is loaded, it'll be displayed in
the GUI figure.  Hovering over the GUI figure will reveal a toolbar
with a '+' icon which can be used to highlight a ROI and display
info about localizations from the ROI.

REQUIRES:
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

```
properties
   % Single Molecule Data structure. (see smi_core.SingleMoleculeData)
   SMD

   % Super-resolution reconstruction image. (YxXx3 float)
   SRImage

   % Figure containing the interactive GUI.
   GUIFigure

   % Axes containing the SR image.
   ImageAxes

   % ROI in SMD coordinates. ([YStart, XStart, YEnd, XEnd])
   ROI
```
