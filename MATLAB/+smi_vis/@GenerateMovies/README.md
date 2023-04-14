### +smi_vis/@GenerateMovies

GenerateMovies contains methods useful for generating movies,
particularly for single-particle tracking data.

REQUIRES:
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

---
    
```
properties:
   % Structure of parameters enforced in the generated movie.
   % See prepDefaults() for a description of the parameter options and
   % playMovie() for usage.
   Params struct = struct([])

   % Raw data displayed under trajectories. (YSizexXSizex3xNFrames)
   RawData

   % Single Molecule Fitting structure with pixel size and framerate.
   SMF = smi_core.SingleMoleculeFitting;

   % Optional SMD containing points to mark in the movie.
   SMD = smi_core.SingleMoleculeData.createSMD();

   % Tracking Results structure for the trajectories.
   TR = smi_core.TrackingResults.createTR();
```

---

methods:
- **[addAxesTicks](addAxesTicks.m)**:
  prepares/displays tick labels on the axes.
- **[addTimeStamp](addTimeStamp.m)**:
  places a time stamp in a corner of the PlotAxes
- **[cropRawData](cropRawData.m)**:
  crops the input data to the region specified in Params
- **[defineCropROI](defineCropROI.m)**:
  defines the ROI surrounding the data for an auto-cropping
- **[generateMovie](generateMovie.m)**:
  generates a movie of 2D trajectories over time
- **[gui](gui.m)**:
  prepares a movie generation GUI for the GenerateMovies class
- **[makeFrame](makeFrame.m)**:
  plots a single frame for the trajectory movie
- **[playMovie](playMovie.m)**:
  prepares a movie of the trajectories in TR
- **[plotTrajectories](plotTrajectories.m)**:
  plots trajectories in the specified axes
- **[prepAxes](prepAxes.m)**:
  prepAxes prepares some features (e.g., axis limits) of the movie axes.
- **[prepDefaults](prepDefaults.m)**:
  prepares a default 'Params' structure
- **[prepRawData](prepRawData.m)**:
  crops and rescales the movie data for display purposes.
- **[saveMovie](saveMovie.m)**:
  generates and saves a movie of 2D trajectories
- **[saveRawDataMovie](saveRawDataMovie.m)**:
  generates and saves a movie of raw data
- **[setVitalParams](setVitalParams.m)**:
  ensures vital parameters in obj.Params are set (i.e.,
  values that'll let the user create a movie without crashing the code)
