### +smi_vis/@GenerateMovies

GenerateMovies contains methods useful for generating movies,
particularly for single-particle tracking data.

REQUIRES:
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
    
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

Methods:
- **addAxesTicks**: prepares/displays tick labels on the axes.
- **addTimeStamp**: places a time stamp in a corner of the PlotAxes
- **cropRawData**: crops the input data to the region specified in Params
- **defineCropROI**: defines the ROI surrounding the data for an auto-cropping
- **generateMovie**: generates a movie of 2D trajectories over time
- **gui**: prepares a movie generation GUI for the GenerateMovies class
- **makeFrame**: plots a single frame for the trajectory movie
- **playMovie**: prepares a movie of the trajectories in TR
- **plotTrajectories**: plots trajectories in the specified axes
- **prepAxes**: prepAxes prepares some features (e.g., axis limits) of the
  movie axes.
- **prepDefaults**: prepares a default 'Params' structure
- **prepRawData**: crops and rescales the movie data for display purposes.
- **saveMovie**: generates and saves a movie of 2D trajectories
- **saveRawDataMovie**: generates and saves a movie of raw data
- **setVitalParams**: ensures vital parameters in obj.Params are set (i.e.,
  values that'll let the user create a movie without crashing the code)
