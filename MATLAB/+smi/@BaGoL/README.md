### +smi/@BaGoL

BaGoL Implements a Bayesian Grouping of Localizations (BaGoL)

Single molecule localization based super-resolution data can contain
repeat localizations from the same emitter. These localizations can be
grouped to give better localization precision.  BaGoL explores the
possible number of emitters and their positions that can explain the
observed localizations and uncertainties (data). An 'emitter' is a
blinking/binding point source that generates a single 'localization'
each time there is a blinking/binding event. Localizations must be
frame connected.

The core algorithm uses Reversible Jump Markov Chain Monte Carlo to
add, remove and move emitters and to explore the allocation of
localizations to emitters. The localization precisions are assumed to be
accurate. The prior distribution is parameterized by either a Poisson 
or Gamma distribution function which can be either given as an input or 
learned via the hierarchical Bayes approach.

The primary BaGoL outputs are a 'Posterior Image' and MAPN coordinates
that can also be used to generate a 'MAPN image'. The Posterior Image
shows a probability distribution of emitter locations that is a weighted
average over the number of emitters, emitter locations, and
allocation of localizations to emitters. The Maximum a Posteriori
Number of emitters (MAPN) result uses only the information from
the most likely number of emitters within the chain.  MAPN emitter 
coordinates and their uncertainties are returned and images from these 
coordinates are generated similarly to a traditional SR reconstruction.

This class implements the pre/post-processing steps needed to analyze
a spatially extended data set typical of super-resolution experiments.
This includes breaking data into subregions and collating the results
from the subregions.

REQUIRES:
- MATLAB 2016 or higher versions.
- Statistics and Machine learning toolbox.

CITATION: "High-Precision Estimation of Emitter Positions using Bayesian
          Grouping of Localizations", Mohamadreza Fazel, Michael J. Wester,
          David J. Schodt, Sebastian Restrepo Cruz, Sebastian Strauss,
          Florian Schueder, Thomas Schlichthaerle, Jennifer M. Gillette,
          Diane S. Lidke, Bernd Rieger, Ralf Jungmann and Keith A. Lidke,
          Nature Communications, **13**(7152), November 22, 2022, 1--11,
          (DOI: 10.1038/s41467-022-34894-2).

---

USAGE:
```
  B=BaGoL()       % create object
  B.SMD=....      % set properties
  B.analyze_all() % run complete analysis
```
The class also has several methods for visualizing and saving results.
See 'doc BaGoL' for the complete list of methods.

```
PROPERTIES:
  SMD:        A structure containing the fields:
      X:          Vector of X localization positions (Pixel)(Nx1)
      Y:          Vector of Y localization positions (Pixel)(Nx1)
      Z:          Vector of Z localization positions (Pixel)(Nx1)(Optional)
      X_SE:       Vector of X localization standard error (Pixel)(Nx1)
      Y_SE:       Vector of Y localization standard error (Pixel)(Nx1)
      Z_SE:       Vector of Z localization standard error (Pixel)(Nx1)(Optional)
      FrameNum:   Vector of localization frame numbers (Nx1)
      PixelSize:  Camera pixel size (nm/Pixel) 
  Xi:         Loc./emitter [lambda] (Poisson) or [k theta] (Gamma).
              When learning Xi this is used to initialize a chain. If
              two initial values are given it uses a gamma and otherwise
              a poisson prior.
  Alpha_Xi:   Shape parameter of Xi gamma hyper prior
  Alpha_Xi:   Scale parameter of Xi gamma hyper prio
  ROIsize:    ROI size for RJMCMC (nm) (Default=200)
  Overlap:    Allowed overlap between subregions (nm)(Default=50)
  Drift:      Expected magnitude of drift (nm/frame)(Default=0)
  SE_Adjust:  Adjustement of localization precisions (nm) (Default=0)
  N_Burnin:   Number of samples in burn in of RJMCMC chain (Default=2000)
  N_Trials:   Number of samples in RJMCMC chain post burn in (Default=3000)
  P_Jumps:    Proposal probabilities for RJMCMC Jumps
              [Move, Allocate, Add, Remove]
              sum(P_Jumps) must equal 1.
              (Default = [0.25, 0.25, 0.25, 0.25])
  PixelSize:  The pixel size for output posterior images (nm) (Default=1)
  PImageFlag: Generate Posterior image. 0 or 1. (Default=0)
  HierarchFlag: Use hierarchical Bayse to learn Xi. 0 or 1. (Default=0)
  NSamples:   Number of RJMCMC samples before sampling Xi (Default=10)
  PImageSize: Size of the output posterior images (nm)
  ChainFlag:  Save RJMCMC chain. 0 or 1. (Default=0)
  XStart:     X starting coordinate of output posterior images (nm)
  YStart:     Y starting coordinate of output posterior images (nm)

  ClusterSMD: An array of SMD structures for each cluster
  MAPN:       SMD output structure containing fields:
      X:      Vector of X emitter positions (nm)(Kx1)
      Y:      Vector of Y emitter positions (nm)(Kx1)
      Z:      Vector of Z emitter positions (nm)(Kx1)
      X_SE:   Vector of X localization standard error (nm)(Kx1)
      Y_SE:   Vector of Y localization standard error (nm)(Kx1)
      Z_SE:   Vector of Z localization standard error (nm)(Kx1)
      Nmean:  Mean number of localizations per emitter (Kx1)
  PImage:     Posterior Image
  Chain:      Cell array of RJMCMC chains for pre-clusters (See BaGoL_RJMCMC)
  XiChain:    Chain of Xi (locs per emitters dist. either Poisson or gamma)
  SaveName:   Final results are saved under this name (Default: BaGoL)
```

---

methods:
- **[analyze_all](BaGoL.m)** Implements complete BaGoL analysis of SR dataset

- **[BaGoL_RJMCMC](BaGoL_RJMCMC.m)**:
  This the core BaGoL algorithm. It uses Reversible Jump Markov Chain 
  Monte Carlo to add, remove and move emitters, and to explore the 
  classification of localizations to emitters.
- **[BaGoL_RJMCMC_Hierarchical](BaGoL_RJMCMC_Hierarchical.m)**:
  BaGoL's core RJMCMC algorithm that takes
  NSamples samples from the posterior using RJMCMC and return the last 
  sample to be used in the hierarchical Bayes approach
- **[SEfilter](SEfilter.m)**:
  filters out NaNs and localizations with zero precisions
- **[align_template](align_template.m)**:
  Aligns a set of points to a template
- **[assignROIs](assignROIs.m)**:
  applies SE_Adjust and filters zero precisions or nans in the data
- **[dispIm](dispIm.m)**:
  GUI to display BaGoL output images with some useful tools
- **[errPlot](errPlot.m)**:
  Plots the input precisions, where circles represent the errors
- **[frameConnect](frameConnect.m)**:
  Connects coordinates from a blinking event across consecutive frames
- **[genCluster](genCluster.m)**:
  Simulates localization data from several types of emitter patterns
- **[genMAPN](genMAPN.m)**:
  Generates the MAPN coordinates from the most repeated model in the chain
- **[genMAPNIm](genMAPNIm.m)**:
  Produces a Gaussian blob image from either SMD or MAPN
- **[genPosterior](genPosterior.m)**:
  Updates Posterior Image using RJMCMC chain for a subregion (ROI)
- **[genROIs](genROIs.m)**:
  takes the input coordinates and splits them into smaller regions
- **[genSRMAPNOverlay](genSRMAPNOverlay.m)**:
  generates a multicolor overlay containing circles with radii
  proportional to the average localization precision
- **[hierBaGoL_analysis](hierBaGoL_analysis.m)**:
  This function is adapted from EGFR_dSTORM.m in the BaGoL distribution
- **[hierBaGoL_run](hierBaGoL_run.m)**:
  runs one or more BaGoL analyses  (see
  [MATLAB/examples/hierBaGoL_wrapper.m](../../examples/hierBaGoL_wrapper.m)).
- **[importLLSMD](importLLSMD.m)**:
  Imports a Lidke Lab SMD and converts to BaGoL format
- **[initPostIm](initPostIm.m)**:
  Find the Posterior Image size and position
- **[loadPICASSOh5](loadPICASSOh5.m)**:
  Loads a dataset saved in H5-format by the PICASSO software
- **[makeIm](makeIm.m)**:
  Produces a Gaussian blob image from the list of input coordinates
- **[plotMAPN](plotMAPN.m)**:
  Plots the MAPN positions together with input localizations
- **[plotNND](plotNND.m)**:
  Makes and saves the NND histogram of MAPN coordinates
- **[plotNND_PDF](plotNND_PDF.m)**:
  makes and saves the NND PDF histogram of MAPN coordinates
- **[precluster](precluster.m)**:
  Finds the independent clusters for analysis with RJMCMC
- **[removeOverlap](removeOverlap.m)**:
  Finds coordinates located in overlapping regions
- **[sampleGam](sampleGam.m)**:
  infers parameters of gamma prior on number of localizations per emitter
- **[samplePoiss](samplePoiss.m)**:
  infers parameter of Poisson prior on number of localizations per emitter
- **[saveBaGoL](saveBaGoL.m)**:
  Saves plots of NND, precisions, Xi, SR, MAPN, Posterior, Overlay images
- **[saveMAPN](saveMAPN.m)**:
  Saves the MAPN coordinates
- **[scaleIm](scaleIm.m)**:
  Scales and clips the image intensity to improve image contrast
- **[scalebar](scalebar.m)**:
  Add scalebar to image
