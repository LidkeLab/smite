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

```
USAGE:
  B=BaGoL()       % create object
  B.SMD=....      % set properties
  B.analyze_all() % run complete analysis

The class also has several methods for visualizing and saving results.
See 'doc BaGoL' for the complete list of methods.

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

REQUIRES:
  MATLAB 2016 or higher versions.
  Statistics and Machine learning toolbox.
```

CITATION: "High-Precision Estimation of Emitter Positions using Bayesian
          Grouping of Localizations", Mohamadreza Fazel, Michael J. Wester,
          David J. Schodt, Sebastian Restrepo Cruz, Sebastian Strauss,
          Florian Schueder, Thomas Schlichthaerle, Jennifer M. Gillette,
          Diane S. Lidke, Bernd Rieger, Ralf Jungmann and Keith A. Lidke,
          Nature Communications, **13**(7152), November 22, 2022, 1--11,
          (DOI: 10.1038/s41467-022-34894-2).
