### +smi_sim/@SimSMLM
    
SimSMLM is a Single Molecule Localization Microscopy data generating
class. This class can produce Siemen's star shaped data or k-tets
with particles
distributed uniformly throughout the frames. The blinking events are
produced based on the given rate parameters (K_OnToOff, K_OffToOn and
K_OnToBleach). The data has the same uniform background noise for the
whole sequence. The output is corrupted with Poisson noise.
```    
Typical data flows are
   produce noisy coordinates:
      SMD_True -> SMD_Labeled -> SMD_Model -> SMD_Data
   produce noisy image stacks:
      SMD_True -> SMD_Labeled -> SMD_Model -> Model -> Data
where
  SMD_True      true locations of localizations
  SMD_Labeled   obj.LabelingEfficiency applied to SMD_True localizations,
                removing localizations that are not labeled
  SMD_Model     blinks generated for SMD_Labeled localizations
  SMD_Data      SMD_Model with positional and intensity noise added
  Model         Gaussian blob image stack produced from SMD_Model
  Data          Model image stack to which Poisson noise has been applied
```
Model and Data are image stacks (n x n x f), where n is the linear size
of the image in pixels and f is the total number of frames to be
generated (f = obj.NDatasets * obj.NFrames).

EITHER, generate an SMD structure with positional and intensity noise.
(SMD_Data <- genNoisySMD) OR ALTERNATIVELY, generate the blobs without
Poisson noise (Model) and then add it in (Data) [see genImageStack].

---

Examples using this class can be found in MATLAB/examples:
- Example_Clustering.m
- Example_StatisticsClustering.m
- Example_simSMLM.m
For example,
```
   %% Example of generating sythetic SMLM data

   %Create sim object

   S=smi_sim.SimSMLM()
   NWings=20
   S.NDatasets=20
   S.SZ = 64;
   S.simStar(NWings)

   % Generate Images 
   [Model,Data]=S.genImageStack();

   % Generate Noisy Coordinates
   [SMD_Noisy]=S.genNoisySMD()
   figure;scatter(SMD_Noisy.X,SMD_Noisy.Y)
```

---

```
properties
   SZ=256            % Linear size of image (pixels)
   Rho=30            % Fluorophore Density (fluorophore/pixel)
   NDatasets=1       % Number of datasets
   NFrames=1000      % Number of frames per dataset
   ZoomFactor=20     % It can be either smaller or larger than one- CHANGE
   K_OnToOff=1       % Fluorophore turns Off from On state
                     %    (default:1 frames^-1)
   K_OffToOn=0.0005  % Fluorophore return to On state from Off state
                     %    (default:0.0005 frames^-1)
   K_OnToBleach=0.2  % Fluorophore bleached out (default:1/5 frames^-1)
   EmissionRate=1000 % Emission rate (Intensity) of photons (photons/frame)
   Bg=5              % Background Count Rate (counts/pixel)
   PSFSigma=1.3      % Point Spread Function Sigma size (Pixels)
   LabelingEfficiency=1 % Fluorophore labeling efficiency [range: 0 - 1]
   NOnEvents         % Number of on events per localization
   % A string which determine if the particle starts on or starts randomly
   % on or off.  It can be either 'on' or 'Equib'.
   StartState='Equib'
   % SparseFlag, if true, turns on using sparse matrices in genBlinks,
   % allowing it to accomodate larger examples.  However, sparse matrices
   % are also slower to manipulate than regular ones especially for
   % smaller examples, so shouldn't be used unless needed, thus the
   % default is false.
   SparseFlag=false
   Verbose = 1      % Verbosity level

   % Generic note: SMD_* below are SMD structures with various fields
   % filled in as appropriate at that stage.
   SMD_True         % True coordinates produced by sim*
   SMD_Labeled      % True labeled coordinates produced by applyLabelEffic
   SMD_Model        % Coordinates with blinks produced by genBlinks

   % Note SMD_Model.ConnectID: SMD_Labeled -> SMD_Model (see genBlinks)
   % ConnectID produces an indexing array that associates each
   % localization in SMD_Model with the localization that produced it in
   % SMD_Labeled.  That is,
   %
   %    obj.SMD_Labeled.X(obj.SMD_Model.ConnectID) == obj.SMD_Model.X
```

---

methods:
- **[applyLabelEffic](applyLabelEffic.m)**:
  applies labeling efficiency to an existing set of fluorophores
- **[genBlinks](genBlinks.m)**:
  generates the blinking time traces for a single particle over the given
  number of the frames considering the input parameters and returns the
  SMD_Model structure
- **[genImageStack](genImageStack.m)**:
  generates image stacks without and with Poisson noise
- **[genModel](genModel.m)**:
  generates SMD_Labeled and SMD_Model from SMD_True
- **[genNoisyData](genNoisyData.m)**:
  generates a noisy image
- **[genNoisySMD](genNoisySMD.m)**:
  takes the SMD_Model struct as an input, generates the noisy SMD struct and
  return SMD_Data
- **[kDisc](kDisc.m)**:
  produces a disc of k points
- **[kTet](kTet.m)**:
  produces a circle of k equally spaced points starting at a random place
- **[simkDiscs](simkDiscs.m)**:
  generates 2D k-discs in the simulation region
- **[simStar](simStar.m)**:
  simulates the Siemen star
- **[simkTets](simkTets.m)**:
  generates 2D k-tets in the simulation region

unit test:
- **[unitTest](unitTest.m)**:
  tests the functionality of this class
