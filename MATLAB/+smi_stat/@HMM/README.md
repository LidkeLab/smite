### +smi_stat/@HMM

@HMM contains methods related to/useful for hidden Markov model
(HMM) of single-particle tracking (SPT).  In particular, this class
is designed to perform HMM and related analyses on the output
Tracking Results (TR) structures produced by smi.SPT.

NOTE: This class ALWAYS uses camera units (pixels, frames)
throughout the analysis.  When UnitFlag = 1, units are not
converted to physical units (micrometers, seconds) until the
very end when outputs/plots are being produced.

See reference:
Nitta, C. F., Green, E. W., Jhamba, E. D., Keth, J. M.,
Ortiz-Caraveo, I., Grattan, R. M., Schodt, D. J., Gibson, A. C.,
Rajput, A., Lidke, K. A., Steinkamp, M. P., Wilson, B. S.,
& Lidke, D. S. (2021). EGFR transactivates RON to drive oncogenic
crosstalk. eLife
https://doi.org/10.7554/eLife.63678

REQUIRES:
- Optimization Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox
    
properties:
```
   % Separation between fluorophores on a dimer (pixels)
   DimerSeparation(1, 1) {mustBeFloat(DimerSeparation)} = 0.5;
   
   % Typical domain size for the free, dimer, domain model (pixels)
   DomainSeparation(1, 1) {mustBeFloat(DomainSeparation)} = 2;
   
   % Max. separation for dimer candidates (pre-processing) (pixels)
   MaxSeparation(1, 1) {mustBeFloat(MaxSeparation)} = 5;

   % Registration error inflation (pixels)
   RegErrorInflation(1, 1) {mustBeFloat(RegErrorInflation)} = 0.0;
   
   % Identifier for one of the pre-built models. (Default = 'DF')
   % OPTIONS:
   %   'DF': dimer or free
   %   'DDF': dimer, domain, or free
   ModelSpecifier {mustBeText(ModelSpecifier)} = 'DF';
   
   % Handles to the state PDFs used in the HMM (cell array)
   PDFHandles(:, 1) cell
   
   % Initial guess of rate parameters (NRatesx1 float)
   RateParametersGuess(:, 1) {mustBeFloat(RateParametersGuess)} = ...
       0.01 * [1; 1];
   
   % Array of TR structures corresponding to dimer candidates. (Nx2)
   % This is organized as a NCandidatex2 structure, with each column
   % being a dimer candidate.
   TRArray(:, 2) struct
   
   % Data channel names added to certain outputs. (cell array of char)
   ChannelNames cell = {'Channel 1'; 'Channel 2'};
   
   % Model state names added to certain outputs. (cell array of char)
   StateNames cell = {'Dimer'; 'Free'};
   
   % Label for save directory to indicate a condition. (char/string)
   ConditionLabel {mustBeText(ConditionLabel)} = '';
   
   % Indicate results should be saved. (Default = true)
   % NOTE: This is used when running obj.performFullAnalysis().
   SaveResults logical = true;
   
   % Indicate outputs should be in physical units. (Default = false)
   UnitFlag logical = false;
   
   % Indicate movies should be generated and saved. (Default = false)
   % NOTE: This is used when running obj.performFullAnalysis().
   GenerateMovies logical = false;
   
   % Set of movie parameters (see GenerateMovies)
   MovieParams struct = smi_vis.GenerateMovies.prepDefaults();
   
   % Set of plot parameters (see DisplayParams in plotDimerPairInfo())
   PlotParams struct = struct();
   
   % Indicate plots should be generated and saved. (Default = true)
   % GeneratePlots(1) indicates the histogram of dimer durations
   %   should be generated.
   % GeneratePlots(2) indicates the summaray plots should be
   %   generated.
   % NOTE: This is used when running obj.performFullAnalysis() inside
   %       of obj.saveResults().
   GeneratePlots = [true; true];

   % Structure of parameters (see smi_core.SingleMoleculeFitting)
   % NOTE: This can be passed as an array of two SMFs, one for each
   %       channel, which can be used when information about each
   %       channel is needed (e.g., in obj.createAllMovies()).
   %       Alternatively, this can be passed as an array of SMFs
   %       matching the size of obj.TRArray, allowing for things like
   %       loading raw data for each pair in TRArraay.
   SMF = smi_core.SingleMoleculeFitting;
   
   % Top level directory for saving results.
   % NOTE: If left empty, obj.saveResults() will try to use
   %       obj.SMF.Data.ResultsDir.  If that is also empty, a default
   %       will be set to pwd().
   SaveDir {mustBeText(SaveDir)} = '';
   
   % Verbosity level of obj.performFullAnalysis(). (Default = 1)
   Verbose {mustBeInteger(Verbose)} = 1;
```
Methods:
- **computeDimerDurations**:
  computes the length of dimer events in StateSequence
- **computeLogLikelihood**: computes the LogLikelihood of a sequence in the HMM
- **computePairSeparation**: computes separation between traj. pairs in TRArray
- **computeViterbiPath**: estimates state sequence for a HMM using Viterbi alg
- **createAllMovies**: creates dimer movies for all pairs in TRArray
- **createDimerMovie**:
  creates a movie of trajectories superimposed on raw data
- **createSummaryPlot**: creates a multi-panel summary plot of the HMM analysis
- **estimateRateParameters**: estimates rate parameters from dimer candidates
- **findDimerCandidates**: finds dimer candidate pairs between TR1 and TR2
- **findDimerCandidatesFromFiles**: creates a TRArray from the provided files
- **generateEmissionMatrix**:
  computes emission matrix for an observed separation
- **generateEmissionPDFs**: creates array of emission density function handles
- **generateTransitionMatrix**:
  generates the transition matrix of a Markov model
- **isolateCandidateTRArray**: isolates dimer candidate pairs in TRArray
- **performFullAnalysis**: performs all analyses on the data in obj.TRArray
- **plotDimerPairInfo**:
  creates various plots related to dimer pair trajectories
- **plotSepDistribs**: makes histograms of the distributions of separations
- **saveResults**: saves useful results of the Hidden Markov model analysis.

For a detailed example, see
[Example_HMM.m](../../examples/Example_HMM.m).
