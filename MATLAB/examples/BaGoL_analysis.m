%% Bayesian Grouping of Localizations (BaGoL)
%
%  This function is adapted from EGFR_dSTORM.m in the BaGoL distribution.
%
%  BaGoL is run for a part of the region of the data with a broad gamma   
%  prior to find the distribution of the number of localizations per  
%  emitter.  In the second run of BaGoL, the entire region is processed 
%  using the found distribution as a prior.

% Requirements and Setup:
%   1. MATLAB 2016 or higher versions
%   2. Statistics and Machine Learning Toolbox
%   3. BaGoL class
%   4. Set MATLAB directory to BaGoL directory.
%
% Description of how to run...
%   1. Set the parameters in the following section
%   2. Set a cutoff for filtering too bright localizations 
%   3. Select the region to be processed.
%   4. Set the parameter 'DataROI' for the region of data to be processed.
%
% Results include:
%   Saved Results:
%     Histogram of nearest neighbor distances between emitters.
%     Histogram of mean number of binding/blinking events per emitter.
%     Histogram of localization precisions for emitters.
%     SR-image of the selected region.
%     SR-image of the selected region after NND filtering.
%     MAPN image of the selected region.
%     Posterior image of the selected region.
%   Output available on work space:
%     MAPN: Clusters information are stored in this propertiy:
%     MAPN.X: X-Centers (nm)
%     MAPN.Y: Y-Centers (nm)
%     MAPN.X_SE: X-Centers precisions (nm)
%     MAPN.Y_SE: Y-Centers precisions (nm)
%     MAPN.AlphaX: X-Drifts of clusters (nm/frame)
%     MAPN.AlphaY: Y-Drifts of clusters (nm/frame)
%     MAPN.AlphaX_SE: X-Drift precisions (nm/frame)
%     MAPN.AlphaY_SE: Y-Drift precisions (nm/frame)
%     MAPN.Nmean: Mean number of binding events per docking strand

function BGL = BaGoL_analysis(FileNameIn, DataDir, SaveDir, BaGoLParams)
%
% INPUTS:
%    FileNameIn    name of file containing coordinate (SMR/SMD) structure
%    DataDir       directory in which FileNameIn is located
%    SaveDir       directory in which saved results are put
%    BaGoLParams   structure with the following parameters:
%       ImageSize         Image size (pixel)
%       PixelSize         Pixel size (nm)
%       OutputPixelSize   Pixel size for posterior images (nm)
%       IntensityCutoff   Intensity cutoff
%       NNR               At least NNN localizations within NNR (nm)
%       NNN               Number of Nearest Neighbors
%       NNNmax            No more than NNN locs within NNR (nm)
%       SE_Adjust         Precision inflation applied to SE (nm)
%       ClusterDrift      Expected magnitude of drift (nm/frame)
%       ROIsz             ROI size for RJMCMC (nm)
%       OverLap           Size of overlapping region (nm)
%       PreCluster        Pre-clustering parameter (nm)
%       Lambda            Loc./emitter parameters for [lambda] (Poisson) or
%                         [k theta] (Gamma) prior
%       DataROI           1st pass [Xmin, Xmax, Ymin, Ymax] (pixel)
%       N_Burnin1         Length of Burn-in chain      (1st pass)
%       N_Trials1         Length of post-burn-in chain (1st pass)
%       N_Burnin2         Length of Burn-in chain      (2nd pass)
%       N_Trials2         Length of post-burn-in chain (2nd pass)
%       Skip1stPass       If true, skip the 1st pass and directly use the
%                         values for Lambda defined above for the prior
%       Make2ndPassExpPrior If true, adjust the 1st pass Lambda posterior to an
%                         exponential to used as the prior for the 2nd pass
%
% NOTES:
%    If Make2ndPassExpPrior is true, then the posterior Lambda from the 1st
%    pass % is converted into an exponential prior for the 2nd pass via:
%       [k2_prior, theta2_prior] = [1, k1_posterior * theta1_posterior]
%    noting that k is the shape parameter (1 produces an exponential) and theta
%    is the scale parameter for a gamma distribution, where the product
%    k * theta is the mean of the distribution for the number of localizations
%    per emitter, so the idea is to perserve this mean value while retaining an
%    exponential shape.  This is important in dSTORM.
%
%    For very sparse datasets, turn off filtering (NNR = inf, NN = 0,
%    NNNmax = inf).
%
%    SE_Adjust adds to X_SE and Y_SE, so inflates the precision.  For DNA_PAINT
%    data, SE_Adjust = 1--2 nm, while for dSTORM, slightly bigger values should
%    be used.
%
%    k and theta are the shape and scale parameters for the Gamma probability
%    distribution function.
%
%    DataROI is defined when running BaGoL's first pass over only part of the
%    image.  If DataROI is empty, use the whole image.
%
% OUTPUTS:
%    BGL           BaGoL object containing the results of the analysis
%    ...           various plots, images and saved results detailed above

% Created by
%    Mohamadreza Fazel (2019) and Michael J. Wester (2021), Lidke Lab

%% Important Parameters
PixelSize = BaGoLParams.PixelSize;       %Pixel size (nm)
OutputPixelSize = BaGoLParams.OutputPixelSize;
IntensityCutoff = BaGoLParams.IntensityCutoff; %Intensity cutoff
NNR = BaGoLParams.NNR;%There should be at least NNN localizations within NNR(nm)
NNN = BaGoLParams.NNN;                   %Number of Nearest Neighbors
SE_Adjust = BaGoLParams.SE_Adjust;       %Precision inflation (nm)
ClusterDrift = BaGoLParams.ClusterDrift; %No cluster drift
ROIsz = BaGoLParams.ROIsz;               %ROI size (nm)
OverLap = BaGoLParams.OverLap;           %Size of overlapping region (nm)
PreCluster = BaGoLParams.PreCluster;     %Pre-clustering parameter (nm)
SZ = BaGoLParams.ImageSize;              %Image size (pixel)
N_Burnin1 = BaGoLParams.N_Burnin1;       %Length of Burn-in chain      (1st)
N_Trials1 = BaGoLParams.N_Trials1;       %Length of post-burn-in chain (1st)
N_Burnin2 = BaGoLParams.N_Burnin2;       %Length of Burn-in chain      (2nd)
N_Trials2 = BaGoLParams.N_Trials2;       %Length of post-burn-in chain (2nd)

% --------- Initialize BaGoL

%% Load data
load(fullfile(DataDir,FileNameIn));
% The above data is assumed to be an SMR structure, however, as smite comes in,
% the above may be an SMD structure, so deal with that.  The two structures are
% the same for the fields (Photons, X, Y, X_SE, Y_SE, DatasetNum, FrameNum)
% used below, but different for Nframes/NFrames.
if ~exist('SMR', 'var') && exist('SMD', 'var')
   SMR = SMD;
   SMR.Nframes = SMD.NFrames;
end
% BaGoL seems to work better with non-negative coordinates.
minX = min(SMR.X);   % pixels
minY = min(SMR.Y);   % pixels
if minX < 0
   SMR.X = SMR.X - minX;
end
if minY < 0
   SMR.Y = SMR.Y - minY;
end
% Eliminate trailing _Results* from the FileName for saving results.
FileName = regexprep(FileNameIn, '_Results.*$', '');

% Save the BaGoL _ResultsStruct.mat file in SaveDir and the rest of the BaGoL
% outputs in SaveDirLong.  This arrangement is chosen so that Results_BaGoL
% holds only uniquely named files/directories for the situation where several
% _ResultsStruct.mat files reside in the same (higher level) directory,
% therefore the results of multiple BaGoL runs will share common space in
% Results_BaGoL.
if ~isfolder(SaveDir)
   mkdir(SaveDir); 
end
SaveDirLong = fullfile(SaveDir, FileName);
if ~isfolder(SaveDirLong)
   mkdir(SaveDirLong);
end

% Remove bright localizations that are likely to be more than one emitters 
IndP = SMR.Photons < IntensityCutoff;

% First run: Estimating Lambda using a wide gamma distribution
Lambda = BaGoLParams.Lambda; %[gamma,etha] parameters for gamma prior.
% Make the first run on a smaller subregion (commented out for now).
%DataROI = [80 120 120 160]; %Region to find Lambda (pixel) [XStart XEnd YStart YEnd]
if ~isempty(BaGoLParams.DataROI)
   Ind = SMR.X > BaGoLParams.DataROI(1) & SMR.X < BaGoLParams.DataROI(2) & ...
         SMR.Y > BaGoLParams.DataROI(3) & SMR.Y < BaGoLParams.DataROI(4) & ...
         IndP;
else
   Ind = IndP;
end
n_Ind = sum(Ind);
fprintf('[1] Localizations kept = %d\n', n_Ind);
if n_Ind == 0
   error('No localizations kept!');
end

nPre = numel(Ind);
SMD = [];
SMD.X = PixelSize*SMR.X(Ind);
SMD.Y = PixelSize*SMR.Y(Ind);
SMD.Z = [];
SMD.X_SE = PixelSize*SMR.X_SE(Ind)+SE_Adjust;
SMD.Y_SE = PixelSize*SMR.Y_SE(Ind)+SE_Adjust;
SMD.Z_SE = [];
SMD.FrameNum = SMR.Nframes*single((SMR.DatasetNum(Ind)-1))+single(SMR.FrameNum(Ind));

%Setting the class properties
BGL = smi.BaGoL;
%BGL = BaGoL;
BGL.SMD = SMD;
BGL.ROIsize = ROIsz; %size of the subregions to be processed
BGL.Overlap = OverLap; %Overlapping region size between the adjacent regions.
BGL.Lambda = Lambda; %Parameters for prior distribution (gamma in this case)
BGL.Cutoff = PreCluster; %Parameter for hierarchical clustering (nm)
BGL.NNN = NNN; %Localizations with less than NNN neighbors within NNR distance
               %are considered outliers
BGL.NNR = NNR; %Localizations with less than NNN neighbors within NNR distance
               %are considered outliers
BGL.N_Burnin = N_Burnin1; %Length of Burn-in chain
BGL.N_Trials = N_Trials1; %Length of post-burn-in chain
BGL.ChainFlag = 0; %Save the chain if 1
BGL.Drift = ClusterDrift; %Expected magnitude of drift (nm/frame)
BGL.PixelSize = OutputPixelSize; %Pixel size for the posterior image

% ---------- 1st Pass of BaGoL (to estimate Lambda)

if ~BaGoLParams.Skip1stPass
   %Analyzing the data
   BGL.analyze_all()

   %Histogram of the number of blinking per emitter used to calculate the
   %average number of blinkings per emitter, which is used in the next section.
   % This is the posterior from the first pass and the prior for the second
   % pass of the algorithm.
   NMean = BGL.MAPN.Nmean;
   if isempty(NMean)
      error('BaGoL MAPN.Nmean is empty---cannot proceed further!');
   end
   [ParamPoiss, ParamGam] = BGL.fitLambda(NMean);
   print(gcf, fullfile(SaveDirLong, 'Lambda_Hist_PriorEst'), '-dpng');
   close

   % Estimate lambda prior for the second run using the results from the first
   % run
   %PdP = fitdist(NMean,'gamma');
   PdP = ParamGam;

   if BaGoLParams.Make2ndPassExpPrior
      PdPsave = PdP;
      PdP.a = 1;
      PdP.b = PdPsave.a * PdPsave.b;
      fprintf('Lambda adjust: [%.3f, %.3f] -> [%.3f, %.3f]\n', ...
              PdPsave.a, PdPsave.b, PdP.a, PdP.b);
   end
end

% ---------- 2nd Pass of BaGoL

% Adjusting the class properties for the second run, using a more  
%restrictive gamma prior with the found average number of blinkings in the  
%previous section and saving the results.

%DataROI = [110 158 110 158]; %Region to analyze (pixel) [XStart XEnd YStart YEnd]
%ImSize = SZ;% (DataROI(2)-DataROI(1))*PixelSize;

ImSize = SZ*PixelSize;

%Populating SMD with localizations from the selected region.
%FileList = dir(fullfile(DataDir,'*.mat'));
    
%if ~isempty(BaGoLParams.DataROI)
%   Ind = SMR.X > BaGoLParams.DataROI(1) & SMR.X < BaGoLParams.DataROI(2) & ...
%         SMR.Y > BaGoLParams.DataROI(3) & SMR.Y < BaGoLParams.DataROI(4) & ...
%         IndP;
%else
   Ind = IndP;
%end
n_Ind = sum(Ind);
fprintf('[2] Localizations kept = %d\n', n_Ind);
if n_Ind == 0
   error('No localizations kept!');
end

SMD = [];
SMD.X = PixelSize*SMR.X(Ind);
SMD.Y = PixelSize*SMR.Y(Ind);
SMD.Z = [];
SMD.X_SE = PixelSize*SMR.X_SE(Ind)+SE_Adjust;
SMD.Y_SE = PixelSize*SMR.Y_SE(Ind)+SE_Adjust;
SMD.Z_SE = [];
SMD.FrameNum = SMR.Nframes*single((SMR.DatasetNum(Ind)-1))+single(SMR.FrameNum(Ind));

%Ajusting the class properties
BGL.SMD = SMD;

if ~BaGoLParams.Skip1stPass
   BGL.Lambda = [PdP.a, PdP.b]; %Mean number of blinking estimated from the
                                %former section, used for Poisson dist.
end

BGL.N_Burnin = N_Burnin2; %Length of Burn-in chain
BGL.N_Trials = N_Trials2; %Length of post-burn-in chain
BGL.PImageSize = ImSize; %size of the posterior image
BGL.PixelSize = OutputPixelSize; %Pixel size for the posterior image
BGL.PImageFlag = 1; %Producing the posterior image
BGL.XStart = 0;
BGL.YStart = 0;

%Analyzing the data
tic;
BGL.analyze_all()
T = toc;
nPost = numel(BGL.MAPN.X);
fprintf('%s: It took %g seconds to process the selected region.\n',FileName,T)
fprintf('%s pre- to post-BaGoL localizations: %d -> %d\n', ...
        FileName, nPre, nPost);

% ---------- Save Results and Plots

% Restore original displacements.
if minX < 0
   BGL.SMD.X  = BGL.SMD.X  + minX * PixelSize;
   BGL.MAPN.X = BGL.MAPN.X + minX * PixelSize;
   for i = 1 : numel(BGL.ClusterSMD)
      BGL.ClusterSMD(i).X = BGL.ClusterSMD(i).X + minX * PixelSize;
   end
end
if minY < 0
   BGL.SMD.Y  = BGL.SMD.Y  + minY * PixelSize;
   BGL.MAPN.Y = BGL.MAPN.Y + minY * PixelSize;
   for i = 1 : numel(BGL.ClusterSMD)
      BGL.ClusterSMD(i).Y = BGL.ClusterSMD(i).Y + minY * PixelSize;
   end
end

save(fullfile(SaveDir, ...
              sprintf('BaGoL_Results_%s_ResultsStruct', FileName)), 'BGL');
ScaleBarLength = 1000;   % nm
BGL.saveBaGoL(ScaleBarLength, SaveDirLong);
BGL.plotMAPN(3, SaveDirLong);
BGL.plotNND_PDF(SaveDirLong)
BGL.genSRMAPNOverlay(BGL.SMD, BGL.MAPN, SZ*PixelSize, SZ*PixelSize, ...
                     'rescale', SaveDirLong, 0, 0, 1);
%BGL.errPlot(BGL.MAPN);
%saveas(gcf, fullfile(SaveDirLong, 'MAPN_SE'), 'fig');
close all

end
