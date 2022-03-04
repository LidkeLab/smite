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
%
% Description of how to run...
%   1. Set the parameters in the following section via BaGoLHier_wrapper.
%   2. Results are placed in ResultsDir.  The main results (the .mat file
%      containing the BGL object which includes the input SMD structure,
%      output MAPN structure, the posterior image [PImage], the Xi chain, etc.,
%         BaGoL_Results_*_ResultsStruct.mat
%      ) are placed at the top level, while various plots are placed in the
%      subdirectory identified by the dataset that was processed.
%
% Results include:
%   Saved Results:
%     BaGoL_X-SE.png:            Histogram of X-localization precisions after
%                                grouping. 
%     BaGoL_Y-SE.png:            Histogram of Y-Localization precisions after
%                                grouping.
%     LocsScatter-MAPN.fig:      Plot of time color-coded localizations and
%                                MAPN-coordinates.
%     MAPN.mat:                  Structure containing the MAPN-coordinates of
%                                emitters.
%     MAPN-Im.png:               MAPN image which is the image of localizations
%                                from the most likely model. 
%     NND.png:                   Histogram of nearest neighbor distances from
%                                MAPN-coordinates. 
%     NNDScaledData.png:         PDFs of nearest neighbor distances + random
%                                at the same density scaled by 99% of the data.
%     NNDScaledRandom.png:       PDFs of nearest neighbor distances + random
%                                at the same density scaled by 99% of the
%                                random PDF.
%     Overlay_SR_Map_circle.png: Overlay of the SR and MAPN coordinates where 
%                                every coordinate is represented by a circle  
%                                located at the given location and a radius 
%                                of double of the given precision.  Input
%                                localizations (data) are shown by green
%                                circles and found emitter locations are shown
%                                by magenta circles.
%     Overlay_SR_Map.png:        Overlay of grayscale SR-image and color MAPN
%                                image.
%     Overlay_SR_Post.png:       Overlay of grayscale SR-image and color
%                                posterior image. 
%     Post-Im.png:               Posterior image or histogram image of the
%                                chain (weighted average over all models).
%     SR_Im.png:                 Traditional super-resolution image. 
%     Xi.png:                    Distribution of localizations per emitter.
%     XiChain.png                Plot of the Xi chain after burnin.
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
%    FileNameIn    name of file containing coordinate SMD structure
%    DataDir       directory in which FileNameIn is located
%    SaveDir       directory in which saved results are put
%    BaGoLParams   structure with the following parameters:
%       ImageSize         Image size (pixel)
%       PixelSize         Pixel size (nm)
%       OutputPixelSize   Pixel size for posterior images (nm)
%       IntensityCutoff   Intensity cutoff
%       SE_Adjust         Precision inflation applied to SE (nm)
%       ClusterDrift      Expected magnitude of drift (nm/frame)
%       ROIsz             ROI size for RJMCMC (nm)
%       OverLap           Size of overlapping region (nm)
%       Xi                Loc./emitter parameters for [lambda] (Poisson) or
%                         [k theta] (Gamma) prior
%       DataROI           [Xmin, Xmax, Ymin, Ymax] (pixel)
%       N_Burnin          Length of Burn-in chain
%       N_Trials          Length of post-burn-in chain
%       NSamples          Number of samples (in N_Burnin and N_Trials) before
%                         sampling Xi
%       Y_Adjust          Apply coordinate transform for y-coordinates if
%                         non-empty in the form Y = Y_Adjust - Y (pixels)
%       N_NN              Minimum # of nearest neighbors to survive filtering
%
%    SE_Adjust adds to X_SE and Y_SE, so inflates the precision.  For DNA_PAINT
%    data, SE_Adjust = 1--2 nm, while for dSTORM, slightly bigger values should
%    be used.
%
%    k and theta are the shape and scale parameters for the Gamma probability
%    distribution function.
%
%    DataROI is defined when running BaGoL's over only part of the
%    image.  If DataROI is empty, use the whole image.
%
% OUTPUTS:
%    BGL           BaGoL object containing the results of the analysis
%    ...           various plots, images and saved results detailed above

% Created by
%    Mohamadreza Fazel (2019) and Michael J. Wester (2022), Lidke Lab

%% Important Parameters
PixelSize = BaGoLParams.PixelSize;       %Pixel size (nm)
OutputPixelSize = BaGoLParams.OutputPixelSize;
IntensityCutoff = BaGoLParams.IntensityCutoff; %Intensity cutoff
SE_Adjust = BaGoLParams.SE_Adjust;       %Precision inflation (nm)
ClusterDrift = BaGoLParams.ClusterDrift; %Emitter drift velocities (nm/frame)
ROIsz = BaGoLParams.ROIsz;               %ROI size (nm)
OverLap = BaGoLParams.OverLap;           %Size of overlapping region (nm)
SZ = BaGoLParams.ImageSize;              %Image size (pixel)
N_Burnin = BaGoLParams.N_Burnin;         %Length of Burn-in chain
N_Trials = BaGoLParams.N_Trials;         %Length of post-burn-in chain
NSamples = BaGoLParams.NSamples;         %Number of samples before sampling Xi
DataROI  = BaGoLParams.DataROI;          %Optional analysis ROI (nm)
Y_Adjust = BaGoLParams.Y_Adjust;         %LL vs UL origin transform
ImSize = SZ*PixelSize;                   %Image size (nm)
N_NN = BaGoLParams.N_NN;                 %Min # for nearest neighbor filter

% --------- Initialize BaGoL

%% Load data
load(fullfile(DataDir,FileNameIn));
% The above data is assumed to be an SMD structure.  Note that X/Y and
% X_SE/Y_SE units are pixels, later to be transformed into nm for the BaGoL
% analysis.
if ~exist('SMD', 'var')
   if exist('SMR', 'var')
      SMD = SMR;
   else
      error('SMD structure expected!');
   end
end
% If Y_Adjust is non-empty (y-coordinate origin is in the upper left and y
% increases downward), adjust the Y values so their origin is in the lower left
% and y increases upward.
%if ~isempty(Y_Adjust)
%   SMD.Y = Y_Adjust - SMD.Y; 
%end
% BaGoL seems to work better with non-negative coordinates.
minX = min(SMD.X);   % pixels
minY = min(SMD.Y);   % pixels
if minX < 0
   SMD.X = SMD.X - minX;
   %if ~isempty(DataROI)
   %   DataROI(1 : 2) = DataROI(1 : 2) - minX;
   %end
end
if minY < 0
   SMD.Y = SMD.Y - minY;
   %if ~isempty(DataROI)
   %   DataROI(3 : 4) = DataROI(3 : 4) - minY;
   %end
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
% This should now be done earlier in preparing the dataset.
%if isfield('SMD', 'Photons')
%   IndP = SMD.Photons < IntensityCutoff;
%else
    IndP = SMD.FrameNum > 0;   % This should be all the localizations
%end

Xi = BaGoLParams.Xi; %[k, theta] parameters for gamma prior.
% Make the run on a smaller subregion.
%DataROI = [80 120 120 160];%Region to find Xi (pixel) [XStart XEnd YStart YEnd]
if ~isempty(DataROI)
   if ~isempty(Y_Adjust)
      SMD.Y = Y_Adjust - SMD.Y;
   end
   Ind = SMD.X >= DataROI(1) & SMD.X <= DataROI(2) & ...
         SMD.Y >= DataROI(3) & SMD.Y <= DataROI(4) & IndP;
   if ~isempty(Y_Adjust)
      SMD.Y = Y_Adjust - SMD.Y;
   end
   ImSize = (DataROI(2) - DataROI(1))*PixelSize; 
   XStart = DataROI(1)*PixelSize;
   YStart = DataROI(3)*PixelSize;
else
   Ind = IndP;
end
n_Ind = sum(Ind);
fprintf('DataROI localizations kept = %d out of %d\n', n_Ind, numel(Ind));
if n_Ind == 0
   error('No localizations kept!');
end

figure; hold on; plot(SMD.X, SMD.Y, 'k.'); title(FileName); hold off;
saveas(gcf, fullfile(SaveDirLong, 'FULL'), 'png');
SMD.X = PixelSize*SMD.X(Ind);
SMD.Y = PixelSize*SMD.Y(Ind);
SMD.Z = [];
%SMD.X_SE = PixelSize*SMD.X_SE(Ind)+SE_Adjust;
%SMD.Y_SE = PixelSize*SMD.Y_SE(Ind)+SE_Adjust;
SMD.X_SE = PixelSize*SMD.X_SE(Ind);
SMD.Y_SE = PixelSize*SMD.Y_SE(Ind);
SMD.Z_SE = [];
if isfield('SMD', 'NFrames')
   SMD.FrameNum = SMD.NFrames*single((SMD.DatasetNum(Ind)-1))+single(SMD.FrameNum(Ind));
end

%% Localizations are filtered based on NND within 3 median of localization sigma
if N_NN > 0
   Prec_Median = median([SMD.X_SE;SMD.Y_SE]);
   [~,D]=knnsearch([SMD.X,SMD.Y],[SMD.X,SMD.Y],'K',length(SMD.X));
   D(:,1)=[];

   ID = D < 3*Prec_Median;
   N = sum(ID,2);
   Ind = N >= N_NN;

   SMD.X = SMD.X(Ind);
   SMD.Y = SMD.Y(Ind);
   SMD.X_SE = SMD.X_SE(Ind);
   SMD.Y_SE = SMD.Y_SE(Ind);
   SMD.FrameNum = SMD.FrameNum(Ind);
   fprintf('NN localizations kept = %d out of %d\n', sum(Ind), numel(Ind));
   if n_Ind == 0
      error('No localizations kept!');
   end
end

figure; hold on; plot(SMD.X/PixelSize, SMD.Y/PixelSize, 'k.'); title(FileName); hold off;
saveas(gcf, fullfile(SaveDirLong, 'ROI'), 'png');
% Setting the class properties
BGL = smi.BaGoL;
BGL.SMD = SMD;
BGL.HierarchFlag = 1;      %Use a Hierarchial Prior to learn Xi if 1
BGL.NSamples = NSamples;   %Number of samples before sampling Xi
BGL.SE_Adjust = SE_Adjust; %Localization precision adjustment (nm)
BGL.PImageFlag = 1;        %Generate Posterior image if 1
BGL.PImageSize = ImSize;   %Size of the output posterior images
BGL.ROIsize = ROIsz;       %Size of the subregions to be processed
BGL.Overlap = OverLap;     %Overlapping region size between adjacent regions
BGL.Xi = Xi;   %Parameters for prior distribution (gamma in this case)
BGL.N_Burnin = N_Burnin;   %Length of Burn-in chain
BGL.N_Trials = N_Trials;   %Length of post-burn-in chain
BGL.ChainFlag = 0;         %Save the chain if 1
BGL.Drift = ClusterDrift;  %Expected magnitude of drift (nm/frame)
BGL.PixelSize = OutputPixelSize; %Pixel size for the posterior image
if ~isempty(DataROI)
   BGL.XStart = XStart;
   BGL.YStart = YStart;
end

% ---------- Run BaGoL

%Analyzing the data
BGL.analyze_all()

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

fprintf('Saving BGL ...\n');
try
   save(fullfile(SaveDir, ...
                 sprintf('BaGoL_Results_%s_ResultsStruct', FileName)), 'BGL');
end
ScaleBarLength = 1000;   % nm
fprintf('saveBaGoL ...\n');
try
   BGL.saveBaGoL(ScaleBarLength, SaveDirLong, 1);
end
fprintf('plotMAPN ...\n');
try
   BGL.plotMAPN(SaveDirLong, 'on');
end
fprintf('plotNND_PDF ...\n');
try
   BGL.plotNND_PDF(SaveDirLong)
end
fprintf('genSRMAPNOverlay ...\n');
try
   MAPN = BGL.MAPN;
   MAPN.X_SE = max(1, MAPN.X_SE);
   MAPN.Y_SE = max(1, MAPN.Y_SE);
   RadiusScale = 2;
   if isempty(DataROI)
      BGL.genSRMAPNOverlay(BGL.SMD, MAPN, ImSize, ImSize, 'rescale', ...
                           SaveDirLong, 0, 0, RadiusScale, ScaleBarLength);
   else
      BGL.genSRMAPNOverlay(BGL.SMD, MAPN, ImSize, ImSize, 'rescale',    ...
                           SaveDirLong, XStart, YStart, RadiusScale,        ...
                           ScaleBarLength);
   end
catch ME
   fprintf('### PROBLEM with genSRMAPNOverlay ###\n');
   fprintf('%s\n', ME.identifier);
   fprintf('%s\n', ME.message);
end
%BGL.errPlot(BGL.MAPN);
%saveas(gcf, fullfile(SaveDirLong, 'MAPN_SE'), 'fig');
try
   if numel(Xi) == 1
      plot(BGL.XiChain(:, 1), 'k.');
   else
      plot(BGL.XiChain(:, 1) .* BGL.XiChain(:, 2), 'k.');
   end
   saveas(gcf, fullfile(SaveDirLong, 'XiChain'), 'png');
end
close all
L = BGL.XiChain;
L = L(N_Burnin/NSamples + 1 : end, :);
l = L(:, 1) .* L(:, 2);
m = mean(l);
v = var(l);
theta = v / m;
k = m / theta;
fid = fopen(fullfile(SaveDir, 'prior.txt'), 'w');
fprintf(fid, '(k, theta) = (%f, %f)\n', k, theta);
fclose(fid);

end
