function BGL = hierBaGoL_analysis(SMD, FileNameIn, SaveDir, BaGoLParams)
%Bayesian Grouping of Localizations (BaGoL)
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
%   1. Set the parameters in the following section via hierBaGoL_wrapper.
%   2. Results are placed in ResultsDir.  The main results (the .mat file
%      containing the BGL object which includes the input SMD structure,
%      output MAPN structure, the posterior image [PImage], the Xi chain, etc.,
%         BaGoL_Results_*_ResultsStruct.mat
%      and the .mat file which contains the MAPN coordinates, etc.,
%         MAPN_*.mat
%      ) are placed at the top level, while various plots are placed in the
%      subdirectory identified by the dataset that was processed.
%
% Results include:
%   Saved Results:
%     See above as well.
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
%     MAPN_NmeanHist.png:        Histogram of the number of localizations per
%                                emitter found for the MAPN results.
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
%     MAPN: Clusters information are stored in this property:
%     MAPN.X: X-Centers (nm)
%     MAPN.Y: Y-Centers (nm)
%     MAPN.X_SE: X-Centers precisions (nm)
%     MAPN.Y_SE: Y-Centers precisions (nm)
%     MAPN.AlphaX: X-Drifts of clusters (nm/frame)
%     MAPN.AlphaY: Y-Drifts of clusters (nm/frame)
%     MAPN.AlphaX_SE: X-Drift precisions (nm/frame)
%     MAPN.AlphaY_SE: Y-Drift precisions (nm/frame)
%     MAPN.Nmean: Mean number of binding events per docking strand or
%                 localizations per emitter
%
% INPUTS:
%    SMD           Single Molecule Data structure containg X, Y, X_SE, Y_SE
%    FileNameIn    name of file containing coordinate SMD structure
%    SaveDir       directory in which saved results are put
%    BaGoLParams   structure with the following parameters:
%       ImageSize         Image size (pixel)
%       PixelSize         Pixel size (nm)
%       OutputPixelSize   Pixel size for posterior images (nm)
%       SE_Adjust         Precision inflation applied to SE (nm)
%       ClusterDrift      Expected magnitude of drift (nm/frame)
%       ROIsz             ROI size for RJMCMC (nm)
%       OverLap           Size of overlapping region (nm)
%       Cutoff            Pre-clustering cutoff (nm) [default = ROIsz]
%       Xi                Loc./emitter parameters for [lambda] (Poisson) or
%                         [k theta] (Gamma) prior
%       DataROI           [Xmin, Xmax, Ymin, Ymax] (pixel)
%       N_Burnin          Length of Burn-in chain
%       N_Trials          Length of post-burn-in chain
%       NSamples          Number of samples (in N_Burnin and N_Trials) before
%                         sampling Xi
%       Y_Adjust          Apply coordinate transform for y-coordinates if
%                         non-empty in the form Y = Y_Adjust - Y (pixels)
%       
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
SZ        = BaGoLParams.ImageSize;   % Image size (pixel)
PixelSize = BaGoLParams.PixelSize;   % Pixel size (nm)
DataROI   = BaGoLParams.DataROI;     % Optional analysis ROI (nm)
Y_Adjust  = BaGoLParams.Y_Adjust;    % LL vs UL origin transform
ImSize    = SZ*PixelSize;            % Image size (nm)
Xi        = BaGoLParams.Xi;          % [k, theta] parameters for gamma prior
XStart    = 0;
YStart    = 0;

% --------- Initialize BaGoL

%% Load data
%load(fullfile(DataDir, FileNameIn));
% The above data is assumed to be an SMD structure.  Note that X/Y and
% X_SE/Y_SE units are pixels, later to be transformed into nm for the BaGoL
% analysis.
if ~exist('SMD', 'var')
   error('SMD structure expected!');
end
% If Y_Adjust is non-empty (y-coordinate origin is in the upper left and y
% increases downward), adjust the Y values so their origin is in the lower left
% and y increases upward.
%if ~isempty(Y_Adjust)
%   SMD.Y = Y_Adjust - SMD.Y; 
%end

% Eliminate trailing _Results* from the FileName for saving results.
FileName = regexprep(FileNameIn, '\.mat$', '');
FileName = regexprep(FileName, '_Results$', '');
FileName = regexprep(FileName, '_ResultsStruct', '');
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

%% Filter localizations
% BaGoL works better with non-negative coordinates.
Verbose = 2;
SMD = smi_helpers.Filters.filterNonNeg(SMD, Verbose);

% Make the run on a smaller subregion.
%DataROI = [80 120 120 160];%Region to find Xi (pixel) [XStart XEnd YStart YEnd]

%% Examine a ROI if desired
if ~isempty(DataROI)
   if ~isempty(Y_Adjust)
      SMD.Y = Y_Adjust - SMD.Y;
   end
   Ind = SMD.X >= DataROI(1) & SMD.X <= DataROI(2) & ...
         SMD.Y >= DataROI(3) & SMD.Y <= DataROI(4);
   if ~isempty(Y_Adjust)
      SMD.Y = Y_Adjust - SMD.Y;
   end
   % Convert to nm as BaGoL is expecting nm.
   ImSize = (DataROI(2) - DataROI(1))*PixelSize; 
   XStart = DataROI(1)*PixelSize;
   YStart = DataROI(3)*PixelSize;

   n_Ind = sum(Ind);
   fprintf('DataROI localizations kept = %d out of %d\n', n_Ind, numel(Ind));
   if n_Ind == 0
      error('No localizations kept!');
   end
else
   % This should be all the localizations not previously filtered out.
   Ind = SMD.FrameNum > 0;
end

% FULL plot (full set of SMD localizations).
figure;
hold on;
plot(SMD.X, SZ - SMD.Y, 'k.');
title('ALL localizations');
xlabel('x (pixels)');
ylabel('y (pixels)');
hold off;
saveas(gcf, fullfile(SaveDirLong, 'FULL'), 'png');

% Convert from pixels to nm.
SMD.X = PixelSize*SMD.X(Ind);
SMD.Y = PixelSize*SMD.Y(Ind);
SMD.Z = [];
SMD.X_SE = PixelSize*SMD.X_SE(Ind);
SMD.Y_SE = PixelSize*SMD.Y_SE(Ind);
SMD.Z_SE = [];
if isfield('SMD', 'NFrames')
   SMD.FrameNum = ...
      SMD.NFrames*single((SMD.DatasetNum(Ind)-1))+single(SMD.FrameNum(Ind));
end

% ROI plot (localizations only in ROI).
figure;
hold on;
plot(SMD.X/PixelSize, SZ - SMD.Y/PixelSize, 'k.');
title('ROI localizations');
xlabel('x (pixels)');
ylabel('y (pixels)');
hold off;
saveas(gcf, fullfile(SaveDirLong, 'ROI'), 'png');

%% Set the class properties
BGL = smi.BaGoL;
BGL.SMD = SMD;
   % Use a Hierarchial Prior to learn Xi if 1
BGL.HierarchFlag = 1;
   % Save the chain if 1
BGL.ChainFlag = 0;
   % Number of samples before sampling Xi
BGL.NSamples = BaGoLParams.NSamples;
   % Generate Posterior image if 1
BGL.PImageFlag = 1;
   % Size of the output posterior images
BGL.PImageSize = ImSize;
   % Size of the subregions to be processed
BGL.ROIsize = BaGoLParams.ROIsz;
   % Overlapping region size between adjacent regions
BGL.Overlap = BaGoLParams.OverLap;
   % Parameters for prior distribution (gamma in this case)
BGL.Cutoff = BaGoLParams.Cutoff;
   % Pre-clustering cutoff (nm) [default = ROIsz]
BGL.Xi = Xi;
   % Length of Burn-in chain
BGL.N_Burnin = BaGoLParams.N_Burnin;
   % Length of post-burn-in chain
BGL.N_Trials = BaGoLParams.N_Trials;
   % Expected magnitude of drift (nm/frame)
BGL.Drift = BaGoLParams.ClusterDrift;
   % Pixel size for the posterior image
BGL.PixelSize = BaGoLParams.OutputPixelSize;
   % Localization precision adjustment (nm)
BGL.SE_Adjust = BaGoLParams.SE_Adjust;
if ~isempty(DataROI)
   BGL.XStart = XStart;
   BGL.YStart = YStart;
end

% ---------- Run BaGoL

% Analyzing the data
BGL.analyze_all()

% ---------- Save Results and Plots

% This file can be huge for many localizations, so only produce it if the
% number of input localizations is not too large.  Most of the space is taken
% up by the Chain.
if numel(SMD.X) <= 100000
   fprintf('Saving BGL ...\n');
   try
%     BGL.Chain  = [];
%     BGL.PImage = [];
%     BGL.SMD    = [];
      save(fullfile(SaveDir, ...
              sprintf('BaGoL_Results_%s_ResultsStruct', FileName)), 'BGL');
   catch ME
      fprintf('### PROBLEM with saving BGL ###\n');
      fprintf('%s\n', ME.identifier);
      fprintf('%s\n', ME.message);
   end
end

% This file contains just the MAPN coordinates, so much smaller than the
% ResultsStruct file and should always be saved.  Also, save a few extra
% parameters that are handy for displaying ROIs.
fprintf('saveMAPN ...\n');
try
   MAPN = BGL.MAPN;
   PImageSize = BGL.PImageSize;
   PixelSizeSAVE = PixelSize;
   PixelSize = BaGoLParams.OutputPixelSize;
   save(fullfile(SaveDir, sprintf('MAPN_%s', FileName)), 'MAPN', ...
        'XStart', 'YStart', 'PImageSize', 'PixelSize');
   PixelSize = PixelSizeSAVE;
catch ME
   fprintf('### PROBLEM with saveMAPN ###\n');
   fprintf('%s\n', ME.identifier);
   fprintf('%s\n', ME.message);
end

% In general, the number of pixels that make up the scale bar is (from
% scalebar.m)
%    barlength = round(Length / PixelSize);
% where Length is the desired scale bar length in um (SMITE) or nm (BaGoL),
% and the input camera PixelSize is in corresponding units (um/pixel or
% nm/pixel).
%
% For SMITE Gaussian images, the SMD coordinates (pixels) are magnified by
% SRImageZoom (see +smi_vis/@GenerateImages/gaussianImage.m where
%    SMDin.X = SMD.X * SRImageZoom, etc.
% ), so the true scale bar length in original coordinates is thus
%    ScalebarLength / SRImageZoom
% which by default (see +smi_vis/@GenerateImages/gaussianImage.m:
%    ScalebarLength   scalebar length (um) [default: 10 um]
% and +smi/@SMLM/SMLM.m:
%    SRImageZoom  = 20 % magnification factor for SR     images generated
%    SRCircImZoom = 25 % magnification factor for circle images generated
% ) is 10 um / 20 = 0.5 um = 500 nm.
% For circle images, the default is 10 um / 25 = 0.4 um = 400 nm.
%
% Now for BaGoL Gaussian and circle plots, the input SMD structure from SMITE
% provided to BaGoL has pixel coordinates.  These are converted into nm via
% SMD.PixelSize (see "Convert from pixels to nm" above where there are lines
% like
%    SMD.X = PixelSize * SMD.X(Ind);
% in which PixelSize = BaGoLParams.PixelSize, that is, the input pixel size,
% meaning the camera number of nm per pixel).  To produce the Gaussian/circle
% plots, it is necessary to convert nm into the OutputPixelSize, which
% confusingly is called PixelSize in BaGoL (+smi/@BaGoL/hierBaGoL_analysis.m):
%    BGL.PixelSize = BaGoLParams.OutputPixelSize;
% For example, in +smi/@BaGoL/genSRMAPNOverlay.m, there are lines like
%    X = ((SMD.X - Xstart) / PixelSize);
% Therefore, the true scale bar length for these BaGoL plots will be simply
% ScaleBarLength (as the conversion to output pixels is done properly in
% scalebar.m).  ScaleBarLength defaults to 0.5 um = 500 nm (see line below).
ScaleBarLength = 500;   % nm
fprintf('ScaleBarLength = %g nm, OutputPixelSize = %g nm\n', ...
        ScaleBarLength, BaGoLParams.OutputPixelSize);

fprintf('saveBaGoL ...\n');
try
   BGL.saveBaGoL(ScaleBarLength, SaveDirLong, 1);
catch ME
   fprintf('### PROBLEM with saveBaGoL ###\n');
   fprintf('%s\n', ME.identifier);
   fprintf('%s\n', ME.message);
end

fprintf('plotMAPN ...\n');
try
   BGL.plotMAPN(SaveDirLong, 'on');
catch ME
   fprintf('### PROBLEM with plotMAPN ###\n');
   fprintf('%s\n', ME.identifier);
   fprintf('%s\n', ME.message);
end

fprintf('plotNND_PDF ...\n');
try
   BGL.plotNND_PDF(SaveDirLong)
catch ME
   fprintf('### PROBLEM with plotNND_PDF ###\n');
   fprintf('%s\n', ME.identifier);
   fprintf('%s\n', ME.message);
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

fprintf('XiChain ...\n');
try
   if numel(Xi) == 1
      plot(BGL.XiChain(:, 1), 'k.');
   else
      plot(BGL.XiChain(:, 1) .* BGL.XiChain(:, 2), 'k.');
   end
   hold on
   title('Xi chain');
   xlabel(sprintf('RJMCMC jumps sampled every %d', BaGoLParams.NSamples));
   ylabel('\lambda (localizations per emitter)');
   hold off
   saveas(gcf, fullfile(SaveDirLong, 'XiChain'), 'png');
catch ME
   fprintf('### PROBLEM with XiChain ###\n');
   fprintf('%s\n', ME.identifier);
   fprintf('%s\n', ME.message);
end

if numel(Xi) == 2
   fprintf('kChain ...\n');
   try
      plot(BGL.XiChain(:, 1), 'k.');
      hold on
      title('k chain');
      xlabel(sprintf('RJMCMC jumps sampled every %d', BaGoLParams.NSamples));
      ylabel('k in gamma(k, \theta)');
      hold off
      saveas(gcf, fullfile(SaveDirLong, 'kChain'), 'png');
   catch ME
      fprintf('### PROBLEM with kChain ###\n');
      fprintf('%s\n', ME.identifier);
      fprintf('%s\n', ME.message);
   end
end

if numel(Xi) == 2
   fprintf('thetaChain ...\n');
   try
      plot(BGL.XiChain(:, 2), 'k.');
      hold on
      title('\theta chain');
      xlabel(sprintf('RJMCMC jumps sampled every %d', BaGoLParams.NSamples));
      ylabel('\theta in gamma(k, \theta)');
      hold off
      saveas(gcf, fullfile(SaveDirLong, 'thetaChain'), 'png');
   catch ME
      fprintf('### PROBLEM with thetaChain ###\n');
      fprintf('%s\n', ME.identifier);
      fprintf('%s\n', ME.message);
   end
end

fprintf('MAPN_NmeanHistogram ...\n');
try
   h = histogram(MAPN.Nmean);
   hold on
   lo = BGL.N_Burnin/BGL.NSamples + 1;
   AB = median(BGL.XiChain(lo : end, :));
   X = h.BinLimits(1) : 0.1 : h.BinLimits(2); 
   plot(X, sum(h.Data(:)) * h.BinWidth * gampdf(X, AB(1), AB(2)));
%  plot(X, sum(BGL.XiChain(lo : end, 1) .* BGL.XiChain(lo : end, 2)) .* ...
%              gampdf(X, AB(1), AB(2)));
   title('MAPN localizations per emitter');
   xlabel('localizations per emitter');
   ylabel('frequency');
   legend('MAPN.Nmean', 'gampdf(XiChain)', 'Location', 'NorthEast');
   hold off
   saveas(gcf, fullfile(SaveDirLong, 'MAPN_NmeanHist'), 'png');
catch ME
   fprintf('### PROBLEM with MAPN_NmeanHist ###\n');
   fprintf('%s\n', ME.identifier);
   fprintf('%s\n', ME.message);
end

% Visualize Precluster Results
fprintf('Preclusters ...\n');
try
   figure;
   hold on
   for i = 1 : numel(BGL.ClusterSMD)
      plot(BGL.ClusterSMD(i).X, BGL.ClusterSMD(i).Y, '.');
   end
   axis equal
   hold off
   saveas(gcf, fullfile(SaveDirLong, 'Preclusters'), 'png');
catch ME
   fprintf('### PROBLEM with Preclusters ###\n');
   fprintf('%s\n', ME.identifier);
   fprintf('%s\n', ME.message);
end

close all

fprintf('Producing prior.txt ...\n');
try
   L = BGL.XiChain;
   L = L(BGL.N_Burnin/BGL.NSamples + 1 : end, :);
   l = L(:, 1) .* L(:, 2);
   m = mean(l);
   v = var(l);
   theta = v / m;
   k = m / theta;
   fid = fopen(fullfile(SaveDirLong, 'prior.txt'), 'w');
   fprintf(fid, '(k, theta) = (%f, %f)\n', k, theta);
   fclose(fid);
catch ME
   fprintf('### PROBLEM with producing prior.txt ###\n');
   fprintf('%s\n', ME.identifier);
   fprintf('%s\n', ME.message);
end

end
