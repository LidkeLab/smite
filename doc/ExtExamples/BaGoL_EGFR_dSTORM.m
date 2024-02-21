%% Bayesian Grouping of Localizations (BaGoL) Example for dSTORM EGFR
%  BaGoL is run for a region of the data using hierarchical  
%  Bayes to find the distribution of the number of localizations per  
%  emitter (Xi) from the data itself. In the second run of BaGoL, the entire 
%  region is processed using the found distribution as an input.
%
% Requirements and Setup:
%   1. Windows 64 bit OS
%   2. MATLAB 2017b or higher versions
%   3. Statistics and Machine Learning Toolbox
%   4. BaGoL class
%   5. Set MATLAB directory to BaGoL directory.
%
% Description of how to run...
%   1. Set the parameters in the following section
%   2. Set a cutoff for filtering too bright localizations 
%   3. Select a region of data to learn Xi.
%   4. Set the parameter 'DataROI' for the entire data set to be processed.
%
% Results include:
%   Saved Results:
%     SR_Im.png:                 Traditional super-resolution image. 
%     Post-Im.png:               Posterior image or histogram image of the chain
%                                (weighted average over all models).
%     MAPN-Im.png:               MAPN image which is the image of localizations from the
%                                most likely model. 
%     Overlay_SR_Map.png:        Overlay of grayscale SR-image and color MAPN image.
%     Overlay_SR_Post.png:       Overlay of grayscale SR-image and color posterior image. 
%     Overlay_SR_Map_circle.png: Overlay of the SR & MAPN coordinates where 
%                                every coordinate is represented by a circle  
% 		                         located at the given location and a radius 
%                                of double of the given precision.
%     Xi.png:                    Number of localizations per emitter dist.
%     NND.png:                   Histogram of nearest neighbor distances from
%                                MAPN-coordinates. 
%     BaGoL_X-SE.png:            Histogram of X-localization precisions after grouping. 
%     BaGoL_Y-SE.png:            Histogram of Y-Localization precisions after grouping.
%     LocsScatter-MAPN.fig:      Plot of time color-coded localizations and
%                                MAPN-coordinates.
%     MAPN.mat:                  Structure containing the MAPN-coordinates of emitters.
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
%     MAPN.Nmean: Mean number of binding events per docking strand

warning('OFF', 'stats:kmeans:FailedToConvergeRep');

%% Important Parameters

PixelSize = 97;   % (nm)
IntensityCutoff = 4000; %Intensity cutoff (double of the main intensity peak)
PrecCorrect = 0.75; %Precision inflation. Added to X_SE, Y_SE (nm)
ClusterDrift = 0; %No cluster drift
ROIsz = 200; %ROI size (nm)
OverLap = 15; %Size of overlapping region (nm)
SaveDir = fullfile('Results_EGFR'); %Saving Directory
if ~isdir(SaveDir)
    mkdir(SaveDir); 
end
%% Load data
DataDir = fullfile('Data');
FileName = 'SMR_dSTORM_EGFR';
load(fullfile(DataDir,FileName))

%% Remove bright localizations that are likely to be more than one emitters 
IndP = SMR.Photons < IntensityCutoff;

%% First run: Estimating Xi using hierarchical Bayes
Xi = [1 6]; %[gamma,etha] parameters for gamma prior.
DataROI = [110 134 110 134]; %Region to find Xi (pixel) [XStart XEnd YStart YEnd]

Ind = SMR.X >= DataROI(1) & SMR.X < DataROI(2) & SMR.Y > DataROI(3) & SMR.Y < DataROI(4) & IndP;
SMD = [];
SMD.X = PixelSize*SMR.X(Ind); %Converting units from pixel to nm
SMD.Y = PixelSize*SMR.Y(Ind); %Converting units from pixel to nm
SMD.Z = [];
SMD.X_SE = PixelSize*SMR.X_SE(Ind)+PrecCorrect; %Converting units from pixel to nm
SMD.Y_SE = PixelSize*SMR.Y_SE(Ind)+PrecCorrect; %Converting units from pixel to nm
SMD.Z_SE = [];
SMD.FrameNum = SMR.Nframes*single((SMR.DatasetNum(Ind)-1))+single(SMR.FrameNum(Ind));

%Setting the class properties
EGF = smi.BaGoL;
EGF.SMD = SMD;
EGF.ROIsize = ROIsz; %size of the subregions to be processed
EGF.Overlap = OverLap; %Overlapping region size between the adjacent regions.
EGF.Xi = Xi; %Parameters for prior distribution (gamma in this case)
EGF.N_Burnin = 2000; %Length of Burn-in chain
EGF.N_Trials = 3000; %Length of post-burn-in chain
EGF.ChainFlag = 0; %Save the chain
EGF.Drift = ClusterDrift; %No cluster drift;
EGF.HierarchFlag = 1; %Flag indicating to learn Xi
EGF.Cutoff = 50;
EGF.PixelSize = 2;

%Analyzing the data
EGF.analyze_all();

%% Adjusting the class properties for the second run, using results from the
% first run to fix \xi for the second run

DataROI = [110 158 110 158]; %Region to analyze (pixel) [XStart XEnd YStart YEnd]
ImSize = (DataROI(2)-DataROI(1))*PixelSize;
XStart = DataROI(1)*PixelSize;
YStart = DataROI(3)*PixelSize;
PixSZ = 2; %Pixel size of the output images (nm)

%Populating SMD with localizations from the selected region.

Ind = SMR.X >= DataROI(1) & SMR.X < DataROI(2) & SMR.Y > DataROI(3) & SMR.Y < DataROI(4) & IndP;
SMD = [];
SMD.X = PixelSize*SMR.X(Ind); %Converting units from pixel to nm
SMD.Y = PixelSize*SMR.Y(Ind); %Converting units from pixel to nm
SMD.Z = [];
SMD.X_SE = PixelSize*SMR.X_SE(Ind)+PrecCorrect; %Converting units from pixel to nm
SMD.Y_SE = PixelSize*SMR.Y_SE(Ind)+PrecCorrect; %Converting units from pixel to nm
SMD.Z_SE = [];
SMD.FrameNum = SMR.Nframes*single((SMR.DatasetNum(Ind)-1))+single(SMR.FrameNum(Ind));

%Ajusting the class properties
EGF.SMD = SMD;
EGF.Xi = [1,mean(prod(EGF.XiChain'))]; %Mean number of blinking estimated from the former section
EGF.N_Burnin = 2000; %Length of Burn-in chain
EGF.N_Trials = 3000; %Length of post-burn-in chain
EGF.PImageSize = ImSize; %size of the posterior image
EGF.PixelSize = PixSZ; %Pixel size for the posterior image
EGF.PImageFlag = 1; %Producing the posterior image
EGF.XStart = XStart;
EGF.YStart = YStart;
EGF.HierarchFlag = 0;

%Analyzing the data
tic;
EGF.analyze_all();
T = toc();
fprintf('It took %g seconds to process the selected region of dSTORM EGFR data.\n',T)

%Generating output plots and images.
ScaleBar = 1000; %length of scale bars (nm)
EGF.saveBaGoL(ScaleBar,SaveDir);
EGF.plotMAPN(SaveDir)

RadiusScale = 2;
BaGoL.genSRMAPNOverlay(EGF.SMD, EGF.MAPN, ImSize, ImSize, PixSZ, SaveDir, ...
                       XStart, YStart, RadiusScale, ScaleBar);

%NND from EGFR vs NND of randomly distributed data.
[~,D]=knnsearch([EGF.MAPN.X,EGF.MAPN.Y],[EGF.MAPN.X,EGF.MAPN.Y],'k',2);
Dis = D(:,2);
figure('Visible','off');histogram(Dis(Dis<4*median(Dis)),'normalization','pdf','BinEdges',0:2:250)
Rho = length(EGF.MAPN.X)/(ImSize*ImSize);
X = 0:0.1:250;
Y = 2*pi*Rho*X.*exp(-pi*Rho*X.^2);
hold;plot(X,Y,'r','linewidth',1.5)
xlabel('NND(nm)');ylabel('PDF')
legend('Found-NND','Random-NND')
xlim([0 250]);ylim([0 0.050])
set(gca,'FontSize',15)
print(gcf,fullfile(SaveDir,'NND-Hist+Random'),'-dpng')
