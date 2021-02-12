classdef BaGoL < handle
%BaGoL Implements a Bayesian Grouping of Localizations (BaGoL)
%
%Single molecule localization based super-resolution data can contain
%repeat localizations from the same emitter. These localizations can be
%grouped to give better localization precision.  BaGoL explores the
%possible number of emitters and their positions that can explain the
%observed localizations and uncertainties. An 'emitter' is a
%blinking/binding point source that generates a single 'localization'
%each time there is a blinking/binding event. Localizations must be
%frame connected.
%
%The core algorithm uses Reversible Jump Markov Chain Monte Carlo to
%add, remove and move emitters and to explore the classification of
%localizations to emitters.  Prior information on the distribution of
%localizations per emitter is required and the localization precisions
%are assumed to be accurate. The prior distribution is
%parameterized by either a Poisson or Gamma distribution function.
%
%The primary BaGoL outputs are a 'Posterior Image' and MAPN coordinates
%that can also be used to generate a 'MAPN image'. The Posterior Image
%shows a probablity distribution of emitter locations that is a weighted
%average over the number of emitters, emitter locations, and
%classification of localizations to emitters. The Maximum a Posteriori
%Number of emitters (MAPN) result uses only the information from
%the most likely number of emitters.  MAPN emitter coordinates and their
%uncertainties are returned and images from these coordinates are
%generated similarly to a traditional SR reconstruction.
%
%This class implements the pre/post-processing steps needed to analyze
%a spatially extended data set typical of super-resolution experiments.
%This includes breaking data into subregions and collating the results
%from the subregions.
%
%
% USAGE:
%   B=BaGoL()       %create object
%   B.SMD=....      %set properties
%   B.analyze_all() %run complete analysis
%
%The class also has several methods for visualizing and saving results.
%See 'doc BaGoL' for the complete list of methods.
%
% PROPERTIES:
%   SMD:        A structure containing the fields:
%       X:          Vector of X localization positions (Pixel)(Nx1)
%       Y:          Vector of Y localization positions (Pixel)(Nx1)
%       Z:          Vector of Z localization positions (Pixel)(Nx1)(Optional)
%       X_SE:       Vector of X localization standard error (Pixel)(Nx1)
%       Y_SE:       Vector of Y localization standard error (Pixel)(Nx1)
%       Z_SE:       Vector of Z localization standard error (Pixel)(Nx1)(Optional)
%       FrameNum:   Vector of localization frame numbers (Nx1)
%       PixelSize:  Camera pixel size (nm/Pixel) ??is it microns?
%   Lambda:     Loc./emitter params [lambda] (Poisson) or [k theta] (Gamma)
%   Cutoff:     Max distance between localizations in analyzed cluster (nm) (Default=20)
%   NNN:        Minimum Number of Nearest Neighbors (Default=2)
%   NNR:        Nearest Neighbor Radius (nm)(Default = Inf)
%   ROIsize:    ROI size for RJMCMC (nm)(Default=200)
%   Overlap:    Allowed overlap between subregions (nm)(Default=50)
%   Drift:      Expected magnitude of drift (nm/frame)(Default=0)
%   SE_Adjust:  Adjustement of localization precisions (nm) (Default=0)
%   N_Burnin:   Number of jumps in burn in of RJMCMC chain (Default=2000)
%   N_Trials:   Number of jumps in RJMCMC chain post burn in (Default=3000)
%   P_Jumps:    Proposal probabilities for RJMCMC Jumps
%               [Move, Allocate, Add, Remove]
%               sum(P_Jumps) must equal 1.
%               (Default = [0.25, 0.25, 0.25, 0.25])
%   PixelSize:  The pixel size for output posterior images (nm) (Default=1)
%   PImageFlag: Generate Posterior image. 0 or 1. (Default=0)
%   PImageSize: Size of the output posterior images (nm)
%   ChainFlag:  Save RJMCMC chain. 0 or 1. (Default=0)
%   XStart:     X starting coordinate of output posterior images (nm)
%   YStart:     Y starting coordinate of output posterior images (nm)
%
%   ClusterSMD: An array of SMD structures for each cluster
%   MAPN:       SMD output structure containing fields:
%       X:      Vector of X emitter positions (nm)(Kx1)
%       Y:      Vector of Y emitter positions (nm)(Kx1)
%       Z:      Vector of Z emitter positions (nm)(Kx1)
%       X_SE:   Vector of X localization standard error (nm)(Kx1)
%       Y_SE:   Vector of Y localization standard error (nm)(Kx1)
%       Z_SE:   Vector of Z localization standard error (nm)(Kx1)
%       Nmean:  Mean number of localizations per emitter (Kx1)
%   PImage:     Posterior Image
%   Chain:      Cell array of RJMCMC chains for pre-clusters (See BaGoL_RJMCMC)
%
% REQUIRES:
%   MATLAB 2016 or higher versions.
%   Statistics and Machine learning toolbox.
%
% CITATION: "Sub-Nanometer Precision using Bayesian Grouping of Localizations"
%           Mohamadreza Fazel, Michael J. Wester, Sebastian Restrepo Cruz,
%           Sebastian Strauss, Florian Schueder, Jennifer M. Gillette,
%           Diane S. Lidke, Bernd Rieger, Ralf Jungmann, Keith A. Lidke
%

% Created by:
%   Mohamadreza Fazel, Lidke Lab 2019
%
    properties
        SMD %Structure containing localization coordinates and uncertainties
        ClusterSMD %SMD array corresponding to clusters   
        MAPN %MAPN output coordinates and uncertainties.
        PImage; %Posterior image
        PImageFlag = 1; %Generate Posterior image. 0 or 1. (Default=1)
        PImageSize = []; %Size of Posterior image (nm)
        PixelSize = 2; %Pixel size for posterior images (nm) (Default=2)
        XStart = []; %X starting coordinate of posterior images (nm)
        YStart = []; %Y starting coordinate of posterior images (nm)
        P_Jumps=[0.25 0.25 0.25 0.25]; %Proposal probabilities for RJMCMC Jumps 
        N_Burnin=2000; %Number of jumps in burn in of RJMCMC chain (Default=2000)
        N_Trials=3000; %Number of jumps in RJMCMC chain post burn in (Default=3000)
        ROIsize=200; %ROI size for RJMCMC (nm)(Default=200)
        Overlap=50; %Allowed overlap between subregions (nm)(Default=50)  
        NNR = Inf; %Nearest Neighbor Radius (nm)(Default = Inf)
        NNN = 2; %Minimum Number of Nearest Neighbors (Default=2)
        Drift = 0; %Expected magnitude of drift velocity (nm/frame)(Default=0) 
        SE_Adjust = 0; %Localization precision adjustment (nm) (Default=0)
        Cutoff=20; %Max distance between localizations in analyzed cluster (nm) (Default=20)
        Lambda; %Loc./emitter params [lambda] (Poisson) or [k, theta] (Gamma)
        ChainFlag = 0; %Save RJCMC chain. 0 or 1. (Default=0)
        Chain = {}; %Cell array of RJMCMC chains for clusters
    end
    
    methods
        
        function obj=analyze_all(obj)
        %analyze_all Implements complete BaGoL analysis of SR dataset
        %
        %This  method performs all steps necessary to analyze a
        %super-resolution data set with BaGoL. Internally, it intializes
        %output structures, generates ROIs and then clusters for RJMCMC,
        %then loops over clusters and collates the output. 
        %
        %This method sequentially calls: 
        %   ROIs = obj.genROIs()
        %   obj.precluster(...)
        %       BaGoL_RJMCMC(...)
        %       genPosterior(...)
        %       genMAPN(...)
        %   removeOverlap(obj,ROIs,MeanX,MeanY,ii)
        %   genMAPN(...); %% why is this done again??
        %   genROIs(...) %% why is this done again??
        %   genPosterior(...)
        %   plotNND(...)
        %
        % USAGE:
        %   B.analyze_all
            
       % plotNND(obj,SaveDir)
            %Save results first, then ROIs are generated. Next, the outliers
            %are filtered and the preclusters are identified. Then, each
            %precluster is processed usin the RJMCMC algorithm. Finally, the
            %chain from RJMCMC is used to generate the posterior image and the
            %MAPN coordinates.
            
            %auto calculate posterior size
            if isempty(obj.XStart)
               obj.initPostIm(); 
            end
            
            Ind = obj.SMD.X_SE == 0 | obj.SMD.Y_SE == 0;
            obj.SMD.X(Ind) = [];
            obj.SMD.Y(Ind) = [];
            obj.SMD.X_SE(Ind) = [];
            obj.SMD.Y_SE(Ind) = [];
            obj.SMD.FrameNum(Ind) = [];
            obj.XStart = obj.XStart;
            obj.YStart = obj.YStart;
            if isempty(obj.SMD.FrameNum)
                obj.SMD.FrameNum = zeros(size(obj.SMD.X));
            else
                obj.SMD.FrameNum = single(obj.SMD.FrameNum);
            end
            ROIs = obj.genROIs();
            if obj.PImageFlag == 1
                if isempty(obj.PImageSize) || ~isscalar(obj.PImageSize)
                    error('PImageSize must be given as a scalar');
                end
                SZ = ceil(obj.PImageSize/obj.PixelSize);
                PostIm = zeros(SZ,'single');
            end
            obj.precluster(ROIs);
            if obj.ChainFlag == 1
                obj.Chain = cell(length(obj.ClusterSMD),1);
            end
            obj.MAPN.X = [];
            obj.MAPN.Y = [];
            obj.MAPN.Z = [];
            obj.MAPN.X_SE = [];
            obj.MAPN.Y_SE = [];
            obj.MAPN.Z_SE = [];
            obj.MAPN.AlphaX = [];
            obj.MAPN.AlphaY = [];
            obj.MAPN.AlphaX_SE = [];
            obj.MAPN.AlphaY_SE = [];
            obj.MAPN.Nmean = [];
            obj.MAPN.N = 0;
            ClustNumHeirar = length(obj.ClusterSMD);
            if isempty(obj.SMD.FrameNum)
               MaxAlpha = 0; 
            elseif max(obj.SMD.FrameNum~=0)
               MaxAlpha = obj.Drift; 
            else
               MaxAlpha = 0; 
            end
            warning('OFF', 'stats:kmeans:FailedToConvergeRep');
            for nn = 1:ClustNumHeirar
                if nn/50 == ceil(nn/50)
                    fprintf('Cluster: %g out of %g\n',nn,ClustNumHeirar);
                end
                
                AnimFlag = 0;
                [TChain]=BaGoL.BaGoL_RJMCMC(obj.ClusterSMD(nn),obj.Lambda,MaxAlpha,obj.P_Jumps,obj.N_Trials,obj.N_Burnin,AnimFlag);
                
                if obj.PImageFlag == 1
                    PostIm = obj.genPosterior(PostIm,SZ,TChain,ROIs,nn);
                end
                obj.genMAPN(TChain,ROIs,nn);
                if obj.ChainFlag == 1
                    obj.Chain{nn} = TChain; 
                end
            end
            warning('ON', 'stats:kmeans:FailedToConvergeRep');
            if obj.PImageFlag == 1
                obj.PImage = PostIm;
            end
            
            if obj.ChainFlag == 1
                for mm = length(obj.ClusterSMD):-1:1
                     if isempty(obj.Chain{mm})
                         obj.Chain(mm) = []; 
                         obj.ClusterSMD(mm) = [];
                     end
                end
            end
        end
    end
    
    methods (Static) 
       [SMD,SMC]=genCluster(StuctType,Scale,Ndist,PSF,MeanPhoton,Prec_Cutoff,DriftVec,PlotFlag); 
       [SMD,SMD_combined]=frameConnect(SMDin,LOS,MaxDistance,MaxFrameGap,FitType);
       saveMAPN(Directory,FileType,MAPN)
       errPlot(SMD);
       SMD=loadPICASSOh5(DataDir,FileName)
       [Chain]=BaGoL_RJMCMC(SMD,Lambda,MaxAlpha,PMove,NChain,NBurnin,DEBUG)
       [SRIm,MapIm]=makeIm(SMD,MAPN,SZ,PixSize,XStart,YStart)
       [Alpha,XShift,YShift,Aligned,Chain] = align_template(Temp,Input,Start,Cutoff,NChain,PlotFlag)
       [ImageOut,Image] = scalebar(Image,PixelSize,Length,Location)
       [Im]=scaleIm(Im,Percentile)
       [ParamPoiss,ParamGam]=fitLambda(NMean,Percent)
       dispIm()
       Im=genSRMAPNOverlay(SMD,MAPN,XSize,YSize,PixelSize,SaveDir,XStart,YStart,RadiusScale,ScaleBarLength)
    end
    
end
