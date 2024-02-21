classdef BaGoL < handle
% BaGoL Implements a Bayesian Grouping of Localizations (BaGoL)
%
% Single molecule localization based super-resolution data can contain
% repeat localizations from the same emitter. These localizations can be
% grouped to give better localization precision.  BaGoL explores the
% possible number of emitters and their positions that can explain the
% observed localizations and uncertainties (data). An 'emitter' is a
% blinking/binding point source that generates a single 'localization'
% each time there is a blinking/binding event. Localizations must be
% frame connected.
%
% The core algorithm uses Reversible Jump Markov Chain Monte Carlo to
% add, remove and move emitters and to explore the allocation of
% localizations to emitters. The localization precisions are assumed to be
% accurate. The prior distribution is parameterized by either a Poisson 
% or Gamma distribution function which can be either given as an input or 
% learned via the hierarchical Bayes approach.
%
% The primary BaGoL outputs are a 'Posterior Image' and MAPN coordinates
% that can also be used to generate a 'MAPN image'. The Posterior Image
% shows a probability distribution of emitter locations that is a weighted
% average over the number of emitters, emitter locations, and
% allocation of localizations to emitters. The Maximum a Posteriori
% Number of emitters (MAPN) result uses only the information from
% the most likely number of emitters within the chain.  MAPN emitter 
% coordinates and their uncertainties are returned and images from these 
% coordinates are generated similarly to a traditional SR reconstruction.
%
% This class implements the pre/post-processing steps needed to analyze
% a spatially extended data set typical of super-resolution experiments.
% This includes breaking data into subregions and collating the results
% from the subregions.
%
%
% USAGE:
%   B=BaGoL()       % create object
%   B.SMD=....      % set properties
%   B.analyze_all() % run complete analysis
%
% The class also has several methods for visualizing and saving results.
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
%       PixelSize:  Camera pixel size (nm/Pixel) 
%   Xi:         Loc./emitter [lambda] (Poisson) or [k theta] (Gamma).
%               When learning Xi this is used to initialize a chain. If
%               two initial values are given it uses a gamma and otherwise
%               a poisson prior.
%   Alpha_Xi:   Shape parameter of Xi gamma hyper prior
%   Alpha_Xi:   Scale parameter of Xi gamma hyper prio
%   ROIsize:    ROI size for RJMCMC (nm) (Default=200)
%   Overlap:    Allowed overlap between subregions (nm)(Default=50)
%   Drift:      Expected magnitude of drift (nm/frame)(Default=0)
%   SE_Adjust:  Adjustement of localization precisions (nm) (Default=0)
%   N_Burnin:   Number of samples in burn in of RJMCMC chain (Default=2000)
%   N_Trials:   Number of samples in RJMCMC chain post burn in (Default=3000)
%   P_Jumps:    Proposal probabilities for RJMCMC Jumps
%               [Move, Allocate, Add, Remove]
%               sum(P_Jumps) must equal 1.
%               (Default = [0.25, 0.25, 0.25, 0.25])
%   PixelSize:  The pixel size for output posterior images (nm) (Default=1)
%   PImageFlag: Generate Posterior image. 0 or 1. (Default=0)
%   HierarchFlag: Use hierarchical Bayse to learn Xi. 0 or 1. (Default=0)
%   NSamples:   Number of RJMCMC samples before sampling Xi (Default=10)
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
%   XiChain:    Chain of Xi (locs per emitters dist. either Poisson or gamma)
%   SaveName:   Final results are saved under this name (Default: BaGoL)
%
% REQUIRES:
%   MATLAB 2016 or higher versions.
%   Statistics and Machine learning toolbox.
%
% CITATION: "High-Precision Estimation of Emitter Positions using Bayesian
%           Grouping of Localizations", Mohamadreza Fazel, Michael J. Wester,
%           David J. Schodt, Sebastian Restrepo Cruz, Sebastian Strauss,
%           Florian Schueder, Thomas Schlichthaerle, Jennifer M. Gillette,
%           Diane S. Lidke, Bernd Rieger, Ralf Jungmann and Keith A. Lidke,
%           Nature Communications, 13(7152), November 22, 2022, 1--11,
%           (DOI: 10.1038/s41467-022-34894-2).

% Created by:
%   Mohamadreza Fazel, Lidke Lab 2019
%
    properties
        SMD % Structure containing localization coordinates and uncertainties
        ClusterSMD % SMD array corresponding to subregions   
        MAPN % MAPN output coordinates and uncertainties.
        PImage; % Posterior image
        PImageFlag = 1; % Generate Posterior image. 0 or 1. (Default=1)
        PImageSize = []; % Size of Posterior image (nm)
        PixelSize = 2; % Pixel size for posterior images (nm) (Default=2)
        XStart = []; % X starting coordinate of posterior images (nm)
        YStart = []; % Y starting coordinate of posterior images (nm)
        P_Jumps=[0.25 0.25 0.25 0.25]; % Proposal probabilities for RJMCMC Jumps 
        N_Burnin=2000; % Number of samples in burn in of RJMCMC chain (Default=2000)
        N_Trials=3000; % Number of samples in RJMCMC chain post burn in (Default=3000)
        NSamples = 10; % Number of samples in other params before taking a Xi sample
        ROIsize=200; % ROI size for RJMCMC (nm)(Default=200)
        Overlap=20; % Allowed overlap between subregions (nm)(Default=20)  
        Drift = 0; % Expected magnitude of drift velocity (nm/frame)(Default=0) 
        SE_Adjust = 0; % Localization precision adjustment (nm) (Default=0)
        Xi=20; % Loc./emitter params [lambda] (Poisson) or [k, theta] (Gamma) (Default=20)
        Alpha_Xi = 1; %Shape parameter of Xi hyper prior (Default = 1)
        Beta_Xi = 50; %Scale parameter of Xi hyper prior (Default = 50)
        ChainFlag = 0; % Save RJCMC chain. 0 or 1. (Default=0)
        HierarchFlag = 0; % Use a Hierarchial Prior to learn XI
        Chain = {}; % Cell array of RJMCMC chains for clusters
        XiChain; % Chain of samples for Xi. [lambda]->Poisson, [k, theta]-> gamma
        SaveName = 'BaGoL'; % Final results are saved under this name
        Cutoff; %preclustering parameter used in hierarchical clustering (nm) (Default = ROIsize)
    end
    
    methods
        
        function obj=analyze_all(obj)
        % analyze_all Implements complete BaGoL analysis of SR dataset
        %
        % This  method performs all steps necessary to analyze a
        % super-resolution data set with BaGoL. Internally, it initializes
        % output structures, generates ROIs for RJMCMC, then loops over 
        % ROIs and collates the output. 
        %
        % This method sequentially calls: 
        %   ROIs = obj.genROIs()
        %   obj.assignROIs(...)
        %       BaGoL_RJMCMC(...)
        %       genPosterior(...)
        %       genMAPN(...)
        %
        % USAGE:
        %   B.analyze_all
        %
            
            %auto calculate posterior image size
            if isempty(obj.XStart)
               obj.initPostIm(); 
            end
            
            Ind = obj.SMD.X_SE == 0 | obj.SMD.Y_SE == 0;
            obj.SMD.X(Ind) = [];
            obj.SMD.Y(Ind) = [];
            obj.SMD.X_SE(Ind) = [];
            obj.SMD.Y_SE(Ind) = [];
            obj.SMD.FrameNum(Ind) = [];
%             obj.XStart = obj.XStart;
%             obj.YStart = obj.YStart;
            if isempty(obj.Cutoff)
                obj.Cutoff = obj.ROIsize;
            end
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
            %obj.assignROIs(ROIs);
            obj.precluster(ROIs);

            % Visualize Precluster Results
            figure
            hold on
            for ii = 1:length(obj.ClusterSMD)
                plot(obj.ClusterSMD(ii).X,obj.ClusterSMD(ii).Y,'.')
            end
            axis equal
          
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
            
            ClustNumHeirar = length(obj.ClusterSMD);
            if isempty(obj.SMD.FrameNum)
               MaxAlpha = 0; 
            elseif max(obj.SMD.FrameNum~=0)
               MaxAlpha = obj.Drift; 
            else
               MaxAlpha = 0; 
            end
            warning('OFF', 'stats:kmeans:FailedToConvergeRep');
            
            if obj.HierarchFlag == 0
                %using a calibrated Xi distribution (HierarchFlag=0)
                if obj.ChainFlag == 1
                    obj.Chain = cell(length(obj.ClusterSMD),1);
                end
                for nn = 1:ClustNumHeirar
                    if nn/10 == floor(nn/10)
                        fprintf('Subregion: %g out of %g \n',nn,ClustNumHeirar);
                    end
                    
                    AnimFlag = 0;
                    [TChain]=smi.BaGoL.BaGoL_RJMCMC(obj.ClusterSMD(nn),obj.Xi,MaxAlpha,obj.P_Jumps,obj.N_Trials,obj.N_Burnin,AnimFlag);

                    if obj.PImageFlag == 1
                        PostIm = obj.genPosterior(PostIm,SZ,TChain,ROIs,nn);
                    end
                    obj.genMAPN(TChain,ROIs,nn);
                    if obj.ChainFlag == 1
                        obj.Chain{nn} = TChain; 
                    end
                end
            else
                %learning Xi (HierarchFlag=1)
                tXi = obj.Xi;
                obj.Chain = cell(length(obj.ClusterSMD),1);
                NChain = floor(obj.N_Trials/obj.NSamples);
                NIter = floor((obj.N_Trials+obj.N_Burnin)/obj.NSamples);
                NBurn = floor(obj.N_Burnin/obj.NSamples);
                NPoints = zeros(length(obj.ClusterSMD),1);
                tK = zeros(length(obj.ClusterSMD),1);
                for nn = 1:length(obj.ClusterSMD)
                    obj.Chain{nn}.Chain(NChain).N = [];
                    obj.Chain{nn}.Chain(NChain).X = [];
                    obj.Chain{nn}.Chain(NChain).Y = [];
                    obj.Chain{nn}.Chain(NChain).AlphaX = [];
                    obj.Chain{nn}.Chain(NChain).AlphaY = [];
                    obj.Chain{nn}.Chain(NChain).ID = [];
                    NPoints(nn)=length(obj.ClusterSMD(nn).X);
                end
                %tmp is a temporary structure to store some parameters
                %before saving them in the returned chain. The returned chain 
                %is sparsed by only saving the last values in the tmp struct. 
                tmp(ClustNumHeirar,1).MuX = 0;
                tmp(ClustNumHeirar,1).MuY = 0;
                tmp(ClustNumHeirar,1).AlphaX = 0;
                tmp(ClustNumHeirar,1).AlphaY = 0;
                
                PDFgrids = cell(ClustNumHeirar,1);
                for nn = 1:ClustNumHeirar
                    DX = 1;
                    X_min = min(obj.ClusterSMD(nn).X-3*obj.ClusterSMD(nn).X_SE);
                    X_max = max(obj.ClusterSMD(nn).X+3*obj.ClusterSMD(nn).X_SE);
                    X_range = X_min:DX:X_max;
                    Y_min = min(obj.ClusterSMD(nn).Y-3*obj.ClusterSMD(nn).Y_SE);
                    Y_max = max(obj.ClusterSMD(nn).Y+3*obj.ClusterSMD(nn).Y_SE);
                    Y_range = Y_min:DX:Y_max;

                    [Xg,Yg] = meshgrid(X_range,Y_range);
                    PDFgrid = zeros(size(Xg),'single');
                    for pp = 1:length(obj.ClusterSMD(nn).X)
                        PDFgrid = PDFgrid + single(normpdf(Xg,obj.ClusterSMD(nn).X(pp),obj.ClusterSMD(nn).X_SE(pp))...
                            .*normpdf(Yg,obj.ClusterSMD(nn).Y(pp),obj.ClusterSMD(nn).Y_SE(pp))); 
                    end
                    PDFgrids{nn}.PDFgrid = PDFgrid;
                end
                
                %If the chain is initialized with a single paramter use a
                %Poisson prior otherwise gamma prior (default is Poisson).
                IsPoiss = length(obj.Xi)==1;
                if IsPoiss
                    obj.XiChain = zeros(NIter,1);
                else
                    obj.XiChain = zeros(NIter,2);
                end
                
                for mm = 1:NIter
                    if mm/100 == floor(mm/100)
                        fprintf('Iteration %g out of %g over all subregions\n',mm,NIter);
                    end
                    for nn=1:ClustNumHeirar
                        if mm == 1
                            %For the first sample, the locations are intitialized to random values
                            [K,MuX,MuY,AlphaX,AlphaY,ID]=smi.BaGoL.BaGoL_RJMCMC_Hierarchical(obj.ClusterSMD(nn) ...
                                ,PDFgrids{nn}.PDFgrid,MaxAlpha,obj.P_Jumps,obj.NSamples,tXi);
                        else
                            %For jumps other than the first one, the previous values in the chain are used
                            [K,MuX,MuY,AlphaX,AlphaY,ID]=smi.BaGoL.BaGoL_RJMCMC_Hierarchical(obj.ClusterSMD(nn) ...
                                ,PDFgrids{nn}.PDFgrid,MaxAlpha,obj.P_Jumps,obj.NSamples,tXi,tmp(nn).MuX,...
                                tmp(nn).MuY,tmp(nn).AlphaX,tmp(nn).AlphaY);
                        end
                        %tmp is a temporary structure to save some params
                        tmp(nn).MuX = MuX;
                        tmp(nn).MuY = MuY;
                        tmp(nn).AlphaX = AlphaX;
                        tmp(nn).AlphaY = AlphaY;
                        %Emitters inside overlapping regions are removed before 
                        %sampling Xi.
                        Ind = removeOverlap(obj,ROIs,MuX,MuY,nn);
                        tK(nn) = length(Ind);
                        NPoints(nn) = sum(ismember(ID,Ind));
                        if mm > NBurn 
                            obj.Chain{nn}.Chain(mm-NBurn).N = K;
                            obj.Chain{nn}.Chain(mm-NBurn).X = MuX;
                            obj.Chain{nn}.Chain(mm-NBurn).Y = MuY;
                            obj.Chain{nn}.Chain(mm-NBurn).AlphaX = AlphaX;
                            obj.Chain{nn}.Chain(mm-NBurn).AlphaY = AlphaY;
                            obj.Chain{nn}.Chain(mm-NBurn).ID = ID;
                        end
                    end
                    %Sampling Xi for either Poisson or gamma 
                    if IsPoiss
                        tXi = smi.BaGoL.samplePoiss(NPoints,tK,tXi,obj.Alpha_Xi,obj.Beta_Xi);
                        obj.XiChain(mm,:) = tXi;
                    else
                        tXi = smi.BaGoL.sampleGam(NPoints,tK,tXi,obj.Alpha_Xi,obj.Beta_Xi);
                        obj.XiChain(mm,:) = tXi;
                    end
                    
                end
                
                %Generating MAPN coordinates and posterior image
                for nn = 1:ClustNumHeirar
                     TChain = obj.Chain{nn}.Chain;
                     obj.genMAPN(TChain,ROIs,nn);
                     if obj.PImageFlag == 1
                         PostIm = obj.genPosterior(PostIm,SZ,TChain,ROIs,nn);
                     end
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
       [Chain]=BaGoL_RJMCMC(SMD,Xi,MaxAlpha,PMove,NChain,NBurnin,DEBUG,MuX,MuY,AlphaX,AlphaY)
       [K,X,Y,AlphaX,AlphaY,ID]=BaGoL_RJMCMC_Hierarchical(SMD,PDFgrid,AlphaSTD,P_Jumps,NSamples,Xi,MuX,MuY,AlphaX,AlphaY)
       [SRIm,MapIm]=makeIm(SMD,MAPN,SZ,PixSize,XStart,YStart)
       Xi = sampleGam(NPoints,K,Xi,Alpha,Beta);
       Xi = samplePoiss(NPoints,K,Xi,Alpha,Beta);
       [Alpha,XShift,YShift,Aligned,Chain] = align_template (Temp,Input,Start,Cutoff,NChain,PlotFlag)
       [ImageOut,Image] = scalebar (Image,PixelSize,Length,Location)
       [Im]=scaleIm (Im,Percentile)
       [ParamPoiss,ParamGam]=fitLambda(NMean,Percent)
       dispIm()
       [Coords,Ind] = SEfilter(ROIs)
       Im=genSRMAPNOverlay(SMD,MAPN,XSize,YSize,PixelSize,SaveDir,XStart,YStart,RadiusScale,ScaleBarLength)
       BGL = hierBaGoL_analysis(SMD, FileNameIn, SaveDir, BaGoLParams)
       BGL = hierBaGoL_run(Files, DataROI, Results_BaGoL, BaGoLParams, ROIs)
    end
    
end
