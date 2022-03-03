function [SMD,SMC]=genCluster(StuctType,Scale,Ndist,PSF,MeanPhoton,Prec_Cutoff,DriftVec,PlotFlag) 
%genCluster Simulates localization data from several types of emitter patterns 
% [SMD,SMC]=BaGoL.genCluster(StuctType,Scale,Ndist,PSF,MeanPhoton,Prec_Cutoff,DriftVec,PlotFlag) 
%
%  Generates synthetic localization data in the form used by BaGoL.
%   
%  Localization data is generated from emitter positions arranged in an 
%  MxN grid, an N-mer arranged in circle, or a user input set of emitter 
%  positions. Uncertainties of the positions are calculated from the 
%  number of photons collected from each localization and converting to a 
%  localization precision using PSF_Sigma/sqrt(N_Photons).  The number of 
%  photons for each localization is drawn from an exponential distribution 
%  parameterized by the mean number of photons. The number of localizations
%  per emitter is taken from a Poisson distribution with a given mean
%  (lambda). The number of frames is the same as the length of drift vector.
%
% INPUTS:
%    StuctType:    'MxNgrid', 'Nmer', 'user'. M,N must be integers
%    Scale:         Spatial scale (nm). Options are: 
%       'MxNGrid':      Smallest nearest neighbor distance in grid
%       'Nmer':         Radius of circular N-mer
%       'user':         SMD input with X,Y fields for emitter position
%    Ndist:         Mean number of localizations per emitter taken from a Poisson
%    PSF:           PSF_Sigma (nm)
%    MeanPhotons:   Mean Photons from an exponential distribution
%    Prec_Cutoff:   Precision threshold for generated localizations (nm)
%    DriftVec:      Linear drift, D is dimensionality (nm) (DxFrame)(Default = 0)
%    PlotFlag:      0 or 1. Plots generated localizations. (Default = 0)
%
% OUTPUTS:
%    SMD:           SMD (Single Molecule Data) structure with fields:
%       SMD.X         x-coordinate (nm) (Kx1)
%       SMD.X_SE      x-precision  (nm) (Kx1)
%       SMD.Y         y-coordinate (nm) (Kx1)
%       SMD.Y_SE      y-precision  (nm) (Kx1)
%       SMD.Z         z-coordinate (nm) (Kx1) (OPTIONAL)
%       SMD.Z_SE      z-precision  (nm) (Kx1) (OPTIONAL)
%       SMD.FrameNum: time stamp of localization (used for drift)
%       SMD.ID        Emitter index  
%    SMC:           SMD (Single Molecule Data) structure with fields:
%       SMC.X         Emitter True X Postions (nm)(Nx1)
%       SMC.Y         Emitter True Y Postions (nm)(Nx1)
%       SMC.Z         Emitter True Z Postions (nm)(Nx1)
%

% Created by:
%    Mohamadreza Fazel (Lidkelab 2019)   

if isscalar(Ndist)
    Ndist = Ndist*(1+expcdf((PSF/Prec_Cutoff)^2,MeanPhoton));
end
if ~exist('DriftVec','var')
   DriftVec = zeros(1,1000); 
end
if ~exist('PlotFlag','var')
   PlotFlag = 0;  
end

STR = StuctType(end-2:end);
MaxFrame = length(DriftVec);
SMD.X = [];
SMD.Y = [];
SMD.Z = [];
SMD.X_SE = [];
SMD.Y_SE = [];
SMD.Z_SE = [];
SMD.FrameNum = [];
SMD.ID = [];
SMC.X = [];
SMC.Y = [];
if strcmp(STR,'rid')
    N = str2double(StuctType(1))*str2double(StuctType(3));
    XC = Scale:Scale:str2double(StuctType(1))*Scale;
    YC = Scale:Scale:str2double(StuctType(3))*Scale;
    SMC.X = repmat(XC,1,length(YC));
    for mm = 1:length(YC)
        SMC.Y = cat(2,SMC.Y,YC(mm)*ones(1,length(XC)));
    end
    SMC.Z = [];
    for nn = 1:N
        NThisEmitter = poissrnd(Ndist);
        if NThisEmitter == 0
            NThisEmitter = poissrnd(Ndist);
        end
        FrameNum = randperm(MaxFrame,NThisEmitter);
        for mm = 1:length(FrameNum)
            Photons = exprnd(MeanPhoton); 
            Prec = PSF/sqrt(Photons);
            SMD.X = cat(1,SMD.X,SMC.X(nn)+Prec*randn()+DriftVec(FrameNum(mm),1));
            SMD.Y = cat(1,SMD.Y,SMC.Y(nn)+Prec*randn()+DriftVec(FrameNum(mm),2));
            SMD.X_SE = cat(1,SMD.X_SE,Prec);
            SMD.Y_SE = cat(1,SMD.Y_SE,Prec);
            SMD.FrameNum = cat(1,SMD.FrameNum,FrameNum(mm));
            SMD.ID = cat(1,SMD.ID,nn);
        end
    end
    Ind = SMD.X_SE < 40 & SMD.Y_SE < 40;
    SMD.X = SMD.X(Ind);
    SMD.Y = SMD.Y(Ind);
    SMD.X_SE = SMD.X_SE(Ind);
    SMD.Y_SE = SMD.Y_SE(Ind);
    SMD.FrameNum = SMD.FrameNum(Ind);
    SMD.ID = SMD.ID(Ind);
elseif strcmp(STR,'mer') || strcmp(STR,'Mer')
    N = str2double(StuctType(1:end-3));
    Theta = linspace(0,2*pi,N+1)';
    Theta(end) = [];
    Theta = Theta + (Theta(2)-Theta(1))/2;
    R = Scale*ones(N,1);
    [XC,YC] = pol2cart(Theta,R);
    SMC.X = XC + 2*Scale;
    SMC.Y = YC + 2*Scale;
    SMC.Z = [];
    for nn = 1:N
        NThisEmitter = poissrnd(Ndist);
        if NThisEmitter == 0
            NThisEmitter = poissrnd(Ndist);
        end
        FrameNum = randperm(MaxFrame,NThisEmitter);
        for mm = 1:length(FrameNum)
            Photons = exprnd(MeanPhoton); 
            Prec = PSF/sqrt(Photons);
            SMD.X = cat(1,SMD.X,SMC.X(nn)+Prec*randn()+DriftVec(FrameNum(mm),1));
            SMD.Y = cat(1,SMD.Y,SMC.Y(nn)+Prec*randn()+DriftVec(FrameNum(mm),2));
            SMD.X_SE = cat(1,SMD.X_SE,Prec);
            SMD.Y_SE = cat(1,SMD.Y_SE,Prec);
            SMD.FrameNum = cat(1,SMD.FrameNum,FrameNum(mm)); 
            SMD.ID = cat(1,SMD.ID,nn);
        end
    end
    Ind = SMD.X_SE < 40 & SMD.Y_SE < 40;
    SMD.X = SMD.X(Ind);
    SMD.Y = SMD.Y(Ind);
    SMD.X_SE = SMD.X_SE(Ind);
    SMD.Y_SE = SMD.Y_SE(Ind);
    SMD.FrameNum = SMD.FrameNum(Ind);
    SMD.ID = SMD.ID(Ind);
elseif strcmp(STR,'ser') 
    SMC.X = Scale.X;
    SMC.Y = Scale.Y;
    SMC.Z = [];
    for nn = 1:length(SMC.X)
        NThisEmitter = poissrnd(Ndist);
        if NThisEmitter == 0
            NThisEmitter = poissrnd(Ndist);
        end
        FrameNum = randperm(MaxFrame,NThisEmitter);
        for mm = 1:length(FrameNum)
            Photons = exprnd(MeanPhoton); 
            Prec = PSF/sqrt(Photons);
            SMD.X = cat(1,SMD.X,SMC.X(nn)+Prec*randn())+DriftVec(FrameNum(mm));
            SMD.Y = cat(1,SMD.Y,SMC.Y(nn)+Prec*randn())+DriftVec(FrameNum(mm));
            SMD.X_SE = cat(1,SMD.X_SE,Prec);
            SMD.Y_SE = cat(1,SMD.Y_SE,Prec);
            SMD.FrameNum = cat(1,SMD.FrameNum,FrameNum(mm)); 
            SMD.ID = cat(1,SMD.ID,nn);
        end
    end
    Ind = SMD.X_SE < 40 & SMD.Y_SE < 40;
    SMD.X = SMD.X(Ind);
    SMD.Y = SMD.Y(Ind);
    SMD.X_SE = SMD.X_SE(Ind);
    SMD.Y_SE = SMD.Y_SE(Ind);
    SMD.FrameNum = SMD.FrameNum(Ind);
    SMD.ID = SMD.ID(Ind);
end

IndR = SMD.X_SE > Prec_Cutoff | SMD.Y_SE > Prec_Cutoff;
SMD.X(IndR) = [];
SMD.Y(IndR) = [];
SMD.X_SE(IndR) = [];
SMD.Y_SE(IndR) = [];
SMD.FrameNum(IndR) = [];
SMD.ID(IndR) = [];

if PlotFlag
    Ang = 0:0.05:2*pi+0.05;
    for nn = 1:length(SMD.X)
        Prec = 2*sqrt((SMD.X_SE(nn)^2+SMD.Y_SE(nn)^2)/2);
        XCir = Prec*cos(Ang)+SMD.X(nn);
        YCir = Prec*sin(Ang)+SMD.Y(nn);
        if nn == 1
            figure; 
        end
        plot(XCir,YCir,'color',[0 0.4470 0.7410])
        axis equal
        if nn == 1
            hold; 
        end
    end
    plot(SMC.X,SMC.Y,'ok','linewidth',1.5)
end

end


