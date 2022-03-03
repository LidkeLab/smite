function [SMD,SMD_combined]=frameConnect(SMD,LOS,MaxDistance,MaxFrameGap,FitType)
%frameConnect Connects coordinates from a blinking event across consecutive frames. 
% [SMD,SMD_combined]=BaGoL.frameConnect(SMD,LOS,MaxDistance,MaxFrameGap,FitType)
%
% INPUTS:
%    SMD:            Single Molecule Data (SMD). It is a structure including
%                    all localization properties such as X,Y,Photons,...
%    LOS:            Level of Significance (Default = 0.01)
%    MaxDistance:    Maximum distance between two localizations for being
%                    connected (Pixels) (Default = 1)
%    MaxFrameGap:    Maximum Frame Gap between two localizations for being
%                    connected (Frame) (Default = 4)
%    FitType:        'Gaussian2D', 'Gaussian3D', 'SinglePSF', or 'DoublePSF',
%                    (Default = 'Gaussian2D')
%
% OUTPUTS:
%    SMD:            Single Molecule Data; this function creates one more field
%                    to this structure: SMD.ConnectID
%    SMD_combined:   Single Molecule Data of connected localization.
%                    It includes only these fields (X,Y,X_SE,Y_SE,NCombined,
%                    FrameNum,Photons,Bg)
%
% REQUIRES:
%     Statistics Toolbox
%     mex file of c_code: c_FrameConnect on MATLAB's path or current directory
%

% Created by:
%    Hanieh Mazloom-Farsibaf (LidkeLab 2017)

%check input
if nargin<1
    error('You have to enter input Single Molecule Data (SMD) structure')
end
if nargin<2
    LOS=0.01;
end
if nargin <3
    MaxDistance=1;
end
if nargin<4
    MaxFrameGap=4;
end
if nargin<5
    FitType='GaussianBasic';
end
if isempty(SMD)
    error('Input Structure should include at least X,Y,I,Bg for all localizations')
end

%set the parameters for c_FrameConnect
LOS=single(LOS);

%Label the localizations based on frame Connection function
%0 if they didn't pass Threshold and ~=0 if they passed. Its values
%are an integer number for each cluster of localizations combined as
%a higher precesion.
ConnectID=zeros(size(SMD.X,1),1);
SMD.ConnectID=ConnectID;
MaxConID=0;
if ~isfield(SMD,'DatasetNum')
    SMD.DatasetNum = ones(size(SMD.X));
end
if ~isfield(SMD,'ThreshFlag')
    SMD.ThreshFlag = zeros(size(SMD.X));
end
if ~isfield(SMD,'Bg')
    SMD.Bg = zeros(size(SMD.X));
end
if ~isfield(SMD,'LogL')
    SMD.LogL = ones(size(SMD.X));
end

%create output for c_FrameConnect
Coordout_all=[];  % coordinate of combined localization (x,y,or z)
CoordSEout_all=[];
Fout_all=[];     % last frame number which emitter is combined
Nout_all=[];     % number of localizations combined
P_Aveout_all=[];
P_ave_SEout_all=[];
P_Addout_all=[];

% start Frame Connection loop for each sequence with a bunch of frames
SMD_combined.DatasetNum=[];
for ii=1:max(SMD.DatasetNum) %loop over data sets (FileName is the number of sequence.)
    
    mask = (SMD.DatasetNum==ii & SMD.ThreshFlag==0); % extract one sequence from whole saved image
    
    if sum(mask)==0;continue;end
    
    X=SMD.X(mask);
    Y=SMD.Y(mask);
    Photons=single(SMD.Photons(mask));
    Bg=SMD.Bg(mask);
    X_SE=SMD.X_SE(mask);
    Y_SE=SMD.Y_SE(mask);
    LogL=single(SMD.LogL(mask));
    FrameNum=uint32(SMD.FrameNum(mask));
    Bg=single(Bg);
    
    switch FitType
        case {'GaussianBasic'}
            Coord_in=single(cat(2,X,Y));
            Coord_SE_in=single(cat(2,X_SE,Y_SE));
            P_Ave_in=[];
            P_Ave_SE_in=[];
            P_Add_in=single(cat(2,Photons,Bg,LogL));
        case {'GaussianZ'}
            Z=SMD.Z(mask);
            Z_SE=SMD.Z_SE(mask);
            Coord_in=single(cat(2,X,Y,Z));
            Coord_SE_in=single(cat(2,X_SE,Y_SE,Z_SE));
            P_Ave_in=[];
            P_Ave_SE_in=[];
            P_Add_in=single(cat(2,Photons,Bg,LogL));
        case {'GaussianSigma'}
            Coord_in=single(cat(2,X,Y));
            Coord_SE_in=single(cat(2,X_SE,Y_SE));
            P_Ave_in=single(SMD.PSFSigma(mask));
            P_Ave_SE_in=single(SMD.PSFSigma_SE(mask).^2);
            P_Add_in=single(cat(2,Photons,Bg,LogL));
        case {'GaussianSigmaXY'}
            Sx=SMD.PSFSigmaX(mask);
            Sx_SE=SMD.PSFSigmaX_SE(mask);
            Sy=SMD.PSFSigmaY(mask);
            Sy_SE=SMD.PSFSigmaY_SE(mask);
            Coord_in=single(cat(2,X,Y));
            Coord_SE_in=single(cat(2,X_SE,Y_SE));
            P_Ave_in=single(cat(2,Sx,Sy));
            P_Ave_SE_in=single(cat(2,Sx_SE.^2,Sy_SE.^2));
            P_Add_in=single(cat(2,Photons,Bg,LogL));
         case {'SinglePSF'}
            Z=SMD.Z(mask);
            Z_SE=SMD.Z_SE(mask);
            Coord_in=single(cat(2,X,Y,Z));
            Coord_SE_in=single(cat(2,X_SE,Y_SE,Z_SE));
            P_Ave_in=[];
            P_Ave_SE_in=[];
            P_Add_in=single(cat(2,Photons,Bg,LogL));
        otherwise
            error('No FitType Given');
    end
    
    % test frame connection hypothesis histogram (using c-function)
    [Coordout, CoordSEout, Nout, Fout,P_Aveout,P_ave_SEout,P_Addout,ConnectID]=c_FrameConnect_BaGoL(LOS,Coord_in,Coord_SE_in,FrameNum,P_Ave_in,P_Ave_SE_in,P_Add_in,MaxDistance,MaxFrameGap,MaxConID);
    
    SMD.ConnectID(mask)=ConnectID;
    Coordout_all=cat(1,Coordout_all,Coordout);
    CoordSEout_all=cat(1,CoordSEout_all,CoordSEout);
    Nout_all=cat(1,Nout_all,Nout);
    Fout_all=cat(1,Fout_all,Fout);
    if P_Aveout~=0
        P_Aveout_all=cat(1,P_Aveout_all,P_Aveout);
        P_ave_SEout_all=cat(1,P_ave_SEout_all,P_ave_SEout);
    end
    P_Addout_all=cat(1,P_Addout_all,P_Addout);
    MaxConID=max(ConnectID);
    SMD_combined.DatasetNum=cat(1,SMD_combined.DatasetNum,ii*ones(length(Nout), 1));
end

SMD_combined.X=Coordout_all(:,1);
SMD_combined.Y=Coordout_all(:,2);
SMD_combined.X_SE=(CoordSEout_all(:,1));
SMD_combined.Y_SE=(CoordSEout_all(:,2));
SMD_combined.NCombined=Nout_all;
SMD_combined.FrameNum=Fout_all;
SMD_combined.Photons=P_Addout_all(:,1);
SMD_combined.Bg=P_Addout_all(:,2);
SMD_combined.LogL=P_Addout_all(:,3);

if (strcmp(FitType,'GaussianZ') ||strcmp(FitType,'SinglePSF') )
    SMD_combined.Z=Coordout_all(:,3);
    SMD_combined.Z_SE=sqrt(CoordSEout_all(:,3));
end


if (strcmp(FitType,'GaussianSigma'))
    SMD_combined.PSFSigma=P_Aveout_all;
    SMD_combined.PSFSigma_SE=P_ave_SEout_all;
end
if (strcmp(FitType,'GaussianSigmaXY'))
    SMD_combined.PSFSigmaX=P_Aveout_all(:,1);
    SMD_combined.PSFSigmaY=P_Aveout_all(:,2);
    SMD_combined.PSFSigmaX_SE=P_ave_SEout_all(:,1);
    SMD_combined.PSFSigmaY_SE=P_ave_SEout_all(:,2);
end

end
