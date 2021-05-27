function genBlinks(obj,StartState)

%This function generates the blinking time traces for a single particle 
%over the given number of the frames considering the input parameters and
%returns the SMD_Model structure. 

 % INPUTS:
 
 % obj: The object of the SimSMLM() class.
 
 % [obj.SMD_Labeled]: This is a structure with the following three fields:
 
 % obj.SMD_Labeled.X: The X-positions of particles generated randomly.
 % (Number of the generated particles x 1)(Pixels)
    
 % obj.SMD_Labeled.Y: The Y-positions of particles generated randomly
 % (Number of the generated particles x 1),(Pixels)
    
 % obj.SMD_Labeled.Z: The Z-positions of particles generated randomly
 % (Number of the generated particles x 1), (um)

 % obj.SparseFlag: If true, turns on using sparse matrices in genBlinks,
 % allowing it to accomodate larger examples.  However, sparse matrices
 % are also slower to manipulate than regular ones especially for
 % smaller examples, so shouldn't be used unless needed, thus the
 % default is false.
 
 % StartState: A string which determine if the particle starts on or
 % starts randomly on or off. It can be either 'on' or 'Equib'.
 
 % OUTPUTS:
 
 % [obj.SMD_Model]: This is a structure with the following fields:
 
 % obj.SMD_Model.X: The X-positions of the particles seen on the frames.
 % (Number of the seen particles x 1),(Pixels)
    
 % obj.SMD_Model.Y: The Y-positions of the particles seen on the frames.
 % (Number of the seen particles x 1),(Pixels)
    
 % obj.SMD_Model.Z: The Z-positions of particles (um)
    
 % obj.SMD_Model.Photons: The intensity of the particles.
 % (Number of the seen particles x 1),(Photon counts)
    
 % obj.SMD_Model.FrameNum:The frames that the particles have been detected.
 % (Number of the seen particles x 1)

 % obj.SMD_Model.NFrames:The number of frames that the particles have been detected.

 % obj.SMD_Model.DatasetNum:The dataset that corresponds to the FrameNum.
 % (Number of the seen particles x 1)

 % obj.SMD_Model.NDatasets: The number of datasets in which to organize the frames.
 
 % obj.SMD_Model.PSFSigma: Point Spread Function Sigma size (Pixels)
 
 % obj.SMD_Model.Bg: Background Count Rate (counts/pixel), this will be empty.

% First making empty vectors that will be filled later.
NLabels = numel(obj.SMD_Labeled.X);
Photons=[];
X=[];
Y=[];
Z=[];
FrameNum=[];
PSFSigma=[];
Bg=[];
TotalNFrames = obj.NDatasets*obj.NFrames;
NOnEvents = zeros(NLabels, 1, 'uint8');   % number of on events

% obj.SMD_Model.ConnectID: SMD_Labeled -> SMD_Model
% ConnectID produces an indexing array that associates each localization in
% SMD_Model with the localization that produced it in SMD_Labeled.  That is,
%
%    obj.SMD_Labeled.X(obj.SMD_Model.ConnectID) == obj.SMD_Model.X
%
% where SMD_Model.ConnectID has the same dimensions as SMD_Model.X/Y/etc.  For
% example,
%
%    obj.SMD_Model.X(k) == obj.SMD_Labeled.X(obj.SMD_Model.ConnectID(k))
%
% for k = 1 : numel(obj.SMD_Model.X).
ConnectID = [];

if obj.SparseFlag
   IntArray = sparse(NLabels, TotalNFrames);
else
   IntArray = zeros(NLabels, TotalNFrames);
end

%The following loop iterates over each particle to generate the blinking
%events for them.

for mm=1:NLabels
    [Temp, NTime] = blinks(obj.K_OnToOff, obj.K_OffToOn, obj.K_OnToBleach, ...
                           TotalNFrames, StartState);
    NOnEvents(mm) = NTime;
    
    %Blinks() makes the blinking events. It takes the following inputs:
    
    %K_OnToOff: Fluorophore turns Off from On state (default:1 frames^-1)
 
    %K_OffToOn: Fluorophore return to On state from Off state (default:0.0005 frames^-1)
 
    %K_OnToBleach: Fluorophore bleaches (default:1/5 frames^-1)
 
    %TotalNFrames: Total number of frames (pixels x pixels)
    
    %StartState: A string which determine if the particle starts on or
    %starts randomly on or off. It can be either 'on' or 'Equib'.
    
    IntArray(mm,:)=Temp';
    %Finding the frames where the particle was ON. Note that the
    %particle might not be on at all and this would be an empty
    %array. In this case, we won't have any FrameNum, Photons or
    %found positions for this particle.
    FrameNumIndiv = find(Temp~=0);
    if ~isempty(FrameNumIndiv) 
        FrameNum = cat(1,FrameNum,FrameNumIndiv);
        Indiv = obj.EmissionRate*Temp(FrameNumIndiv);
        Photons = cat(1,Photons,Indiv);
        %Indiv(:,1)=obj.LabelCoords(mm,1);
        Indiv(:,1)=obj.SMD_Labeled.X(mm);
        X = cat(1,X,Indiv);
        %Indiv(:,1)=obj.LabelCoords(mm,2);
        Indiv(:,1)=obj.SMD_Labeled.Y(mm);
        Y = cat(1,Y,Indiv);
        nF = numel(FrameNumIndiv);
        ConnectID = [ConnectID; repmat(mm, nF, 1)];
    end
end
obj.SMD_Model = obj.SMD_Labeled;
obj.SMD_Model.X = X;
obj.SMD_Model.Y = Y;
if isscalar (obj.PSFSigma) 
    obj.SMD_Model.Z = [];
    obj.SMD_Model.PSFSigma = obj.PSFSigma*ones([length(Photons),1]);
end
obj.SMD_Model.Photons    = Photons;
obj.SMD_Model.Bg         = 0;
obj.SMD_Model.NDatasets  = obj.NDatasets;
obj.SMD_Model.NFrames    = obj.NFrames;
AbsoluteFrameNum = FrameNum;
obj.SMD_Model.DatasetNum = zeros(size(AbsoluteFrameNum));
obj.SMD_Model.FrameNum   = zeros(size(AbsoluteFrameNum));
% Convert absolute frame numbers to per dataset frame numbers.
lo = 1;
for i = 1 : obj.NDatasets
   hi = lo + obj.NFrames - 1;
   indx = find(lo <= AbsoluteFrameNum & AbsoluteFrameNum <= hi);
   obj.SMD_Model.DatasetNum(indx) = i;
   obj.SMD_Model.FrameNum(indx) = AbsoluteFrameNum(indx) - lo + 1;
   lo = lo + obj.NFrames;
end
obj.SMD_Model.ConnectID = ConnectID;
obj.NOnEvents = NOnEvents;
    
    %Nested function to generate blinking events.
    function [IvsT,NTime]=...
       blinks(K_OnToOff,K_OffToOn,K_OnToBleach,NFrames,StartState)
    %Blinks() generates blinking time trace for a single
    %particle over the given number of the frames considering
    %the parameters K_OffToOn, K_OnToOff and K_OnToBleach.
    
    NTime=0; %Number of the times that the particle goes on.
    
    %Based on the input 'StartState' the particle can start in the
    %on-state or it can be started randomly in either on-state or
    %off-state.
    switch StartState
        case 'On'
            T=0; %The time when the particle goes on. T=0 implies
            %that the particle is on from the beginning.
        case 'Equib'
            %Find start state:
            State=rand < (K_OffToOn/(K_OffToOn+K_OnToOff));
            %Randomly determine if the particle starts off or on.
            if State %starts on
                T=0;
            else %starts from off
                T=exprnd(1/K_OffToOn); %The random time when the particle goes on.
            end
    end
    
    %The following while-loop gives the times when the particle
    %goes on, goes off and photobleach.
    
    while T<NFrames
        NTime=NTime+1;
        % TRate (Time rate) is an array of 3 columns, where the first column gives
        % the times when the particle goes on, the second column gives the
        % time when the particle goes off and third column gives the time
        % when particle photobleaches.
        TRate(NTime,1)=T;
        D=exprnd(1/(K_OnToOff+K_OnToBleach)); %Generate blink duratrion
        TRate(NTime,2)=min(T+D,NFrames); %The On-time plus the duration gives the off-time.
        if rand() > (K_OnToOff/(K_OnToOff+K_OnToBleach)) %fluorophore bleaches
           TRate(NTime,3)=rand();
           break;
        end
        %if this condition is met.
        T=T+D+exprnd(1/K_OffToOn); %Advance to the next blinking event.
    end
    
    %Turn blinking events to I vs T
    IvsT=zeros(NFrames,1);
    for nn=1:NTime
        StartT=floor(TRate(nn,1))+1; %index for start frame
        EndT=min(NFrames,floor(TRate(nn,2))+1); %index for start frame
        if StartT==EndT %Blinking happens within one frame
            IvsT(StartT)=TRate(nn,2)-TRate(nn,1);
        else
            %This for-loop goes over the frames where the particle is on.
            for ii=StartT:EndT
                if ii==StartT
                    IvsT(ii)=StartT-TRate(nn,1);
                elseif ii==EndT
                    IvsT(ii)=1-(EndT-TRate(nn,2));
                else
                    IvsT(ii)=1;
                end
            end
        end
    end

    end % Blinks

end % genBlinks
