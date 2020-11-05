function [SMD_Model] = genBlinks(obj,SMD_True,K_OnToOff,K_OffToOn,K_OnToBleach,NFrames,StartState)
%Making empty vectors that will be filled later.
Photons=[];
X=[];
Y=[];
Z=[];
FrameNum=[];
Bg=[];


% The following loop iterates over each particle to generate the blinking
% events for them.

for mm=1:obj.NLabels
    %genBlinks() makes the blinking events.
    Temp=Blinks(K_OnToOff,K_OffToOn,K_OnToBleach,NFrames,StartState);
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
        Indiv(:,1)=obj.LabelCoords(mm,1);
        X = cat(1,X,Indiv);
        Indiv(:,1)=obj.LabelCoords(mm,2);
        Y = cat(1,Y,Indiv);
    end
    SMD_Model.X = X;
    SMD_Model.Y = Y;
    SMD_Model.Z = [];
    SMD_Model.FrameNum = FrameNum;
    SMD_Model.Photons = Photons;
    SMD_Model.Bg = 0;

    
    % NoiseIm = Bg*ones(SZ); Noise factor will be included later in Data. 
end
    
    %Nested function to generate blinking events.
    function IvsT=Blinks(K_OnToOff,K_OffToOn,K_OnToBleach,NFrames,StartState)
    %genBlinks() generates blinking time trace for a single
    %particle over the given number of the frames considering
    %the parameters K_OffToOn, K_OnToOff and K_OnToBleach.
    
    NPairs=0; %Number of the times that the particle goes on.
    
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
        NPairs=NPairs+1;
        % OnOffPairs is an array of 3 columns, where the first column gives
        % the times when the particle goes on, the second column gives the
        % time when the particle goes off and third column gives the time
        % when particle photobleaches.
        OnOffPairs(NPairs,1)=T;
        D=exprnd(1/K_OnToOff); %Generate blink duratrion
        OnOffPairs(NPairs,2)=min(T+D,NFrames); %The On-time plus the duration gives the off-time.
        OnOffPairs(NPairs,3)= rand > (K_OnToOff/(K_OnToOff+K_OnToBleach)); %fluorophore bleaches
        %if this is condition met.
        T=T+D+exprnd(1/K_OffToOn); %Advance to the next blinking event.
    end
    
    %Turn blinking events to I vs T
    IvsT=zeros(NFrames,1);
    for nn=1:NPairs
        StartT=floor(OnOffPairs(nn,1))+1; %index for start frame
        EndT=min(NFrames,floor(OnOffPairs(nn,2))+1); %index for start frame
        if StartT==EndT %Blinking happens within one frame
            IvsT(StartT)=OnOffPairs(nn,2)-OnOffPairs(nn,1);
        else
            %This for-loop goes over the frames where the particle is on.
            for ii=StartT:EndT
                if ii==StartT
                    IvsT(ii)=StartT-OnOffPairs(nn,1);
                elseif ii==EndT
                    IvsT(ii)=1-(EndT-OnOffPairs(nn,2));
                else
                    IvsT(ii)=1;
                end
            end
        end
    end

    end % Blinks
end % genBlinks

%Calling SimSMLM.gaussBlobImage() to generate the blobs.
%[Model] = SimSMLM.gaussBlobImage(SZ,NFrames,SMD_Model,Bg,0,0);
