classdef SimSMLM < handle
    
    % Sim_SMLM is a Single Molecule Localization Microscopy data generating
    % Class. This class can produce Siemen's star shaped data with particles
    % distributed uniformly throughout the frames. The blinking events are
    % produced based on the given rate parameters (K_OnToOff, K_OffToOn and
    % K_OnToBleach). The data has the same uniform background noise for the
    % whole sequence. The output is corrupted with the Poisson noise.
    
    
    properties
        SZ               % Linear size of image (pixels)
        Rho              % Fluorophore Density (flourophore/pixel)
        NFrames          % Number of frames
        ZoomFactor       % It can be either smaller or larger than one
        K_OnToOff        % Fluorophore turns Off from On state (1 frames^-1)
        K_OffToOn        % Fluorophore return to On state from Off state (0.0005 frames^-1)
        K_OnToBleach     % Fluorophore bleached out (1/5 frames^-1)
        EmissionRate     % Emission rate (Intensity) of photons (photons/frame)
        Bg               % Background Count Rate (counts/pixel)
        PSFSigma         % Point Spread Function Sigma size (Pixels).
        
    end
    
    methods (Static)
        
        function [Data,SMD,Model] = SimStar(obj,NWings,StartState)
            
            % This function simulates the Siemen star and returns the frames
            % with particles distributed uniformaly on the wings of the star.
            
            % INPUTS:
            
            % NWings: The number of wings of the Siemen's star.
            
            % StartState:  A string which determine if the particle starts
            % on or starts randomly on or off. It can be either 'On' or 'Equib'.
            
            % OUTPUTS:
            
            % [Data]: The frames of the data (Row x Column x NFrames).
            
            % [SMD]: This is a structure with the following fields:
            % SMD.XTrue: The X-positions of particles generated randomly.
            % (Number of the generated particles x 1)(Pixels)
            
            % SMD.YTrue: The Y-positions of particles generated randomly
            % (Number of the generated particles x 1),(Pixels)
            
            % SMD.ZTrue: The Z-positions of particles generated randomly
            % (Number of the generated particles x 1), (um)
            
            % SMD.X: The X-positions of the particles seen on the frames.
            % (Number of the seen particles x 1),(Pixels)
            
            % SMD.Y: The Y-positions of the particles seen on the frames.
            % (Number of the seen particles x 1),(Pixels)
            
            % SMD.Z: The Z-positions of particles (um)
            
            % SMD.Photons: The intensity of the particles.
            % (Number of the seen particles x 1),(Photon counts)
            
            % SMD.FrameNum:The frames that the particles have been detected.
            % (Number of the seen particles x 1)
            
            % SMD.Bg: The background noise (1x1)(Photon counts)
            
            % [Model]: The frames of the model (Row x Column x NFrames).
            
            
            R = obj.SZ/3; %The length of the wing is a third of the size of the frame
            
            % To facilitate the calculations, we will work in both polar and
            % Cartesian coordinates. The random particles will be generated
            % throughout a cycle and then we remove those that are not
            % inside any wing.
            
            Alpha = -pi:2*pi/NWings:pi; %The angle that each wing starts and ends.
            Nn=poissrnd(floor(pi*R^2*obj.Rho)); % Nn are the number of emitters.
            X = 2*(rand(Nn,1)-0.5);
            Y = 2*(rand(Nn,1)-0.5);
            IndZero=find(X.^2 + Y.^2>1);
            X(IndZero)=0;
            Y(IndZero)=0;
            [Theta,~] = cart2pol(X,Y);
            
            for ij = 1:2:NWings
                % i and j represent the ith particle in the jth frame.
                BInd = find(Alpha(ij)<Theta & Theta<Alpha(ij+1));
                X(BInd)=0;
                Y(BInd)=0;
            end
            
            % Delete those particles from the list that are not inside any wing.
            DInd = find(X==0);
            X(DInd)=[];
            Y(DInd)=[];
            LabelCoords(:,1)=X*R+obj.SZ/2;
            LabelCoords(:,2)=Y*R+obj.SZ/2;
            LabelCoords = LabelCoords*obj.ZoomFactor;
            obj.SZ = obj.SZ*obj.ZoomFactor;
            NLabels=length(LabelCoords); % Number of the generated particles.
            IntArray=zeros(NLabels,obj.NFrames); % This is a 2D-array with
            % the size of(Number of the particles)x(Number of the frames)
            % to store the trace of the blinking events. The elements of
            % this array can be either zero or one. One signifies an on event.
            % Let say the element on the ith row and jth column is one. This
            % signifies the ith particle is ON in the jth frame.
            Photons=[]; %Making emepty vectors that will be filled later.
            X=[];
            Y=[];
            Z=[];
            FrameNum = [];
            % The following loop iterates over each particle to generate
            % the blinking events for them.
            for mm=1:NLabels
                %genBlinks() makes the blinking events.
                Temp=genBlinks(obj,StartState);
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
                    Indiv(:,1)=LabelCoords(mm,1);
                    X = cat(1,X,Indiv);
                    Indiv(:,1)=LabelCoords(mm,2);
                    Y = cat(1,Y,Indiv);
                end
            end
            %Saving the generated data in the structure SMD.
            SMD.XTrue = LabelCoords(:,1)+1;
            SMD.YTrue = LabelCoords(:,2)+1;
            if isscalar(obj.PSFSigma)
                SMD.ZTrue = [];
            end
            SMD.FrameNum = FrameNum;
            SMD.Photons = Photons;
            SMD.X = X+1;
            SMD.Y = Y+1;
            if isscalar(obj.PSFSigma)
                SMD.Z = [];
                SMD.PSFSigma = obj.PSFSigma*ones([length(Photons),1]);
            end
            
            SMD.Bg = 0;
            NoiseIm = obj.Bg*ones(obj.SZ);
            
            %Calling SMA_Sim.gaussBlobImage() to generate the blobs.
            [Model] = smi_sim.gaussBlobImage(obj.SZ,obj.NFrames,SMD,obj.Bg,0,0);
            %Add background and noise
            [Model] = Model+obj.Bg;
            [Data] = poissrnd(single(Model));
            NoiseIm = zeros(size(Data(:,:,1)));
            [Data] = Data+randn(size(Data)).*repmat(NoiseIm,[1 1 obj.NFrames]);
            SMD.Bg = obj.Bg;
            
            %Nested function to generate blinking events.
            function IvsT=genBlinks(obj,StartState)
                
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
                        State=rand < (obj.K_OffToOn/(obj.K_OffToOn+obj.K_OnToOff));
                        %Randomly determine if the particle starts off or on.
                        if State %starts on
                            T=0;
                        else %starts from off
                            T=exprnd(1/obj.K_OffToOn); %The random time when the particle goes on.
                        end
                end
                
                %The following while-loop gives the times when the particle
                %goes on, goes off and photobleach.
                
                while T<obj.NFrames
                    NPairs=NPairs+1;
                    % OnOffPairs is an array of 3 columns, where the first column gives
                    % the times when the particle goes on, the second column gives the
                    % time when the particle goes off and third column gives the time
                    % when particle photobleaches.
                    OnOffPairs(NPairs,1)=T;
                    D=exprnd(1/obj.K_OnToOff); %Generate blink duratrion
                    OnOffPairs(NPairs,2)=min(T+D,obj.NFrames); %The On-time plus the duration gives the off-time.
                    OnOffPairs(NPairs,3)= rand > (obj.K_OnToOff/(obj.K_OnToOff+obj.K_OnToBleach)); %fluorophore bleaches
                    %if this condition met.
                    T=T+D+exprnd(1/obj.K_OffToOn); %Advance to the next blinking event.
                end
                
                %Turn blinking events to I vs T
                IvsT=zeros(obj.NFrames,1);
                for nn=1:NPairs
                    StartT=floor(OnOffPairs(nn,1))+1; %index for start frame
                    EndT=min(obj.NFrames,floor(OnOffPairs(nn,2))+1); %index for start frame
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
                
            end
            
        end
    end
    
end

