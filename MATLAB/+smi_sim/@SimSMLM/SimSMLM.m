classdef SimSMLM < handle
    
    % SimSMLM is a Single Molecule Localization Microscopy data generating
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
    
    methods 
        
        function [SMD_True] = SimStar(obj,NWings)
            
            % This function simulates the Siemen star and returns the frames
            % with particles distributed uniformaly on the wings of the star.
            
            % INPUTS:
            
            % NWings: The number of wings of the Siemen's star.
            
            % OUTPUTS:
            
            % [SMD]: This is a structure with the following fields:
            
            % SMD_True.X: The X-positions of particles generated randomly.
            % (Number of the generated particles x 1)(Pixels)
            
            % SMD_True.Y: The Y-positions of particles generated randomly
            % (Number of the generated particles x 1),(Pixels)
            
            % SMD_True.Z: The Z-positions of particles generated randomly
            % (Number of the generated particles x 1), (um)
            
            % [Model]: The frames of the model (Row x Column x NFrames).
            
            % SMD_Model.X: The X-positions of the particles seen on the frames.
            % (Number of the seen particles x 1),(Pixels)
            
            % SMD_Model.Y: The Y-positions of the particles seen on the frames.
            % (Number of the seen particles x 1),(Pixels)
            
            % SMD_Model.Z: The Z-positions of particles (um)
            
            % SMD_Model.Photons: The intensity of the particles.
            % (Number of the seen particles x 1),(Photon counts)
            
            % SMD_Model.FrameNum:The frames that the particles have been detected.
            % (Number of the seen particles x 1)
            
            % SMD_Model.Bg: The background noise (1x1)(Photon counts)
            
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
            
            %Saving the generated data in the structure SMD.
            SMD_True.X = LabelCoords(:,1);
            SMD_True.Y = LabelCoords(:,2);
            if isscalar(obj.PSFSigma)
                SMD_True.Z = [];
            end
        end
        % Call the genBlinks () function to generate the model
        [SMD_Model] = genBlinks(SMD_True,K_OnToOff,K_OffToOn,K_OnToBleach,NFrames,StartState)
    end
    
end

