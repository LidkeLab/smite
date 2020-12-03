function [SMD_True] = simStar(obj,NWings)
    
    % This function simulates the Siemen star and returns the frames
    % with particles distributed uniformaly on the wings of the star.
    
    % INPUT:
    
    % obj: The object of the SimSMLM() class.
    
    % NWings: The number of wings of the Siemen's star.
    
    % OUTPUTS:
    
    % [SMD_True]: This is a structure with the following fields:
    
    % SMD_True.X: The X-positions of particles generated randomly.
    % (Number of the generated particles x 1)(Pixels)
    
    % SMD_True.Y: The Y-positions of particles generated randomly
    % (Number of the generated particles x 1),(Pixels)
    
    % SMD_True.Z: The Z-positions of particles generated randomly
    % (Number of the generated particles x 1), (um)
    
    
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
    DInd = find(Y==0);
    X(DInd)=[];
    Y(DInd)=[];
    LabelCoords(:,1)=X.*R+obj.SZ/2;
    LabelCoords(:,2)=Y.*R+obj.SZ/2;
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
    SMD_True.X = LabelCoords(:,1)+1;
    SMD_True.Y = LabelCoords(:,2)+1;
    if isscalar(obj.PSFSigma)
        SMD_True.Z = [];
    end

    obj.LabelCoords = LabelCoords;
    obj.NLabels = NLabels;
end
