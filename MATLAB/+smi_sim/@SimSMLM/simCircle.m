function simCircle(obj, Radius)
    % This function simulates a Circle for SMLM simulation.
    %
    % INPUT:
    % obj: The object of the SimSMLM() class.
    % Radius: Radius of the circle in normalized units (default: 0.5).
    
    % Default parameters if not provided
    if nargin < 2
        Radius = 0.5; % Default radius
    end
    
    % Reference scaling factor
    R = obj.SZ / 3; 
    
    % Generate random points in a normalized [-1, 1] space
    Nn = poissrnd(floor(pi * R^2 * obj.Rho)); % Number of emitters
    X = 2 * (rand(Nn, 1) - 0.5); % Random X positions
    Y = 2 * (rand(Nn, 1) - 0.5); % Random Y positions
    
    % Initialize mask for the circle
    Mask = (X.^2 + Y.^2 <= Radius^2 & X.^2 + Y.^2 >= (Radius - 0.05)^2);
    
    % Filter points within the circle's edge
    X = X(Mask);
    Y = Y(Mask);
    
    % Rescale to canvas size
    LabelCoords(:, 1) = X .* R + obj.SZ / 2;
    LabelCoords(:, 2) = Y .* R + obj.SZ / 2;
    
    % Save the generated data into the SMD structure
    obj.SMD_True = smi_core.SingleMoleculeData.createSMD();
    obj.SMD_True.X = LabelCoords(:, 1);
    obj.SMD_True.Y = LabelCoords(:, 2);
    
    if isscalar(obj.PSFSigma)
        obj.SMD_True.Z = [];
    end
    
    obj.SMD_True.XSize = obj.SZ;
    obj.SMD_True.YSize = obj.SZ;

    % Generate the simulation model
    obj.genModel();
end
% simCircle(obj, 0.4); % Simulate a circle with radius 0.6