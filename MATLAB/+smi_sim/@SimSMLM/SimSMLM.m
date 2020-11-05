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

    properties(SetAccess = protected)
        LabelCoords
        NLabels
    end
    
    methods 
        
        [SMD_True] = simStar(obj,NWings)
        % Call the genBlinks() function to generate the model
        [SMD_Model] = genBlinks(obj,SMD_True,K_OnToOff,K_OffToOn,K_OnToBleach,NFrames,StartState)

    end 

    
end

