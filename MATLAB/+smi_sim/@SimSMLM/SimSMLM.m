classdef SimSMLM < handle
    
    % SimSMLM is a Single Molecule Localization Microscopy data generating
    % Class. This class can produce Siemen's star shaped data with particles
    % distributed uniformly throughout the frames. The blinking events are
    % produced based on the given rate parameters (K_OnToOff, K_OffToOn and
    % K_OnToBleach). The data has the same uniform background noise for the
    % whole sequence. The output is corrupted with the Poisson noise.
    
    % Typical data flows are
    %    produce noisy coordinates:
    %       SMD_True -> SMD_Labeled -> SMD_Model -> SMD_Model_Noisy
    %    produce noisy image stacks
    %       SMD_True -> SMD_Labeled -> SMD_Model -> Model -> Data
    % where
    %   SMD_True           true locations of localizations
    %   SMD_Labeled        obj.LabelingEfficiency applied to SMD_True
    %                      localizations, removing localizations that are not
    %                      labeled
    %   SMD_Model          blinks generated for SMD_True_Labeled localizations
    %   SMD_Model_Noisy    SMD_Model with positional and intensity noise added
    %   Model              Gaussian blob image stack produced from SMD_Model
    %   Data               Model image stack to which Poisson noise has been
    %                      applied
    %
    % Model and Data are image stacks (n x n x f), where n is the linear size
    % of the image in pixels and f is the total number of frames to be
    % generated (f = obj.NDatasets * obj.NFrames).

    properties
        SZ=256            % Linear size of image (pixels)
        Rho=30            % Fluorophore Density (fluorophore/pixel)
        NDatasets=1       % Number of datasets
        NFrames=1000      % Number of frames per dataset
        ZoomFactor=20     % It can be either smaller or larger than one- CHANGE
        K_OnToOff=1       % Fluorophore turns Off from On state
                          %    (default:1 frames^-1)
        K_OffToOn=0.0005  % Fluorophore return to On state from Off state
                          %    (default:0.0005 frames^-1)
        K_OnToBleach=0.2  % Fluorophore bleached out (default:1/5 frames^-1)
        EmissionRate=1000 % Emission rate (Intensity) of photons (photons/frame)
        Bg=5              % Background Count Rate (counts/pixel)
        PSFSigma=1.3      % Point Spread Function Sigma size (Pixels)
        LabelingEfficiency=1 % Fluorophore labeling efficiency [range: 0 - 1]
        % A string which determine if the particle starts on or starts randomly
        % on or off.  It can be either 'on' or 'Equib'.
        StartState='Equib'
        Verbose = 1      % Verbosity level

        % Generic note: SMD_* below are SMD structures with various fields
        % filled in as appropriate at that stage.
        SMD_True         % True coordinates produced by sim*
        SMD_Labeled      % True labeled coordinates produced by applyLabelEffic
        SMD_Model        % Coordinates with blinks produced by genBlinks
    end

    methods 
         
        simkTets(obj, kk, radius_kTet)
        applyLabelEffic(obj)
        
    end 

    properties(SetAccess = protected)

    end

    methods(Static)

        SMD_True = kTet(k, center, radius, startAngle)
        unitTest()

    end

end

