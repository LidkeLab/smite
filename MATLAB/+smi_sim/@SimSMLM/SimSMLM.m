classdef SimSMLM < handle
    
    % SimSMLM is a Single Molecule Localization Microscopy data generating
    % Class. This class can produce Siemen's star shaped data with particles
    % distributed uniformly throughout the frames. The blinking events are
    % produced based on the given rate parameters (K_OnToOff, K_OffToOn and
    % K_OnToBleach). The data has the same uniform background noise for the
    % whole sequence. The output is corrupted with the Poisson noise.
    
    % Typical data flows are
%    produce noisy coordinates:
%       SMD_True -> SMD_True_Labeled -> SMD_Model -> SMD_Model_Noisy
%    produce noisy image stacks
%       SMD_True -> SMD_True_Labeled -> SMD_Model -> Model -> Data


    properties
        SZ=256               % Linear size of image (pixels)
        Rho=30              % Fluorophore Density (fluorophore/pixel)
        NDatasets=1      % Number of datasets
        NFrames=1000     % Number of frames per dataset
        ZoomFactor=20    % It can be either smaller or larger than one- change
        K_OnToOff=1      % Fluorophore turns Off from On state
                         %    (default:1 frames^-1)
        K_OffToOn=.0005  % Fluorophore return to On state from Off state
                         %    (default:0.0005 frames^-1)
        K_OnToBleach=.2   % Fluorophore bleached out (default:1/5 frames^-1)
        EmissionRate=1000 % Emission rate (Intensity) of photons (photons/frame)
        Bg=5             % Background Count Rate (counts/pixel)
        PSFSigma=1.3      % Point Spread Function Sigma size (Pixels)
        LabelingEfficiency = 1
        % Fluorophore labeling efficiency [range: 0 - 1]
        
        StartState='Equib'   % A string which determine if the particle starts on
        % or starts randomly on or off.  It can be either 'on'
        % or 'Equib'.
        Verbose = 1      % Verbosity level
        
        SMD_True        
        SMD_Labeled
        SMD_Model
                
        NWings           % The number of wings of the Siemen's star
        OrderkTet        % Order of kTet (the value of k)
        RadiuskTet       % Radius of kTet
        
        
    end

    %properties(SetAccess = protected)
    %    LabelCoords
    %    NLabels
    %end
    
    methods
        
        function [SMD_True, SMD_Model, SMD_Data] = SiemenStar(obj,NWings)
        
            if nargin()<2
                NWings=16;
            end
            
        if nargout == 1
            [SMD_True] = simStar(obj,NWings);
        end
        
        if nargout == 2
             [SMD_True] = simStar(obj,NWings);
             [SMD_Model] = genBlinks(obj,SMD_True,obj.StartState);
        end
         
        if nargout == 3
             [SMD_True] = simStar(obj,NWings);
             [SMD_Model] = genBlinks(obj,SMD_True,obj.StartState);
             [SMD_Data] = genNoisySMD(obj,SMD_Model);
        end
        
        end
    
    end
    
    methods 
         
        [SMD_True] = kTets(obj, kk, radius_kTet)
        [SMD_True_Labeled] = applyLabelEffic(obj, SMD_True)
        
    end 

    methods(Static)

        SMD_True = kTet(k, center, radius, startAngle)
        unitTest()
        SimDimers()

    end

end

