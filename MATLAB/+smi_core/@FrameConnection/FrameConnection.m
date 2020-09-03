classdef FrameConnection < handle
    %FrameConnection contains methods for performing frame connection.
    %
    % REQUIRES:
    %   c_FrameConnect.mex*
    
    % Created by: 
    %   David J. Schodt (Lidke Lab 2020) 
    %       based on a code written Hanieh Mazloom-Farsibaf
    
    
    properties
        LoS = 0.01; % level of significance
        MaxSeparation = 1; % Maximum allowable separation between two 
                           % localizations to still be considered 
                           % candidates for connection.
                           % (Pixels)(Default = 1 pixel)
        MaxFrameGap = 5; % Maximum allowable frame gap between two 
                         % localizations to still be considered 
                         % candidates for connection.
                         % (Pixels)(Default = 5 frames)
        FitType = 'XYNB'; % Fit type used to fit localizations 
                          % (see GaussMLE class for details)
                          % (Character array)(Default = 'XYNB')
        SMD; % Single Molecule Data structure (see 
             % SingleMoleculeData class for details) set by the user.  This
             % is the original, non frame-connected SMD structure.
    end
    
    properties (SetAccess = 'protected')
        SMDCombined; % Single Molecule Data structure (see 
                     % SingleMoleculeData class for details) whose member
                     % localizations have been frame-connected.
                     % SMDCombined is the intended "output" to be prepared
                     % by methods in this class.
    end
    
    methods
        function [obj] = FrameConnection(SMD, SMF)
            %FrameConnection is a constructor for the FrameConnection class
            % This constructor can be used to (optionally) set the input
            % SMD structure and/or "unwrap" an SMF structure to set the 
            % available class properties.
        end
    end
    
    methods (Static)
    end
    
    
end