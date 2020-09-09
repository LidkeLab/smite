classdef FrameConnection < handle
    %FrameConnection contains methods for performing frame connection.
    % This class contains methods to perform frame-connection and to do
    % associated tasks.  More specifically, the main usage of this class
    % is to combine a time series of localizations arising from the same
    % emitter into a single localization with great precision.
    %
    % EXAMPLE USAGE:
    %   FC = smi_core.FrameConnection(SMD, SMF);
    %   FC.performFrameConnection();
    %   The primary "outputs" are then FC.SMDCombined and FC.SMD .
    %
    % REQUIRES:
    %   c_FrameConnect.mex*
    %
    % PROPERTY LIST:
    %   LoS
    %   MaxSeparation
    %   MaxFrameGap
    %   FitType
    %   SMD
    %   SMDCombined
    %
    % METHOD LIST:
    %   findConnected()
    %   performFrameConnection()
    %   unitTest()
    
    % Created by:
    %   David J. Schodt (Lidke Lab 2020)
    
    
    properties
        % Level of Significance
        % (Scalar)(Type double)(Default = 0.01)
        LoS(1, 1) double = 0.01;
        
        % Maximum allowable separation between two localizations to still
        % be considered candidates for connection.
        % (Scalar)(Type double)(Pixels)(Default = 1 pixel)
        MaxSeparation(1, 1) double = 1;
        
        % Maximum allowable frame gap between two localizations to still be
        % considered candidates for connection.
        % (Pixels)(Type uint32)(Default = 5 frames)
        MaxFrameGap(1, 1) uint32 = 5;
        
        % Fit type used to fit localizations (see GaussMLE class for
        % details)
        % (Type char)(Default = 'XYNB')
        FitType char = 'XYNB';
        
        % Single Molecule Data structure (see SingleMoleculeData class for
        % details) set by the user.  This is the, non frame-connected SMD
        % structure.  After frame-connection, an additional field
        % (SMD.ConnectID) which indicates which localizations were
        % connected together.
        SMD struct;
    end
    
    properties (SetAccess = 'protected')
        % Single Molecule Data structure (see SingleMoleculeData class for
        % details) whose member localizations have been frame-connected.
        % SMDCombined is the intended "output" to be prepared by methods in
        % this class.
        SMDCombined struct;
    end
    
    methods
        function [obj] = FrameConnection(SMD, SMF)
            %FrameConnection is a constructor for the FrameConnection class
            % This constructor can be used to (optionally) set the input
            % SMD structure and/or "unwrap" an SMF structure to set the
            % available class properties.
            
            % Set the input SMD structure as a class property (if provided)
            if (exist('SMD', 'var') && ~isempty(SMD))
                obj.SMD = SMD;
            end
            
            % "Unwrap" the SMF structure to see if it contains any of the
            % class properties that we'll need.
            if (exist('SMF', 'var') && ~isempty(SMF))
                FrameConnectFields = {'MaxSeparation', ...
                    'MaxFrameGap', ...
                    'LoS'};
                FittingFields = {'FitType'};
                for ff = 1:numel(FrameConnectFields)
                    % Check if this field exists in SMF.FrameConnect, and
                    % if it does, copy it into the appropriate class
                    % property.
                    if isfield(SMF.FrameConnection, FrameConnectFields{ff})
                        obj.(FrameConnectFields{ff}) = ...
                            SMF.FrameConnection.(FrameConnectFields{ff});
                    end
                end
                for ff = 1:numel(FittingFields)
                    % Check if this field exists in SMF.Fitting, and if it
                    % does, copy it into the appropriate class property.
                    if isfield(SMF.Fitting, FittingFields{ff})
                        obj.(FittingFields{ff}) = ...
                            SMF.Fitting.(FittingFields{ff});
                    end
                end
            end
        end
        
        performFrameConnection(obj)
    end
    
    methods (Static)
        [SMDIndex] = findConnected(SMR, SMD, ID);
        [Success] = unitTest();
    end
    
    
end