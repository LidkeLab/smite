classdef FrameConnection < handle
    % FrameConnection performs frame-connection on data in an SMD structure
    % This class contains methods to perform frame-connection and to do
    % associated tasks.  More specifically, the main usage of this class
    % is to combine a time series of localizations arising from the same
    % emitter into a single localization with precision greater than any
    % one of the localizations in the time series.
    %
    % EXAMPLE USAGE:
    %   FC = smi_core.FrameConnection(SMD, SMF);
    %   FC.performFrameConnection();
    %   The primary "outputs" are then FC.SMDCombined and FC.SMD .
    %
    % REQUIRES:
    %   c_FrameConnect.mex*
    
    % Created by:
    %   David J. Schodt (Lidke Lab 2020)
    
    
    properties
        FitType char = 'XYNB'; % (Default = 'XYNB') see GaussMLE class
        LoS(1, 1) double = 0.01; % (Default = 0.01), Level of Significance
        MaxFrameGap(1, 1) uint32 = 5; % (Frames)(Default = 5)
        MaxSeparation(1, 1) double = 1; % (Pixels)(Default = 1)
        SMD struct; % see SingleMoleculeData class
    end
    
    properties (SetAccess = 'protected')
        SMDCombined struct; % frame-connected SMD structure.
    end
    
    methods
        function [obj] = FrameConnection(SMD, SMF)
            % FrameConnection prepares the FrameConnection class for use.
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