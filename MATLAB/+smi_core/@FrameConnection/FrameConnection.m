classdef FrameConnection < handle
    %FrameConnection performs frame-connection on data in an SMD structure
    %
    % This class contains methods to perform frame-connection and to do
    % associated tasks.  More specifically, the main usage of this class
    % is to combine a time series of localizations arising from the same
    % emitter into a single localization with precision greater than any
    % one of the localizations in the time series.
    %
    % EXAMPLE USAGE:
    %   FC = smi_core.FrameConnection(SMD, SMF);
    %   [SMDCombined, SMD] = FC.performFrameConnection();
    %   Alternatively, you can use this class as a "function" as follows:
    %   [~, SMD, SMDCombined] = smi_core.FrameConnection(SMD, SMF, 1);
    %
    % REQUIRES:
    %   smi_c_FrameConnection.mex*
    
    % Created by:
    %   David J. Schodt (Lidke Lab 2020)
    
    
    properties
        BoxSize(1, 1) double = 7; % (Pixels)(Default = 7) see GaussMLE
        FitType = 'XYNB'; % (Default = 'XYNB') see GaussMLE class
        LoS(1, 1) double = 0.01; % (Default = 0.01), Level of Significance
        MaxFrameGap(1, 1) uint32 = 5; % (Frames)(Default = 5)
        MaxSeparation(1, 1) double = 1; % (Pixels)(Default = 1)
        NParams(1, 1) double = 4; % (Default = 4) see GaussMLE
        SMD % see SingleMoleculeData class
    end
    
    properties (SetAccess = 'protected')
        SMDCombined % frame-connected SMD structure.
    end
    
    methods
        function [obj, SMDCombined, SMD] = FrameConnection(SMD, SMF, ...
                AutoRun)
            %FrameConnection prepares the FrameConnection class for use.
            % This constructor can be used to (optionally) set the input
            % SMD structure and/or "unwrap" an SMF structure to set the
            % available class properties.
            
            % Set a default for the AutoRun flag, which specifies whether
            % or not we should attempt to perform frame-connection
            % immediately (i.e., perform frame-connection in this
            % constructor).
            if (~exist('AutoRun', 'var') || isempty(AutoRun))
                AutoRun = 0;
            end            
            
            % Set the input SMD structure as a class property (if provided)
            if (exist('SMD', 'var') && ~isempty(SMD))
                obj.SMD = SMD;
                AllFieldsSet = 1;
            else
                AllFieldsSet = 0;
            end
            
            % Set class properties based on the input SMF structure (if it
            % was provided).
            if (exist('SMF', 'var') && ~isempty(SMF))
                % Pad the input SMF structure to ensure it contains all
                % fields defined in SingleMoleculeFitting.createSMF().
                SMF = smi_core.SingleMoleculeFitting.padSMF(SMF);
                
                % Set the desired SMF fields to class properties.
                obj.BoxSize = SMF.BoxFinding.BoxSize;
               	obj.FitType = SMF.Fitting.FitType;
                obj.NParams = SMF.Fitting.NParams;
                obj.LoS = SMF.FrameConnection.LoS;
                obj.MaxFrameGap = SMF.FrameConnection.MaxFrameGap;
                obj.MaxSeparation = SMF.FrameConnection.MaxSeparation;
            else
                AllFieldsSet = 0;
            end
            
            % If all required fields are set, run the frame-connection
            % process.
            if (AutoRun && AllFieldsSet)
                [SMDCombined, SMD] = obj.performFrameConnection();
            else
                if (nargout() > 1)
                    warning(['Constructor outputs SMD and ', ...
                        'SMDCombined were requested but ', ...
                        'frame-connection was not performed.'])
                end
                SMD = obj.SMD;
                SMDCombined = SMD;
            end
        end
        
        [SMDCombined, SMD, OutputMessage] = performFrameConnection(obj)
        
    end
    
    methods (Static)
        [SMDIndex] = findConnected(SMR, SMD, ID);
        [Success] = unitTest();
    end
    
    
end
