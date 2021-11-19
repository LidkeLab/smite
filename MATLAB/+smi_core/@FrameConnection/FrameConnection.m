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
    %   [~, SMDCombined, SMD] = smi_core.FrameConnection(SMD, SMF, 1);
    %
    % REQUIRES:
    %   smi_c_FrameConnection.mex*
    
    % Created by:
    %   David J. Schodt (Lidke Lab 2020)
    
    
    properties
        SMF % see SingleMoleculeFitting class
        SMD % see SingleMoleculeData class
        Verbose = 1; % (Default = 1) Verbosity level of main workflow
    end
    
    properties (Hidden)
        InternalParams = struct([]); % set of parameters used in LAP-FC
    end
    
    properties (SetAccess = 'protected')
        SMDCombined % frame-connected SMD structure
    end
    
    methods
        function [obj, SMDCombined, SMD] = FrameConnection(SMD, SMF, ...
                Verbose, AutoRun)
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
            
            % Set the verbosity level if provided.
            if (exist('Verbose', 'var') && ~isempty(Verbose))
                obj.Verbose = Verbose;
            end
            if ((obj.Verbose>2) && (nargin()>0))
                fprintf(['\tFrameConnection constructor: ', ...
                    'Setting class properties based on constructor ', ...
                    'inputs...\n'])
            end
            
            % Set the input SMD structure as a class property (if provided)
            if (exist('SMD', 'var') && ~isempty(SMD))
                if (obj.Verbose > 2)
                    fprintf(['\tFrameConnection constructor: ', ...
                        'Storing input SMD as a class property...\n'])
                end
                obj.SMD = SMD;
                AllFieldsSet = 1;
            else
                AllFieldsSet = 0;
            end
            
            % Store the SMF as a class property (if it was provided).
            if (exist('SMF', 'var') && ~isempty(SMF))
                % Pad the input SMF structure to ensure it contains all
                % fields defined in SingleMoleculeFitting.createSMF().
                obj.SMF = smi_core.SingleMoleculeFitting.padSMF(SMF);
            else
                obj.SMF = smi_core.SingleMoleculeFitting;
            end
            
            % If all required fields are set, run the frame-connection
            % process.
            if (AutoRun && AllFieldsSet)
                if (obj.Verbose > 2)
                    fprintf(['\tFrameConnection constructor: ', ...
                        'Auto-running obj.performFrameConnection()...\n'])
                end
                [SMDCombined, SMD] = obj.performFrameConnection();
            else
                if ((obj.Verbose>0) && (nargout()>1))
                    warning(['Constructor outputs SMD and ', ...
                        'SMDCombined were requested but ', ...
                        'frame-connection was not performed.'])
                end
                SMD = obj.SMD;
                SMDCombined = SMD;
            end
        end
        
        [SMDCombined, SMD] = performFrameConnection(obj)
        
    end
    
    methods (Static)
        [SMDIndex] = findConnected(SMR, SMD, ID);
        [SMDCombined, SMD] = hypothesisTestFC(SMD, SMF, Verbose);
        [SMD, InternalParams] = lapFC(SMD, SMF, Verbose, InternalParams);
        [SMD] = classicalFC(SMD, SMF, Verbose);
        [SMD] = revisedClassicalFC(SMD, SMF, Verbose);
        [ClusterData] = organizeClusterData(SMD);
        [KOn, KOff, KBleach, PMiss, NEmitters] = ...
            estimateRateParameters(ClusterData, Verbose);
        [InitialDensity] = estimateLocalDensity(ClusterData, ...
            NNearestClusters, KOn, KOff, KBleach, PMiss);
        [CostMatrix] = createCostMatrix(ClusterData, ...
            KOn, KOff, KBleach, PMiss, InitialDensity, MaxFrameGap, ...
            EndFrame, NonLinkMarker)
        [ConnectID, MaxConnectID] = ...
            linkClusters(ConnectID, MaxConnectID, UpdateIndices, Link12);
        [SMDCombined] = combineLocalizations(SMD, SMF);
        [Success] = unitTest();
    end
    
    
end