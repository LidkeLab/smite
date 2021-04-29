function [TR, SMD] = performFullAnalysis(obj)
%performFullAnalysis fits and tracks data pointed to by obj.SMF
% This method is the main run method for the smi.SPT class, meaning that it
% will load raw data, perform gain/offset correction, fit localizations to
% the data, create trajectories from the localizations, and then save the
% results.
%
% OUTPUTS:
%   TR: Tracking Results structure (see smi_core.TrackingResults)
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Prepare an SMLM class (we'll use this to load data, perform gain/offset
% correction, and fit the data).
if obj.TryLowPValueLocs
    obj.SMFCopy = smi_core.SingleMoleculeFitting.reloadSMF(obj.SMF);
    obj.SMF.Thresholding.MinPValue = 0;
end
SMLM = smi.SMLM(obj.SMF);
SMLM.Verbose = obj.Verbose;

% Load data, perform gain/offset correction, and fit the data.
SMLM.analyzeAll()
obj.SMF = SMLM.SMF;
obj.SMD = SMLM.SMD;
obj.SMDPreThresh = SMLM.SMDPreThresh;

% Generate trajectories from the localizations in obj.SMD.
% NOTE: On the first pass, the diffusion constant in obj.SMF.Tracking.D
%       will be used for all cost matrices.
obj.DiffusionConstant = [];
obj.generateTrajectories()

% If needed, estimate diffusion constants and recursively track with the
% new values.
if obj.UseTrackByTrackD
    OnesArray = ones(numel(obj.SMD.FrameNum), 1);
    for rr = 1:obj.NRecursions
        % Send an update to the command window.
        if (obj.Verbose > 1)
            fprintf(['\tsmi.spt.performFullAnalysis(): ', ...
                'tracking recursion iteration %i\n'], rr)
        end
        
        % Estimate the diffusion constants from the previous tracking 
        % results.
        obj.DiffusionEstimator.TR = obj.TR;
        DiffusionStruct = ...
            obj.DiffusionEstimator.estimateDiffusionConstant();
        DiffusionConstantCurrent = ...
            [DiffusionStruct(1).DiffusionConstant, ...
            DiffusionStruct(1).DiffusionConstantSE];
        obj.DiffusionConstant = OnesArray ...
            * [DiffusionStruct(2).DiffusionConstant, ...
            DiffusionStruct(2).DiffusionConstantSE];            
        for ii = 1:numel(obj.TR)
            obj.DiffusionConstant(obj.TR(ii).IndSMD, 1) = ...
                DiffusionConstantCurrent(ii, 1);
            obj.DiffusionConstant(obj.TR(ii).IndSMD, 2) = ...
                DiffusionConstantCurrent(ii, 2);
        end
        
        % Filter the diffusion constants, setting those which are invalid 
        % to the previously estimated ensemble value. For the corresponding
        % SEs, we'll set the bad values to an SE of inf (since we don't
        % actually know the value, and having a low SE might prevent those
        % localizations from being connected to others).
        BadValueBool = ((obj.DiffusionConstant(:, 1)<0) ...
            | isnan(obj.DiffusionConstant(:, 1)) ...
            | isinf(obj.DiffusionConstant(:, 1)));
        obj.DiffusionConstant(BadValueBool, 1) = ...
            DiffusionStruct(2).DiffusionConstant;
        obj.DiffusionConstant(BadValueBool, 2) = inf;
        
        % Re-track the data.
        obj.SMD.ConnectID = [];
        obj.generateTrajectories();
    end
end

% Remove short trajectories from the TR structure.
% NOTE: I'm leaving everything in SMD.  It might be nice to also threshold
%       short trajectories in SMD, but for now I'll leave it this way.
obj.TR = smi_core.TrackingResults.threshTrajLength(obj.TR, ...
    obj.SMF.Tracking.MinTrackLength);

% Make copies of TR and SMD for the outputs.
if nargout
    TR = obj.TR;
    SMD = obj.SMD;
end

% Save the results.
obj.saveResults()


end