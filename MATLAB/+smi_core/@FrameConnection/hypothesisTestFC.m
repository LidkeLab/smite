function [SMDCombined, SMD] = hypothesisTestFC(SMD, SMF, Verbose)
%hypothesisTestFC connects localizations in 'SMD' by the hyp. test method.
% This method solves the frame-connection problem using the hypothesis
% testing method for connecting localizations.  That is, a p-value is
% computed and compared to the level of significance to test the null
% hypothesis that the tested localizations arose from the same emitter.
%
% NOTE: This method will add an additional field to SMD called
%       "ConnectID".  SMD.ConnectID is an integer array indicating
%       which localizations were connected during the frame connection
%       process.  For example, if
%       (SMD.ConnectID(nn) == SMD.ConnectID(mm)), the localizations
%       in SMD identified by the indices nn and mm were connected during
%       frame connection.  The exact value of the field "ConnectID" is
%       itself arbitrary and carries no meaning further than associating
%       localizations. This field is directly related to
%       SMDCombined.ConnectID as follows:
%           For a given ConnectID, say nn, the indices in arrays of SMD
%           that were combined to generate a field in SMDCombined can be
%           found as IndicesSMD = find(SMD.ConnectID == nn) (alternatively,
%           IndicesSMD = smi_core.FrameConnection.findConnected(...
%               SMDCombined, SMD, nn) )
%
% INPUTS:
%   SMD: SingleMoleculeData structure with the localizations that we wish
%        to frame-connect.
%   SMF: SingleMoleculeFitting structure defining relevant parameters.
%   Verbose: Integer specifying the verbosity level. (Default = 1)
%
% OUTPUTS:
%   SMDCombined: SMDCombined contains the "frame-connected" localizations,
%                i.e., the result of performing frame-connection on SMD
%   SMD: SMD but with the field 'ConnectID' populated.
%
% CITATION:
%   Wester, M. J., Schodt, D. J., Mazloom-Farsibaf, H., Fazel, M., 
%   Pallikkuth, S., and Lidke, K. A. (2021). 
%   Robust, Fiducial-Free Drift Correction for Super-resolution Imaging.
%   bioRxiv , 2021.03.26.437196doi:10.1101/2021.03.26.437196

% Created by:
%   David J. Schodt (Lidke Lab, 2020)
%       based on a code written Hanieh Mazloom-Farsibaf


% Set defaults if needed.
if (~exist('SMF', 'var') || isempty(SMF))
    SMF = smi_core.SingleMoleculeFitting;
end
if (~exist('Verbose', 'var') || isempty(Verbose))
    Verbose = 1;
end    
    
% Initialize some temporary arrays (these are arrays like X, Y, ... which
% are concatenated inside of the main for loop below, and later stored in
% the output 'SMDCombined').
NCombined = [];
ConnectID = [];
X = [];
Y = [];
X_SE = [];
Y_SE = [];
Z = [];
Z_SE = [];
FrameNum = [];
Photons = [];
Bg = [];
LogLikelihood = [];
DatasetNum = [];
PSFSigma = [];
PSFSigma_SE = [];
PSFSigmaX = [];
PSFSigmaY = [];
PSFSigmaX_SE = [];
PSFSigmaY_SE = [];

% Loop through each dataset and perform the frame connection process.
% NOTE: Since numel(unique(SMD.DatasetNum)) is typically small (e.g.,
%       < 100) I don't anticipate any benefit to making local copies of any
%       arrays from SMD (as opposed to accessing them from SMD
%       inside the loop, as I'm doing here).
InputExtras = [];
InputExtrasSE = [];
SMD.ConnectID = zeros(numel(SMD.X), 1, 'uint32');
MaxConnectID = max(SMD.ConnectID);
UniqueDatasetNum = unique(SMD.DatasetNum);
if iscolumn(UniqueDatasetNum)
    UniqueDatasetNum = UniqueDatasetNum.';
end
for nn = UniqueDatasetNum
    % Provide a Command Window update if needed.
    if (Verbose > 2)
        fprintf(['\tFrameConnection.hypothesisTestFC(): ', ...
            'Performing frame connection for dataset %i...\n'], ...
            nn)
    end
    
    % Isolate all valid localizations in the nn-th dataset and typecast
    % certain arrays for smi_c_FrameConnection.mex* .
    CurrentBool = ((SMD.DatasetNum==nn) & (SMD.ThreshFlag==0));
    if ~any(CurrentBool)
        % All localizations in this dataset were thresholded out.
        if (Verbose > 2)
            fprintf(['\t\tAll dataset %i localizations ', ...
                'were thresholded.\n'], nn)
        end
        continue
    end
    InputFrameNum = uint32(SMD.FrameNum(CurrentBool));
    InputCoords = single([SMD.X(CurrentBool), ...
        SMD.Y(CurrentBool)]);
    InputCoordsSE = single([SMD.X_SE(CurrentBool), ...
        SMD.Y_SE(CurrentBool)]);
    InputPhotonsBgLogL = single([SMD.Photons(CurrentBool), ...
        SMD.Bg(CurrentBool), ...
        SMD.LogLikelihood(CurrentBool)]);
    
    % Determine if we need to append additional arrays (e.g., for 3D
    % fitting we'll need to append 'Z' associated arrays).
    switch SMF.Fitting.FitType
        case 'XYNBS'
            InputExtras = single(SMD.PSFSigma(CurrentBool));
            InputExtrasSE = single(SMD.PSFSigma_SE(CurrentBool));
        case 'XYNBSXSY'
            InputExtras = single([SMD.PSFSigmaX(CurrentBool), ...
                SMD.PSFSigmaY(CurrentBool)]);
            InputExtrasSE = single([SMD.PSFSigmaX_SE(CurrentBool), ...
                SMD.PSFSigmaY_SE(CurrentBool)]);
        case 'XYZNB'
            InputCoords = [InputCoords, ...
                single(SMD.Z(CurrentBool))];
            InputCoordsSE = [InputCoordsSE, ...
                single(SMD.Z_SE(CurrentBool))];
    end
    
    % Perform the frame-connection using smi_c_FrameConnection.mex*.
    if (Verbose > 2)
        fprintf('\t\tCalling smi_c_FrameConnection.mex*...\n')
    end
    [OutputCoords, OutputCoordsSE, NConnected, OutputFrames, ...
        OutputExtras, OutputExtrasSE, OutputPhotonsBgLogL, ...
        OutputConnectID, OutputConnectIDCombined] = ...
        smi_c_FrameConnection(SMF.FrameConnection.LoS, ...
        InputCoords, InputCoordsSE, InputFrameNum, ...
        InputExtras, InputExtrasSE, InputPhotonsBgLogL, ...
        SMF.FrameConnection.MaxSeparation, ...
        SMF.FrameConnection.MaxFrameGap, MaxConnectID);
    
    % Update 'SMD' to contain the current connection information.
    if (Verbose > 2)
        fprintf(['\t\t%i localizations combined to ', ...
            '%i localizations in dataset %i.\n'], ...
            numel(OutputConnectID), numel(OutputConnectIDCombined), nn)
    end
    SMD.ConnectID(CurrentBool) = OutputConnectID;
    MaxConnectID = max(OutputConnectID);
    
    % Store the outputs from smi_c_FrameConnection.mex* in temporary arrays
    % (these arrays will later be copied into the 'SMDCombined' output).
    X = [X; OutputCoords(:, 1)];
    Y = [Y; OutputCoords(:, 2)];
    X_SE = [X_SE; OutputCoordsSE(:, 1)];
    Y_SE = [Y_SE; OutputCoordsSE(:, 2)];
    NCombined = [NCombined; NConnected];
    FrameNum = [FrameNum; OutputFrames];
    Photons = [Photons; OutputPhotonsBgLogL(:, 1)];
    Bg = [Bg; OutputPhotonsBgLogL(:, 2)];
    LogLikelihood = [LogLikelihood; OutputPhotonsBgLogL(:, 3)];
    ConnectID = [ConnectID; OutputConnectIDCombined];
    DatasetNum = [DatasetNum; nn*ones(numel(OutputFrames), 1, 'uint32')];
    
    % Store some FitType dependent outputs if needed.
    switch SMF.Fitting.FitType
        case 'XYNBS'
            PSFSigma = [PSFSigma; OutputExtras];
            PSFSigma_SE = [PSFSigma_SE; OutputExtrasSE];
        case 'XYNBSXSY'
            PSFSigmaX = [PSFSigmaX; OutputExtras(:, 1)];
            PSFSigmaY = [PSFSigmaY; OutputExtras(:, 2)];
            PSFSigmaX_SE = [PSFSigmaX_SE; OutputExtrasSE(:, 1)];
            PSFSigmaY_SE = [PSFSigmaY_SE; OutputExtrasSE(:, 2)];
        case 'XYZNB'
            Z = [Z; OutputCoords(:, 3)];
            Z_SE = [Z_SE; OutputCoordsSE(:, 3)];
    end
end

% Store the temporary arrays from the main loop above into 'SMDCombined'
% (as well as some other things we'd like to carry along from SMD).
if (Verbose > 2)
    fprintf('\t\tPreparing output SMD structures...\n')
end
SMDCombined = smi_core.SingleMoleculeData.createSMD();
SMDCombined.NDims = SMD.NDims;
SMDCombined.NDatasets = SMD.NDatasets;
SMDCombined.NFrames = SMD.NFrames;
SMDCombined.XSize = SMD.XSize;
SMDCombined.YSize = SMD.YSize;
SMDCombined.NCombined = NCombined;
SMDCombined.ConnectID = ConnectID;
SMDCombined.X = X;
SMDCombined.Y = Y;
SMDCombined.X_SE = X_SE;
SMDCombined.Y_SE = Y_SE;
SMDCombined.Z = Z;
SMDCombined.Z_SE = Z_SE;
SMDCombined.FrameNum = FrameNum;
SMDCombined.Photons = Photons;
SMDCombined.Bg = Bg;
SMDCombined.LogLikelihood = LogLikelihood;
SMDCombined.PValue = smi_core.GaussMLE.pValue(...
    SMF.Fitting.NParams, SMF.BoxFinding.BoxSize, LogLikelihood);
SMDCombined.DatasetNum = DatasetNum;
SMDCombined.PSFSigma = PSFSigma;
SMDCombined.PSFSigma_SE = PSFSigma_SE;
SMDCombined.PSFSigmaX = PSFSigmaX;
SMDCombined.PSFSigmaY = PSFSigmaY;
SMDCombined.PSFSigmaX_SE = PSFSigmaX_SE;
SMDCombined.PSFSigmaY_SE = PSFSigmaY_SE;

% Add zeros to the ThreshFlag of SMDCombined (we should never be keeping
% localizations in SMDCombined which have non-zero ThreshFlag).
SMDCombined.ThreshFlag = zeros(numel(SMDCombined.FrameNum), 1);

% Provide a final message to summarize the results.
if (Verbose > 2)
    fprintf(['\tFrameConnection.hypothesisTestFC(): ', ...
        '%i localizations combined to %i localizations.\n'], ...
        numel(SMD.FrameNum), numel(SMDCombined.FrameNum))
end


end