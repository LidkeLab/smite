function [SMDCombined, SMD, OutputMessage] = performFrameConnection(obj)
%performFrameConnection is the "run" method of the FrameConnection class.
%
% This method is intended to be used as the main "run" method of the
% FrameConnection class, meaning that most users will set the properties of
% the FrameConnection class, call this method, gather desired properties,
% and then move on to other analyses.
%
% NOTE: This method will add an additional field to obj.SMD called
%       "ConnectID".  obj.SMD.ConnectID is an integer array indicating
%       which localizations were connected during the frame connection
%       process.  For example, if
%       (obj.SMD.ConnectID(nn) == obj.SMD.ConnectID(mm)), the localizations
%       in SMD identified by the indices nn and mm were connected during
%       frame connection.  The exact value of the field "ConnectID" is
%       itself arbitrary and carries no meaning further than associating
%       localizations. This field is directly related to
%       obj.SMDCombined.ConnectID as follows:
%           For a given ConnectID, say nn, the indices in arrays of SMD
%           that were combined to generate a field in SMDCombined can be
%           found as IndicesSMD = find(SMD.ConnectID == nn) (alternatively,
%           IndicesSMD = smi_core.FrameConnection.findConnected(...
%               SMDCombined, SMD, nn) )
%
% INPUTS:
%   obj: An instance of the class smi_core/FrameConnection with all fields
%        populated with meaningful entries.
%
% OUTPUTS:
%   SMDCombined: SMDCombined contains the "frame-connected" localizations,
%                i.e., the result of performing frame-connection on obj.SMD
%   SMD: obj.SMD but with the field 'ConnectID' populated (see
%        smi_core.FrameConnection.findConnected() for a careful description
%        of 'ConnectID'.
%   OutputMessage: A summary message that describes how many localizations
%                  in SMD were collapsed into localizations in SMDCombined.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)
%       based on a code written Hanieh Mazloom-Farsibaf


% Make a local copy of obj.SMD (accessing obj.SMD.(...) inside a loop can
% sometimes be slower than accessing SMD.(...)).
SMD = obj.SMD; % temporary local variable

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
% NOTE: Since numel(unique(obj.SMD.DatasetNum)) is typically small (e.g.,
%       < 100) I don't anticipate any benefit to making local copies of any
%       arrays from SMD (as opposed to accessing them from SMD
%       inside the loop, as I'm doing here).
InputExtras = [];
InputExtrasSE = [];
SMD.ConnectID = zeros(numel(SMD.X), 1, 'uint32');
MaxConnectID = max(SMD.ConnectID);
for nn = unique(SMD.DatasetNum)
    % Isolate all valid localizations in the nn-th dataset and typecast
    % certain arrays for c_FrameConnect.mex* .
    CurrentBool = ((SMD.DatasetNum==nn) & (SMD.ThreshFlag==0));
    if ~any(CurrentBool)
        % All localizations in this dataset were thresholded out.
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
    switch obj.FitType
        case 'XYNBS'
            InputExtras = single(SMD.PSFSigma(CurrentBool));
            InputExtrasSE = single(SMD.PSFSigma_SE(CurrentBool).^2);
        case 'XYNBSXSY'
            InputExtras = single([SMD.PSFSigmaX(CurrentBool), ...
                SMD.PSFSigmaY(CurrentBool)]);
            InputExtrasSE = single([SMD.PSFSigmaX_SE(CurrentBool), ...
                SMD.PSFSigmaY_SE(CurrentBool)].^2);
        case 'XYZNB'
            InputCoords = [InputCoords, ...
                single(SMD.Z(CurrentBool))];
            InputCoordsSE = [InputCoordsSE, ...
                single(SMD.Z_SE(CurrentBool))];
    end
    
    % Perform the frame-connection using c_FrameConnect.mex*.
    [OutputCoords, OutputCoordsSE, NConnected, OutputFrames, ...
        OutputExtras, OutputExtrasSE, OutputPhotonsBgLogL, ...
        OutputConnectID, OutputConnectIDCombined] = ...
        c_FrameConnect(obj.LoS, ...
        InputCoords, InputCoordsSE, InputFrameNum, ...
        InputExtras, InputExtrasSE, InputPhotonsBgLogL, ...
        obj.MaxSeparation, obj.MaxFrameGap, MaxConnectID);
    
    % Update 'SMD' to contain the current connection information.
    SMD.ConnectID(CurrentBool) = OutputConnectID;
    MaxConnectID = max(OutputConnectID);
    
    % Store the outputs from c_FrameConnect.mex* in temporary arrays (these
    % arrays will later be copied into the 'SMDCombined' output).
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
    switch obj.FitType
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
SMDCombined = smi_core.SingleMoleculeData.createSMD();
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

% Add a new field (IndSMD) to SMDCombined that specifies which indices of 
% SMD were used to generate each entry.
%NLocTotal = numel(SMDCombined.ConnectID);
%IndSMD = cell(NLocTotal, 1);
%for ii = 1:NLocTotal
%    IndSMD{ii} = uint32(obj.findConnected(SMDCombined, SMD, ...
%        SMDCombined.ConnectID(ii)));
%end
%SMDCombined.IndSMD = IndSMD;

% Store the updated SMD and SMDCombined in obj.
obj.SMDCombined = SMDCombined;
obj.SMD = SMD;

% A helpful message enumerating how many localizations were collapsed here.
OutputMessage = sprintf('Frame connecting: %d -> %d localizations\n', ...
        numel(SMD.X), numel(SMDCombined.X));
    

end
