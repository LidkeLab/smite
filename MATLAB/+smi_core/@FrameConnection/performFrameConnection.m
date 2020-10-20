function [SMDCombined, SMD] = performFrameConnection(obj)
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

% Created by:
%   David J. Schodt (Lidke Lab, 2020)
%       based on a code written Hanieh Mazloom-Farsibaf


% Initialize several new fields in SMD and SMDCombined.
SMD = obj.SMD; % temporary local variable
SMD.ConnectID = zeros(numel(SMD.X), 1, 'uint32');
SMDCombined = SMD;

% Delete several fields from the initialized SMDCombined to ensure they
% don't interfere with the analysis.
SMDCombined.NCombined = [];
SMDCombined.ConnectID = [];
SMDCombined.X = [];
SMDCombined.Y = [];
SMDCombined.X_SE = [];
SMDCombined.Y_SE = [];
SMDCombined.Z = [];
SMDCombined.Z_SE = [];
SMDCombined.FrameNum = [];
SMDCombined.Photons = [];
SMDCombined.Bg = [];
SMDCombined.LogLikelihood = [];
SMDCombined.DatasetNum = [];
SMDCombined.PSFSigma = [];
SMDCombined.PSFSigma_SE = [];
SMDCombined.PSFSigmaX = [];
SMDCombined.PSFSigmaY = [];
SMDCombined.PSFSigmaX_SE = [];
SMDCombined.PSFSigmaY_SE = [];

% Loop through each dataset and perform the frame connection process.
% NOTE: Since numel(unique(obj.SMD.DatasetNum)) is typically small (e.g.,
%       < 100) I don't anticipate any benefit to making local copies of any
%       arrays from SMD (as opposed to accessing them from SMD
%       inside the loop, as I'm doing here).
InputExtras = []; % extra inputs sent to c_FrameConnect.mex*
InputExtrasSE = [];
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
    
    % Store the outputs from c_FrameConnect.mex* in their appropriate place
    % in the class objects.
    SMDCombined.X = [SMDCombined.X; OutputCoords(:, 1)];
    SMDCombined.Y = [SMDCombined.Y; OutputCoords(:, 2)];
    SMDCombined.X_SE = [SMDCombined.X_SE; OutputCoordsSE(:, 1)];
    SMDCombined.Y_SE = [SMDCombined.Y_SE; OutputCoordsSE(:, 2)];
    SMDCombined.NCombined = [SMDCombined.NCombined; NConnected];
    SMDCombined.FrameNum = [SMDCombined.FrameNum; OutputFrames];
    SMDCombined.Photons = [SMDCombined.Photons; ...
        OutputPhotonsBgLogL(:, 1)];
    SMDCombined.Bg = [SMDCombined.Bg; ...
        OutputPhotonsBgLogL(:, 2)];
    SMDCombined.LogLikelihood = [SMDCombined.LogLikelihood; ...
        OutputPhotonsBgLogL(:, 3)];
    SMDCombined.ConnectID = [SMDCombined.ConnectID; ...
        OutputConnectIDCombined];
    SMDCombined.DatasetNum = [SMDCombined.DatasetNum; ...
        nn*ones(numel(OutputFrames), 1, 'uint32')];
    SMD.ConnectID(CurrentBool) = OutputConnectID;
    MaxConnectID = max(OutputConnectID);
    
    % Store some FitType dependent outputs if needed.
    if strcmpi(obj.FitType, 'XYNBS')
        SMDCombined.PSFSigma = [SMDCombined.PSFSigma; ...
            OutputExtras];
        SMDCombined.PSFSigma_SE = [SMDCombined.PSFSigma_SE; ...
            OutputExtrasSE];
    elseif strcmpi(obj.FitType, 'XYNBSXSY')
        SMDCombined.PSFSigmaX = [SMDCombined.PSFSigmaX; ...
            OutputExtras(:, 1)];
        SMDCombined.PSFSigmaY = [SMDCombined.PSFSigmaY; ...
            OutputExtras(:, 2)];
        SMDCombined.PSFSigmaX_SE = [SMDCombined.PSFSigmaX_SE; ...
            OutputExtrasSE(:, 1)];
        SMDCombined.PSFSigmaY_SE = [SMDCombined.PSFSigmaY_SE; ...
            OutputExtrasSE(:, 2)];
    elseif strcmpi(obj.FitType, 'XYZNB')
        SMDCombined.Z = [SMDCombined.Z; OutputCoords(:, 3)];
        SMDCombined.Z_SE = [SMDCombined.Z_SE; OutputCoordsSE(:, 3)];
    end
end

% Add zeros to the ThreshFlag of SMDCombined (we should never be keeping
% localizations in SMDCombined which have non-zero ThreshFlag).
SMDCombined.ThreshFlag = zeros(numel(SMDCombined.FrameNum), 1);

% Add a new field (IndSMD) to SMDCombined that specifies which indices of 
% SMD were used to generate each entry.
NLocTotal = numel(SMDCombined.ConnectID);
IndSMD = cell(NLocTotal, 1);
for ii = 1:NLocTotal
    IndSMD{ii} = uint32(obj.findConnected(SMDCombined, SMD, ...
        SMDCombined.ConnectID(ii)));
end
SMDCombined.IndSMD = IndSMD;

% Store the updated SMD and SMDCombined in obj.
obj.SMDCombined = SMDCombined;
obj.SMD = SMD;

fprintf('Frame connecting: %d -> %d localizations\n', ...
        numel(SMD), numel(SMDCombined.X));

end
