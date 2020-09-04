function performFrameConnection(obj)
%performFrameConnection is the "run" method of the FrameConnection class
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
%       localizations.
%
% INPUTS:
%   obj: An instance of the class smi_core/FrameConnection with all fields
%        populated with meaningful entries.

% Created by:
%   David J. Schodt (Lidke Lab 2020)


% Initialize the field "ConnectID", which will be used to string together
% trajectories.
obj.SMD.ConnectID = zeros(numel(obj.SMD.X), 1, 'uint32');
MaxConnectID = 0;

% Loop through each dataset and perform the frame connection process.
% NOTE: Since numel(unique(obj.SMD.DatasetNum)) is typically small (e.g.,
%       < 100) I don't anticipate any benefit to making local copies of any
%       arrays from obj.SMD (as opposed to accessing them from obj.SMD
%       inside the loop, as I'm doing here).
InputExtras = []; % extra inputs sent to c_FrameConnect.mex*
InputExtrasSE = [];
obj.SMDCombined = smi_core.SingleMoleculeData.createSMD();
obj.SMDCombined.NCombined = [];
obj.SMDCombined.CombinedID = [];
for nn = unique(obj.SMD.DatasetNum)
    % Isolate all valid localizations in the nn-th dataset and typecast
    % certain arrays for c_FrameConnect.mex* .
    CurrentBool = ((obj.SMD.DatasetNum==nn) & (obj.SMD.ThreshFlag==0));
    if ~any(CurrentBool)
        % All localizations in this dataset were thresholded out.
        continue
    end
    InputFrameNum = uint32(obj.SMD.FrameNum(CurrentBool));
    InputCoords = single([obj.SMD.X(CurrentBool), ...
        obj.SMD.Y(CurrentBool)]);
    InputCoordsSE = single([obj.SMD.X_SE(CurrentBool), ...
        obj.SMD.Y_SE(CurrentBool)]);
    InputPhotonsBgLogL = single([obj.SMD.Photons(CurrentBool), ...
        obj.SMD.Bg(CurrentBool), ...
        obj.SMD.LogL(CurrentBool)]);
    
    % Determine if we need to append additional arrays (e.g., for 3D
    % fitting we'll need to append 'Z' associated arrays).
    switch obj.FitType
        case 'XYNBS'
            InputExtras = single(obj.SMD.PSFSigma(CurrentBool));
            InputExtrasSE = single(obj.SMD.PSFSigma_SE(CurrentBool).^2);
        case 'XYNBSXSY'
            InputExtras = single([obj.SMD.PSFSigmaX(CurrentBool), ...
                obj.SMD.PSFSigmaY(CurrentBool)]);
            InputExtrasSE = single([obj.SMD.PSFSigmaX_SE(CurrentBool), ...
                obj.SMD.PSFSigmaY_SE(CurrentBool)].^2);
        case 'XYZNB'
            InputCoords = [InputCoords, ...
                single(obj.SMD.Z(CurrentBool))];
            InputCoordsSE = [InputCoordsSE, ...
                single(obj.SMD.Z_SE(CurrentBool))];
    end
    
    % Perform the frame-connection using c_FrameConnect.mex*.
    [OutputCoords, OutputCoordsSE, NConnected, OutputFrames, ...
        OutputExtras, OutputExtrasSE, OutputPhotonsBgLogL, ...
        OutputConnectID, OutputCombinedID] = c_FrameConnect(obj.LoS, ...
        InputCoords, InputCoordsSE, InputFrameNum, ...
        InputExtras, InputExtrasSE, InputPhotonsBgLogL, ...
        obj.MaxSeparation, obj.MaxFrameGap, MaxConnectID);
    
    % Store the outputs from c_FrameConnect.mex* in their appropriate place
    % in the class objects.
    obj.SMDCombined.X = [obj.SMDCombined.X; OutputCoords(:, 1)];
    obj.SMDCombined.Y = [obj.SMDCombined.Y; OutputCoords(:, 2)];
    obj.SMDCombined.X_SE = [obj.SMDCombined.X_SE; OutputCoordsSE(:, 1)];
    obj.SMDCombined.Y_SE = [obj.SMDCombined.Y_SE; OutputCoordsSE(:, 2)];
    obj.SMDCombined.NCombined = [obj.SMDCombined.NCombined; NConnected];
    obj.SMDCombined.FrameNum = [obj.SMDCombined.FrameNum; OutputFrames];
    obj.SMDCombined.Photons = [obj.SMDCombined.Photons; ...
        OutputPhotonsBgLogL(:, 1)];
    obj.SMDCombined.Bg = [obj.SMDCombined.Bg; ...
        OutputPhotonsBgLogL(:, 2)];
    obj.SMDCombined.LogL = [obj.SMDCombined.LogL; ...
        OutputPhotonsBgLogL(:, 3)];
    obj.SMDCombined.CombinedID = [obj.SMDCombined.CombinedID; ...
        OutputCombinedID];
    obj.SMDCombined.DatasetNum = [obj.SMDCombined.DatasetNum; ...
        nn*ones(numel(OutputFrames), 1, 'uint32')];
    obj.SMD.ConnectID(CurrentBool) = OutputConnectID;
    MaxConnectID = max(OutputConnectID);
    
    % Store some FitType dependent outputs if needed.
    if strcmpi(obj.FitType, 'XYNBS')
        SMDCombined.PSFSigma = [obj.SMDCombined.PSFSigma; ...
            OutputExtras];
        SMDCombined.PSFSigma_SE = [obj.SMDCombined.PSFSigma_SE; ...
            OutputExtrasSE];
    elseif strcmpi(obj.FitType, 'XYNBSXSY')
        SMDCombined.PSFSigmaX = [obj.SMDCombined.PSFSigmaX; ...
            OutputExtras(:, 1)];
        SMDCombined.PSFSigmaY = [obj.SMDCombined.PSFSigmaY; ...
            OutputExtras(:, 2)];
        SMDCombined.PSFSigmaX_SE = [obj.SMDCombined.PSFSigmaX_SE; ...
            OutputExtrasSE(:, 1)];
        SMDCombined.PSFSigmaY_SE = [obj.SMDCombined.PSFSigmaY_SE; ...
            OutputExtrasSE(:, 2)];
    elseif strcmpi(obj.FitType, 'XYZNB')
        SMDCombined.Z = [SMDCombined.Z; OutputCoords(:, 3)];
        SMDCombined.Z_SE = [SMDCombined.Z_SE; OutputCoordsSE(:, 3)];
    end
end


end