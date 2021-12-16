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
%   David J. Schodt (Lidke Lab, 2021)


% Make sure obj.SMD contains temporal data.  If it doesn't, exit here
% (since there's nothing to frame-connect).
if (isempty(obj.SMD.FrameNum) || (numel(obj.SMD.FrameNum)==1))
    SMD = obj.SMD;
    SMDCombined = SMD;
    obj.SMDCombined = obj.SMDCombined;
    return
end

% Perform frame-connection using the requested method.
if (obj.Verbose > 0)
    fprintf(['FrameConnection.performFrameConnection(): ', ...
        'Performing frame connection...\n'])
end
switch lower(obj.SMF.FrameConnection.Method)
    case 'hypothesis test'
        % As it's written now,  obj.hypothesisTestFC() is fully
        % self-contained, so we should return after running it.
        [SMDCombined, SMD] = obj.hypothesisTestFC(obj.SMD, obj.SMF, ...
            obj.Verbose);
        obj.SMDCombined = SMDCombined;
        obj.SMD = SMD;
        return
    case 'lap-fc'
        [obj.SMD, obj.InternalParams] = ...
            obj.lapFC(obj.SMD, obj.SMF, obj.Verbose, obj.InternalParams);
    case 'classical'
        obj.SMD = obj.classicalFC(obj.SMD, obj.SMF, obj.Verbose);
    case 'revised classical'
        obj.SMD = obj.revisedClassicalFC(obj.SMD, obj.SMF, obj.Verbose);
end
SMDCombined = obj.combineLocalizations(obj.SMD, obj.SMF);
obj.SMDCombined = SMDCombined;
if (obj.Verbose > 0)
    fprintf(['FrameConnection.performFrameConnection(): ', ...
        'Frame connection complete.\n'])
    fprintf('FrameConnection: %d -> %d localizations\n', ...
        numel(SMD.X), numel(SMDCombined.X));
end


end
