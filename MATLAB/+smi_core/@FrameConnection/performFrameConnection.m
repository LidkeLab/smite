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
ConnectID = zeros(numel(SMD.X), 1, 'uint64');
MaxConnectID = 0;

% Initialize some arrays that we may populate inside the for loop below.
% Most of these are the frame-connected fields from SMD, e.g., ConnectedXY
% is an array containing the frame-connected SMD.X and SMD.Y.
ConnectedXY = [];
ConnectedSE = [];
ConnectedFrameNum = [];
ConnectedPSFSigma = [];
NCombined = [];
CombinedID = [];


end