function [Success] = unitTest()
%unitTest tests vital functionality of the class smi_core.FrameConnection
% This method tests various methods in the class smi_core.FrameConnection
% class to ensure that the vital functionality is working as intended.
% This is done by simulating a populated SMD structure with known
% associations between entries, performing the frame-connection procedure
% using smi_core.FrameConnection methods, and then checking the result.
%
% INPUTS:
%
% OUTPUTS: 
%   Success: An array of flags specifying the success of various methods in
%            smi_core.FrameConnection, where 1 means success and 0 means
%            failure (of a specific method). (Boolean array)
%               Success(1): FrameConnection.FrameConnection
%               Success(2): FrameConnection.performFrameConnection()

% Created by: 
%   David J. Schodt (Lidke Lab, 2020)


% Initialize the Success output.
Success = zeros(2, 1, 'logical');

% Seed the random number generator so that simulated SMD is predictable.
% NOTE: If this is changed, there will almost certainly be entries of
%       Success that are 0.
rng(1234)

% Simulate several emitter positions.
% NOTE: Most of these parameters were chosen arbitrarily, HOWEVER, due to
%       the "checksum" style of checking I do below, these should not be
%       changed (at least not without changing the checks below).
FrameSize = 64;
NEmitters = 7;
NFrames = 13;
PhotonSum = 1e3;
Background = 5; 
SigmaNoise = 0.2;
LogL = -40; % this was the average in a (different) simulation I did
Coordinates = zeros(NEmitters, 2, NFrames);
Coordinates(:, :, 1) = FrameSize * rand(NEmitters, 2);
Coordinates(:, :, 2:NFrames) = Coordinates(:, :, 1) ...
    + SigmaNoise*randn(NEmitters, 2, NFrames-1);

% Construct the simulated SMD structure.
OnesArray =  ones(NEmitters*NFrames, 1);
SMDSim = smi_core.SingleMoleculeData.createSMD();
SMDSim.X = reshape(Coordinates(:, 1, :), NEmitters*NFrames, 1);
SMDSim.Y = reshape(Coordinates(:, 2, :), NEmitters*NFrames, 1);
SMDSim.X_SE = SigmaNoise * OnesArray;
SMDSim.Y_SE = SMDSim.X_SE;
SMDSim.Photons = PhotonSum * OnesArray;
SMDSim.Bg = Background * OnesArray;
SMDSim.LogL = LogL * OnesArray; % -40 is ~ what I got in another simulation
SMDSim.FrameNum = repelem((1:NFrames).', NEmitters);
SMDSim.DatasetNum = 1;
SMDSim.ThreshFlag = 0 * OnesArray;

% Prepare the frame-connection class and ensure the constructor is working
% as intended.
SMF = smi_core.SingleMoleculeFitting.createSMF();
SMF.FrameConnection.LoS = 0.01023;
SMF.FrameConnection.MaxSeparation = 1.023;
SMF.FrameConnection.MaxFrameGap = 4.0321;
FC = smi_core.FrameConnection(SMDSim, SMF);
FCParams = SMF.FrameConnection;
Success(1) = (~isempty(FC.SMD) ...
    && (double(FC.LoS)==double(FCParams.LoS)) ...
    && (double(FC.MaxSeparation)==double(FCParams.MaxSeparation)) ...
    && (uint32(FC.MaxFrameGap)==uint32(FCParams.MaxFrameGap)));

% Perform the frame-connection procedure and perform various checks to see
% if it worked correctly.
% NOTE: I determined the proper equalities of the sums by testing this
%       first...
FC.performFrameConnection();
ExpectedSE = round(sqrt(1 / (NFrames/(SigmaNoise^2))), 4);
Success(2) = ((round(sum(double(FC.SMDCombined.X)), 4)==215.4026) ...
    && (round(sum(double(FC.SMDCombined.Y)), 4)==313.1330) ...
    && all(round(double(FC.SMDCombined.X_SE), 4)==ExpectedSE) ... ...
    && all(round(double(FC.SMDCombined.Y_SE), 4)==ExpectedSE) ...
    && all(round(double(FC.SMDCombined.Photons))==PhotonSum*NFrames) ...
    && all(round(double(FC.SMDCombined.LogL))==LogL*NFrames) ...
    && all(round(double(FC.SMDCombined.Bg))==Background*NFrames) ...
    && all(uint32(FC.SMDCombined.DatasetNum)==uint32(1)) ...
    && all(uint32(FC.SMDCombined.FrameNum)==uint32(NFrames)) ...
    && all(uint32(FC.SMDCombined.ConnectID)==uint32(1:NEmitters).') ...
    && all(uint32(FC.SMDCombined.NCombined)==uint32(NFrames)));
    
   
end