function [Success] = unitTest()
%unitTest tests vital functionality of the class smi_core.FrameConnection
%
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
%               Success(1): FrameConnection() (constructor)
%               Success(2): performFrameConnection(), 'XYNB'
%               Success(3): performFrameConnection(), 'XYNBS'
%               Success(4): performFrameConnection(), 'XYNBSXSY'
%               Success(5): performFrameConnection(), 'XYZNB'
%               Success(6): findConnected()

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Initialize the Success output.
Success = zeros(6, 1, 'logical');

% Seed the random number generator so that the simulated SMD is predictable
% NOTE: If this is changed, there will almost certainly be entries of
%       Success that are 0.
rng(1234)

% Simulate several emitter positions.
% NOTE: Most of these parameters were chosen arbitrarily, HOWEVER, there
%       may be some checks below that might be affected by changing these
%       parameters (i.e., causing some elements of Success to be 0 even if
%       things worked correctly).
FrameSize = 64;
NEmitters = 7;
NFrames = 13;
PhotonSum = 1e3;
Background = 5;
SigmaNoise = 0.2;
PSFSigma = 1.3;
PSFSigmaXY = [1.2; 1.4];
PSFSigmaNoise = 0.123;
LogLikelihood = -40;
Coordinates = zeros(NEmitters, 2, NFrames);
Coordinates(:, :, 1) = FrameSize * rand(NEmitters, 2);
Coordinates(:, :, 2:NFrames) = Coordinates(:, :, 1) ...
    + SigmaNoise*randn(NEmitters, 2, NFrames-1);

% Construct the simulated SMD structure.
OnesArray =  ones(NEmitters*NFrames, 1);
SMDSim = smi_core.SingleMoleculeData.createSMD();
SMDSim.X = reshape(Coordinates(:, 1, :), NEmitters*NFrames, 1);
SMDSim.Y = reshape(Coordinates(:, 2, :), NEmitters*NFrames, 1);
SMDSim.Z = SigmaNoise * randn(NEmitters*NFrames, 1);
SMDSim.X_SE = SigmaNoise * OnesArray;
SMDSim.Y_SE = SMDSim.X_SE;
SMDSim.Z_SE = SMDSim.X_SE;
SMDSim.Photons = PhotonSum * OnesArray;
SMDSim.Bg = Background * OnesArray;
SMDSim.LogLikelihood = LogLikelihood * OnesArray;
SMDSim.PSFSigma = PSFSigma*OnesArray ...
    + PSFSigmaNoise*randn(NEmitters*NFrames, 1);
SMDSim.PSFSigma_SE = PSFSigmaNoise * OnesArray;
SMDSim.PSFSigmaX = PSFSigmaXY(1)*OnesArray ...
    + PSFSigmaNoise*randn(NEmitters*NFrames, 1);
SMDSim.PSFSigmaY = PSFSigmaXY(2)*OnesArray ...
    + PSFSigmaNoise*randn(NEmitters*NFrames, 1);
SMDSim.PSFSigmaX_SE = PSFSigmaNoise * OnesArray;
SMDSim.PSFSigmaY_SE = SMDSim.PSFSigmaX_SE;
SMDSim.FrameNum = repelem((1:NFrames).', NEmitters);
SMDSim.DatasetNum = OnesArray;
SMDSim.ThreshFlag = 0 * OnesArray;

% Prepare the frame-connection class and ensure the constructor is working
% as intended.
SMF = smi_core.SingleMoleculeFitting;
SMF.FrameConnection.FitType = 'XYNB';
SMF.FrameConnection.LoS = 0.01023;
SMF.FrameConnection.MaxSeparation = 1.023;
SMF.FrameConnection.MaxFrameGap = 4;
FC = smi_core.FrameConnection(SMDSim, SMF, 3);
FCParams = SMF.FrameConnection;
Success(1) = (~isempty(FC.SMD) ...
    && strcmp(FC.SMF.Fitting.FitType, FCParams.FitType) ...
    && (double(FC.SMF.FrameConnection.LoS)==double(FCParams.LoS)) ...
    && (double(FC.SMF.FrameConnection.MaxSeparation) ...
        ==double(FCParams.MaxSeparation)) ...
    && (uint32(FC.SMF.FrameConnection.MaxFrameGap) ...
        ==uint32(FCParams.MaxFrameGap)));

% Compute some expected arrays for the simulation.
% NOTE: We should also compute the expected values of PSFSigma,
%       PSFSigma_SE, ..., but I don't know how those work at the time of
%       writing this code!
ExpectedVariance = 1 / (NFrames/(SigmaNoise^2));
ExpectedXYZ_SE = sqrt(ExpectedVariance);
ExpectedX = sum(reshape(SMDSim.X, [], NFrames)/(SigmaNoise^2), 2) ...
    * ExpectedVariance;
ExpectedY = sum(reshape(SMDSim.Y, [], NFrames)/(SigmaNoise^2), 2) ...
    * ExpectedVariance;
ExpectedZ = sum(reshape(SMDSim.Z, [], NFrames)/(SigmaNoise^2), 2) ...
    * ExpectedVariance;
ExpectedPhotons = PhotonSum * NFrames;
ExpectedBg = Background * NFrames;
ExpectedLogLikelihood = LogLikelihood * NFrames;

% Perform the frame-connection procedure for the fit type 'XYNB' and
% perform various checks to see if it worked correctly.
[SMDCombined] = FC.performFrameConnection();
Success(2) = baseSuccess(SMDCombined, ...
    ExpectedX, ExpectedY, ExpectedXYZ_SE, ...
    ExpectedPhotons, ExpectedBg, ExpectedLogLikelihood, ...
    NFrames, NEmitters);

% Perform the frame-connection procedure for the fit type 'XYNBS' and
% perform various checks to see if it worked correctly.
% NOTE: The precision of my floating point comparisons below were chosen
%       arbitrarily.
FC.SMF.Fitting.FitType = 'XYNBS';
[SMDCombined] = FC.performFrameConnection();
BaseSuccess = baseSuccess(SMDCombined, ...
    ExpectedX, ExpectedY, ExpectedXYZ_SE, ...
    ExpectedPhotons, ExpectedBg, ExpectedLogLikelihood, ...
    NFrames, NEmitters);
Success(3) = (BaseSuccess ...
    && (abs(sum(SMDCombined.PSFSigma)-9.2274) <= 1e-4) ...
    && all(abs(SMDCombined.PSFSigma_SE-0.0042) <= 1e-4));

% Perform the frame-connection procedure for the fit type 'XYNBSXSY' and
% perform various checks to see if it worked correctly.
FC.SMF.Fitting.FitType = 'XYNBSXSY';
[SMDCombined] = FC.performFrameConnection();
BaseSuccess = baseSuccess(SMDCombined, ...
    ExpectedX, ExpectedY, ExpectedXYZ_SE, ...
    ExpectedPhotons, ExpectedBg, ExpectedLogLikelihood, ...
    NFrames, NEmitters);
Success(4) = (BaseSuccess ...
    && (abs(sum(SMDCombined.PSFSigmaX)-8.5437) <= 1e-4) ...
    && (abs(sum(SMDCombined.PSFSigmaY)-9.5851) <= 1e-4) ...
    && all(abs(SMDCombined.PSFSigmaX_SE-0.0042) <= 1e-4) ...
    && all(abs(SMDCombined.PSFSigmaY_SE-0.0042) <= 1e-4));

% Perform the frame-connection procedure for the fit type 'XYZNB' and
% perform various checks to see if it worked correctly.
FC.SMF.Fitting.FitType = 'XYZNB';
[SMDCombined] = FC.performFrameConnection();
BaseSuccess = baseSuccess(SMDCombined, ...
    ExpectedX, ExpectedY, ExpectedXYZ_SE, ...
    ExpectedPhotons, ExpectedBg, ExpectedLogLikelihood, ...
    NFrames, NEmitters);
Success(5) = (BaseSuccess ...
    && all(abs(SMDCombined.Z-ExpectedZ) < 1e-4) ...
    && all(abs(SMDCombined.Z_SE-ExpectedXYZ_SE) < 1e-4));

% Test FrameConnection.findConnected() on some new, incomplete 'SMD' type
% structures.  This is done by randomly permuting initial (pre frame
% connection) SMD.ConnectID indices, using findConnected() to determine the
% new permuted indices, and then summing over the permuted indices found
% for each element of SMDCombined.ConnectID to check the result.
% NOTE: In the computation of 'Success' below, I've used ismember instead
%       of == just in case the output dimensions of findConnected() are
%       changed someday.
SMDCombinedInc.ConnectID = (1:NEmitters).';
ConnectIDSMD = repelem((1:NEmitters).', NFrames, 1);
SMDInc.ConnectID = ConnectIDSMD(randperm(NEmitters * NFrames));
IndSMDSum = zeros(NEmitters, 1);
for ii = 1:NEmitters
    IndSMDSum(ii) = sum(FC.findConnected(SMDCombinedInc, SMDInc, ...
        SMDCombinedInc.ConnectID(ii)));
end
Success(6) = all(ismember(double(IndSMDSum), ...
    [573; 585; 693; 550; 605; 487; 693]));


    function [BaseSuccess] = baseSuccess(SMDCombined, ...
            ExpectedX, ExpectedY, ExpectedXYZ_SE, ...
            ExpectedPhotons, ExpectedBg, ExpectedLogLikelihood, ...
            NFrames, NEmitters)
        % This function computes the baseline success of the
        % frame-connection, meaning that it checks several fields in the
        % SMDCombined structure that will not change based on the fit-type.
        % NOTE: The precision of my floating point comparisons below were
        %       chosen arbitrarily.
        
        BaseSuccess = ...
            (all(abs(SMDCombined.X-ExpectedX) < 1e-4) ...
            && all(abs(SMDCombined.Y-ExpectedY) < 1e-4) ...
            && all(abs(SMDCombined.X_SE-ExpectedXYZ_SE) < 1e-4) ...
            && all(abs(SMDCombined.Y_SE-ExpectedXYZ_SE) < 1e-4) ...
            && all(abs(SMDCombined.Photons-ExpectedPhotons) < 1e-4) ...
            && all(abs(SMDCombined.LogLikelihood...
                -ExpectedLogLikelihood) < 1e-4) ...
            && all(abs(SMDCombined.Bg-ExpectedBg) < 1e-4) ...
            && all(uint32(SMDCombined.DatasetNum)==uint32(1)) ...
            && all(uint32(SMDCombined.FrameNum)==uint32(NFrames)) ...
            && all(uint32(SMDCombined.ConnectID)==uint32(1:NEmitters).')...
            && all(uint32(SMDCombined.NCombined)==uint32(NFrames)));
        
    end


end