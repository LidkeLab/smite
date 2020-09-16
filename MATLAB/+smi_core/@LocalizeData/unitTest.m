function [Success] = unitTest()
%unitTest tests vital functionality of the LocalizeData class.
% This method performs various tests to ensure that vital functionality of
% the smi_core.LocalizeData class is working as intended.
%
% NOTE: Failure of this unit test may in fact be caused by failure of some
%       other class methods OUTSIDE of the LocalizeData class.  A more
%       direct/better test of this class could not be conceptualized at the
%       time of writing.
%
% INPUTS:
%
% OUTPUTS:
%   Success: An array of boolean flags to indicate success of various tests
%            performed.
%            Success(1): genLocalizations - SMD generated successfully.
%            Success(2): genLocalizations - ScaledData generated
%                        successfully

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Initialize the Success output.
Success = zeros(2, 1, 'logical');

% Seed the random number generator so that we always get the same results.
rng(1234)

% Generate some simulated raw data (I grabbed some numbers from a random
% sCMOS calibration file and tried to roughly match those numbers for this
% simulation).
% NOTE: Adding structure to the raw data probably doesn't add anything to
%       the unitTest(), but it seems worth it in the case that there is a
%       failure in the unitTest(): it might be easier to track down if
%       there is something meaningful to look at.
FrameSizeFull = 256; % don't change this! other numbers assume = 256
NFrames = 10;
Background = 10; % e- (treated as photons)
CameraGain = 2 + 0.2*randn(FrameSizeFull, 'single'); % ADU / e-
CameraOffset = 100 + 0.5*randn(FrameSizeFull, 'single'); % ADU
CameraReadNoise = 0*(3 + randn(FrameSizeFull, 'single')).^2; % ADU^2
SMD = smi_core.SingleMoleculeData.createSMD();
SMD.X = repmat(128 + 64*[0; 1; 1; -1; -1], [NFrames, 1]);
SMD.Y = repmat(128 + 64*[0; 1; -1; 1; -1], [NFrames, 1]);
SMD.Photons = 1e3 * ones(5*NFrames, 1);
SMD.PSFSigma = 1.3;
SMD.FrameNum = repelem((1:NFrames).', 5);
SMD.Bg = zeros(5*NFrames, 1);
[~, Data] = smi_sim.GaussBlobs.gaussBlobImage(FrameSizeFull, NFrames, ...
    SMD, Background);

% Add read noise to the simulated data and then convert to ADU.
% NOTE: gaussBlobImage() can also add read noise (in photons) but I wanted 
%       to keep ReadNoiseVariance in units of ADU^2 for unit consistency.
ReadNoiseVariancePhotons = CameraReadNoise ./ (CameraGain.^2);
DataWithReadNoise = Data ...
    + sqrt(ReadNoiseVariancePhotons).*randn(FrameSizeFull);
DataWithReadNoise(DataWithReadNoise < 0) = 0;
RawData = CameraGain.*DataWithReadNoise + CameraOffset;

% Generate an SMF structure.
SMF = smi_core.SingleMoleculeFitting.createSMF();
SMF.Data.CameraType = 'sCMOS';
SMF.Data.CameraGain = CameraGain;
SMF.Data.CameraOffset = CameraOffset;
SMF.Data.CameraReadNoise = CameraReadNoise;
SMF.BoxFinding.BoxSize = 10;
SMF.Fitting.PSFSigma = PSFSigma;

% Attempt to generate localizations from the simulated data.
LD = smi_core.LocalizeData(RawData, SMF);
[SMD, ScaledData] = LD.genLocalizations();

% Check that the SMD makes sense (an exact check won't be done though).
NEmitters = 5;
Success(1) = (all(SMD.FrameNum==repelem((1:NFrames).', NEmitters)) ...
    && (numel(SMD.X)==NEmitters*NFrames) ...
    && ~any(SMD.ThreshFlag));

% Check that ScaledData makes sense.
Success(2) = all(abs(ScaledData(:)-DataWithReadNoise(:)) < 0.1);


end