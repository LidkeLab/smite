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
%            Success(2): LocalizeData - constructor working as intended

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Initialize the Success output.
Success = zeros(2, 1, 'logical');

% Seed the random number generator so that we always get the same results.
rng(1234)

% Generate some simulated data in units of photons.
FrameSizeFull = 256; % don't change this! other numbers assume = 256
NFrames = 10;
Background = 10;
SMD = smi_core.SingleMoleculeData.createSMD();
SMD.X = repmat(128 + 64*[0; 1; 1; -1; -1], [NFrames, 1]);
SMD.Y = repmat(128 + 64*[0; 1; -1; 1; -1], [NFrames, 1]);
SMD.Photons = 1e3 * ones(5*NFrames, 1);
SMD.PSFSigma = 1.3;
SMD.FrameNum = repelem((1:NFrames).', 5);
SMD.Bg = zeros(5*NFrames, 1);
[~, ScaledData] = smi_sim.GaussBlobs.gaussBlobImage(...
    FrameSizeFull, NFrames, SMD, Background);

% Generate an SMF structure.
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.CameraType = 'SCMOS';
SMF.BoxFinding.BoxSize = 10;
SMF.Fitting.PSFSigma = SMD.PSFSigma;

% Attempt to generate localizations from the simulated data.
LD = smi_core.LocalizeData(ScaledData, SMF);
[SMD, SMDPreThresh] = LD.genLocalizations();

% Check that SMD and SMDPreThresh make sense and are consistent with each
% other (this isn't meant to be an exact check of individual fields).
NEmitters = 5;
Success(1) = ...
    (all(SMDPreThresh.FrameNum==repelem((1:NFrames).', NEmitters)) ...
    && (numel(SMDPreThresh.X)==NEmitters*NFrames) ...
    && (numel(SMD.X)==sum(SMDPreThresh.ThreshFlag==0)));

% Check that the constructor is setting class properties as intended.
% NOTE: I'm just choosing numbers that'll stick out from the defaults.
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.CameraType = 'SCMOS';
SMF.BoxFinding.BoxSize = 8;
SMF.BoxFinding.BoxOverlap = 3;
SMF.BoxFinding.MinPhotons = 123;
SMF.Fitting.FitType = 'XYNBS';
SMF.Fitting.PSFSigma = 1.23;
SMF.Fitting.Iterations = 31;
SMF.Fitting.ZFitStruct.Ax = 1.2;
SMF.Fitting.ZFitStruct.Ay = 2.1;
SMF.Fitting.ZFitStruct.Bx = 1.2;
SMF.Fitting.ZFitStruct.By = 2.1;
SMF.Fitting.ZFitStruct.Gamma = 1.2;
SMF.Fitting.ZFitStruct.D = 2.1;
SMF.Thresholding.MaxXY_SE = 0.123;
SMF.Thresholding.MaxZ_SE = 0.0321;
SMF.Thresholding.MinPValue = 0.0123;
SMF.Thresholding.MinPSFSigma = 1.23;
SMF.Thresholding.MaxPSFSigma = 2.13;
SMF.Thresholding.MinPhotons = 213;
SMF.Thresholding.MaxBg = 12.3;
MinMaxTrue.X_SE = [0, SMF.Thresholding.MaxXY_SE];
MinMaxTrue.Y_SE = [0, SMF.Thresholding.MaxXY_SE];
MinMaxTrue.Z_SE = [0, SMF.Thresholding.MaxZ_SE];
MinMaxTrue.PValue = [SMF.Thresholding.MinPValue, 1];
MinMaxTrue.PSFSigma = [SMF.Thresholding.MinPSFSigma, ...
    SMF.Thresholding.MaxPSFSigma];
MinMaxTrue.Photons = [SMF.Thresholding.MinPhotons, inf];
MinMaxTrue.Bg = [0, SMF.Thresholding.MaxBg];
LD = smi_core.LocalizeData(ScaledData, SMF);
Success(2) = (strcmp(LD.CameraType, SMF.Data.CameraType) ...
    && (LD.BoxSize==SMF.BoxFinding.BoxSize) ...
    && (LD.BoxOverlap==SMF.BoxFinding.BoxOverlap) ...
    && strcmp(LD.FitType, SMF.Fitting.FitType) ...
    && (LD.PSFSigma==SMF.Fitting.PSFSigma) ...
    && (LD.MinPhotons==SMF.BoxFinding.MinPhotons) ...
    && (LD.Iterations==SMF.Fitting.Iterations) ...
    && (LD.ZFitStruct.Ax==SMF.Fitting.ZFitStruct.Ax) ...
    && (LD.ZFitStruct.Ay==SMF.Fitting.ZFitStruct.Ay) ...
    && (LD.ZFitStruct.Bx==SMF.Fitting.ZFitStruct.Bx) ...
    && (LD.ZFitStruct.Gamma==SMF.Fitting.ZFitStruct.Gamma) ...
    && (LD.ZFitStruct.D==SMF.Fitting.ZFitStruct.D) ...
    && all(LD.MinMax.X_SE==MinMaxTrue.X_SE) ...
    && all(LD.MinMax.Y_SE==MinMaxTrue.Y_SE) ...
    && all(LD.MinMax.Z_SE==MinMaxTrue.Z_SE) ...
    && all(LD.MinMax.PValue==MinMaxTrue.PValue) ...
    && all(LD.MinMax.PSFSigma==MinMaxTrue.PSFSigma) ...
    && all(LD.MinMax.Photons==MinMaxTrue.Photons) ...
    && all(LD.MinMax.Bg==MinMaxTrue.Bg));


end