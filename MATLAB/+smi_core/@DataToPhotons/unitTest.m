function [Success] = unitTest()
%unitTest checks vital functionality of the DataToPhotons class.
% This method will perform various tests to determine if the vital
% methods of the DataToPhotons class are working as intended.
%
% INPUTS:
%
% OUTPUTS: 
%   Success: A boolean array whose elements indicate success/failure of
%            various methods of the DataToPhotons class. (Boolean)
%            Success(1): method convertToPhotons(), RawData is gain and
%                        offset corrected succesfully for non-scalar gain
%                        and offset.
%            Success(2): method convertToPhotons(), the input read-noise 
%                        (given as a variance) is succesfully converted to 
%                        units of photons^2 for non-scalar gain and offset.
%            Success(3): Same as Success(1) for scalar gain.
%            Success(4): Same as Success(2) for scalar offset.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


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
SMF = smi_core.SingleMoleculeFitting.createSMF();
SMF.Data.CameraGain = 2 + 0.2*randn(FrameSizeFull, 'single'); % ADU / e-
SMF.Data.CameraOffset = 100 + 0.5*randn(FrameSizeFull, 'single'); % ADU
SMF.Data.CameraReadNoise = ...
    (3 +  1*randn(FrameSizeFull, 'single')).^2; % ADU^2
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
ReadNoiseVariancePhotons = SMF.Data.CameraReadNoise ...
    ./ (SMF.Data.CameraGain.^2);
DataWithReadNoise = Data ...
    + sqrt(ReadNoiseVariancePhotons).*randn(FrameSizeFull);
DataWithReadNoise(DataWithReadNoise < 0) = 0;
RawDataFull = SMF.Data.CameraGain.*DataWithReadNoise ...
    + SMF.Data.CameraOffset;

% Divide RawDataFull into 5 quadrants, one centered on each of the blobs
% in the simulation (this is just to test that I'm indexing into arrays
% like CameraGain correctly in convertToPhotons()).
IndicesTopLeft = [1:128; 1:128].'; % [rows, columns]
IndicesTopRight = [1:128; 129:256].';
IndicesBottomRight = [129:256; 129:256].';
IndicesBottomLeft = [129:256; 1:128].';
IndicesCenter = (128-64) + IndicesTopLeft;
RawDataTopLeft = RawDataFull(IndicesTopLeft(:, 1), ...
    IndicesTopLeft(:, 2), :);
RawDataTopRight = RawDataFull(IndicesTopRight(:, 1), ...
    IndicesTopRight(:, 2), :);
RawDataBottomRight = RawDataFull(IndicesBottomRight(:, 1), ...
    IndicesBottomRight(:, 2), :);
RawDataBottomLeft = RawDataFull(IndicesBottomLeft(:, 1), ...
    IndicesBottomLeft(:, 2), :);
RawDataCenter = RawDataFull(IndicesCenter(:, 1), ...
    IndicesCenter(:, 2), :);
ROITopLeft = [min(IndicesTopLeft), max(IndicesTopLeft)];
ROITopRight = [min(IndicesTopRight), max(IndicesTopRight)];
ROIBottomRight = [min(IndicesBottomRight), max(IndicesBottomRight)];
ROIBottomLeft = [min(IndicesBottomLeft), max(IndicesBottomLeft)];

% Test the gain/offset corrections in convertToPhotons() in the corner test
% quadrants.
Success = zeros(2, 1, 'logical');
CalibrationROI = [1, 1, FrameSizeFull, FrameSizeFull];
CorrectedData = zeros(FrameSizeFull, FrameSizeFull, NFrames, 'single');
CorrectedNoise = zeros(FrameSizeFull, 'single');
[CorrectedData(IndicesTopLeft(:, 1), IndicesTopLeft(:, 2), :), ...
    CorrectedNoise(IndicesTopLeft(:, 1), IndicesTopLeft(:, 2))] = ...
    smi_core.DataToPhotons.convertToPhotons(...
    RawDataTopLeft, SMF, ROITopLeft, CalibrationROI);
[CorrectedData(IndicesTopRight(:, 1), IndicesTopRight(:, 2), :), ...
    CorrectedNoise(IndicesTopRight(:, 1), IndicesTopRight(:, 2))] = ...
    smi_core.DataToPhotons.convertToPhotons(...
    RawDataTopRight, SMF, ROITopRight, CalibrationROI);
[CorrectedData(IndicesBottomRight(:, 1), IndicesBottomRight(:, 2), :), ...
    CorrectedNoise(IndicesBottomRight(:, 1), IndicesBottomRight(:, 2))] = ...
    smi_core.DataToPhotons.convertToPhotons(...
    RawDataBottomRight, SMF, ROIBottomRight, CalibrationROI);
[CorrectedData(IndicesBottomLeft(:, 1), IndicesBottomLeft(:, 2), :), ...
    CorrectedNoise(IndicesBottomLeft(:, 1), IndicesBottomLeft(:, 2))] = ...
    smi_core.DataToPhotons.convertToPhotons(...
    RawDataBottomLeft, SMF, ROIBottomLeft, CalibrationROI);
Success(1) = all(abs(CorrectedData(:)-DataWithReadNoise(:)) < 0.1);
Success(2) = all(abs(CorrectedNoise(:)-ReadNoiseVariancePhotons(:)) < 0.1);

% Test the gain/offset corrections in convertToPhotons() in the center test
% quadrant, avoiding setting the input ROI to test the default indexing.
[CorrectedDataCenter, CorrectedNoiseCenter] = ...
    smi_core.DataToPhotons.convertToPhotons(...
    RawDataCenter, SMF, [], CalibrationROI);
TrueDataCenter = DataWithReadNoise(IndicesCenter(:, 1), ...
    IndicesCenter(:, 2), :);
TrueNoiseCenter = ReadNoiseVariancePhotons(IndicesCenter(:, 1), ...
    IndicesCenter(:, 2));
Success(1) = (Success(1) ...
    & all(abs(CorrectedDataCenter(:)-TrueDataCenter(:)) < 0.1));
Success(2) = (Success(2) ...
    & all(abs(CorrectedNoiseCenter(:)-TrueNoiseCenter(:)) < 0.1));

% Repeat the above tests for scalar gain and offset (the quadrants no
% longer matter in this case so we can just test the entire image).
% Add read noise to the simulated data and then convert to ADU.
% NOTE: gaussBlobImage() can also add read noise (in photons) but I wanted 
%       to keep ReadNoiseVariance in units of ADU^2 for unit consistency.
SMF.Data.CameraGain = mean(SMF.Data.CameraGain(:));
SMF.Data.CameraOffset = mean(SMF.Data.CameraOffset(:));
SMF.Data.CameraReadNoise = mean(SMF.Data.CameraReadNoise(:));
ScalarReadNoisePhotons = SMF.Data.CameraReadNoise / SMF.Data.CameraGain.^2;
DataWithScalarReadNoise = Data ...
    + sqrt(ScalarReadNoisePhotons)*randn(FrameSizeFull);
DataWithScalarReadNoise(DataWithScalarReadNoise < 0) = 0;
RawData = SMF.Data.CameraGain.*DataWithScalarReadNoise ...
    + SMF.Data.CameraOffset;
[CorrectedData, CorrectedNoise] = ...
    smi_core.DataToPhotons.convertToPhotons(RawData, SMF);
Success(3) = all(abs(CorrectedData(:)-DataWithScalarReadNoise(:)) < 0.1);
Success(4) = (abs(CorrectedNoise-ScalarReadNoisePhotons) < 0.1);


end

