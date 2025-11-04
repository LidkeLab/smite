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
%            Success(1): Test of the class constructor.
%            Success(2): method convertToPhotons(), RawData is gain and
%                        offset corrected succesfully for non-scalar gain
%                        and offset.
%            Success(3): method convertToPhotons(), the input read-noise
%                        (given as a variance) is succesfully converted to
%                        units of photons^2 for non-scalar gain and offset.
%            Success(4): Same as Success(1) for scalar gain.
%            Success(5): Same as Success(2) for scalar offset.
%            Success(6): Test sCMOS camera with scalar calibration (Orca
%                        Fusion scenario) - scalars auto-expand to 2D arrays.

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
CameraGain = 2 + 0.2*randn(FrameSizeFull, 'single'); % ADU / e-
CameraOffset = 100 + 0.5*randn(FrameSizeFull, 'single'); % ADU
CameraReadNoise = (3 + randn(FrameSizeFull, 'single')).^2; % ADU^2
SMD = smi_core.SingleMoleculeData.createSMD();
SMD.X = repmat(128 + 64*[0; 1; 1; -1; -1], [NFrames, 1]);
SMD.Y = repmat(128 + 64*[0; 1; -1; 1; -1], [NFrames, 1]);
SMD.Photons = 1e3 * ones(5*NFrames, 1);
SMD.PSFSigma = 1.3;
SMD.FrameNum = repelem((1:NFrames).', 5);
SMD.Bg = zeros(5*NFrames, 1);
SMD.YSize = FrameSizeFull;
SMD.XSize = FrameSizeFull;
SMD.NFrames = NFrames;
[~, Data] = smi_sim.GaussBlobs.gaussBlobImage(SMD);

% Add read noise to the simulated data and then convert to ADU.
% NOTE: gaussBlobImage() can also add read noise (in photons) but I wanted
%       to keep ReadNoiseVariance in units of ADU^2 for unit consistency.
ReadNoiseVariancePhotons = CameraReadNoise ./ (CameraGain.^2);
DataWithReadNoise = Data ...
    + sqrt(ReadNoiseVariancePhotons).*randn(FrameSizeFull);
DataWithReadNoise(DataWithReadNoise < 0) = 0;
RawDataFull = CameraGain.*DataWithReadNoise + CameraOffset;

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
CalibrationROI = [1, 1, FrameSizeFull, FrameSizeFull];

% Create an instance of the DataToPhotons class and test the constructor.
% If this runs without throwing an exception, I'll just assume everything 
% worked.
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.CameraGain = CameraGain;
SMF.Data.CameraOffset = CameraOffset;
SMF.Data.CameraReadNoise = CameraReadNoise;
Success = zeros(1, 6, 'logical');
try
    DTP = smi_core.DataToPhotons(SMF, RawDataBottomRight, ...
        ROIBottomRight, CalibrationROI, true);
    DTP = smi_core.DataToPhotons(SMF, RawDataBottomRight, ...
        ROIBottomRight, CalibrationROI, false);
    Success(1) = true;
catch MException
    warning('error message: %s', ...
        MException.identifier, MException.message)
end

% Test the gain/offset corrections in convertToPhotons() in the corner test
% quadrants.
CorrectedData = zeros(FrameSizeFull, FrameSizeFull, NFrames, 'single');
CorrectedNoise = zeros(FrameSizeFull, 'single');
DTP.RawData = RawDataTopLeft;
DTP.RawDataROI = ROITopLeft;
[CorrectedData(IndicesTopLeft(:, 1), IndicesTopLeft(:, 2), :), ...
    CorrectedNoise(IndicesTopLeft(:, 1), IndicesTopLeft(:, 2))] = ...
    DTP.convertData();
DTP.RawData = RawDataTopRight;
DTP.RawDataROI = ROITopRight;
[CorrectedData(IndicesTopRight(:, 1), IndicesTopRight(:, 2), :), ...
    CorrectedNoise(IndicesTopRight(:, 1), IndicesTopRight(:, 2))] = ...
    DTP.convertData();
DTP.RawData = RawDataBottomRight;
DTP.RawDataROI = ROIBottomRight;
[CorrectedData(IndicesBottomRight(:, 1), IndicesBottomRight(:, 2), :), ...
    CorrectedNoise(IndicesBottomRight(:, 1), IndicesBottomRight(:, 2))] = ...
    DTP.convertData();
DTP.RawData = RawDataBottomLeft;
DTP.RawDataROI = ROIBottomLeft;
[CorrectedData(IndicesBottomLeft(:, 1), IndicesBottomLeft(:, 2), :), ...
    CorrectedNoise(IndicesBottomLeft(:, 1), IndicesBottomLeft(:, 2))] = ...
    DTP.convertData();
Success(2) = all(abs(CorrectedData(:)-DataWithReadNoise(:)) < 0.1);
Success(3) = all(abs(CorrectedNoise(:)-ReadNoiseVariancePhotons(:)) < 0.1);

% Test the gain/offset corrections in the center test quadrant, avoiding
% setting the input ROI to test the default indexing.
DTP.RawData = RawDataCenter;
DTP.RawDataROI = [];
[CorrectedDataCenter, CorrectedNoiseCenter] = DTP.convertData();
TrueDataCenter = DataWithReadNoise(IndicesCenter(:, 1), ...
    IndicesCenter(:, 2), :);
TrueNoiseCenter = ReadNoiseVariancePhotons(IndicesCenter(:, 1), ...
    IndicesCenter(:, 2));
Success(2) = (Success(2) ...
    & all(abs(CorrectedDataCenter(:)-TrueDataCenter(:)) < 0.1));
Success(3) = (Success(3) ...
    & all(abs(CorrectedNoiseCenter(:)-TrueNoiseCenter(:)) < 0.1));

% Repeat the above tests for scalar gain and offset (the quadrants no
% longer matter in this case so we can just test the entire image).
% Add read noise to the simulated data and then convert to ADU.
% NOTE: gaussBlobImage() can also add read noise (in photons) but I wanted
%       to keep ReadNoiseVariance in units of ADU^2 for unit consistency.
CameraGain = mean(CameraGain(:));
CameraOffset = mean(CameraOffset(:));
CameraReadNoise = mean(CameraReadNoise(:));
ScalarReadNoisePhotons = CameraReadNoise / CameraGain.^2;
DataWithScalarReadNoise = Data ...
    + sqrt(ScalarReadNoisePhotons)*randn(FrameSizeFull);
DataWithScalarReadNoise(DataWithScalarReadNoise < 0) = 0;
RawData = CameraGain.*DataWithScalarReadNoise + CameraOffset;
DTP.RawData = RawData;
DTP.CameraGain = CameraGain;
DTP.CameraOffset = CameraOffset;
DTP.CameraReadNoise = CameraReadNoise;
[CorrectedData, CorrectedNoise] = DTP.convertData();
Success(4) = all(abs(CorrectedData(:)-DataWithScalarReadNoise(:)) < 0.1);
% CorrectedNoise may be expanded to 2D array by scalar expansion logic
if isscalar(CorrectedNoise)
    Success(5) = (abs(CorrectedNoise-ScalarReadNoisePhotons) < 0.1);
else
    % All elements should equal the scalar value
    Success(5) = all(abs(CorrectedNoise(:)-ScalarReadNoisePhotons) < 0.1);
end

% Test sCMOS camera with scalar calibration (modern uniform sensors like
% Orca Fusion). This tests the auto-expansion of scalar values to 2D arrays.
SMF_sCMOS = smi_core.SingleMoleculeFitting;
SMF_sCMOS.Data.CameraType = 'SCMOS';
SMF_sCMOS.Data.CameraGain = CameraGain;  % Use the scalar value already defined
SMF_sCMOS.Data.CameraOffset = CameraOffset;
SMF_sCMOS.Data.CameraReadNoise = CameraReadNoise;
SMF_sCMOS.Data.CalibrationFilePath = '';  % Empty to trigger scalar path
try
    DTP_sCMOS = smi_core.DataToPhotons(SMF_sCMOS, RawData, [], [], 0, false);
    [CorrectedData_sCMOS, CorrectedNoise_sCMOS] = DTP_sCMOS.convertData();
    % Verify the conversion worked correctly
    % Note: The scalar expansion happens internally in convertToPhotons,
    % so DTP_sCMOS.CameraGain remains scalar, but the output arrays are
    % properly sized and the conversion is correct.
    Success(6) = all(abs(CorrectedData_sCMOS(:)-DataWithScalarReadNoise(:)) < 0.1) ...
        && (size(CorrectedNoise_sCMOS, 1) == FrameSizeFull) ...
        && (size(CorrectedNoise_sCMOS, 2) == FrameSizeFull);
catch MException
    warning('sCMOS scalar test failed: %s - %s', ...
        MException.identifier, MException.message)
    Success(6) = false;
end


end
