function [Success] = unitTest()
%unitTest tests vital functionality of smi_core.ChannelRegistration
%
% This method tests various methods in the class
% smi_core.ChannelRegistration class to ensure that the vital functionality
% is working as intended.  This is done by simulating some data, applying
% an artificial "warping" of the field of view, and then trying to undo
% this effect.
%
% OUTPUTS:
%   Success: An array of flags specifying the success of various methods in
%            smi_core.ChannelRegistration, where 1 means success and 0
%            means failure (of a specific method/functionality). 
%            (Boolean array)
%               Success(1): ChannelRegistration() (constructor)
%               Success(2): findTransform()
%               Success(3): estimateRegistrationError()
%               Success(4): visualizeCoordTransform()
%               Success(5): visualizeImageTransform()
%               Success(6): visualizeRegistrationError()
%               Success(7): visualizeRegistrationResults()

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


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
FrameSize = 128;
NEmittersPerRow = 10;
NEmitters = NEmittersPerRow ^ 2;
PhotonSum = 1e4;
Background = 5;
PSFSigma = 1.3;
SigmaNoise = PSFSigma / sqrt(PhotonSum);
PSFSigmaNoise = 0.123;
LogLikelihood = -40;
GridSpacing = FrameSize / (NEmittersPerRow+1);
[XMesh, YMesh] = meshgrid(GridSpacing:GridSpacing:(FrameSize-GridSpacing));
Coordinates = [XMesh(:), YMesh(:)] + SigmaNoise*randn(NEmitters, 2);

% Construct a "fixed" (reference) SMD structure.
OnesArray =  ones(NEmitters, 1);
SMDFixed = smi_core.SingleMoleculeData.createSMD();
SMDFixed.X = reshape(Coordinates(:, 1), NEmitters, 1);
SMDFixed.Y = reshape(Coordinates(:, 2), NEmitters, 1);
SMDFixed.Z = [];
SMDFixed.X_SE = SigmaNoise * OnesArray;
SMDFixed.Y_SE = SMDFixed.X_SE;
SMDFixed.Z_SE = [];
SMDFixed.Photons = PhotonSum * OnesArray;
SMDFixed.Bg = Background * OnesArray;
SMDFixed.LogLikelihood = LogLikelihood * OnesArray;
SMDFixed.PSFSigma = PSFSigma*OnesArray ...
    + PSFSigmaNoise*randn(NEmitters, 1);
SMDFixed.PSFSigma_SE = PSFSigmaNoise * OnesArray;
SMDFixed.FrameNum = ones(NEmitters, 1);
SMDFixed.DatasetNum = OnesArray;
SMDFixed.ThreshFlag = 0 * OnesArray;

% Construct the artificially shifted "moving" SMD structure by applying an
% affine transform to it.
% NOTE: These numbers were just made up at the time of writing, and I
%       didn't really check through the math carefully!
RotationAngle = 0.0123 * (pi/180);
Shift = [1; 2];
Shear = [0.01; 0.005];
Scale = [1.005; 1.001];
RotationMatrix = [cos(RotationAngle), -sin(RotationAngle); ...
    sin(RotationAngle), cos(RotationAngle)];
ShearMatrix = [1, Shear(1); ...
    Shear(2), 1];
ScaleMatrix = [Scale(1), 0; ...
    0, Scale(2)];
WarpingMatrix = RotationMatrix * ShearMatrix * ScaleMatrix;
Coords = [SMDFixed.X, SMDFixed.Y].';
WarpedCoords = Coords.';
for ii = 1:size(Coords, 2)
    WarpedCoords3D = WarpingMatrix*Coords(:, ii) + Shift;
    WarpedCoords(ii, :) = WarpedCoords3D(1:2).';
end
SMDMoving = SMDFixed;
SMDMoving.X = WarpedCoords(:, 1);
SMDMoving.Y = WarpedCoords(:, 2);

% Simulate some raw data.
% NOTE: Warping the first image to make the moving image would be faster,
%       but I'm not doing that because this is slightly easier to code...
[~, RawDataFixed] = smi_sim.GaussBlobs.gaussBlobImage(...
    FrameSize, 1, SMDFixed, Background);
[~, RawDataMoving] = smi_sim.GaussBlobs.gaussBlobImage(...
    FrameSize, 1, SMDMoving, Background);

% Find a transform to align these channels (testing various components
% along the way).
SMF = smi_core.SingleMoleculeFitting;
SMF.BoxFinding.MinPhotons = 100;
SMF.Fitting.FitType = 'XYNBS';
SMF.Fitting.PSFSigma = PSFSigma;
try
    ChannelReg = smi_core.ChannelRegistration([], [], SMF);
    Success(1) = true;
catch MException
    warning('ChannelRegistration constructor failed with error: %s-%s',...
        MException.identifier, MException.message)
end
ChannelReg.NNeighborPoints = 12;
ChannelReg.TransformationType = 'lwm';
ChannelReg.ManualCull = false;
ChannelReg.SeparationThreshold = GridSpacing / 2;
ChannelReg.FiducialImages = cat(3, RawDataFixed, RawDataMoving);
try
    ChannelReg.findTransform();
    Success(2) = true;
catch MException
    warning(['ChannelRegistration.findTransform() ', ...
        'failed with error: %s-%s'], ...
        MException.identifier, MException.message)
end

% Test estimateRegistrationError() and make sure the error is "reasonable"
% (BASED ON THE HARD-CODED THRESHOLD BELOW).
ErrorThreshold = 0.1;
try
    SquaredDifference = ChannelReg.estimateRegistrationError(...
        ChannelReg.RegistrationTransform{2}, ...
        ChannelReg.Coordinates{2}(:, :, 2), ...
        ChannelReg.Coordinates{2}(:, :, 1));
    Success(3) = (sqrt(mean(SquaredDifference)) <= ErrorThreshold);
catch MException
    warning(['ChannelRegistration.findTransform() ', ...
        'failed with error: %s-%s'], ...
        MException.identifier, MException.message)
end

% Test several visualization tools.
PlotFigure = figure('Visible', 'off');
try
    ChannelReg.visualizeCoordTransform(PlotFigure, ...
        ChannelReg.RegistrationTransform{2}, FrameSize * [1, 1]);
    Success(4) = true;
catch MException
    warning(['ChannelReg.visualizeCoordTransform() ', ...
        'failed with error: %s-%s'], ...
        MException.identifier, MException.message)
end
clf(PlotFigure);
PlotAxes = axes(PlotFigure);
try
    ChannelReg.visualizeImageTransform(PlotAxes, ...
        ChannelReg.RegistrationTransform{2}, FrameSize * [1, 1]);
    Success(5) = true;
catch MException
    warning(['ChannelReg.visualizeImageTransform() ', ...
        'failed with error: %s-%s'], ...
        MException.identifier, MException.message)
end
cla(PlotAxes);
try
ChannelReg.visualizeRegistrationError(PlotAxes, ...
    ChannelReg.RegistrationTransform{2}, ...
    ChannelReg.Coordinates{2}(:, :, 2), ...
    ChannelReg.Coordinates{2}(:, :, 1));
    Success(6) = true;
catch MException
    warning(['ChannelReg.visualizeRegistrationError() ', ...
        'failed with error: %s-%s'], ...
        MException.identifier, MException.message)
end
cla(PlotAxes);
try
    ChannelReg.visualizeRegistrationResults(PlotAxes, ...
        ChannelReg.RegistrationTransform{2}, ...
        ChannelReg.Coordinates{2}(:, :, 2), ...
        ChannelReg.Coordinates{2}(:, :, 1), ...
        ChannelReg.FiducialImages(:, :, 2), ...
        ChannelReg.FiducialImages(:, :, 1));
    Success(7) = true;
catch MException
    warning(['ChannelReg.visualizeRegistrationResults() ', ...
        'failed with error: %s-%s'], ...
        MException.identifier, MException.message)
end
close(PlotFigure);


end