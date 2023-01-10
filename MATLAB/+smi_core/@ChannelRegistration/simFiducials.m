function [RawData, Files] = simFiducials(ROISize, NRows, SaveDir)
%simFiducials simulates two channel fiducials to test channel registration.
% This method simulates a grid of test points with an affine transform
% shifting points from the reference channel.
%
% INPUTS:
%   ROISize: Size of the simulated fiducials. (Scalar)(Default=128)
%   NRows: Number of rows and columns of fiducial points.
%          (Scalar)(Default=8)
%   SaveDir: Optional save directory.  If provided, fiducials will be saved
%            in two separate files in the directory as 
%            "fiducial1<time>.mat" and "fiducial2<time>.mat", where
%            "<time>" is a timestamp generated upon saving.
%
% OUTPUTS:
%   RawData: Simulated two-channel noise-free fiducial data. Channel 1 is
%            the reference channel and channel 2 is the moving channel.
%            (ROISize x ROISize x 2)
%   Files: If "SaveDir" was given, this will be a cell array of the
%          filenames where the fiducials were saved.  Files{1} points to
%          the reference, Files{2} to the moving channel, and Files{3} to a
%          file containing a single image with both channels stored
%          side-by-side.

% Created by:
%   David J. Schodt (Lidke Lab, 2023)


% Set defaults.
if (~exist('ROISize', 'var') || isempty(ROISize))
    ROISize = 128;
end
if (~exist('NRows', 'var') || isempty(NRows))
    NRows = 8;
end

% Simulate several emitter positions.
% NOTE: Most of these parameters were chosen arbitrarily.
rng(1234)
GridSpacing = ROISize / (NRows+1);
[XMesh, YMesh] = meshgrid(GridSpacing:GridSpacing:(ROISize-GridSpacing));
Coordinates = [XMesh(:), YMesh(:)];

% Construct a "fixed" (reference) SMD structure.
NEmitters = NRows ^ 2;
PhotonSum = 1e4;
Background = 5;
PSFSigma = 1.3;
OnesArray =  ones(NEmitters, 1);
SMDFixed = smi_core.SingleMoleculeData.createSMD();
SMDFixed.X = reshape(Coordinates(:, 1), NEmitters, 1);
SMDFixed.Y = reshape(Coordinates(:, 2), NEmitters, 1);
SMDFixed.Photons = PhotonSum * OnesArray;
SMDFixed.Bg = Background * OnesArray;
SMDFixed.PSFSigma = PSFSigma;
SMDFixed.FrameNum = ones(NEmitters, 1);
SMDFixed.NFrames = 1;
SMDFixed.DatasetNum = OnesArray;
SMDFixed.ThreshFlag = 0 * OnesArray;

% Construct the artificially shifted "moving" SMD structure by applying an
% affine transform to it.
% NOTE: These numbers were just made up at the time of writing to give
%       something that looked reasonable.
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
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.DataROI = [1, 1, ROISize, ROISize];
[~, RawData] = smi_sim.GaussBlobs.gaussBlobImage(...
    SMDFixed, SMF, Background);
[~, RawData(:, :, 2)] = smi_sim.GaussBlobs.gaussBlobImage(...
    SMDMoving, SMF, Background);

% If needed, save the fiducials as .mat files.
if exist('SaveDir', 'var')
    if ~isfolder(SaveDir)
        mkdir(SaveDir)
    end

    % Save fiducials in separate files.
    TimeString = smi_helpers.genTimeString();
    sequence = RawData(:, :, 1);
    Files{1, 1} = fullfile(SaveDir, ...
        sprintf('fiducial1_%s.mat', TimeString));
    save(Files{1}, 'sequence')
    sequence = RawData(:, :, 2);
    Files{2, 1} = fullfile(SaveDir, ...
        sprintf('fiducial2_%s.mat', TimeString));
    save(Files{2}, 'sequence')

    % Save fiducials in a single file, side-by-side.
    sequence =  [RawData(:, :, 1), RawData(:, :,  2)];
    Files{3, 1} = fullfile(SaveDir, ...
        sprintf('fiducial12_%s.mat', TimeString));
    save(Files{3}, 'sequence')
end


end
