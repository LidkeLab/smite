% This script demonstrates some of the functionality in the
% ChannelRegistration class.

% NOTE: This class also has a basic GUI which can do (most) of what is done
%       below. You can start using the GUI with 
%           ChannelReg = smi_core.ChannelRegistration;
%           ChannelReg.gui();
%       The GUI cannot be used to transform coordinates/images, so you'll
%       have to do that in script form (as demonstrated below) after
%       computing the transform with the GUI.

%% Simulate some fiducials to mimic format of real fiducial data.
% NOTE: This will save the simulated fiducials in 'FiducialFileDir' as .mat
%       files.
SMITEPath = fileparts(which('setupSMITE.m'));
FiducialFileDir = fullfile(SMITEPath, ...
    'examples', 'example_data', 'channel_registration');
[~, FiducialFiles] = smi_core.ChannelRegistration.simFiducials(128, 8, ...
    FiducialFileDir);
[~, FileNames, Extension] = fileparts(FiducialFiles);

%% Two separate fiducial files, one image per file (2 channels total).
% Define the path to the fiducial files.  The first entry is for the
% reference fiducial and the second is for the moving fiducial.
FiducialFileNames = {sprintf('%s.mat', FileNames{1}); 
    sprintf('%s.mat', FileNames{2})};

% Prepare an SMF structure for fitting the fiducials.
% NOTE: If using ChannelReg.AutoscaleFiducials = true, you can often get
%       away with only changing SMF.Fitting.PSFSigma.
SMF = smi_core.SingleMoleculeFitting;
SMF.Fitting.PSFSigma = 1.3;

% Prepare the channel registration class.
Verbose = 1;
ChannelReg = smi_core.ChannelRegistration(...
    FiducialFileDir, FiducialFileNames, SMF, Verbose);
ChannelReg.AutoscaleFiducials = true;
ChannelReg.ManualCull = true;

% Compute a locally weighted mean transform from the fiducials.
% NOTE: The transform from fiducial 2 to fiducial 1 will be stored in 
%       ChannelReg.RegistrationTransform{2}.
% NOTE: This will produce an interactive plot showing the fiducials and
%       overlain localizations, in which you can cull points that you don't 
%       want to use. If you don't want this to happen, you can set 
%       ChannelReg.ManualCull = false (not recommended).
ChannelReg.TransformationType = 'lwm';
ChannelReg.TransformationBasis = 'coordinates';
ChannelReg.NNeighborPoints = 12;
ChannelReg.findTransform();

% Estimate the resulting registration error in two ways: (1) RMSE between
% fiducial coordinates before and after the transform, and (2) RMSE between
% fiducial coordinates before and after transform using a leave-one-out
% (LOO) analysis.
% NOTE: For the LOO analysis, we have to provide some additional inputs to
%       estimate the RMSE.  These inputs will change depending on the value
%       of ChannelReg.TransformationType, e.g., for 'lwm', we have to input
%       NNeighborPoints, but for 'polynomial', we have to input the
%       PolynomialDegree.
FixedCoordinates = ChannelReg.Coordinates{2}(:, :, 1);
MovingCoordinates = ChannelReg.Coordinates{2}(:, :, 2);
RMSE = sqrt(mean(ChannelReg.estimateRegistrationError(...
    ChannelReg.RegistrationTransform{2}, ...
    MovingCoordinates, FixedCoordinates)));
RMSELOO = sqrt(mean(ChannelReg.estimateRegErrorLOO(...
    ChannelReg.TransformationType, {ChannelReg.NNeighborPoints}, ...
    MovingCoordinates, FixedCoordinates)));

% Visualize the performance of the channel registration.
FixedImages = ChannelReg.FiducialImages(:, :, 1);
MovingImages = ChannelReg.FiducialImages(:, :, 2);
PlotFigure = figure();
ChannelReg.visualizeRegistrationResults(PlotFigure, ...
    ChannelReg.RegistrationTransform{2}, ...
    MovingCoordinates, FixedCoordinates, ...
    MovingImages, FixedImages);

% Visualize the registration error.
PlotFigure = figure();
PlotAxes = axes(PlotFigure);
ChannelReg.visualizeRegistrationError(PlotAxes, ...
    ChannelReg.RegistrationTransform{2}, ...
    MovingCoordinates, FixedCoordinates);

% Visualize the transform magnitude and gradient (this isn't usually useful
% unless something went very wrong, in which case it might be obvious in
% this plot!).
FiducialSize = diff(ChannelReg.FiducialROI([1, 2; 3, 4])) + 1;
PlotFigure = figure();
ChannelReg.visualizeCoordTransform(PlotFigure, ...
    ChannelReg.RegistrationTransform{2}, FiducialSize);

% Visualize the transforms effect on images.
% WARNING: This one hurts my eyes a bit... Also, it's not too useful unless
%          the transform is very dramatic.
PlotFigure = figure();
PlotAxes = axes(PlotFigure);
ChannelReg.visualizeImageTransform(PlotAxes, ...
    ChannelReg.RegistrationTransform{2}, FiducialSize);

% Apply the transform to an SMD structure.
% NOTE: The important code is the call to ChannelReg.transformSMD().  The
%       rest is just me preparing an arbitrary SMD structure.
SMD = smi_core.SingleMoleculeData.createSMD();
SMD.X = linspace(1, FiducialSize(2), 100).';
SMD.Y = linspace(1, FiducialSize(1), 100).';
SMDTransformed = ChannelReg.transformSMD(...
    ChannelReg.RegistrationTransform{2}, ...
    SMD);

% Apply the transform to an image (not meant to be viewed, just showing how
% it can be done!).
% NOTE: You should never transform your raw images before analyzing!  This
%       should only be used for qualitative purposes, e.g., to make a movie
%       with SPT data.
NImages = 100;
TestImages = smi_sim.GaussBlobs.genRandomBlobImage(FiducialSize(1), NImages);
TransformedImages = ChannelReg.transformImages(...
    ChannelReg.RegistrationTransform{2}, ...
    TestImages);

% Save the transform.
ChannelReg.exportTransform(FiducialFileDir)

% Open the GUI, which can do all of the above actions (except transform an
% SMD/transform images).
ChannelReg.gui()

%% Single fiducial image with two channels side by side.
% Define the path to the fiducial files.
FiducialFileNames = {sprintf('%s.mat', FileNames{3})};

% Prepare an SMF structure for fitting the fiducials.
% NOTE: If using ChannelReg.AutoscaleFiducials = true, you can often get
%       away with only changing SMF.Fitting.PSFSigma.
SMF = smi_core.SingleMoleculeFitting;
SMF.Fitting.PSFSigma = 1.3;

% Prepare the channel registration class.
Verbose = 1;
ChannelReg = smi_core.ChannelRegistration(...
    FiducialFileDir, FiducialFileNames, SMF, Verbose);
ChannelReg.AutoscaleFiducials = true;
ChannelReg.ManualCull = true;

% Specify the fiducial file formatting. In this case, I'm specifying 
% [1, 2], which means that the fiducial image will be evenly split along
% its columns, with the right half considered the "moving" channel. See
% ChannelRegistration.convertSplitFormatToROIs() for more options.
ChannelReg.SplitFormat = [1, 2];

% Compute a locally weighted mean transform from the fiducials.
% NOTE: The transform from fiducial 2 to fiducial 1 will be stored in 
%       ChannelReg.RegistrationTransform{2}.
% NOTE: This will produce an interactive plot showing the fiducials and
%       overlain localizations, in which you can cull points that you don't 
%       want to use. If you don't want this to happen, you can set 
%       ChannelReg.ManualCull = false (not recommended).
ChannelReg.TransformationType = 'lwm';
ChannelReg.TransformationBasis = 'coordinates';
ChannelReg.NNeighborPoints = 12;
ChannelReg.findTransform();
