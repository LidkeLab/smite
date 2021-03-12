% This script demonstrates some of the functionality in the
% ChannelRegistration class.

% NOTE: This class also has a basic GUI which can do (most) of what is done
%       below. You can start using the GUI with 
%           ChannelReg = smi_core.ChannelRegistration;
%           ChannelReg.gui();
%       The GUI cannot be used to transform coordinates/images, so you'll
%       have to do that in script form (as demonstrated below) after
%       computing the transform with the GUI.

%% Two separate fiducial files, one image per file (2 channels total).
% Define the path to the fiducial files.
SMITEPath = fileparts(which('setupSMITE.m'));
FiducialFileDir = fullfile(SMITEPath, ...
    'examples', 'example_data', 'channel_registration');
FiducialFileNames = {'Fiducial1.mat', 'Fiducial2.mat'};

% Prepare an SMF structure for fitting the fiducials.
% NOTE: If using ChannelReg.AutoscaleFiducials = true, you can often get
%       away with only changing SMF.Fitting.PSFSigma.
SMF = smi_core.SingleMoleculeFitting;
SMF.Fitting.PSFSigma = 1.3;

% Prepare the channel registration class.
ChannelReg = smi_core.ChannelRegistration(...
    FiducialFileDir, FiducialFileNames, SMF);

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

% Apply the transform to an SMD structure.
FiducialSize = diff(ChannelReg.FiducialROI([1, 2; 3, 4])) + 1;
SMD = smi_core.SingleMoleculeData.createSMD();
SMD.X = linspace(1, FiducialSize(2), 100).';
SMD.Y = linspace(1, FiducialSize(1), 100).';
SMDTransformed = ChannelReg.transformSMD(...
    ChannelReg.RegistrationTransform{2}, ...
    SMD);

% Apply the transform to an image (not meant to be viewed, just showing how
% it can be done!).
NImages = 100;
TestImages = randi(123, [FiducialSize, NImages]);
TransformedImages = ChannelReg.transformImages(...
    ChannelReg.RegistrationTransform{2}, ...
    TestImages);

% Visualize the performance of the channel registration.
PlotFigure = figure();
ChannelReg.visualizeRegistrationResults(PlotFigure, ...
    ChannelReg.RegistrationTransform{2}, ...
    ChannelReg.Coordinates{2}(:, :, 2), ...
    ChannelReg.Coordinates{2}(:, :, 1), ...
    ChannelReg.FiducialImages(:, :, 2), ...
    ChannelReg.FiducialImages(:, :, 1));

% Visualize the registration error.
PlotFigure = figure();
PlotAxes = axes(PlotFigure);
ChannelReg.visualizeRegistrationError(PlotAxes, ...
    ChannelReg.RegistrationTransform{2}, ...
    ChannelReg.Coordinates{2}(:, :, 2), ...
    ChannelReg.Coordinates{2}(:, :, 1));

% Visualize the transform magnitude and gradient (this isn't usually useful
% unless something went very wrong, in which case it might be obvious in
% this plot!).
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

% % Save the transform.
% ChannelReg.exportTransform([], pwd())

%% Single fiducial image with two channels side by side.
% Define the path to the fiducial files.
SMITEPath = fileparts(which('setupSMITE.m'));
FiducialFileDir = fullfile(SMITEPath, ...
    'examples', 'example_data', 'channel_registration');
FiducialFileNames = 'Fiducials2Channel.mat';

% Prepare an SMF structure for fitting the fiducials.
% NOTE: If using ChannelReg.AutoscaleFiducials = true, you can often get
%       away with only changing SMF.Fitting.PSFSigma.
SMF = smi_core.SingleMoleculeFitting;
SMF.Fitting.PSFSigma = 1.3;

% Prepare the channel registration class.
ChannelReg = smi_core.ChannelRegistration(...
    FiducialFileDir, FiducialFileNames, SMF);

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
