FiducialFileDir = 'C:\Users\David\Documents\MATLAB\channel_reg';
FiducialFileNames{1} = 'channel1_fiducial.mat';
FiducialFileNames{2} = 'channel2_fiducial.mat';
FiducialFileNames{3} = 'channel3_fiducial.mat';
SMF = smi_core.SingleMoleculeFitting;
SMF.Fitting.FitType = 'XYNB';
SMF.Fitting.PSFSigma = 1.3;
ChannelReg = smi_core.ChannelRegistration(...
    FiducialFileDir, FiducialFileNames, SMF);
ChannelReg.NNeighborPoints = 12;
ChannelReg.TransformationType = 'lwm';
ChannelReg.findTransform();
figure()
ChannelReg.visualizeCoordTransform([], ...
    ChannelReg.RegistrationTransform{2}, [256, 256]);
figure()
ChannelReg.visualizeImageTransform([], ...
    ChannelReg.RegistrationTransform{2}, [256, 256])

Coords1 = ChannelReg.Coordinates{2}(:, :, 1);
Coords2 = ChannelReg.Coordinates{2}(:, :, 2);
RegistrationTransform = ChannelReg.RegistrationTransform{2};
MSE = ChannelReg.estimateRegistrationError(...
    RegistrationTransform, Coords2, Coords1);
figure()
ChannelReg.visualizeRegistrationError([], ...
    RegistrationTransform, Coords2, Coords1);
figure()
ChannelReg.visualizeRegistrationResults([], ...
    ChannelReg.RegistrationTransform{2}, ...
    ChannelReg.Coordinates{2}(:, :, 2), ...
    ChannelReg.Coordinates{2}(:, :, 1), ...
    ChannelReg.FiducialImages(:, :, 2), ...
    ChannelReg.FiducialImages(:, :, 1));
% ChannelReg.TransformationBasis = 'images';
% ChannelReg.TransformationType = 'affine';
% ChannelReg.findTransform();

%%
ChannelReg = smi_core.ChannelRegistration();
% ChannelReg.TransformationBasis = 'images';
ChannelReg.gui()

%%
FiducialFileDir = 'C:\Users\David\Documents\MATLAB\channel_reg';
FiducialFileNames{1} = 'channel1_fiducial.mat';
SMF = smi_core.SingleMoleculeFitting;
SMF.Fitting.FitType = 'XYNB';
SMF.Fitting.PSFSigma = 1.3;
ChannelReg = smi_core.ChannelRegistration(...
    FiducialFileDir, FiducialFileNames, SMF);
ChannelReg.NNeighborPoints = 12;
ChannelReg.TransformationType = 'lwm';
ChannelReg.findTransform();

