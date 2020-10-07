function Success = unitTest()
%unitTest Tests all functionality of smi.SMLM .
%
% OUTPUTS:
%   Success:    Short Description
%
% REQUIRES:
%   Statistics Toolbox
%   Parallel Procesing Toolbox
%   NVidia GPU

% start unitTest
Success = 0;
fprintf(['Running smi.SMLM.unitTest.\n', ...
         'Testing all smi.SMLM functionality ...\n']);

%% Test loading various datatypes.
fprintf('\nTesting loading of various datatypes.\n')

saveName = 'SMLM_testData';
file1 = [saveName, '1.mat'];
file2 = [saveName, '2.mat'];

% Create SMF structure.
SMF = smi_core.SingleMoleculeFitting.createSMF();
SMF.Data.FileDir      = tempdir;
%SMF.DataVariable      = 'Data';
SMF.Data.CameraType   = 'EMCCD';
SMF.Data.CameraGain   = 1;
SMF.Data.CameraOffset = 0;
%SMF.Fitting.FitType   = 'XYNB';
%SMF.ROISearchType     = 'Flat';

% Simulate two small datasets.
fprintf('Simulating data.\n')
[SimData1, ~] = smi_sim.GaussBlobs.gaussBlobImage(100, 10);
[SimData2, ~] = smi_sim.GaussBlobs.gaussBlobImage(100, 10);

% Save datasets as mat files.
fprintf('Saving data as mat files.\n')
Data = SimData1;
save(fullfile(tempdir, file1), 'Data');
Data = SimData2;
save(fullfile(tempdir, file2), 'Data');

% Try running smi.SMLM.  If it fails, delete files before returning error.
fprintf(['Loading and analyzing data saved as mat files.\n', ...
         '(Only doing box finding and fitting.)\n']);
% Update SMF object.
SMF.Data.FileName = {file1, file2};
% Create smi.SMLM object.
%SMLMobj = smi.SMLM('nogui');
SMLMobj = smi.SMLM(SMF);
try
    % Analyze all datasets.
    SMLMobj.analyzeAll();
    clear SMLMobj
catch ME
    delete(fullfile(tempdir, [saveName, '.*']));
    fprintf('Caught following error during smi.SMLM.unitTest:\n')
    disp(ME.identifier);
    disp(ME.message);
end
fprintf('Loading and analyzing data saved as mat file successful.\n');
delete(fullfile(tempdir, [saveName, '.*']));

% Save datasets as ics files.
fprintf('Saving data as ics files.\n')
file1 = [saveName '1.ics'];
file2 = [saveName '2.ics'];
writeim(SimData1,fullfile(tempdir,file1));
writeim(SimData2,fullfile(tempdir,file2));
% Try running smi.SMLM.  If it fails, delete files before returning error,
fprintf(['Loading and analyzing data saved as ics files.\n', ...
         '(Only doing box finding and fitting.)\n']);
% Update SMF object.
SMF.Data.FileName = {file1, file2};
% Create smi.SMLM object.
SMLMobj = smi.SMLM(SMF);
try
    %  Analyze all datasets.
    SMLMobj.analyzeAll();
    clear SMLMobj
catch ME
    delete(fullfile(tempdir,[saveName '.*']));
    fprintf('Caught following error during smi.SMLM.unitTest:\n')
    disp(ME.identifier)
    disp(ME.message);
end
fprintf('Loading and analyzing data saved as ics file successful.\n');
delete(fullfile(tempdir, [saveName, '.*']));

% Save datasets as h5 files.
fprintf('Saving data as h5 files.\n')
h5create(fullfile(tempdir,[saveName '.h5']),'/Data/Channel01/Data0001',size(SimData1));
h5write(fullfile(tempdir,[saveName '.h5']),'/Data/Channel01/Data0001',SimData1);
h5create(fullfile(tempdir,[saveName '.h5']),'/Data/Channel01/Data0002',size(SimData2));
h5write(fullfile(tempdir,[saveName '.h5']),'/Data/Channel01/Data0002',SimData2);
% Try running smi.SMLM.  If it fails, delete files before returning error,
fprintf(['Loading and analyzing data saved as h5 files.\n', ...
         '(Only doing box finding and fitting.)\n']);
% Update SMF object.
SMF.Data.FileName = {[saveName, '.h5']};
%SMF.RawImageSize = [size(SimData1,1),size(SimData1,2)];
% Create smi.SMLM object.
SMLMobj = smi.SMLM(SMF);
try
    % Analyze all datasets.
    SMLMobj.analyzeAll();
    clear SMLMobj
catch ME
    delete(fullfile(tempdir, [saveName '.*']));
    fprintf('Caught following error during smi.SMLM.unitTest:\n')
    disp(ME.identifier)
    disp(ME.message);
end
fprintf('Loading and analyzing data saved as h5 file successful.\n');
delete(fullfile(tempdir,[saveName '*.*']));

%% Simulate and save realistic SMLM data.
fprintf('\nSimulating realistic 2D SMLM data\n');
[SimData1, ~] = smi_sim.GaussBlobs.gaussBlobImage(256, 1000);
[SimData2, ~] = smi_sim.GaussBlobs.gaussBlobImage(256, 1000);
fprintf('\nSaving realistic SMSR data.\n');
saveName = 'SMLM_testData';
h5create(fullfile(tempdir,[saveName '.h5']),'/Data/Channel01/Data0001',size(SimData1));
h5write(fullfile(tempdir,[saveName '.h5']),'/Data/Channel01/Data0001',SimData1);
h5create(fullfile(tempdir,[saveName '.h5']),'/Data/Channel01/Data0002',size(SimData2));
h5write(fullfile(tempdir,[saveName '.h5']),'/Data/Channel01/Data0002',SimData2);

%% Test complete data flow
fprintf('\nTesting 2D analysis.\n');
% Create SMF structure.
SMF = smi_core.SingleMoleculeFitting.createSMF();
% Cannot save in folders inside tempdir, so everything must be saved in tempdir
% directly.
SMF.Data.FileDir          = tempdir;
SMF.Data.FileName         = {[saveName, '.h5']};
SMF.Data.ResultsDir       = tempdir;
SMF.Data.CameraType       = 'EMCCD';
SMF.Data.CameraGain       = 1;
SMF.Data.CameraOffset     = 0;
SMF.Thresholding.On       = true;
SMF.Thresholding.MaxSE_XY = 0.1;
SMF.FrameConnection.On    = true;
SMF.DriftCorrection.On    = true;
%SMF.RawImageSize = [size(SimData1,1),size(SimData1,2)];
% Create smi.SMLM object.
%SMLMobj = smi.SMLM('nogui');
SMLMobj = smi.SMLM(SMF);
% fullFit: fitting -> thresholding -> frame connection -> drift correction
SMLMobj.fullFit();
% generate output plots
fprintf('Generating output plots.\n');
SMLMobj.genPlots();
% generate color overlay
fprintf('Generating color overlay.\n');
%SMLMobj.SMD.FitBoxSize = SMLMobj.BoxSize;   % needed for genBlobOverlay
%SMLMobj.SMR.FitBoxSize = SMLMobj.BoxSize;   % needed for genBlobOverlay
DatasetNum = 1;
SMLMobj.genBlobOverlay(DatasetNum);
% save
SMLMobj.save();
SMLMobj.saveResults();
SMLMobj.exportFileType='Excel';
SMLMobj.exportResults();
SMLMobj.exportFileType='txt';
SMLMobj.exportResults();
% delete object and data
clear SMLMobj
delete(fullfile(tempdir,[saveName '*.*']));
fprintf('Done.\n');
