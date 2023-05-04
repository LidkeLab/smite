function Success = unitTest()
%unitTest Tests all functionality of smi.SMLM .
%
% OUTPUTS:
%   Success:    Flags indicating which tests passed (1 - yes, 0 - no)
%
% REQUIRES:
%   Statistics Toolbox
%   Parallel Procesing Toolbox
%   NVidia GPU

%% ----------------------------------------------------------------------------

% start unitTest
Success(1) = 1;
fprintf(['Running smi.SMLM.unitTest.\n', ...
         'Testing all smi.SMLM functionality.\n']);

% Test loading various datatypes.
fprintf('Testing loading of various datatypes.\n')

saveName = 'SMLM_testData';

SaveDir = fullfile(tempdir, 'smite', 'unitTest', 'SMLM');
if ~isfolder(SaveDir)
   mkdir(fullfile(tempdir, 'smite'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest', 'SMLM'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest', 'SMLM', saveName));
end

file1 = [saveName, '1.mat'];
file2 = [saveName, '2.mat'];

% Create SMF structure.
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir      = SaveDir;
SMF.Data.DataVariable = 'Data';
SMF.Data.CameraType   = 'EMCCD';
SMF.Data.CameraGain   = 1;
SMF.Data.CameraOffset = 0;
%SMF.Fitting.FitType   = 'XYNB';
%SMF.ROISearchType     = 'Flat';

% Simulate two small datasets.
fprintf('Simulating data.\n')
SimData1 = smi_sim.GaussBlobs.genRandomBlobImage(64, 250);
SimData2 = smi_sim.GaussBlobs.genRandomBlobImage(64, 250);

% Save datasets as mat files.
fprintf('Saving data as mat files.\n')
Data = SimData1;
save(fullfile(SaveDir, file1), 'Data');
Data = SimData2;
save(fullfile(SaveDir, file2), 'Data');

% Try running smi.SMLM.  If it fails, delete files before returning error.
fprintf(['Loading and analyzing data saved as mat files.\n', ...
         '   (Only doing box finding and fitting.)\n']);
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
    delete(fullfile(SaveDir, [saveName, '*.*']));
    fprintf('Caught following error during smi.SMLM.unitTest:\n')
    disp(ME.identifier);
    disp(ME.message);
    Success(1) = 0;
end
fprintf('Loading and analyzing data saved as mat file done.\n');
delete(fullfile(SaveDir, [saveName, '*.*']));

%% ----------------------------------------------------------------------------

Success(2) = 1;
% Save datasets as h5 files.
fprintf('\nSaving data as h5 files.\n')
h5create(fullfile(SaveDir,[saveName '.h5']),'/Data/Channel01/Data0001',size(SimData1));
h5write(fullfile(SaveDir,[saveName '.h5']),'/Data/Channel01/Data0001',SimData1);
h5create(fullfile(SaveDir,[saveName '.h5']),'/Data/Channel01/Data0002',size(SimData2));
h5write(fullfile(SaveDir,[saveName '.h5']),'/Data/Channel01/Data0002',SimData2);
% Try running smi.SMLM.  If it fails, delete files before returning error,
fprintf(['Loading and analyzing data saved as h5 files.\n', ...
         '   (Only doing box finding and fitting.)\n']);
% Create SMF structure.
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir      = SaveDir;
SMF.Data.FileName     = {[saveName '.h5']};
SMF.Data.CameraType   = 'EMCCD';
SMF.Data.CameraGain   = 1;
SMF.Data.CameraOffset = 0;
SMF.Data.FileName = {[saveName, '.h5']};
%SMF.RawImageSize = [size(SimData1,1),size(SimData1,2)];
% Create smi.SMLM object.
SMLMobj = smi.SMLM(SMF);
try
    % Analyze all datasets.
    SMLMobj.analyzeAll();
    clear SMLMobj
catch ME
    delete(fullfile(SaveDir, [saveName '.*']));
    fprintf('Caught following error during smi.SMLM.unitTest:\n')
    disp(ME.identifier)
    disp(ME.message);
    Success(2) = 0;
end
fprintf('Loading and analyzing data saved as h5 file done.\n');
delete(fullfile(SaveDir, [saveName, '.*']));

%% ----------------------------------------------------------------------------

Success(3) = 1;
% Simulate and save realistic SMLM data.
fprintf('\nSimulating realistic 2D SMLM data\n');

SimData1 = smi_sim.GaussBlobs.genRandomBlobImage(256, 1000);
SimData2 = smi_sim.GaussBlobs.genRandomBlobImage(256, 1000);
fprintf('Saving realistic SMSR data.\n');
saveName = 'SMLM_testData';
h5create(fullfile(SaveDir,[saveName '.h5']),'/Data/Channel01/Data0001',size(SimData1));
h5write(fullfile(SaveDir,[saveName '.h5']),'/Data/Channel01/Data0001',SimData1);
h5create(fullfile(SaveDir,[saveName '.h5']),'/Data/Channel01/Data0002',size(SimData2));
h5write(fullfile(SaveDir,[saveName '.h5']),'/Data/Channel01/Data0002',SimData2);

% Test complete data flow
fprintf('Testing 2D analysis.\n');
% Create SMF structure.
SMF = smi_core.SingleMoleculeFitting();
% Cannot save in folders inside SaveDir, so everything must be saved in SaveDir
% directly.
SMF.Data.FileDir          = SaveDir;
SMF.Data.FileName         = {[saveName, '.h5']};
SMF.Data.ResultsDir       = SaveDir;
SMF.Data.CameraType       = 'EMCCD';
SMF.Data.CameraGain       = 1;
SMF.Data.CameraOffset     = 0;
SMF.Thresholding.On       = true;
SMF.Thresholding.MaxSE_XY = 0.1;
SMF.FrameConnection.On    = true;
SMF.DriftCorrection.On    = true;
%SMF.RawImageSize = [size(SimData1,1),size(SimData1,2)];
fprintf('-> ResultsDir = %s\n', SMF.Data.ResultsDir);
% Create smi.SMLM object.
%SMLMobj = smi.SMLM('nogui');
SMLMobj = smi.SMLM(SMF);
try
   % fullFit: fitting -> thresholding -> frame connection -> drift correction
   SMLMobj.fullAnalysis();
   fprintf('Files produced in %s:\n', SMF.Data.ResultsDir);
   dir(fullfile(SMF.Data.ResultsDir, [saveName, '*.*']))
%  % generate color overlay
%  fprintf('Generating color overlay.\n');
%  %SMLMobj.SMD.FitBoxSize = SMLMobj.BoxSize;   % needed for genBlobOverlay
%  %SMLMobj.SMR.FitBoxSize = SMLMobj.BoxSize;   % needed for genBlobOverlay
%  DatasetNum = 1;
%  SMLMobj.genBlobOverlay(DatasetNum);
%  % save
%  SMLMobj.save();
%  SMLMobj.saveResults();
%  SMLMobj.exportFileType='Excel';
%  SMLMobj.exportResults();
%  SMLMobj.exportFileType='txt';
%  SMLMobj.exportResults();
catch ME
   delete(fullfile(SaveDir, saveName, '*.*'));
   rmdir(fullfile(SaveDir, saveName));
   delete(fullfile(SaveDir, [saveName, '*.*']));
   fprintf('Caught following error during smi.SMLM.unitTest:\n')
   disp(ME.identifier)
   disp(ME.message);
   Success(3) = 0;
end
% delete object and data
%clear SMLMobj
%delete(fullfile(SaveDir, saveName, '*.*'));
%rmdir(fullfile(SaveDir, saveName));
%delete(fullfile(SaveDir, [saveName, '*.*']));
fprintf('Done.\n');
