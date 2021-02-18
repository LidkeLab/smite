function SimDimers()

%Simulate data with a given receptor density and labeling efficiency, 
%then process it through smite SR fitting and finally BaGoL.

clear all;
close all;
SaveDir = 'C:\Users\sajja\Documents\MATLAB\Dimers'; %Change save directory
kk = 2;                  % dimers
PixSize = 100;           % nm in a pixel
radius_kTet = 50/PixSize; % Radius of the k-tets (pixel)
SE = 20/PixSize;         % Standard error of a localization (pixel)
receptorDensity = 60;    % Receptors/um^2
% Units: (hextets/um^2) / ((nm/um) / (nm/pixel))^2 -> hextets/pixel^2
kTetDensity = (receptorDensity/kk) / (1000 / PixSize)^2;
n_Labels = 1;            % Number of different types of labels (1 or 2)
LabelingEfficiency = 1.0;% Labeling efficiency in the range [0, 1] (0.5-1)

SIM = smi_sim.SimSMLM();
%SIM.Rho = 30;            % Fluorophore Density (flourophore/pixel)
SIM.SZ = 256;            % Linear size of image (pixels)
%SIM.NDatasets = 50;      % Number of datasets
SIM.NFrames = 1000;      % Number of frames per dataset
SIM.ZoomFactor = 5;      % It can be either smaller or larger than one
SIM.K_OnToOff = 1;       % Fluorophore turns Off from On state (frames^-1)
SIM.K_OffToOn = 0.0005;  % Fluorophore return to On state from Off (frames^-1)
SIM.K_OnToBleach = 0.2;  % Fluorophore bleached out (frames^-1)
SIM.EmissionRate = 1000; % Emission rate (Intensity) of photons (photons/frame)
SIM.Bg = 15;             % Background Count Rate (counts/pixel)
SIM.PSFSigma = 1.3;      % Point Spread Function Sigma size (Pixels).
SaveFile = sprintf('recpDensity=%d,LabelingEff=%d.mat', ...
                   receptorDensity, 100 * LabelingEfficiency);
% Generate k-tets in the simulation region (units are pixels).
n_kTets = round(kTetDensity * SIM.SZ^2);
fprintf('Generating %d %d-tets ...\n', n_kTets, kk);
center_kTet = rand(1, 2) * SIM.SZ;
SMD_True = SIM.kTet(kk, center_kTet, radius_kTet);
Label = repmat(1:n_Labels, 1, kk/n_Labels)';
for i = 2 : n_kTets
   center_kTet = rand(1, 2) * SIM.SZ;
   SMD_True_tmp = SIM.kTet(kk, center_kTet, radius_kTet);
   Label_tmp = repmat(1:n_Labels, 1, kk/n_Labels)';
   SMD_True = smi_core.SingleMoleculeData.catSMD(SMD_True, SMD_True_tmp);
   Label = [Label; Label_tmp];
end

% Apply labeling efficiency.
Labeled = ones(size(Label));
Labeled(rand(size(Label)) > LabelingEfficiency) = 0;
SMD_True_Labeled = SMD_True;
SMD_True_Labeled.X(~Labeled) = [];
SMD_True_Labeled.Y(~Labeled) = [];
fprintf('%d%% labeling efficiency: %d -> %d localizations\n', ...
        round(100 * LabelingEfficiency), numel(Labeled), sum(Labeled));
    
% Plot labeled localizations.
figure;
plot(SMD_True_Labeled.X, SMD_True_Labeled.Y, 'k.');
hold on
title('SMD\_True\_Labeled');
xlabel('x (pixels)');
ylabel('y (pixels)');
hold off
SimDir = fullfile(SaveDir, 'Results', 'Results_BaGoL', ...
                  regexprep(SaveFile, '.mat$', ''));
if ~isfolder(SimDir)
   mkdir(SimDir);
end
saveas(gcf, fullfile(SimDir, 'True_Labeled'), 'fig');
close

% Generate blinks (units are pixels).
SMD_Model = SIM.genBlinks(SMD_True_Labeled, 'Equib'); 
% Below needed for generating blob images.
SMD_Model.X_SE = SE * ones(size(SMD_Model.X));
SMD_Model.Y_SE = SE * ones(size(SMD_Model.Y));
fprintf('SMD_Model # of frames = %d\n', numel(SMD_Model.FrameNum));

% Generate an SMD structure with positional and intensity noise.  (This 
% bypasses genNoisyData and SR analysis of the generated Data.)
% SMD_Model_Noisy = SIM.genNoisySMD(SMD_Model); 

% Generate the blobs without Poisson noise (units are pixels).
SMD_Model.FitBoxSize = ceil(4 * 2 * SIM.PSFSigma);
Model = smi_sim.GaussBlobs.gaussBlobImage(SIM.SZ, SIM.NFrames, SMD_Model, ...
                                          0, 0, 0);
SumModel = sum(Model, 3);
dipshow(SumModel, 'lin');
diptruesize(gcf, 200);

% Generate the blobs having Poisson noise.
Data = SIM.genNoisyData(Model);
SumData = sum(Data, 3);
dipshow(SumData, 'lin');
diptruesize(gcf, 200);

sequence = Data;
if ~isfolder(SaveDir)
   mkdir(SaveDir);
end
save(fullfile(SaveDir, SaveFile), 'sequence');
close all
 
% ---------- Perform SR analysis

SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir        = SaveDir;
SMF.Data.FileName       = {SaveFile};
SMF.Data.ResultsDir     = fullfile(SaveDir, 'Results');
SMF.Data.CameraType     = 'EMCCD';
% DIPimage cal_readnoise may be helpful for CameraGain and CameraOffset.
SMF.Data.CameraGain     = 1;
SMF.Data.CameraOffset   = 0;
SMF.Thresholding.On     = true;
SMF.Thresholding.MinPValue = 0.005;
SMF.FrameConnection.On  = true;
SMF.DriftCorrection.On  = true;
SMLMobj = smi.SMLM(SMF);
SMLMobj.fullAnalysis();
%SMLMobj.analyzeAll();

% ---------- Perform BaGoL analysis

BaGoLParams.PixelSize = PixSize;    % (nm)
BaGoLParams.OutputPixelSize = 4;    % pixel size for posterior images (nm)
BaGoLParams.IntensityCutoff = 5000; % Intensity cutoff
BaGoLParams.NNR = inf;              % At least NNN localizations within NNR(nm)
BaGoLParams.NNN = 0;                % Number of Nearest Neighbors
BaGoLParams.DriftErr = 1;           % Precision inflation applied to SE (nm)
BaGoLParams.ClusterDrift = 0;       % Expected magnitude of drift (nm/frame)
BaGoLParams.ROIsz = 500;            % ROI size for RJMCMC (nm)
BaGoLParams.OverLap = 50;           % Size of overlapping region (nm)
BaGoLParams.PreCluster = 5;         % Pre-clustering parameter (nm)
BaGoLParams.Lambda = [1, 1];        % [gamma, eta] parameters for gamma prior
% DataROI can speed up the BaGoL 1st pass in certain situations.  If DataROI is
% empty, the whole image is used, otherwise the specified portion of the image
% is used to compute the prior for lambda for the 2nd pass.
%BaGoLParams.DataROI = [100, 200, 100, 200];   % medium dense
%BaGoLParams.DataROI = [100, 150, 100, 150];   % very dense
BaGoLParams.DataROI = [];           % 1st pass [Xmin, Xmax, Ymin, Ymax] (pixel)
Model_Noisy = smi_sim.GaussBlobs.gaussBlobImage(SIM.SZ, SIM.NFrames, ...
                                                SMD_Model_Noisy, 0, 0, 0);
                                            
SaveDir = fullfile(SMF.Data.ResultsDir, 'Results_BaGoL');
[~, SaveFile] = fileparts(SaveFile);
SaveFile = sprintf('%s_Results.mat', SaveFile);
BGL = BaGoL_analysis(SaveFile, SMF.Data.ResultsDir, SaveDir, BaGoLParams);
BGL.genSRMAPNOverlay(SMD_True_Labeled, BGL.MAPN, 256, 256, PixelSize,
                     SaveDirLong,0,0,2);
fprintf('Done.\n');
end