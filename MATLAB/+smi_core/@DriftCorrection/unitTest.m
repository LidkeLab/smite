function [success, SMD2, SMD3, Statistics2, Statistics3] = unitTest()
% unitTest tests smi_core.DriftCorrection.driftCorrectKNN.
% Synthetic data is created by smlmData.m, which then has drift imposed upon
% it.  This data is then drift corrected by driftCorrectKNN, producing a SMD
% structure that has DriftX and DriftY added.
%
% INPUTS:
%   No inputs needed
%
% OUTPUTS:
%   success           0 (failure) or 1 (success)
%   SMD2, SMD3    data structures for 2D or 3D examples with the following
%                 fields:
%       X:            x coordinates (Nx1) where N is total number of points
%       Y:            y coordinates (Nx1)
%       Z:            z coordinates (Nx1)
%       DatasetNum:   dataset number from which localization originates (Nx1)
%       FrameNum:     frame   number from which localization originates (Nx1)
%       NFrames:      number of frames in each dataset
%       NDatasets:    number of datasets
%       DriftX:       found x drift (Nframes x Ndatasets)
%       DriftY:       found y drift (Nframes x Ndatasets)
%       DriftZ:       found z drift (Nframes x Ndatasets)
%   Statistics2:  statistical information about the algorithm performance
%   Statistics3   including various input parameters (2D or 3D)
%
% REQUIRES:
%   Parallel Processing Toolbox
%   Statistics Toolbox
%   NVidia GPU

% Created by:
%   Farzin Farzam (Keith Lidke Lab 2017) [adapted from driftCorrect2D_unitTest]
%   Michael Wester (Lidke Lab 2017/2018)

success = 0;

close all


SaveDir = fullfile(tempdir, 'smite', 'unitTest', 'DriftCorrection');
if ~isfolder(SaveDir)
   mkdir(fullfile(tempdir, 'smite'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest', 'DriftCorrection'));
end

SIM = smi_sim.SimSMLM();
SIM.SZ = 256;
SIM.Rho = 10;
SIM.NDatasets = 10;
SIM.NFrames = 100;
SIM.simStar(16);
SMDsim = SIM.SMD_Model;

XYSize = SIM.SZ;
n_frames = SIM.NDatasets * SIM.NFrames;
PpFX = 0.04;      % x drift (pixels per frame)
PpFY = 0.07;      % y drift (pixels per frame)

rho = SIM.Rho;
FpD = SIM.NFrames; % number of frames per dataset

%load('SMDsim2D');

fprintf('2D\n');
n_particles = numel(SMDsim.X);
fprintf('Number of emitters = %d, per pixel = %f, per dataset = %f\n', ...
        n_particles, n_particles / XYSize^2, n_particles / (n_frames / FpD));
[n_emitters, n_blinks, n_local, n_datasets] = ...
   blinks(n_particles, SMDsim.X, SMDsim.Y, SMDsim.Y, ...
          SMDsim.FrameNum, SMDsim.DatasetNum);
fprintf( ...
   '# of emitters = %d, blinks = %d, localizations = %d, datasets = %d\n', ...
   n_emitters, n_blinks, n_local, n_datasets);
PixelSizeZUnit = 0.1; % um
P2nm = PixelSizeZUnit * 1000;

% Re-organize so X and Y are the first inputs.
SMDin = SMDsim;
SMDin.X_SE = ones(n_particles, 1);
SMDin.Y_SE = ones(n_particles, 1);
SMDin.Bg = zeros(n_particles, 1);
SMDin.PixelSizeZUnit = PixelSizeZUnit;

% Histogram image of this data.
SRImageZoom = 4;
TrueIm = smi_vis.GenerateImages.histogramImage(SMDin, SRImageZoom); %new
P = prctile(TrueIm(TrueIm > 0), 99.9);
TrueIm(TrueIm > P) = P;
TrueIm = 255 * TrueIm / P;
figure; imagesc(TrueIm); colormap(gray); % what we should get afterwards
saveas(gcf, fullfile(SaveDir, 'TrueIm.png'));

X_True = single(SMDin.X);
Y_True = single(SMDin.Y);
% Creating drift per frame: fast way!
frame_num = SMDin.NFrames*(SMDin.DatasetNum - 1) + SMDin.FrameNum - 1;
SMDin.X = SMDin.X + frame_num*PpFX;
SMDin.Y = SMDin.Y + frame_num*PpFY;
% Boundary condtions.
for i = 1 : numel(SMDin.X)
   if SMDin.X(i) >= XYSize;
%     SMDin.X(i) = XYSize;
   end
   if SMDin.Y(i) >= XYSize;
%     SMDin.Y(i) = XYSize;
   end
end
SMDsave = SMDin;

DriftIm = smi_vis.GenerateImages.histogramImage(SMDin, SRImageZoom); %new
P = prctile(DriftIm(DriftIm > 0), 99.9);
DriftIm(DriftIm > P) = P;
DriftIm = 255 * DriftIm / P;
figure; imagesc(DriftIm); colormap(gray);  %synthetic drift image
saveas(gcf, fullfile(SaveDir, 'DriftIm.png'));
%GaussIm = smi_vis.GenerateImages.gaussianImage(SMDin, SRImageZoom);
%figure; imagesc(GaussIm); colormap(gray);

SMF = smi_core.SingleMoleculeFitting();
SMF.DriftCorrection.BFRegistration = false;
DC = smi_core.DriftCorrection(SMF, SMDin);

clear SMD
[SMD, Statistics] = DC.driftCorrectKNN(SMDin);
%Statistics

% Remove any NaNs.
nans = find(isnan(SMDin.X) | isnan(SMDin.Y) | isnan(SMD.X) | isnan(SMD.Y));
n_nans = numel(nans);
if n_nans > 0
   fprintf('%d NaNs removed!\n', n_nans);
   SMDin.X(nans) = [];
   SMDin.Y(nans) = [];
   SMD.X(nans) = [];
   SMD.Y(nans) = [];
   X_True(nans) = [];
   X_Yrue(nans) = [];
end

% Consistency check.
N = numel(SMD.X);
X_inDC = zeros(N, 1, 'single');
Y_inDC = zeros(N, 1, 'single');
X_unDC = zeros(N, 1, 'single');
Y_unDC = zeros(N, 1, 'single');
for k = 1:N
   i = SMD.FrameNum(k);
   j = SMD.DatasetNum(k);
   X_inDC(k) = SMDin.X(k) - SMD.DriftX(i, j);
   Y_inDC(k) = SMDin.Y(k) - SMD.DriftY(i, j);
   X_unDC(k) = SMD.X(k) + SMD.DriftX(i, j);
   Y_unDC(k) = SMD.Y(k) + SMD.DriftY(i, j);
end
consistency_un = sum(abs(SMDin.X - X_unDC) + abs(SMDin.Y - Y_unDC));
consistency_in = sum(abs(SMD.X - X_inDC) + abs(SMD.Y - Y_inDC));
fprintf('SMDin.X/Y - (SMD.X/Y + SMD.DriftX/Y) = %f nm\n', ...
        consistency_un * P2nm);
fprintf('SMD.X/Y - (SMDin.X/Y - SMD.DriftX/Y) = %f nm\n', ...
        consistency_in * P2nm);

correctedDriftIm = smi_vis.GenerateImages.histogramImage(SMD, SRImageZoom);
%figure; imagesc(DriftIm); colormap(gray);
% Clean up the sum image by setting the 0.1% top intensity pixels to the
% 99.9% intensity value.
P = prctile(correctedDriftIm(correctedDriftIm > 0), 99.9);
correctedDriftIm(correctedDriftIm > P) = P;
correctedDriftIm = 255 * correctedDriftIm / P;
figure; imagesc(correctedDriftIm); colormap(gray);
saveas(gcf, fullfile(SaveDir, 'correctedDriftIm.png'));

% Plot the drift correction as a function of time.
DC_fig = DC.plotDriftCorrection(SMD);
figure(DC_fig);
saveas(gcf, fullfile(SaveDir, 'DC_fig.png'));

% Compute absolute drift in pixels per frame.
x_drift_true = PpFX .* (1 : n_frames);
y_drift_true = PpFY .* (1 : n_frames);

% Compute the RMSE between the pre-drift data and the drift corrected
% post-drift data.
[dist1, rmse1, dist2, rmse2, ~] = ...
   smi_core.DriftCorrection.calcDCRMSE(SMD, X_True, Y_True, [], ...
                                       x_drift_true, y_drift_true, []);
fprintf('average distance between true and DC locations = %f nm\n', dist1);
fprintf('RMSE1            between true and DC locations = %f nm\n', rmse1);
fprintf('average distance between true and DC curves    = %f nm\n', dist2);
fprintf('RMSE2            between true and DC curves    = %f nm\n', rmse2);

% Compare computed vs. true drift.
base = 0;
framenums = [];
for j = 1:SMDin.NDatasets
   framenums([1:SMDin.NFrames] + base) = ...
      arrayfun(@(i) SMDin.NFrames*(j - 1) + i - 1, 1:SMDin.NFrames);
   base = base + SMDin.NFrames;
end
x_drift = mean(SMD.DriftX(:) ./ (framenums(:) + 1));
y_drift = mean(SMD.DriftY(:) ./ (framenums(:) + 1));
%fprintf('average x-drift per frame = %f px (true = %f px)\n', x_drift, PpFX);
%fprintf('average y-drift per frame = %f px (true = %f px)\n', y_drift, PpFY);
fprintf('average x-drift per frame = %f nm (true = %f nm)\n', ...
        x_drift * P2nm, PpFX * P2nm);
fprintf('average y-drift per frame = %f nm (true = %f nm)\n', ...
        y_drift * P2nm, PpFY * P2nm);

SMD2 = SMD;
Statistics2.Sim_rho = rho;
Statistics2.N_particles = n_particles;
Statistics2.N_particles_per_pixel   = n_particles / XYSize^2;
Statistics2.N_particles_per_dataset = n_particles / (n_frames / FpD);
Statistics2.N_NaNs = n_nans;
Statistics2.Consistency_un = consistency_un;
Statistics2.Consistency_in = consistency_in;
Statistics2.Dist1 = dist1;
Statistics2.RMSE1 = rmse1;
Statistics2.Dist2 = dist2;
Statistics2.RMSE2 = rmse2;
Statistics2.DriftX_mean = x_drift;
Statistics2.DriftY_mean = y_drift;
Statistics2.DriftX_True = PpFX;
Statistics2.DriftY_True = PpFY;

% -----------------------------------------------------------------------------

% Perform the same calculation, but now separate the intra-dataset and
% inter-dataset portions.  Simulate the` situation of calling the two portions
% from separate functions.
fprintf('\n2D: separated intra-dataset and inter-dataset\n');
clear SMD
SMDin = SMDsave;
SMDIntra = [];
SMDtmp   = [];
X_TrueTmp = [];
Y_TrueTmp = [];
SMF = smi_core.SingleMoleculeFitting();
SMF.DriftCorrection.BFRegistration = false;
obj.DC = smi_core.DriftCorrection(SMF);
for i = 1 : SMDin.NDatasets
   SMDin_i = SMDin;
   SMDin_i.NDatasets = 1;
   mask = SMDin.DatasetNum == i;
   SMDin_i.X          = single(SMDin.X(mask));
   SMDin_i.Y          = single(SMDin.Y(mask));
   SMDin_i.X_SE       = SMDin.X_SE(mask);
   SMDin_i.Y_SE       = SMDin.Y_SE(mask);
   SMDin_i.DatasetNum = SMDin.DatasetNum(mask);
   SMDin_i.FrameNum   = SMDin.FrameNum(mask);
   SMDin_i.Photons    = SMDin.Photons(mask);
   SMDin_i.Bg         = SMDin.Bg(mask);
   X_TrueTmp = [X_TrueTmp; X_True(mask)];
   Y_TrueTmp = [Y_TrueTmp; Y_True(mask)];
   [SMDIntra_i, StatisticsIntra] = obj.DC.driftCorrectKNNIntra(SMDin_i, i, i);
   SMDtmp   = smi_core.SingleMoleculeData.catSMD(SMDtmp, SMDin_i, false);
   SMDIntra = smi_core.SingleMoleculeData.catSMD(SMDIntra, SMDIntra_i, false);
end
[SMDInter, StatisticsInter] = obj.DC.driftCorrectKNNInter(SMDIntra);
SMD   = SMDInter;
SMDin = SMDtmp;
X_True = X_TrueTmp;
Y_True = Y_TrueTmp;

% Remove any NaNs.
nans = find(isnan(SMDin.X) | isnan(SMDin.Y) | isnan(SMD.X) | isnan(SMD.Y));
n_nans = numel(nans);
if n_nans > 0
   fprintf('%d NaNs removed!\n', n_nans);
   SMDin.X(nans) = [];
   SMDin.Y(nans) = [];
   SMD.X(nans) = [];
   SMD.Y(nans) = [];
   X_True(nans) = [];
   X_Yrue(nans) = [];
end

% Consistency checks.  These should be 0 to a few hundred if all is correct.
%
% SMDin.X/Y, SMD.X/Y are the drifted/drift corrected coordinates, respectively.
% SMD.DriftX/Y are the drift corrections defined such that
%    drifted coordinates - drift correction = drift corrected coordinates
N = numel(SMD.X);
X_inDC = zeros(N, 1, 'single');
Y_inDC = zeros(N, 1, 'single');
X_unDC = zeros(N, 1, 'single');
Y_unDC = zeros(N, 1, 'single');
for k = 1:N
   i = SMD.FrameNum(k);
   j = SMD.DatasetNum(k);
   X_inDC(k) = SMDin.X(k) - SMD.DriftX(i, j);
   Y_inDC(k) = SMDin.Y(k) - SMD.DriftY(i, j);
   X_unDC(k) = SMD.X(k) + SMD.DriftX(i, j);
   Y_unDC(k) = SMD.Y(k) + SMD.DriftY(i, j);
end
% consistency_un =
%    drifted coordinates - (drift corrected coordinates + drift correction)
% consistency_in =
%    drifted corrected coordinates - (drifted coordinates - drift correction)
consistency_un = sum(abs(SMDin.X - X_unDC) + abs(SMDin.Y - Y_unDC));
consistency_in = sum(abs(SMD.X - X_inDC) + abs(SMD.Y - Y_inDC));
fprintf('SMDin.X/Y - (SMD.X/Y + SMD.DriftX/Y) = %f nm\n', ...
        consistency_un * P2nm);
fprintf('SMD.X/Y - (SMDin.X/Y - SMD.DriftX/Y) = %f nm\n', ...
        consistency_in * P2nm);

correctedDriftIm = smi_vis.GenerateImages.histogramImage(SMD, SRImageZoom);
%figure; imagesc(DriftIm); colormap(gray);
% Clean up the sum image by setting the 0.1% top intensity pixels to the
% 99.9% intensity value.
P = prctile(correctedDriftIm(correctedDriftIm > 0), 99.9);
correctedDriftIm(correctedDriftIm > P) = P;
correctedDriftIm = 255 * correctedDriftIm / P;
figure; imagesc(correctedDriftIm); colormap(gray);
saveas(gcf, fullfile(SaveDir, 'correctedDriftIm2.png'));

% Plot the drift correction as a function of time.
DC_fig = DC.plotDriftCorrection(SMD);
figure(DC_fig);
saveas(gcf, fullfile(SaveDir, 'DC_fig2.png'));

% Compute absolute drift in pixels per frame.
x_drift_true = PpFX .* (1 : n_frames);
y_drift_true = PpFY .* (1 : n_frames);

% Compute the RMSE between the pre-drift data and the drift corrected
% post-drift data.
[dist1, rmse1, dist2, rmse2, ~] = ...
   smi_core.DriftCorrection.calcDCRMSE(SMD, X_True, Y_True, [], ...
                                       x_drift_true, y_drift_true, []);
fprintf('average distance between true and DC locations = %f nm\n', dist1);
fprintf('RMSE1            between true and DC locations = %f nm\n', rmse1);
fprintf('average distance between true and DC curves    = %f nm\n', dist2);
fprintf('RMSE2            between true and DC curves    = %f nm\n', rmse2);

% Compare computed vs. true drift.
base = 0;
framenums = [];
for j = 1:SMDin.NDatasets
   framenums([1:SMDin.NFrames] + base) = ...
      arrayfun(@(i) SMDin.NFrames*(j - 1) + i - 1, 1:SMDin.NFrames);
   base = base + SMDin.NFrames;
end
x_drift = mean(SMD.DriftX(:) ./ (framenums(:) + 1));
y_drift = mean(SMD.DriftY(:) ./ (framenums(:) + 1));
%fprintf('average x-drift per frame = %f px (true = %f px)\n', x_drift, PpFX);
%fprintf('average y-drift per frame = %f px (true = %f px)\n', y_drift, PpFY);
fprintf('average x-drift per frame = %f nm (true = %f nm)\n', ...
        x_drift * P2nm, PpFX * P2nm);
fprintf('average y-drift per frame = %f nm (true = %f nm)\n', ...
        y_drift * P2nm, PpFY * P2nm);

success = 1;

end

% =============================================================================

function [n_emitters, n_blinks, n_local, n_datasets] = ...
   blinks(n_particles, X, Y, Z, F, D)
% Find all the frame number sets (fns) for each emitter (of which there are
% n_emitters).
%
% INPUTS:
%    n_particles   total number of particles (localizations) given
%    X, Y, Z       particle coordinates
%    F             absolute frame numbers
%    D             dataset numbers
%
% OUTPUTS:
%    n_emitters    number of distinct locations, each corresponding to an
%                  emitter
%    n_blinks      number of blinking events
%    n_local       number of entries
%    n_datasets    number of distinct datasets
%
% Data from smlmData with no imposed drift looks like:
%
%     X         Y        F
%     2.9100    1.1548   22.0000
%     2.9100    1.1548   23.0000
%     2.9100    1.1548   24.0000
%     2.9100    1.1548  452.0000
%    17.6154    1.1999  477.0000
%    58.0133    1.5341   40.0000
%    58.0133    1.5341   41.0000
%
% This corresponds to 3 emitters (distinct positions), 4 blinking events (2 for
% the first emitter, 1 each for the other two), and 7 localizations (before
% frame connection).

   fns = {};
   i = 1;
   n_emitters = 0;
   while i <= n_particles
      x = X(i);
      y = Y(i);
      z = Z(i);
      f = F(i);

      i = i + 1;
      while i <= n_particles & X(i) == x & Y(i) == y & Z(i) == z
         f = [f, F(i)];
         i = i + 1;
      end
      n_emitters = n_emitters + 1;
      fns{n_emitters} = f;
   end

   % Count the number of localizations.
   n_local = n_particles;

   % Count the number of blinking events.
   n_blinks = sum(cellfun(@(e) numel(find(diff(e) ~= 1)) + 1, fns));

   % Count the number of datasets.
   n_datasets = numel(unique(D));

end
