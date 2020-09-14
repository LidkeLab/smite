function [SMD2, SMD3, Statistics2, Statistics3] = unitTest()
% unitTest tests smi_core.DriftCorrection.driftCorrectKNN.
% Synthetic data is created by smlmData.m, which then has drift imposed upon
% it.  This data is then drift corrected by driftCorrectKNN, producing a SMD
% structure that has DriftX and DriftY added.
%
% INPUTS:
%   No inputs needed
%
% OUTPUTS:
%   SMD           data structure with the following fields:
%       X:            x coordinates (Nx1) where N is total number of points
%       Y:            y coordinates (Nx1)
%       DatasetNum:   dataset number from which localization originates (Nx1)
%       FrameNum:     frame   number from which localization originates (Nx1)
%       NFrames:      number of frames in each dataset
%       NDatasets:    number of datasets
%       DriftX:       found x drift (NFrames x NDatasets)
%       DriftY:       found y drift (NFrames x NDatasets)
%   Statistics:   statistical information about the algorithm performance
%                 including various input parameters
%
% REQUIRES:
%   DIPimage Toolbox
%   Parallel Processing Toolbox
%   Statistics Toolbox
%   NVidia GPU
%
% CITATION:
%   Farzin Farzam (Keith Lidke Lab 2017) [adapted from driftCorrect2D_unitTest]
%   Michael Wester (Lidke Lab 2017/2018)

%shape = 'uniform'; % pattern type
shape = 'star';    % pattern type
yn = 'No';         % smlmData compute drift?
%rho = 100;         % fluorophore density (fluorphore/pixel)
rho = 10;           % fluorophore density (fluorphore/pixel)
XYSize = 64;       % linear size of image (pixels)
n_frames = 1000;   % number of frames
%PpFX = 0.004*rand; % x-drift (pixels per frame)
%PpFY = 0.007*rand; % y-drift (pixels per frame)
PpFX = 0.004;      % x drift (pixels per frame)
PpFY = 0.007;      % y drift (pixels per frame)
FpD = 100;         % number of frames per dataset
K_on =  0.0005;    % off to on rate (1/frames)

[Data, SMDsim] = SMA_Sim.smlmData(shape, yn, 1, K_on, 1, 2000, rho, ...
                                  1.3, 15, XYSize, n_frames, 'Equib');
fprintf('2D\n');
fprintf('rho = %d fluorophores/pixel\n', rho);
n_particles = numel(SMDsim.X);
fprintf('Number of emitters = %d, per pixel = %f, per dataset = %f\n', ...
        n_particles, n_particles / XYSize^2, n_particles / (n_frames / FpD));
[n_local, n_blinks] = blinks(n_particles, SMDsim.X, SMDsim.Y, SMDsim.FrameNum);
fprintf('Number of localizations = %d, blinks = %d\n', n_local, n_blinks);
SMD.XSize = XYSize;
SMD.YSize = XYSize;
SMD.X_SE = zeros(n_particles, 1);
SMD.Y_SE = zeros(n_particles, 1);
SMD.Bg = 20*rand(n_particles, 1);
X = single(SMDsim.X);
Y = single(SMDsim.Y);
Photons  = SMDsim.Photons;
FrameNum = SMDsim.FrameNum;
DatasetNum = ones(size(SMDsim.FrameNum));
for i = 1 : numel(DatasetNum)
   DatasetNum(i) = (FrameNum(i) - 1 - mod(FrameNum(i) - 1, FpD)) / FpD + 1;
   FrameNum(i)   = mod(FrameNum(i) - 1, FpD) + 1;
end
PixelSizeZUnit = 0.1; % um

P2nm = PixelSizeZUnit * 1000;

% Re-organize so X and Y are the first inputs.
SMDin.XSize = XYSize;
SMDin.YSize = XYSize;
SMDin.X = X;
SMDin.Y = Y;
SMDin.DatasetNum = DatasetNum;
SMDin.FrameNum   = FrameNum;
SMDin.NDatasets = n_frames / FpD;
SMDin.NFrames   = n_frames / SMDin.NDatasets;
SMDin.Photons = Photons;
SMDin.X_SE = SMD.X_SE;
SMDin.Y_SE = SMD.Y_SE;
SMDin.Bg = SMD.Bg;
SMDin.PixelSizeZUnit = PixelSizeZUnit;

% Histogram image of this data.
SRImageZoom = 4;
TrueIm = SMA_Vis.histogramImage(SMDin, SRImageZoom); %new
P = prctile(TrueIm(TrueIm > 0), 99.9);
TrueIm(TrueIm > P) = P;
TrueIm = 255 * TrueIm / P;
dipshow(TrueIm); % what we should get afterwards

X_True = single(SMDin.X);
Y_True = single(SMDin.Y);
if strcmp(yn, 'No')
   % Creating drift per frame: fast way!
   frame_num = SMDin.NFrames*(SMDin.DatasetNum - 1) + SMDin.FrameNum - 1;
   SMDin.X = SMDin.X + frame_num*PpFX;
   SMDin.Y = SMDin.Y + frame_num*PpFY;
   % Boundary condtions.
   for i = 1 : numel(SMDin.X)
      if SMDin.X(i) >= XYSize;
         SMDin.X(i) = XYSize;
      end
      if SMDin.Y(i) >= XYSize;
         SMDin.Y(i) = XYSize;
      end
   end
end
SMDsave = SMDin;

DriftIm = SMA_Vis.histogramImage(SMDin, SRImageZoom); %new
P = prctile(DriftIm(DriftIm > 0), 99.9);
DriftIm(DriftIm > P) = P;
DriftIm = 255 * DriftIm / P;
dipshow(DriftIm)  %synthetic drift image
%GaussIm = SMA_Vis.gaussianImage(SMDin, SRImageZoom);
%dipshow(GaussIm);

SMF = smi_core.SingleMoleculeFitting.createSMF();
SMF.DriftCorrection.Init_inter = SMDin.NFrames;
DC = smi_core.DriftCorrection(SMF);
%DriftParams.PDegree       = 1;
%DriftParams.TolFun_intra  = 1e-2;
%DriftParams.TolX_intra    = 1e-4;
%DriftParams.TolFun_inter  = 1e-2;
%DriftParams.TolX_inter    = 1e-4;
%DC.DriftParams.Init_inter = SMDin.NFrames;
%DriftParams.Init_inter    = 0;

clear SMDin SMD
SMDin = SMDsave;
[SMD, Statistics] = DC.driftCorrectKNN(SMDin);
%Statistics

% Remove any NaNs.
nans = find(isnan(SMDsave.X) | isnan(SMDsave.Y) | isnan(SMD.X) | isnan(SMD.Y));
n_nans = numel(nans);
if n_nans > 0
   fprintf('%d NaNs removed!\n', n_nans);
   SMDsave.X(nans) = [];
   SMDsave.Y(nans) = [];
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
   X_inDC(k) = SMDsave.X(k) - SMD.DriftX(i, j);
   Y_inDC(k) = SMDsave.Y(k) - SMD.DriftY(i, j);
   X_unDC(k) = SMD.X(k) + SMD.DriftX(i, j);
   Y_unDC(k) = SMD.Y(k) + SMD.DriftY(i, j);
end
consistency_un = sum(abs(SMDin.X - X_unDC) + abs(SMDin.Y - Y_unDC));
consistency_in = sum(abs(SMD.X - X_inDC) + abs(SMD.Y - Y_inDC));
fprintf('SMDin.X/Y - (SMD.X/Y + SMD.DriftX/Y = %f nm\n', ...
        consistency_un * P2nm);
fprintf('SMD.X/Y - (SMDin.X/Y - SMD.DriftX/Y = %f nm\n', ...
        consistency_in * P2nm);

correctedDriftIm = SMA_Vis.histogramImage(SMD, SRImageZoom);
%dipshow(DriftIm)
% Clean up the sum image by setting the 0.1% top intensity pixels to the
% 99.9% intensity value.
P = prctile(correctedDriftIm(correctedDriftIm > 0), 99.9);
correctedDriftIm(correctedDriftIm > P) = P;
correctedDriftIm = 255 * correctedDriftIm / P;
dipshow(correctedDriftIm)

% Plot the drift correction as a function of time.
DC_fig = DC.plotDriftCorrection(SMD);
showm(DC_fig);

% Plot the residual between the pre-drift data and the drift corrected
% post-drift data.
[residual, dist, rmse, nnfig] = DC.calcDCResidual(SMD, X_True, Y_True);
%dipshow(residual);
fprintf('average distance between true and DC locations = %f nm\n', ...
        dist * P2nm);
fprintf('RMSE             between true and DC locations = %f nm\n', ...
        rmse * P2nm);

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
Statistics2.Dist = dist;
Statistics2.RMSE = rmse;
Statistics2.DriftX_mean = x_drift;
Statistics2.DriftY_mean = y_drift;
Statistics2.DriftX_True = PpFX;
Statistics2.DriftY_True = PpFY;

% -----------------------------------------------------------------------------
% 3D

%shape = 'uniform'; % pattern type
shape = 'star';    % pattern type
yn = 'No';         % smlmData compute drift?
rho = 10;          % fluorophore density (fluorphore/pixel)
XYSize = 64;       % linear size of image (pixels)
n_frames = 1000;   % number of frames
%PpFX = 0.004*rand; % x-drift (pixels per frame)
%PpFY = 0.007*rand; % y-drift (pixels per frame)
PpFX = 0.004;      % x drift (pixels per frame)
PpFY = 0.007;      % y drift (pixels per frame)
PpFZ = 0.010;      % z drift (pixels per frame)
FpD = 100;         % number of frames per dataset

%% generating PSF-array
XYSamPerPix = 6;
ZSamPerUnit = 10;
%Create the Sampled PSF using SMA_PSF.scalarPSFZernike
PSFStruct=SMA_PSF.createPSFStruct();
PSFStruct.ZC_Phase(6)=1; %Astigmatism
PSFStruct.PixelSize=.12/XYSamPerPix; %micron from RB setup
PSFStruct.Z=(-0.7:1/ZSamPerUnit:0.7);
PSFStruct.SZ = 64*XYSamPerPix;
PSFStruct.OSZ = 128*XYSamPerPix;
[PSF,PSFStruct_Out]=SMA_PSF.scalarPSFZernike(PSFStruct);
PSF = gather(PSF/sum(PSF(:)));
clear PSFStruct
%% simulate 3D data
PSF_struct.PSF = PSF;
PSF_struct.XYSamPerPix = XYSamPerPix;
PSF_struct.ZSamPerUnit = ZSamPerUnit;

[Data, SMDsim] = SMA_Sim.smlmData(shape, yn ,1, 0.0005, 1, 2000, rho, ...
                                  PSF_struct, 15, XYSize, n_frames, 'Equib');
fprintf('\n3D\n');
fprintf('rho = %d fluorophores/pixel\n', rho);
n_particles = numel(SMDsim.X);
fprintf('Number of particles = %d, per pixel = %f, per dataset = %f\n', ...
        n_particles, n_particles / XYSize^2, n_particles / (n_frames / FpD));
[n_local, n_blinks] = blinks3(n_particles, SMDsim.X, SMDsim.Y, SMDsim.Z, ...
                              SMDsim.FrameNum);
fprintf('Number of localizations = %d, blinks = %d\n', n_local, n_blinks);

SMD.XSize = XYSize;
SMD.YSize = XYSize;
SMD.X_SE = zeros(n_particles, 1);
SMD.Y_SE = zeros(n_particles, 1);
SMD.Z_SE = zeros(n_particles, 1);
SMD.Bg = 20*rand(n_particles, 1);
X = single(SMDsim.X);   % pixels
Y = single(SMDsim.Y);   % pixels
Z = single(SMDsim.Z);   % um
Photons  = SMDsim.Photons;
FrameNum = SMDsim.FrameNum;
DatasetNum = ones(size(SMDsim.FrameNum));
for i = 1 : numel(DatasetNum)
   DatasetNum(i) = (FrameNum(i) - 1 - mod(FrameNum(i) - 1, FpD)) / FpD + 1;
   FrameNum(i)   = mod(FrameNum(i) - 1, FpD) + 1;
end

%re-organize so X and Y are the first inputs
SMDin.XSize = XYSize;
SMDin.YSize = XYSize;
SMDin.PixelSizeZUnit = 0.1; % um per pixel
SMDin.X = X;
SMDin.Y = Y;
SMDin.Z = Z;
SMDin.DatasetNum = DatasetNum;
SMDin.FrameNum   = FrameNum;
SMDin.NDatasets = n_frames / FpD;
SMDin.NFrames   = n_frames / SMDin.NDatasets;
SMDin.Photons = Photons;
SMDin.X_SE = SMD.X_SE;
SMDin.Y_SE = SMD.Y_SE;
SMDin.Z_SE = SMD.Z_SE;
SMDin.Bg = SMD.Bg;

figure();
hold on
cm = colormap(jet);
n_cm = size(cm, 1);
sz = 5;
% Linear color scale.
min_Z = min(SMDin.Z);
max_Z = max(SMDin.Z);
indx = max(1, ceil((n_cm - 1) * (SMDin.Z - min_Z) / (max_Z - min_Z) + 1));
scatter3(SMDin.X * P2nm, SMDin.Y * P2nm, SMDin.Z * 1000, sz, indx);
zlim([-6 * P2nm, 6 * P2nm]);
xlabel('x (pixels)');
ylabel('y (pixels)');
zlabel('z (pixels)');
title('true image');
view([-66, 12])
hold off
%saveas(gcf, '3Dsim_TrueImage', 'png');

X_True = single(SMDin.X);
Y_True = single(SMDin.Y);
Z_True = single(SMDin.Z);
if strcmp(yn, 'No')
   %creating drift per frame: fast way!
   frame_num = SMDin.NFrames*(SMDin.DatasetNum - 1) + SMDin.FrameNum - 1;
   SMDin.X = SMDin.X + frame_num*PpFX;
   SMDin.Y = SMDin.Y + frame_num*PpFY;
   SMDin.Z = SMDin.Z + frame_num*PpFZ*SMDin.PixelSizeZUnit;
   %adjusts the correct xsz and ysz in driftcorrect2D
   for i = 1 : numel(SMDin.X)
      if SMDin.X(i) >= XYSize;
         SMDin.X(i) = XYSize;
      end
      if SMDin.Y(i) >= XYSize;
         SMDin.Y(i) = XYSize;
      end
   end
end
SMDsave = SMDin;

figure();
hold on
cm = colormap(jet);
n_cm = size(cm, 1);
sz = 5;
min_Z = min(SMDsave.Z);
max_Z = max(SMDsave.Z);
indx = max(1, ceil((n_cm - 1) * (SMDsave.Z - min_Z) / (max_Z - min_Z) + 1));
scatter3(SMDsave.X * P2nm, SMDsave.Y * P2nm, SMDsave.Z * 1000, sz, indx);
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('z (nm)');
title('drift image');
view([-66, 12])
hold off
%saveas(gcf, '3Dsim_DriftImage', 'png');

SMF = smi_core.SingleMoleculeFitting.createSMF();
SMF.DriftCorrection.Init_inter = SMDin.NFrames;
DC = smi_core.DriftCorrection(SMF);
%DriftParams.PixelSizeZUnit = SMDin.PixelSizeZUnit;
%DriftParams.PDegree       = 1;
%DriftParams.TolFun_intra  = 1e-2;
%DriftParams.TolX_intra    = 1e-4;
%DriftParams.TolFun_inter  = 1e-2;
%DriftParams.TolX_inter    = 1e-4;
%DC.DriftParams.Init_inter = SMDin.NFrames;
%DriftParams.Init_inter    = 0;

clear SMDin SMD
SMDin = SMDsave;
[SMD, Statistics] = DC.driftCorrectKNN(SMDin);

figure();
hold on
cm = colormap(jet);
n_cm = size(cm, 1);
sz = 5;
min_Z = min(SMD.Z);
max_Z = max(SMD.Z);
indx = max(1, ceil((n_cm - 1) * (SMD.Z - min_Z) / (max_Z - min_Z) + 1));
scatter3(SMD.X * P2nm, SMD.Y * P2nm, SMD.Z * 1000, sz, indx);
xlabel('x (nm)');
ylabel('y (nm)');
zlabel('z (nm)');
title('corrected drift image');
view([-66, 12])
hold off
%saveas(gcf, '3Dsim_correctedDriftImage', 'png');

% Remove any NaNs.
nans = find(isnan(SMDsave.X) | isnan(SMDsave.Y) | isnan(SMD.X) | isnan(SMD.Y));
n_nans = numel(nans);
if n_nans > 0
   fprintf('%d NaNs removed!\n', n_nans);
   SMDsave.X(nans) = [];
   SMDsave.Y(nans) = [];
   SMDsave.Z(nans) = [];
   SMD.X(nans) = [];
   SMD.Y(nans) = [];
   SMD.Z(nans) = [];
   X_True(nans) = [];
   Y_True(nans) = [];
   Z_True(nans) = [];
end

% Consistency check.
N = numel(SMD.X);
X_inDC = zeros(N, 1, 'single');
Y_inDC = zeros(N, 1, 'single');
Z_inDC = zeros(N, 1, 'single');
X_unDC = zeros(N, 1, 'single');
Y_unDC = zeros(N, 1, 'single');
Z_unDC = zeros(N, 1, 'single');
for k = 1:N
   i = SMD.FrameNum(k);
   j = SMD.DatasetNum(k);
   X_inDC(k) = SMDsave.X(k) - SMD.DriftX(i, j);
   Y_inDC(k) = SMDsave.Y(k) - SMD.DriftY(i, j);
   Z_inDC(k) = SMDsave.Z(k) - SMD.DriftZ(i, j);
   X_unDC(k) = SMD.X(k) + SMD.DriftX(i, j);
   Y_unDC(k) = SMD.Y(k) + SMD.DriftY(i, j);
   Z_unDC(k) = SMD.Z(k) + SMD.DriftZ(i, j);
end
consistency_un = sum(abs(SMDin.X - X_unDC) + abs(SMDin.Y - Y_unDC) + ...
                     abs(SMDin.Z - Z_unDC) / SMDin.PixelSizeZUnit);
consistency_in = sum(abs(SMD.X - X_inDC) + abs(SMD.Y - Y_inDC) + ...
                     abs(SMD.Z - Z_inDC) / SMD.PixelSizeZUnit);
fprintf('SMDin.X/Y/Z - (SMD.X/Y/Z + SMD.DriftX/Y/Z = %f nm\n', ...
        consistency_un * P2nm);
fprintf('SMD.X/Y/Z - (SMDin.X/Y/Z - SMD.DriftX/Y/Z = %f nm\n', ...
        consistency_in * P2nm);

% Plot the drift correction as a function of time.
DC_fig = DC.plotDriftCorrection(SMD);
showm(DC_fig);
hold on
view([-66, 12])
hold off
%saveas(gcf, '3Dsim_DC', 'png');

% Plot the residual between the pre-drift data and the drift corrected
% post-drift data.
[residual, dist, rmse, nnfig] = DC.calcDCResidual(SMD, X_True, Y_True);
fprintf('average distance between true and DC locations = %f\n', dist);
fprintf('RMSE             between true and DC locations = %f\n', rmse);

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
z_drift = mean(SMD.DriftZ(:) ./ (framenums(:) + 1) ./ SMD.PixelSizeZUnit);
%fprintf('average x-drift per frame = %f px (true = %f px)\n', x_drift, PpFX);
%fprintf('average y-drift per frame = %f px (true = %f px)\n', y_drift, PpFY);
%fprintf('average z-drift per frame = %f px (true = %f px)\n', z_drift, PpFZ);
fprintf('average x-drift per frame = %f nm (true = %f nm)\n', ...
        x_drift * P2nm, PpFX * P2nm);
fprintf('average y-drift per frame = %f nm (true = %f nm)\n', ...
        y_drift * P2nm, PpFY * P2nm);
fprintf('average z-drift per frame = %f nm (true = %f nm)\n', ...
        z_drift * P2nm, PpFZ * P2nm);

SMD3 = SMD;
Statistics3.Sim_rho = rho;
Statistics3.N_particles = n_particles;
Statistics3.N_particles_per_pixel   = n_particles / XYSize^2;
Statistics3.N_particles_per_dataset = n_particles / (n_frames / FpD);
Statistics3.N_NaNs = n_nans;
Statistics3.Consistency_un = consistency_un;
Statistics3.Consistency_in = consistency_in;
Statistics3.Dist = dist;
Statistics3.RMSE = rmse;
Statistics3.n_local  = n_local;
Statistics3.n_blinks = n_blinks;
Statistics3.DriftX_mean = x_drift;
Statistics3.DriftY_mean = y_drift;
Statistics3.DriftZ_mean = z_drift;
Statistics3.DriftX_True = PpFX;
Statistics3.DriftY_True = PpFY;
Statistics3.DriftZ_True = PpFZ;

end

% =============================================================================

function [n_local, n_blinks] = blinks(n_particles, X, Y, F)
% Find all the frame number sets (fns) for each localization (of which there
% are n_local).

   fns = {};
   i = 1;
   n_local = 0;
   while i <= n_particles
      x = X(i);
      y = Y(i);
      f = F(i);

      i = i + 1;
      while i <= n_particles & X(i) == x & Y(i) == y
         f = [f, F(i)];
         i = i + 1;
      end
      n_local = n_local + 1;
      fns{n_local} = f;
   end

   % Count the number of blinks.
   n_blinks = sum(cellfun(@(e) numel(find(diff(e) ~= 1)) + 1, fns));

end

% -----------------------------------------------------------------------------

function [n_local, n_blinks] = blinks3(n_particles, X, Y, Z, F)
% Find all the frame number sets (fns) for each localization (of which there
% are n_local).

   fns = {};
   i = 1;
   n_local = 0;
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
      n_local = n_local + 1;
      fns{n_local} = f;
   end

   % Count the number of blinks.
   n_blinks = sum(cellfun(@(e) numel(find(diff(e) ~= 1)) + 1, fns));

end
