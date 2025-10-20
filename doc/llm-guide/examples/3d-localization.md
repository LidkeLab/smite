---
title: "3D Localization with Astigmatism"
category: "examples"
level: "advanced"
tags: ["example", "3d", "astigmatism", "z-coordinate", "localization", "psf"]
prerequisites: ["../getting-started/installation.md", "../how-to/localize-molecules.md", "basic-localization.md"]
related: ["../how-to/localize-molecules.md", "../workflows/smlm-analysis.md", "../core-concepts/smd-structure.md"]
summary: "Complete working example of 3D single molecule localization using astigmatism-based PSF engineering with Z coordinate fitting, calibration, and 3D visualization"
estimated_time: "20 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# 3D Localization with Astigmatism

## Purpose

This example demonstrates complete 3D single molecule localization using astigmatism-based PSF engineering. You'll learn to configure the astigmatic PSF model, perform 3D localization including Z coordinates, assess Z precision, and visualize results in 3D. This technique enables axial localization precision of 20-100 nm, extending SMLM into the third dimension.

## Prerequisites

- smite installed and working with GPU support
- Completed [Basic Localization Example](basic-localization.md)
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Basic knowledge of 3D PSF engineering
- 15-20 minutes to run the code

## Overview

Astigmatism-based 3D localization introduces a cylindrical lens in the imaging path, causing the PSF to become elliptical with an orientation that depends on the Z position. As emitters move through focus:
- **Below focus**: PSF elongated in X direction
- **At focus**: PSF circular (symmetric)
- **Above focus**: PSF elongated in Y direction

By fitting the elliptical PSF shape, we can extract the Z coordinate in addition to X and Y.

This example covers:
1. Understanding the astigmatism PSF model and calibration parameters
2. Simulating 3D SMLM data with Z coordinates
3. Configuring SMF for 3D fitting (XYZNB fit type)
4. Performing 3D localization
5. Assessing Z precision and accuracy
6. Visualizing results in 3D
7. Generating side-view projections

Everything runs in memory with no external files required.

## Understanding Astigmatism Calibration

The astigmatic PSF is parameterized by six coefficients that describe how PSF width varies with Z:

**Model equations:**
```
sigma_x(z) = sqrt(Ax * (z - Gamma)^2 + Bx^2 + D)
sigma_y(z) = sqrt(Ay * (z - Gamma)^2 + By^2 - D)
```

Where:
- **Ax, Ay**: Quadratic coefficients for X and Y width vs Z
- **Bx, By**: Minimum widths at focus for X and Y
- **Gamma**: Z offset (defocus position)
- **D**: Astigmatism strength parameter

**Typical calibration values** (from experimental measurements):
- Ax ≈ 0.05 - 0.15 (pixel²/micron²)
- Ay ≈ 0.05 - 0.15 (pixel²/micron²)
- Bx ≈ 1.2 - 1.5 (pixels)
- By ≈ 1.2 - 1.5 (pixels)
- Gamma ≈ 0 microns (at focus)
- D ≈ 0.1 - 0.5 (pixel²)

For this example, we'll use simulated calibration values that produce realistic astigmatic behavior.

## Complete Working Code

Copy and run this complete example:

```matlab
%% 3D Localization with Astigmatism Example
% Demonstrates 3D molecule localization using astigmatic PSF

%% Step 1: Set Up Astigmatism Calibration Parameters
fprintf('=== Configuring 3D Astigmatism System ===\n');

% Astigmatism calibration parameters
% These would normally come from bead calibration measurements
ZFitStruct.Ax = 0.10;           % X width vs Z^2 (pixel^2/micron^2)
ZFitStruct.Ay = 0.10;           % Y width vs Z^2 (pixel^2/micron^2)
ZFitStruct.Bx = 1.3;            % Min X width at focus (pixels)
ZFitStruct.By = 1.3;            % Min Y width at focus (pixels)
ZFitStruct.Gamma = 0.0;         % Z offset / defocus (microns)
ZFitStruct.D = 0.20;            % Astigmatism strength (pixel^2)

% Z range and resolution
Z_min = -0.6;                   % Min Z position (microns)
Z_max = 0.6;                    % Max Z position (microns)
Z_range = Z_max - Z_min;

fprintf('Astigmatism calibration:\n');
fprintf('  Ax = %.3f, Ay = %.3f pixel²/µm²\n', ZFitStruct.Ax, ZFitStruct.Ay);
fprintf('  Bx = %.2f, By = %.2f pixels\n', ZFitStruct.Bx, ZFitStruct.By);
fprintf('  D = %.3f pixel² (astigmatism strength)\n', ZFitStruct.D);
fprintf('  Z range: %.1f to %.1f µm (%.1f µm total)\n', Z_min, Z_max, Z_range);

%% Step 2: Visualize PSF Shape vs Z Position
fprintf('\n=== PSF Shape vs Z Position ===\n');

% Calculate PSF widths across Z range
Z_curve = linspace(Z_min, Z_max, 100);
sigma_x = sqrt(ZFitStruct.Ax * (Z_curve - ZFitStruct.Gamma).^2 + ...
    ZFitStruct.Bx^2 + ZFitStruct.D);
sigma_y = sqrt(ZFitStruct.Ay * (Z_curve - ZFitStruct.Gamma).^2 + ...
    ZFitStruct.By^2 - ZFitStruct.D);

% Plot PSF calibration curve
figure('Name', '3D Astigmatism Calibration', 'Position', [100, 100, 1200, 400]);

subplot(1,3,1);
plot(Z_curve, sigma_x, 'b-', 'LineWidth', 2); hold on;
plot(Z_curve, sigma_y, 'r-', 'LineWidth', 2);
xlabel('Z Position (µm)'); ylabel('PSF Sigma (pixels)');
title('PSF Width vs Z Position');
legend('σ_x', 'σ_y', 'Location', 'best');
grid on;
xlim([Z_min, Z_max]);
fprintf('PSF width range: σ_x = %.2f-%.2f px, σ_y = %.2f-%.2f px\n', ...
    min(sigma_x), max(sigma_x), min(sigma_y), max(sigma_y));

% Plot ellipticity
ellipticity = abs(sigma_x - sigma_y);
subplot(1,3,2);
plot(Z_curve, ellipticity, 'k-', 'LineWidth', 2);
xlabel('Z Position (µm)'); ylabel('|σ_x - σ_y| (pixels)');
title('PSF Ellipticity vs Z');
grid on;
xlim([Z_min, Z_max]);
fprintf('Max ellipticity: %.2f pixels at Z edges\n', max(ellipticity));

% Plot sigma ratio
sigma_ratio = sigma_x ./ sigma_y;
subplot(1,3,3);
plot(Z_curve, sigma_ratio, 'm-', 'LineWidth', 2);
xlabel('Z Position (µm)'); ylabel('σ_x / σ_y');
title('Aspect Ratio vs Z');
grid on;
xlim([Z_min, Z_max]);
yline(1, 'k--', 'Symmetric');
fprintf('Aspect ratio range: %.2f - %.2f\n', min(sigma_ratio), max(sigma_ratio));

%% Step 3: Generate 3D Simulated Data
fprintf('\n=== Generating 3D Simulated Data ===\n');

% Simulation parameters
SZ = 128;                       % Image size (pixels, square)
NFrames = 100;                  % Number of frames
N_emitters_per_frame = 150;     % Emitters per frame
Photons = 1500;                 % Photons per emitter (need more for Z fitting)
Bg = 5;                         % Background (photons/pixel)
BoxSize = 9;                    % Larger boxes for asymmetric PSF

fprintf('Simulation parameters:\n');
fprintf('  Image size: %d × %d pixels\n', SZ, SZ);
fprintf('  Frames: %d\n', NFrames);
fprintf('  Emitters per frame: %d\n', N_emitters_per_frame);
fprintf('  Photons: %d\n', Photons);
fprintf('  Background: %.1f photons/pixel\n', Bg);

% Generate 3D coordinates with uniform Z distribution
N_total = N_emitters_per_frame * NFrames;
X_true = SZ * rand(N_total, 1);
Y_true = SZ * rand(N_total, 1);
Z_true = Z_min + Z_range * rand(N_total, 1);  % Uniform in Z
FrameNum = repelem((1:NFrames)', N_emitters_per_frame);

% Create ground truth SMD
SMD_true = smi_core.SingleMoleculeData.createSMD();
SMD_true.X = X_true;
SMD_true.Y = Y_true;
SMD_true.Z = Z_true;
SMD_true.FrameNum = FrameNum;
SMD_true.Photons = Photons * ones(N_total, 1);
SMD_true.NFrames = NFrames;
SMD_true.NDatasets = 1;

fprintf('Generated %d total emitters\n', N_total);
fprintf('Z distribution: %.2f ± %.2f µm (mean ± std)\n', ...
    mean(Z_true), std(Z_true));

% Generate image stack with astigmatic PSF
% We'll use the smi_sim.GaussBlobs to generate individual blobs
imageStack = zeros(SZ, SZ, NFrames, 'single');

fprintf('Generating astigmatic PSF images...\n');
for nn = 1:NFrames
    if mod(nn, 20) == 0
        fprintf('  Frame %d/%d\n', nn, NFrames);
    end

    % Get emitters in this frame
    idx = (FrameNum == nn);
    x_frame = X_true(idx);
    y_frame = Y_true(idx);
    z_frame = Z_true(idx);

    % Calculate PSF widths for each emitter based on Z
    sigma_x_frame = sqrt(ZFitStruct.Ax * (z_frame - ZFitStruct.Gamma).^2 + ...
        ZFitStruct.Bx^2 + ZFitStruct.D);
    sigma_y_frame = sqrt(ZFitStruct.Ay * (z_frame - ZFitStruct.Gamma).^2 + ...
        ZFitStruct.By^2 - ZFitStruct.D);

    % Generate frame with elliptical Gaussians
    frame = zeros(SZ, SZ, 'single');
    for jj = 1:sum(idx)
        % Create elliptical Gaussian blob
        [X_grid, Y_grid] = meshgrid(1:SZ, 1:SZ);
        dx = X_grid - x_frame(jj);
        dy = Y_grid - y_frame(jj);

        % Elliptical Gaussian
        gauss = Photons * exp(-0.5 * (dx.^2 / sigma_x_frame(jj)^2 + ...
            dy.^2 / sigma_y_frame(jj)^2));
        frame = frame + gauss;
    end

    % Add background and Poisson noise
    frame = frame + Bg;
    frame = poissrnd(frame);
    imageStack(:,:,nn) = single(frame);
end

fprintf('Image stack generation complete\n');

%% Step 4: Configure SMF for 3D Localization
fprintf('\n=== Configuring SMF for 3D Fitting ===\n');

% Create SMF structure
SMF = smi_core.SingleMoleculeFitting();

% Camera parameters
SMF.Data.CameraGain = 1;
SMF.Data.CameraOffset = 0;
SMF.Data.PixelSize = 0.108;     % 108 nm pixels (typical)
SMF.Data.FrameRate = 100;       % 100 Hz

% Box finding parameters
SMF.BoxFinding.BoxSize = BoxSize;
SMF.BoxFinding.MinPhotons = 400;         % Higher threshold for 3D
SMF.BoxFinding.BoxOverlap = 2;

% Fitting parameters - CRITICAL: Use XYZNB fit type
SMF.Fitting.FitType = 'XYZNB';           % 3D fit: X, Y, Z, photons, background
SMF.Fitting.PSFSigma = [1.3, 1.3];       % Initial guess [sigma_x, sigma_y]
SMF.Fitting.Iterations = 25;             % More iterations for 3D
SMF.Fitting.ZFitStruct = ZFitStruct;     % Astigmatism calibration

% Thresholding parameters
SMF.Thresholding.On = true;
SMF.Thresholding.MaxXY_SE = 0.15;        % Max XY precision (pixels)
SMF.Thresholding.MaxZ_SE = 0.05;         % Max Z precision (microns) - CRITICAL
SMF.Thresholding.MinPhotons = 300;
SMF.Thresholding.MinPValue = 0.01;
SMF.Thresholding.MinPSFSigma = 0.8;
SMF.Thresholding.MaxPSFSigma = 2.5;

fprintf('SMF Configuration:\n');
fprintf('  Fit type: %s (3D localization)\n', SMF.Fitting.FitType);
fprintf('  Box size: %d × %d pixels\n', BoxSize, BoxSize);
fprintf('  Initial PSF: [%.1f, %.1f] pixels\n', SMF.Fitting.PSFSigma);
fprintf('  Max Z precision threshold: %.3f µm\n', SMF.Thresholding.MaxZ_SE);
fprintf('  Iterations: %d\n', SMF.Fitting.Iterations);

%% Step 5: Perform 3D Localization
fprintf('\n=== Performing 3D Localization ===\n');

% Create LocalizeData object
LD = smi_core.LocalizeData(imageStack, SMF);
LD.Verbose = 1;

% Run 3D localization
tic;
SMD = LD.genLocalizations();
elapsed_time = toc;

fprintf('3D localization complete in %.2f seconds\n', elapsed_time);
fprintf('Found %d localizations\n', length(SMD.X));

% Check thresholding
if isfield(SMD, 'ThreshFlag')
    passed = sum(SMD.ThreshFlag == 0);
    fprintf('Passed quality filters: %d (%.1f%%)\n', ...
        passed, 100*passed/length(SMD.X));
end

% Verify Z coordinates were fitted
if isfield(SMD, 'Z') && ~isempty(SMD.Z)
    fprintf('Z coordinates successfully fitted!\n');
    fprintf('  Z range: %.3f to %.3f µm\n', min(SMD.Z), max(SMD.Z));
    fprintf('  Z mean: %.3f µm\n', mean(SMD.Z));
else
    error('Z coordinates not found - check FitType and ZFitStruct');
end

%% Step 6: Assess 3D Localization Quality
fprintf('\n=== 3D Localization Quality Assessment ===\n');

% XY precision
median_precision_xy = median(SMD.X_SE);
mean_precision_xy_nm = mean(SMD.X_SE) * SMF.Data.PixelSize * 1000;

fprintf('XY Localization:\n');
fprintf('  Median X precision: %.3f pixels (%.1f nm)\n', ...
    median(SMD.X_SE), median(SMD.X_SE) * SMF.Data.PixelSize * 1000);
fprintf('  Median Y precision: %.3f pixels (%.1f nm)\n', ...
    median(SMD.Y_SE), median(SMD.Y_SE) * SMF.Data.PixelSize * 1000);

% Z precision (in microns)
median_precision_z_um = median(SMD.Z_SE);
median_precision_z_nm = median_precision_z_um * 1000;

fprintf('\nZ Localization:\n');
fprintf('  Median Z precision: %.4f µm (%.1f nm)\n', ...
    median_precision_z_um, median_precision_z_nm);
fprintf('  Mean Z precision: %.4f µm (%.1f nm)\n', ...
    mean(SMD.Z_SE), mean(SMD.Z_SE) * 1000);
fprintf('  Best Z precision: %.4f µm (%.1f nm)\n', ...
    min(SMD.Z_SE), min(SMD.Z_SE) * 1000);
fprintf('  Worst Z precision: %.4f µm (%.1f nm)\n', ...
    max(SMD.Z_SE), max(SMD.Z_SE) * 1000);

% Precision ratio (Z vs XY)
precision_ratio = median_precision_z_nm / mean_precision_xy_nm;
fprintf('\nZ/XY precision ratio: %.1f× (Z is %.1f× worse than XY)\n', ...
    precision_ratio, precision_ratio);

% Photon statistics
fprintf('\nPhoton Statistics:\n');
fprintf('  Photons: %.0f ± %.0f (mean ± std)\n', ...
    mean(SMD.Photons), std(SMD.Photons));
fprintf('  Background: %.1f ± %.1f photons/pixel\n', ...
    mean(SMD.Bg), std(SMD.Bg));

%% Step 7: Compare to Ground Truth
fprintf('\n=== Comparison to Ground Truth ===\n');

% Match localizations to true positions using 3D distance
max_distance = 2;  % pixels in XY
max_z_distance = 0.1;  % microns in Z
matches = 0;
distances_xy = [];
distances_z = [];
distances_3d = [];

for i = 1:length(SMD_true.X)
    % Find nearest localization in XY
    dx = SMD.X - SMD_true.X(i);
    dy = SMD.Y - SMD_true.Y(i);
    dist_xy = sqrt(dx.^2 + dy.^2);

    % Among XY matches, check Z distance
    xy_candidates = find(dist_xy < max_distance);
    if ~isempty(xy_candidates)
        dz = abs(SMD.Z(xy_candidates) - SMD_true.Z(i));
        [min_dz, idx_min] = min(dz);

        if min_dz < max_z_distance
            matches = matches + 1;
            global_idx = xy_candidates(idx_min);
            distances_xy = [distances_xy; dist_xy(global_idx)];
            distances_z = [distances_z; dz(idx_min)];
            distances_3d = [distances_3d; sqrt(dist_xy(global_idx)^2 + ...
                (dz(idx_min)/SMF.Data.PixelSize)^2)];  % 3D dist in pixels
        end
    end
end

detection_rate = 100 * matches / length(SMD_true.X);
fprintf('Detection rate: %.1f%% (%d / %d true emitters)\n', ...
    detection_rate, matches, length(SMD_true.X));

if ~isempty(distances_xy)
    fprintf('\nLocalization Errors:\n');
    fprintf('  XY error: %.3f pixels (%.1f nm) median\n', ...
        median(distances_xy), median(distances_xy) * SMF.Data.PixelSize * 1000);
    fprintf('  Z error: %.4f µm (%.1f nm) median\n', ...
        median(distances_z), median(distances_z) * 1000);
    fprintf('  3D error: %.3f pixels median\n', median(distances_3d));
end

%% Step 8: 3D Visualizations
fprintf('\n=== Creating 3D Visualizations ===\n');

% Main 3D figure
fig3d = figure('Name', '3D Localization Results', ...
    'Position', [100, 100, 1600, 1000], 'Color', 'white');

% Panel 1: XY projection (traditional 2D view)
subplot(2,4,1);
plot(SMD.X, SMD.Y, 'k.', 'MarkerSize', 2);
axis equal; axis([0 SZ 0 SZ]);
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('XY Projection (Traditional 2D View)');
grid on;

% Panel 2: XY colored by Z depth
subplot(2,4,2);
scatter(SMD.X, SMD.Y, 10, SMD.Z, 'filled');
axis equal; axis([0 SZ 0 SZ]);
xlabel('X (pixels)'); ylabel('Y (pixels)');
title('XY Projection Colored by Z Depth');
colormap(gca, jet);
c = colorbar;
c.Label.String = 'Z Position (µm)';
caxis([Z_min, Z_max]);

% Panel 3: XZ side view
subplot(2,4,3);
scatter(SMD.X, SMD.Z, 10, SMD.Z, 'filled');
xlabel('X (pixels)'); ylabel('Z (µm)');
title('XZ Side View');
xlim([0 SZ]); ylim([Z_min, Z_max]);
colormap(gca, jet);
grid on;

% Panel 4: YZ side view
subplot(2,4,4);
scatter(SMD.Y, SMD.Z, 10, SMD.Z, 'filled');
xlabel('Y (pixels)'); ylabel('Z (µm)');
title('YZ Side View');
xlim([0 SZ]); ylim([Z_min, Z_max]);
colormap(gca, jet);
grid on;

% Panel 5: Full 3D scatter plot
subplot(2,4,[5,6]);
scatter3(SMD.X, SMD.Y, SMD.Z, 10, SMD.Z, 'filled');
xlabel('X (pixels)'); ylabel('Y (pixels)'); zlabel('Z (µm)');
title('3D Localization Scatter Plot');
colormap(gca, jet);
colorbar;
view(45, 30);  % 3D viewing angle
grid on;
axis vis3d;

% Panel 6: Z distribution histogram
subplot(2,4,7);
histogram(SMD.Z, 30, 'FaceColor', [0.3, 0.6, 0.9], 'EdgeColor', 'none');
xlabel('Z Position (µm)'); ylabel('Count');
title(sprintf('Z Distribution (range: %.2f µm)', range(SMD.Z)));
xline(mean(SMD.Z), 'r--', 'LineWidth', 2, 'Label', 'Mean');
grid on;

% Panel 7: Z precision histogram
subplot(2,4,8);
z_precision_nm = SMD.Z_SE * 1000;
histogram(z_precision_nm, 30, 'FaceColor', [0.9, 0.4, 0.3], 'EdgeColor', 'none');
xlabel('Z Precision (nm)'); ylabel('Count');
title(sprintf('Z Precision (median: %.1f nm)', median_precision_z_nm));
xline(median_precision_z_nm, 'r--', 'LineWidth', 2, 'Label', 'Median');
grid on;

%% Step 9: Z Precision vs Z Position Analysis
fprintf('\n=== Z Precision vs Z Position ===\n');

figure('Name', 'Z Precision Analysis', 'Position', [150, 150, 1400, 500]);

% Z precision vs Z position
subplot(1,3,1);
scatter(SMD.Z, SMD.Z_SE * 1000, 20, SMD.Photons, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Z Position (µm)'); ylabel('Z Precision (nm)');
title('Z Precision vs Z Position');
colormap(gca, parula);
c = colorbar;
c.Label.String = 'Photons';
grid on;

% Theoretical prediction (CRLB increases near edges)
% Precision roughly proportional to 1/sqrt(dSigma/dZ)
% Calculate derivative of PSF width with Z
Z_theory = linspace(Z_min, Z_max, 100);
dSigmaX_dZ = ZFitStruct.Ax * (Z_theory - ZFitStruct.Gamma) ./ ...
    sqrt(ZFitStruct.Ax * (Z_theory - ZFitStruct.Gamma).^2 + ...
    ZFitStruct.Bx^2 + ZFitStruct.D);
dSigmaY_dZ = ZFitStruct.Ay * (Z_theory - ZFitStruct.Gamma) ./ ...
    sqrt(ZFitStruct.Ay * (Z_theory - ZFitStruct.Gamma).^2 + ...
    ZFitStruct.By^2 - ZFitStruct.D);
dSigma_total = sqrt(dSigmaX_dZ.^2 + dSigmaY_dZ.^2);

subplot(1,3,2);
plot(Z_theory, abs(dSigmaX_dZ), 'b-', 'LineWidth', 2); hold on;
plot(Z_theory, abs(dSigmaY_dZ), 'r-', 'LineWidth', 2);
plot(Z_theory, dSigma_total, 'k-', 'LineWidth', 2);
xlabel('Z Position (µm)'); ylabel('|dσ/dZ| (pixels/µm)');
title('PSF Sensitivity to Z Position');
legend('dσ_x/dZ', 'dσ_y/dZ', 'Combined', 'Location', 'best');
grid on;
fprintf('Peak PSF sensitivity: %.2f pixels/µm at Z = ±%.2f µm\n', ...
    max(dSigma_total), Z_theory(dSigma_total == max(dSigma_total)));

% XY precision vs Z precision
subplot(1,3,3);
xy_precision_nm = SMD.X_SE * SMF.Data.PixelSize * 1000;
scatter(xy_precision_nm, z_precision_nm, 20, SMD.Z, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('XY Precision (nm)'); ylabel('Z Precision (nm)');
title('XY vs Z Precision');
colormap(gca, jet);
c = colorbar;
c.Label.String = 'Z Position (µm)';
grid on;
% Add diagonal line for reference
hold on;
max_prec = max([xy_precision_nm; z_precision_nm]);
plot([0 max_prec], [0 max_prec], 'k--', 'LineWidth', 1.5);
legend('Localizations', 'Equal precision', 'Location', 'northwest');

%% Step 10: Error Analysis vs Z Position
fprintf('\n=== Error Analysis by Z Position ===\n');

if ~isempty(distances_z)
    % Bin errors by Z position
    Z_bins = linspace(Z_min, Z_max, 10);
    Z_bin_centers = (Z_bins(1:end-1) + Z_bins(2:end)) / 2;

    % Match true positions to localizations
    matched_true_Z = zeros(length(distances_z), 1);
    match_idx = 1;
    for i = 1:length(SMD_true.X)
        dx = SMD.X - SMD_true.X(i);
        dy = SMD.Y - SMD_true.Y(i);
        dist_xy = sqrt(dx.^2 + dy.^2);
        xy_candidates = find(dist_xy < max_distance);

        if ~isempty(xy_candidates)
            dz = abs(SMD.Z(xy_candidates) - SMD_true.Z(i));
            [min_dz, ~] = min(dz);

            if min_dz < max_z_distance && match_idx <= length(distances_z)
                matched_true_Z(match_idx) = SMD_true.Z(i);
                match_idx = match_idx + 1;
            end
        end
    end

    mean_error_by_z = zeros(size(Z_bin_centers));
    std_error_by_z = zeros(size(Z_bin_centers));

    for i = 1:length(Z_bin_centers)
        in_bin = matched_true_Z >= Z_bins(i) & matched_true_Z < Z_bins(i+1);
        if sum(in_bin) > 0
            mean_error_by_z(i) = mean(distances_z(in_bin)) * 1000;  % nm
            std_error_by_z(i) = std(distances_z(in_bin)) * 1000;
        end
    end

    figure('Name', 'Z Error vs Position');
    errorbar(Z_bin_centers, mean_error_by_z, std_error_by_z, ...
        'o-', 'LineWidth', 2, 'MarkerSize', 8, 'MarkerFaceColor', 'auto');
    xlabel('True Z Position (µm)');
    ylabel('Z Localization Error (nm)');
    title('Z Error vs Z Position');
    grid on;

    fprintf('Z error by position:\n');
    for i = 1:length(Z_bin_centers)
        if mean_error_by_z(i) > 0
            fprintf('  Z = %.2f µm: %.1f ± %.1f nm\n', ...
                Z_bin_centers(i), mean_error_by_z(i), std_error_by_z(i));
        end
    end
end

%% Step 11: Generate 3D Super-Resolution Image
fprintf('\n=== Generating 3D Super-Resolution Views ===\n');

% Create depth-coded image (color represents Z)
SR_zoom = 20;
SMD.XSize = SZ;
SMD.YSize = SZ;

% Generate histogram for each Z slice
N_z_slices = 5;
Z_slice_edges = linspace(Z_min, Z_max, N_z_slices + 1);

figure('Name', '3D Depth Slices', 'Position', [200, 200, 1400, 800]);

for i = 1:N_z_slices
    % Select localizations in this Z slice
    z_mask = SMD.Z >= Z_slice_edges(i) & SMD.Z < Z_slice_edges(i+1);
    SMD_slice = SMD;
    SMD_slice.X = SMD.X(z_mask);
    SMD_slice.Y = SMD.Y(z_mask);
    SMD_slice.X_SE = SMD.X_SE(z_mask);
    SMD_slice.Y_SE = SMD.Y_SE(z_mask);
    SMD_slice.Photons = SMD.Photons(z_mask);
    SMD_slice.Bg = SMD.Bg(z_mask);

    % Generate Gaussian image
    if ~isempty(SMD_slice.X)
        SR_image = smi_vis.GenerateImages.gaussianImage(SMD_slice, SR_zoom, 0);

        subplot(2, N_z_slices, i);
        imshow(SR_image);
        title(sprintf('Z = %.2f to %.2f µm', Z_slice_edges(i), Z_slice_edges(i+1)));

        % Also show count
        subplot(2, N_z_slices, i + N_z_slices);
        text(0.5, 0.5, sprintf('%d locs', sum(z_mask)), ...
            'HorizontalAlignment', 'center', 'FontSize', 14);
        axis off;
    end
end

%% Step 12: Summary Statistics
fprintf('\n=== 3D Localization Summary ===\n');
fprintf('==========================================\n');
fprintf('Data:\n');
fprintf('  Image size: %d × %d pixels\n', SZ, SZ);
fprintf('  Z range: %.2f to %.2f µm\n', Z_min, Z_max);
fprintf('  Frames: %d\n', NFrames);
fprintf('  True emitters: %d\n', length(SMD_true.X));
fprintf('  Localizations found: %d\n', length(SMD.X));
fprintf('  Detection rate: %.1f%%\n', detection_rate);

fprintf('\nAstigmatism Configuration:\n');
fprintf('  Calibration: Ax=%.2f, Ay=%.2f, Bx=%.2f, By=%.2f, D=%.2f\n', ...
    ZFitStruct.Ax, ZFitStruct.Ay, ZFitStruct.Bx, ZFitStruct.By, ZFitStruct.D);

fprintf('\nLocalization Precision:\n');
fprintf('  XY precision: %.1f nm (median)\n', ...
    median(SMD.X_SE) * SMF.Data.PixelSize * 1000);
fprintf('  Z precision: %.1f nm (median)\n', median_precision_z_nm);
fprintf('  Z/XY ratio: %.1f×\n', precision_ratio);

if ~isempty(distances_z)
    fprintf('\nAccuracy vs Ground Truth:\n');
    fprintf('  XY error: %.1f nm (median)\n', ...
        median(distances_xy) * SMF.Data.PixelSize * 1000);
    fprintf('  Z error: %.1f nm (median)\n', median(distances_z) * 1000);
end

fprintf('\nPerformance:\n');
fprintf('  Localization time: %.2f seconds\n', elapsed_time);
fprintf('  Speed: %.0f localizations/second\n', length(SMD.X) / elapsed_time);
fprintf('==========================================\n');

fprintf('\n3D localization example complete!\n');
```

## What This Example Demonstrates

### Astigmatism PSF Model
- Configures the 6-parameter astigmatic PSF calibration
- Visualizes PSF shape changes across Z positions
- Shows ellipticity and aspect ratio curves
- Demonstrates Z-dependent PSF widths

### 3D Data Generation
- Creates realistic 3D emitter distributions
- Simulates elliptical PSFs with Z-dependent shapes
- Adds Poisson noise to images
- Maintains ground truth for validation

### 3D Localization Pipeline
- Configures XYZNB fit type for 3D fitting
- Sets appropriate box sizes for asymmetric PSFs
- Applies Z-specific quality thresholds
- Extracts Z coordinates with uncertainties

### Quality Assessment
- Compares XY vs Z precision
- Analyzes Z precision across depth range
- Validates against ground truth
- Examines Z-dependent localization errors

### 3D Visualization
- Multiple projection views (XY, XZ, YZ)
- Depth-coded color maps
- 3D scatter plots
- Z-sliced super-resolution images

## Expected Results

When you run this example, you should see:

**Z Precision**: ~30-60 nm median (3-5× worse than XY)
- Better near focus (Z ≈ 0)
- Worse at Z extremes (±0.5 µm)
- Depends on photon count

**Detection Rate**: ~85-95% (slightly lower than 2D due to stricter Z thresholds)

**Z Localization Error**: ~40-80 nm median vs ground truth

**Z/XY Precision Ratio**: ~3-5× (Z is typically 3-5× worse than lateral precision)

## Understanding Z Precision

Z precision depends on:

1. **PSF Asymmetry Gradient**: Precision ∝ 1/(dσ/dZ)
   - Best where PSF shape changes most rapidly with Z
   - Worst near focus (Z ≈ 0) where PSF is symmetric
   - Improves at Z extremes due to higher ellipticity

2. **Photon Count**: Precision ∝ 1/√N
   - More photons = better Z precision
   - Typically need >1000 photons for good 3D localization

3. **Astigmatism Strength**: Parameter D controls Z sensitivity
   - Larger D = more astigmatism = better Z precision
   - But also increases PSF size = worse XY precision
   - Tradeoff between Z and XY performance

4. **Z Position**: Non-uniform precision across depth
   - Best: ±0.3-0.5 µm from focus
   - Worst: At focus (Z ≈ 0)
   - Degrades beyond ±0.6-0.8 µm

## Modifications to Try

### 1. Vary Astigmatism Strength

```matlab
% Stronger astigmatism (better Z, worse XY)
ZFitStruct.D = 0.40;

% Weaker astigmatism (worse Z, better XY)
ZFitStruct.D = 0.10;
```

### 2. Change Z Range

```matlab
% Larger Z range (lower precision at edges)
Z_min = -1.0;
Z_max = 1.0;

% Smaller Z range (better precision)
Z_min = -0.4;
Z_max = 0.4;
```

### 3. Adjust Photon Count

```matlab
% More photons (better Z precision)
Photons = 3000;

% Fewer photons (worse Z precision)
Photons = 800;
% Note: May need to relax thresholds
```

### 4. Modify Z Thresholds

```matlab
% Stricter Z quality (fewer but better localizations)
SMF.Thresholding.MaxZ_SE = 0.03;  % 30 nm max

% More lenient (more localizations, lower quality)
SMF.Thresholding.MaxZ_SE = 0.08;  % 80 nm max
```

### 5. Different Calibration Parameters

```matlab
% Experimental values from real astigmatic system
ZFitStruct.Ax = 0.12;
ZFitStruct.Ay = 0.11;
ZFitStruct.Bx = 1.4;
ZFitStruct.By = 1.35;
ZFitStruct.Gamma = 0.05;  % Slight defocus offset
ZFitStruct.D = 0.25;
```

### 6. Analyze Z-Dependent Structures

```matlab
% Create structure at specific Z
Z_structure = -0.2;  % 200 nm below focus
structure_radius = 20;  % pixels
structure_center = [SZ/2, SZ/2];

% Generate emitters on circular structure
N_structure = 500;
theta = linspace(0, 2*pi, N_structure);
X_structure = structure_center(1) + structure_radius * cos(theta);
Y_structure = structure_center(2) + structure_radius * sin(theta);
Z_structure = Z_structure * ones(N_structure, 1);

% Add to simulation
```

## Real-World Considerations

### Obtaining Calibration Parameters

In practice, ZFitStruct parameters come from bead calibration:

1. Image fluorescent beads at known Z positions (-800 to +800 nm typical)
2. Fit PSF widths (σ_x, σ_y) at each Z
3. Fit calibration curves to model equations
4. Save ZFitStruct for use in analysis

Smite provides tools in `smi_psf.PointSpreadFunction` for calibration.

### Practical Limitations

**Z Range**: Typically limited to ±600-800 nm
- Beyond this, PSF becomes too large/dim
- Z precision degrades significantly

**Photon Requirements**: Need >1000 photons for reliable 3D
- Dimmer emitters may fail Z threshold
- May need to trade detection rate for precision

**Computational Cost**: 3D fitting is slower than 2D
- More fit parameters (5 vs 4)
- More iterations needed
- Larger boxes (asymmetric PSF)

**Field-Dependent Calibration**: PSF may vary across field of view
- May need multiple calibrations
- Or field-dependent correction

## Troubleshooting

**Issue: No Z coordinates in SMD**

Solutions:
- Verify `SMF.Fitting.FitType = 'XYZNB'`
- Check `SMF.Fitting.ZFitStruct` is properly configured
- Ensure all six parameters (Ax, Ay, Bx, By, Gamma, D) are set

**Issue: Poor Z precision (>100 nm)**

Solutions:
- Increase photon count in simulation
- Check astigmatism calibration parameters
- Verify Z thresholds aren't too lenient
- Increase number of fitting iterations

**Issue: Low detection rate (<70%)**

Solutions:
- Relax `MaxZ_SE` threshold
- Lower `MinPhotons` threshold
- Check box size is large enough for elliptical PSF
- Verify Z range matches calibration range

**Issue: Z values out of expected range**

Solutions:
- Check ZFitStruct.Gamma (should be near 0 for centered data)
- Verify calibration parameters are in correct units
- Ensure Z simulation range matches calibration range

**Issue: Z precision worse near focus**

This is expected behavior:
- PSF is most symmetric at Z = 0
- Astigmatism provides least information
- Precision improves away from focus
- Not a bug - inherent to astigmatism method

## See Also

- [Basic Localization Example](basic-localization.md) - 2D localization foundation
- [How to Localize Molecules](../how-to/localize-molecules.md) - Detailed localization guide
- [SMD Structure](../core-concepts/smd-structure.md) - Understanding SMD.Z field
- [smi_psf API](../api-reference/smi-psf.md) - PSF calibration tools
- MATLAB/+smi_psf/@PointSpreadFunction/ - PSF modeling and calibration
- GaussMLE documentation - Details of XYZNB fit type

## References

1. Huang, B., Wang, W., Bates, M., & Zhuang, X. (2008). Three-dimensional super-resolution imaging by stochastic optical reconstruction microscopy. *Science*, 319(5864), 810-813.

2. Smith, C. S., Joseph, N., Rieger, B., & Lidke, K. A. (2010). Fast, single-molecule localization that achieves theoretically minimum uncertainty. *Nature Methods*, 7(5), 373-375.
