---
title: "API Reference: +smi_psf Namespace"
category: "api-reference"
level: "advanced"
tags: ["api", "smi_psf", "psf", "optics", "3d-localization", "zernike", "calibration"]
prerequisites: ["../core-concepts/architecture.md", "../core-concepts/smf-structure.md", "./smi-core.md"]
related: ["../how-to/localize-molecules.md", "../workflows/smlm-analysis.md"]
summary: "Complete API reference for the +smi_psf namespace covering point spread function modeling, Zernike polynomials, phase retrieval, and 3D localization calibration"
estimated_time: "45 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# API Reference: +smi_psf Namespace

## Purpose

The +smi_psf namespace provides advanced tools for point spread function (PSF) modeling, calibration, and optimization in smite. This reference documents the classes and methods for creating theoretical PSF models, performing phase retrieval from experimental data, computing localization precision bounds, and optimizing PSF designs for 3D single molecule localization microscopy.

## Prerequisites

- Understanding of [smite architecture](../core-concepts/architecture.md)
- Familiarity with [optical theory and PSF concepts](https://en.wikipedia.org/wiki/Point_spread_function)
- Knowledge of [Zernike polynomials](https://en.wikipedia.org/wiki/Zernike_polynomials)
- Basic understanding of Fourier optics
- MATLAB GPU programming basics

## Overview

The +smi_psf namespace provides:

**PSF Modeling:**
- `PointSpreadFunction` - Main class for PSF generation and analysis
- Scalar diffraction theory-based PSF computation
- Zernike polynomial-based aberration modeling
- OTF (Optical Transfer Function) rescaling

**3D Localization:**
- Phase retrieval from experimental PSF stacks
- Pupil function estimation
- 3D PSF models for defocus-based localization
- Astigmatic PSF support via phase masks

**Precision Analysis:**
- Cramer-Rao Lower Bound (CRLB) calculations
- Fisher Information Matrix computation
- Theoretical precision limits for X, Y, Z localization

**Optimization:**
- PSF design optimization for 3D imaging
- Zernike coefficient optimization
- Phase mask design (e.g., Prasad zones)

**Mathematical Tools:**
- `Zernike` - Zernike polynomial utilities
- Coordinate transformations (Noll, Wyant indexing)
- Zernike image generation and expansion

## Class Reference

### PointSpreadFunction

**Purpose:** Creates, models, and analyzes point spread functions for single molecule localization microscopy.

**Key Concept:** The PointSpreadFunction class generates PSF models using scalar diffraction theory with the pupil function formalism. It supports arbitrary phase and magnitude pupil functions, enabling modeling of aberrations, specialized phase masks for 3D imaging, and experimental PSF calibration through phase retrieval.

**Class Definition:**
```matlab
classdef PointSpreadFunction < handle
```

**Properties:**

Core optical parameters:
```matlab
Lambda          % Emission wavelength (micrometers)
NA              % Numerical aperture
N               % Index of refraction (immersion medium)
PixelSize       % Camera pixel size (micrometers)
```

PSF data and models:
```matlab
PSF             % Generated PSF stack
Z               % Array of Z positions (micrometers)
PSFStruct       % PSF structure for modeling
PSFModel        % PSF model from PSFStruct
```

Experimental data:
```matlab
DataSet         % PSF data (SZ x SZ x NRepeats x MPlanes)
DataFile        % Path to data file
RawData         % Raw data from file
PSFData         % Cropped PSF data (not averaged)
PSFCenter       % Center of PSF (pixels) [Y, X]
SZ              % Cropped PSF size (pixels)
```

Zernike expansion limits:
```matlab
MaxZCMag        % Max Zernike coefficient for magnitude (default: 81)
MaxZCPhase      % Max Zernike coefficient for phase (default: 81)
```

**Constructor:**
```matlab
PSFObj = smi_psf.PointSpreadFunction()
```

Creates empty PointSpreadFunction object. Properties must be set manually or via loading data.

**Instance Methods:**

`loadData(FileName, DataSetName)` - Load experimental PSF data
```matlab
PSFObj = smi_psf.PointSpreadFunction();
PSFObj.loadData('/path/to/psf_data.mat', 'DataSet');
% Loads DataSet, Z, and optionally PixelSize, NA, N
```

`savePSF(FileName)` - Save PSFStruct and PSF model
```matlab
PSFObj.savePSF('/path/to/psf_model.mat');
% Saves PSFStruct and PSF for later use
```

`cropData(Center)` - Crop data to SZ around center
```matlab
% Interactive center selection
PSFObj.SZ = 64;
PSFObj.cropData();  % Opens interactive window

% Specify center manually
PSFObj.cropData([128, 128]);  % [Y, X] in pixels
```

`phaseRetrieve()` - Retrieve pupil magnitude and phase
```matlab
% After loading and cropping data
PSFObj.Lambda = 0.69;   % micrometers
PSFObj.NA = 1.35;
PSFObj.N = 1.4;
PSFObj.PixelSize = 0.108;
PSFObj.SZ = 64;
PSFObj.MaxZCMag = 22;
PSFObj.MaxZCPhase = 81;

PSFObj.phaseRetrieve();
% Populates PSFStruct and PSFModel
```

**Static Methods:**

#### PSF Structure Creation

`createPSFStruct()` - Create default PSF structure
```matlab
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();

% Default values:
% Z = -2:0.1:2 (micrometers)
% Lambda = 0.69 (micrometers)
% NA = 1.35
% N = 1.4
% PixelSize = 0.1 (micrometers)
% SZ = 64 (pixels)
% OSZ = 256 (oversampling size)
% OTFSigma = [0, 0] (no rescaling)
% ZC_Mag = [1] (uniform magnitude)
% ZC_Phase = [0] (no phase aberration)
```

`createZernikeStruct(SZ, Radius, NMax)` - Create Zernike basis structure
```matlab
% For Zernike expansion up to NMax coefficients
KPixelSize = 1 / (256 * 0.1);
PupilRadius = (1.35 / 0.69) / KPixelSize;
ZStruct = smi_psf.PointSpreadFunction.createZernikeStruct(256, PupilRadius, 81);
% Used internally for phase retrieval and modeling
```

#### PSF Generation

`scalarPSF(PSFStruct)` - Generate PSF from uniform pupil
```matlab
% Basic PSF with no aberrations
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Z = -1:0.1:1;  % Defocus range
[PSF, PSFStruct] = smi_psf.PointSpreadFunction.scalarPSF(PSFStruct);
% PSF: 64 x 64 x 21 stack
```

`scalarPSFPupil(PSFStruct)` - Generate PSF from pupil function
```matlab
% PSFStruct must have Pupil field (OSZ x OSZ x 2)
% Pupil(:,:,1) = magnitude
% Pupil(:,:,2) = phase (radians)
[PSF, PSFStruct] = smi_psf.PointSpreadFunction.scalarPSFPupil(PSFStruct);
```

`scalarPSFZernike(PSFStruct)` - Generate PSF from Zernike coefficients
```matlab
% Create PSF with astigmatism
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Z = -1:0.05:1;

% Zernike coefficients (Noll ordering)
% Index 1: Piston
% Index 2-3: Tip/tilt
% Index 4: Defocus
% Index 5-6: Astigmatism
PSFStruct.ZC_Phase = zeros(15, 1);
PSFStruct.ZC_Phase(5) = 2;  % Add astigmatism

[PSF, PSFStruct] = smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct);
% Creates astigmatic PSF for 3D localization
```

`scalarPSFPrasadZone(PSFStruct, L)` - Generate PSF with Prasad zones
```matlab
% Prasad phase mask for 3D super-resolution
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
L = 4;  % Number of zones
[PSF, PSFStruct] = smi_psf.PointSpreadFunction.scalarPSFPrasadZone(PSFStruct, L);
% Creates engineered PSF with depth-dependent rotation
```

`oversamplePSFPupil(PSFStruct, Sampling)` - Generate oversampled PSF
```matlab
% Higher spatial sampling for accuracy
Sampling = 10;  % 10x oversampling
[PSF, PSFStruct] = smi_psf.PointSpreadFunction.oversamplePSFPupil(PSFStruct, Sampling);
```

#### Phase Retrieval and Calibration

`phaseRetrieval(PSFStruct, Data, MaxZCMag, MaxZCPhase)` - Retrieve pupil from data
```matlab
% Setup PSFStruct with optical parameters
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Lambda = 0.69;
PSFStruct.NA = 1.35;
PSFStruct.N = 1.4;
PSFStruct.PixelSize = 0.108;
PSFStruct.Z = -1:0.1:1;  % Must match data

% Load experimental PSF stack
% Data: 64 x 64 x 21 (matching Z positions)

% Perform phase retrieval
MaxZCMag = 22;    % Magnitude smoothing
MaxZCPhase = 81;  % Phase smoothing
[PSFStruct, PSF] = smi_psf.PointSpreadFunction.phaseRetrieval(...
    PSFStruct, Data, MaxZCMag, MaxZCPhase);

% PSFStruct now contains:
% - Pupil (magnitude and phase)
% - ZC_Mag (Zernike coefficients for magnitude)
% - ZC_Phase (Zernike coefficients for phase)
% - OTFSigma (OTF rescaling parameters)
```

**Algorithm:** Gerchberg-Saxton iterative phase retrieval with Zernike regularization. Alternates between real and Fourier space, constraining pupil with Zernike expansion to ensure smooth solutions.

`phaseRetrievalEM(PSFStruct, Data)` - EM-based phase retrieval
```matlab
% Expectation-maximization approach
[PSFStruct] = smi_psf.PointSpreadFunction.phaseRetrievalEM(PSFStruct, Data);
```

`rescaleOTF(PSFStruct, Data)` - Compute OTF rescaling parameters
```matlab
% After initial PSF generation
PSF = smi_psf.PointSpreadFunction.scalarPSFPupil(PSFStruct);
OTFSigma = smi_psf.PointSpreadFunction.rescaleOTF(PSF, Data);
% Returns [SigmaY, SigmaX] for Gaussian OTF smoothing
```

#### CRLB Calculations

`crlbPSFPupil(PSFStruct, Photons, Bg, PlotFlag)` - Cramer-Rao Lower Bound
```matlab
% Compute theoretical precision limits
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.ZC_Phase(5) = 2;  % Astigmatic PSF
PSFStruct.Z = -1:0.1:1;

Photons = 1000;  % Integrated photons
Bg = 20;         % Background photons/pixel
PlotFlag = 1;    % Show plot

[CRLB, DET] = smi_psf.PointSpreadFunction.crlbPSFPupil(...
    PSFStruct, Photons, Bg, PlotFlag);

% CRLB: NZ x 5 matrix [Y, X, Z, Photons, Bg]
% Each row contains variance for one Z position
fprintf('X precision at focus: %.1f nm\n', sqrt(CRLB(11, 2)) * 1000);
fprintf('Z precision at focus: %.1f nm\n', sqrt(CRLB(11, 3)) * 1000);
```

**Output:**
- `CRLB`: Cramer-Rao lower bounds (variance) for each parameter
- `DET`: Determinant of inverse Fisher Information Matrix (volume in parameter space)

#### PSF Optimization

`optimPSFZernike(PSFStruct, PhaseMask, StartPhase, Photons, Bg)` - Optimize PSF design
```matlab
% Optimize Zernike coefficients for best 3D precision
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Z = -1:0.1:1;

% Mask indicating which coefficients to optimize
% First 3 are piston/tip/tilt (don't help localization)
PhaseMask = [0; 0; 0; ones(12, 1)];  % Optimize coefficients 4-15

Photons = 1000;
Bg = 20;

% Optimize (can take several minutes)
[PSFStruct, PSF, CRLB] = smi_psf.PointSpreadFunction.optimPSFZernike(...
    PSFStruct, PhaseMask, [], Photons, Bg);

% Returns optimized PSF and precision bounds
fprintf('Optimized Z range: %.1f nm\n', ...
    max(sqrt(CRLB(:,3))) * 1000);
```

**Algorithm:** Uses fminsearch to minimize determinant of inverse Fisher Information Matrix over Z range. Finds Zernike coefficients yielding best 3D localization precision.

#### Zernike Polynomial Methods

`zernikeImage(NollCoef, SZ, Radius, R, Theta, Mask)` - Generate Zernike polynomial image
```matlab
% Create single Zernike mode
SZ = 256;
Radius = 100;
[X, Y] = meshgrid(-SZ/2:SZ/2-1, -SZ/2:SZ/2-1);
R = sqrt(X.^2 + Y.^2);
Theta = atan2(Y, X);
Mask = R < Radius;

NollCoef = zeros(15, 1);
NollCoef(5) = 1;  % Astigmatism mode

Image = smi_psf.PointSpreadFunction.zernikeImage(...
    NollCoef, SZ, Radius, R, Theta, Mask);
```

`zernikeSum(NollCoefs, ZStruct)` - Sum Zernike polynomials
```matlab
% Compute weighted sum of Zernike modes
ZStruct = smi_psf.PointSpreadFunction.createZernikeStruct(256, 100, 81);
NollCoefs = randn(15, 1);  % Random aberrations
Image = smi_psf.PointSpreadFunction.zernikeSum(NollCoefs, ZStruct);
```

`zernikeExpansion(Image, ZStruct)` - Decompose into Zernike basis
```matlab
% Project arbitrary phase onto Zernike basis
ZStruct = smi_psf.PointSpreadFunction.createZernikeStruct(256, 100, 81);

% Arbitrary pupil phase
[X, Y] = meshgrid(-128:127, -128:127);
R = sqrt(X.^2 + Y.^2);
Mask = R < 100;
Phase = Mask .* exp(-R.^2 / 1000);  % Gaussian phase

[NollCoef, ImageZ] = smi_psf.PointSpreadFunction.zernikeExpansion(...
    Phase, ZStruct);
% NollCoef: Zernike coefficients
% ImageZ: Reconstructed phase from coefficients
```

#### ROI Generation

`psfROIStack(PSF, XYSamPerPix, ZSamPerUnit, SZ, SMD, NoiseIm)` - Generate PSF ROIs
```matlab
% Create simulated data boxes with noise
% First generate a PSF model
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.PixelSize = 0.108 / 10;  % 10x oversampled
PSFStruct.Z = -1:0.01:1;  % Fine Z sampling
PSF = smi_psf.PointSpreadFunction.scalarPSF(PSFStruct);

% Define localizations
SMD.X = rand(100, 1) * 5 + 5;  % Random X positions
SMD.Y = rand(100, 1) * 5 + 5;  % Random Y positions
SMD.Z = randn(100, 1) * 0.5;   % Random Z positions (microns)
SMD.Photons = 1000 * ones(100, 1);
SMD.Bg = 20;

% Generate ROI stack
XYSamPerPix = 10;  % PSF sampling per pixel
ZSamPerUnit = 100; % PSF samples per micron
SZ = 15;           % Box size (pixels)
NoiseIm = 5 * ones(SZ);  % Read noise (photons)

[Model, Data] = smi_psf.PointSpreadFunction.psfROIStack(...
    PSF, XYSamPerPix, ZSamPerUnit, SZ, SMD, NoiseIm);

% Model: Noiseless PSF boxes (15 x 15 x 100)
% Data: With Poisson and read noise (15 x 15 x 100)
```

**Usage Examples:**

#### Basic PSF Generation

Creating a simple diffraction-limited PSF:

```matlab
% Setup optical parameters
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Lambda = 0.69;    % 690 nm emission
PSFStruct.NA = 1.35;        % Oil immersion objective
PSFStruct.N = 1.4;          % Oil refractive index
PSFStruct.PixelSize = 0.108;  % Camera pixel size (microns)
PSFStruct.SZ = 64;          % PSF size
PSFStruct.Z = -2:0.1:2;     % Z range (microns)

% Generate PSF
[PSF, PSFStruct] = smi_psf.PointSpreadFunction.scalarPSF(PSFStruct);

% Visualize
figure;
for ii = 1:size(PSF, 3)
    imagesc(PSF(:,:,ii));
    axis equal tight;
    colorbar;
    title(sprintf('Z = %.2f μm', PSFStruct.Z(ii)));
    pause(0.1);
end
```

#### Astigmatic PSF for 3D Localization

Creating an astigmatic PSF using Zernike coefficients:

```matlab
% Setup for 3D imaging
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Lambda = 0.69;
PSFStruct.NA = 1.4;
PSFStruct.N = 1.518;  % High-index oil
PSFStruct.PixelSize = 0.108;
PSFStruct.SZ = 64;
PSFStruct.Z = -0.8:0.05:0.8;  % 3D range

% Add astigmatism (Zernike mode 5 in Noll ordering)
PSFStruct.ZC_Mag = [1];  % Uniform magnitude
PSFStruct.ZC_Phase = zeros(15, 1);
PSFStruct.ZC_Phase(5) = 2.0;  % Astigmatism strength

% Generate astigmatic PSF
[PSF, PSFStruct] = smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct);

% Compute precision
Photons = 1000;
Bg = 15;
[CRLB, DET] = smi_psf.PointSpreadFunction.crlbPSFPupil(...
    PSFStruct, Photons, Bg, 1);

% Report precision
ZIdx = find(PSFStruct.Z == 0);  % Focus
fprintf('Lateral precision at focus: %.1f nm\n', ...
    sqrt(mean(CRLB(ZIdx, 1:2))) * 1000);
fprintf('Axial precision at focus: %.1f nm\n', ...
    sqrt(CRLB(ZIdx, 3)) * 1000);

% Save for use in fitting
save('astigmatic_psf_model.mat', 'PSFStruct', 'PSF');
```

#### Phase Retrieval from Experimental Data

Calibrating PSF model from experimental bead measurements:

```matlab
% Step 1: Prepare experimental data
% Measure PSF from fluorescent beads at known Z positions
% Data should be: (Y, X, Z_positions)
% Example: 64 x 64 x 21 stack

% Load experimental data
load('bead_psf_stack.mat', 'Data', 'Z_positions');

% Step 2: Setup PSF structure
PSFObj = smi_psf.PointSpreadFunction();
PSFObj.Lambda = 0.69;
PSFObj.NA = 1.35;
PSFObj.N = 1.4;
PSFObj.PixelSize = 0.108;
PSFObj.Z = Z_positions;  % Must match data
PSFObj.SZ = size(Data, 1);

% Step 3: Normalize and prepare data
Data = Data - min(Data(:));  % Remove offset
for ii = 1:size(Data, 3)
    Data(:,:,ii) = Data(:,:,ii) / sum(sum(Data(:,:,ii)));  % Normalize
end

% Step 4: Perform phase retrieval
PSFObj.MaxZCMag = 22;   % Moderate magnitude smoothing
PSFObj.MaxZCPhase = 81; % Full phase smoothing
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Lambda = PSFObj.Lambda;
PSFStruct.NA = PSFObj.NA;
PSFStruct.N = PSFObj.N;
PSFStruct.PixelSize = PSFObj.PixelSize;
PSFStruct.Z = PSFObj.Z;
PSFStruct.SZ = PSFObj.SZ;

[PSFStruct, PSF] = smi_psf.PointSpreadFunction.phaseRetrieval(...
    PSFStruct, Data, PSFObj.MaxZCMag, PSFObj.MaxZCPhase);

% Step 5: Validate model
figure;
subplot(1,3,1);
imagesc(Data(:,:,11)); title('Data (Z=0)');
axis equal tight; colorbar;

subplot(1,3,2);
imagesc(PSF(:,:,11)); title('Model (Z=0)');
axis equal tight; colorbar;

subplot(1,3,3);
imagesc(Data(:,:,11) - PSF(:,:,11)); title('Residual');
axis equal tight; colorbar;

% Step 6: Examine retrieved pupil
Pupil_Mag = PSFStruct.Pupil(:,:,1);
Pupil_Phase = PSFStruct.Pupil(:,:,2);

figure;
subplot(1,2,1);
imagesc(Pupil_Mag); title('Pupil Magnitude');
axis equal tight; colorbar;

subplot(1,2,2);
imagesc(Pupil_Phase); title('Pupil Phase (rad)');
axis equal tight; colorbar;

% Step 7: Save calibrated PSF
save('calibrated_psf_model.mat', 'PSFStruct', 'PSF');

% Step 8: Report Zernike aberrations
fprintf('\nZernike Phase Coefficients (Noll ordering):\n');
for ii = 1:min(15, length(PSFStruct.ZC_Phase))
    if abs(PSFStruct.ZC_Phase(ii)) > 0.1
        fprintf('  Mode %2d: %+.3f\n', ii, PSFStruct.ZC_Phase(ii));
    end
end
```

#### Optimizing PSF for 3D Imaging

Designing optimal phase mask for 3D super-resolution:

```matlab
% Setup optimization
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Lambda = 0.69;
PSFStruct.NA = 1.4;
PSFStruct.N = 1.518;
PSFStruct.PixelSize = 0.108;
PSFStruct.SZ = 64;
PSFStruct.Z = -0.6:0.05:0.6;  % 1.2 micron range

% Define which Zernike modes to optimize
% Exclude piston, tip, tilt (modes 1-3)
% Optimize defocus, astigmatism, coma, etc. (modes 4-21)
NModes = 21;
PhaseMask = zeros(NModes, 1);
PhaseMask(4:end) = 1;  % Optimize these modes

% Imaging conditions
Photons = 1000;  % Expected photons per molecule
Bg = 20;         % Background photons per pixel

% Run optimization (may take 10-30 minutes)
fprintf('Optimizing PSF design...\n');
[PSFStructOpt, PSFOpt, CRLBOpt] = ...
    smi_psf.PointSpreadFunction.optimPSFZernike(...
        PSFStruct, PhaseMask, [], Photons, Bg);

% Compare to standard astigmatism
PSFStructAstig = PSFStruct;
PSFStructAstig.ZC_Phase = zeros(NModes, 1);
PSFStructAstig.ZC_Phase(5) = 2.0;
[PSFAstig, ~] = smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStructAstig);
[CRLBAstig, ~] = smi_psf.PointSpreadFunction.crlbPSFPupil(...
    PSFStructAstig, Photons, Bg, 0);

% Plot precision comparison
figure;
plot(PSFStruct.Z, sqrt(CRLBAstig(:,3)) * 1000, 'b-', 'LineWidth', 2);
hold on;
plot(PSFStruct.Z, sqrt(CRLBOpt(:,3)) * 1000, 'r-', 'LineWidth', 2);
xlabel('Z Position (μm)');
ylabel('Z Precision (nm)');
legend('Standard Astigmatism', 'Optimized PSF');
title('3D Localization Precision');
grid on;

% Report improvement
MeanPrecAstig = mean(sqrt(CRLBAstig(:,3)));
MeanPrecOpt = mean(sqrt(CRLBOpt(:,3)));
fprintf('\nMean Z precision:\n');
fprintf('  Standard astigmatism: %.1f nm\n', MeanPrecAstig * 1000);
fprintf('  Optimized PSF:        %.1f nm\n', MeanPrecOpt * 1000);
fprintf('  Improvement:          %.1f%%\n', ...
    100 * (1 - MeanPrecOpt/MeanPrecAstig));

% Save optimized PSF
save('optimized_3d_psf.mat', 'PSFStructOpt', 'PSFOpt');
```

#### Computing Localization Precision

Analyzing theoretical precision limits:

```matlab
% Create PSF model
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Lambda = 0.69;
PSFStruct.NA = 1.35;
PSFStruct.PixelSize = 0.108;
PSFStruct.Z = -1:0.1:1;

% Test different photon counts
PhotonLevels = [100, 300, 1000, 3000, 10000];
Bg = 20;

CRLB_All = cell(length(PhotonLevels), 1);

for ii = 1:length(PhotonLevels)
    [CRLB_All{ii}, ~] = smi_psf.PointSpreadFunction.crlbPSFPupil(...
        PSFStruct, PhotonLevels(ii), Bg, 0);
end

% Plot precision vs photons
FocusIdx = find(PSFStruct.Z == 0);
Precision_X = zeros(length(PhotonLevels), 1);
Precision_Y = zeros(length(PhotonLevels), 1);

for ii = 1:length(PhotonLevels)
    Precision_X(ii) = sqrt(CRLB_All{ii}(FocusIdx, 2));
    Precision_Y(ii) = sqrt(CRLB_All{ii}(FocusIdx, 1));
end

figure;
loglog(PhotonLevels, Precision_X * 1000, 'o-', 'LineWidth', 2);
hold on;
loglog(PhotonLevels, Precision_Y * 1000, 's-', 'LineWidth', 2);

% Theoretical sqrt(N) scaling
Theoretical = sqrt(CRLB_All{3}(FocusIdx, 2) * 1000 ./ PhotonLevels) * sqrt(1000);
loglog(PhotonLevels, Theoretical, 'k--', 'LineWidth', 1);

xlabel('Photons');
ylabel('Localization Precision (nm)');
legend('X', 'Y', '1/√N Scaling');
title(sprintf('CRLB vs Photons (Bg = %d)', Bg));
grid on;
```

---

### Zernike Class

**Purpose:** Low-level utilities for Zernike polynomial calculations and coordinate transformations.

**Key Concept:** Zernike polynomials form an orthogonal basis for representing wavefront aberrations on circular pupils. The Zernike class provides conversions between different indexing conventions (Noll, Wyant) and generates polynomial names.

**Class Definition:**
```matlab
classdef Zernike < handle
```

**All methods are static.**

**Indexing Conversion Methods:**

`zNM2Noll(N, M)` - Convert (N, M) radial/azimuthal indices to Noll index
```matlab
% Noll ordering: commonly used in optics
N = 2; M = 2;  % Astigmatism
NollIndex = smi_psf.Zernike.zNM2Noll(N, M);
% Returns: 6
```

`zNoll2NM(NollIndex)` - Convert Noll index to (N, M)
```matlab
NollIndex = 6;
[N, M] = smi_psf.Zernike.zNoll2NM(NollIndex);
% Returns: N=2, M=2
```

`zNM2Wyant(n, m)` - Convert (n, m) to Wyant index
```matlab
% Wyant ordering: alternative convention
n = 2; m = 2;
l = smi_psf.Zernike.zNM2Wyant(n, m);
```

`zWyant2NM(l)` - Convert Wyant index to (n, m)
```matlab
l = 5;
[n, m] = smi_psf.Zernike.zWyant2NM(l);
```

**Utility Methods:**

`zNamesNoll(ll)` - Get Zernike mode names (Noll ordering)
```matlab
NollIndices = [1, 2, 3, 4, 5, 6];
Names = smi_psf.Zernike.zNamesNoll(NollIndices);
% Returns: {'Piston', 'Tip', 'Tilt', 'Defocus', 'Astig 1', 'Astig 2'}
```

`zNamesWyant(ll)` - Get Zernike mode names (Wyant ordering)
```matlab
WyantIndices = [1, 2, 3, 4, 5];
Names = smi_psf.Zernike.zNamesWyant(WyantIndices);
```

`zNZNoll(n)` - Number of Zernike modes up to order n (Noll)
```matlab
n = 5;
nZ = smi_psf.Zernike.zNZNoll(n);
% Returns total modes through 5th radial order
```

`zNZWyant(n)` - Number of Zernike modes up to order n (Wyant)
```matlab
n = 5;
nZ = smi_psf.Zernike.zNZWyant(n);
```

`zProperNollIndex(l_max)` - Ensure proper Noll index
```matlab
% Validates and adjusts index to proper Noll convention
l_max = 20;
l_max_proper = smi_psf.Zernike.zProperNollIndex(l_max);
```

**Usage Examples:**

Understanding Zernike modes:

```matlab
% List first 15 Zernike modes
fprintf('Noll Index  |  N  M  |  Name\n');
fprintf('------------|--------|-------------\n');
for ii = 1:15
    [N, M] = smi_psf.Zernike.zNoll2NM(ii);
    Names = smi_psf.Zernike.zNamesNoll(ii);
    fprintf('    %2d      |  %d  %+d  |  %s\n', ii, N, M, Names{1});
end
```

Converting between conventions:

```matlab
% Working with literature values in different conventions
% Convert from paper using Wyant to Noll for smite
WyantCoefs = [1, 0, 0, 0.5, 1.2, 0];  % From paper

% Convert each coefficient
NollCoefs = zeros(15, 1);
for ii = 1:length(WyantCoefs)
    [n, m] = smi_psf.Zernike.zWyant2NM(ii);
    NollIdx = smi_psf.Zernike.zNM2Noll(n, m);
    NollCoefs(NollIdx) = WyantCoefs(ii);
end

% Now use in PSF generation
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.ZC_Phase = NollCoefs;
[PSF, ~] = smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct);
```

---

## PSF Model Structure Reference

### PSFStruct Fields

Complete specification of a PSF model:

**Optical Parameters:**
```matlab
Lambda          % Emission wavelength (micrometers)
NA              % Numerical aperture
N               % Refractive index of immersion medium
```

**Spatial Parameters:**
```matlab
PixelSize       % Lateral pixel size (micrometers)
SZ              % PSF output size (pixels, even number)
OSZ             % Oversampling size for computation (pixels)
Z               % Vector of Z positions (micrometers)
```

**Pupil Definition:**
```matlab
Pupil           % (OSZ x OSZ x 2) array
                % Pupil(:,:,1) = magnitude
                % Pupil(:,:,2) = phase (radians)
```

**Zernike Representation:**
```matlab
ZC_Mag          % Zernike coefficients for magnitude (Noll ordering)
ZC_Phase        % Zernike coefficients for phase (Noll ordering)
```

**OTF Modification:**
```matlab
OTFSigma        % [SigmaY, SigmaX] Gaussian smoothing (micrometers)
                % Set by rescaleOTF for experimental matching
```

**Example:**
```matlab
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();

% Modify for specific microscope
PSFStruct.Lambda = 0.69;      % 690 nm
PSFStruct.NA = 1.4;           % 1.4 NA objective
PSFStruct.N = 1.518;          % High-index oil
PSFStruct.PixelSize = 0.108;  % Camera pixel
PSFStruct.SZ = 64;            % 64x64 PSF
PSFStruct.Z = -1:0.05:1;      % 2 micron range

% Add astigmatism
PSFStruct.ZC_Phase(5) = 1.5;  % Noll index 5 = astigmatism

[PSF, ~] = smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct);
```

---

## Common Usage Patterns

### Creating a Custom 3D PSF Model

Complete workflow for generating a 3D PSF calibration:

```matlab
% Step 1: Define optical system
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Lambda = 0.69;      % Emission wavelength
PSFStruct.NA = 1.4;           % Objective NA
PSFStruct.N = 1.518;          % Immersion medium
PSFStruct.PixelSize = 0.108;  % Camera pixel size
PSFStruct.SZ = 64;            % PSF size
PSFStruct.Z = -0.8:0.05:0.8;  % 3D range

% Step 2: Design phase mask (astigmatism)
PSFStruct.ZC_Phase = zeros(15, 1);
PSFStruct.ZC_Phase(5) = 2.0;  % Astigmatism strength
PSFStruct.ZC_Phase(6) = 0.5;  % Second astigmatism angle

% Step 3: Generate PSF
[PSF, PSFStruct] = smi_psf.PointSpreadFunction.scalarPSFZernike(PSFStruct);

% Step 4: Validate precision
Photons = 1000;
Bg = 15;
[CRLB, DET] = smi_psf.PointSpreadFunction.crlbPSFPupil(...
    PSFStruct, Photons, Bg, 1);

% Step 5: Save for use in LocalizeData
save('3d_psf_model.mat', 'PSFStruct', 'PSF', 'CRLB');
```

### Experimental PSF Calibration Workflow

From bead measurements to calibrated model:

```matlab
% 1. Load experimental bead data
load('bead_measurements.mat', 'BeadStack', 'Z_um');
% BeadStack: Y x X x Z_positions

% 2. Preprocess
BeadStack = BeadStack - median(BeadStack(:));  % Remove offset
BeadStack(BeadStack < 0) = 0;

% Normalize each plane
for ii = 1:size(BeadStack, 3)
    BeadStack(:,:,ii) = BeadStack(:,:,ii) / sum(sum(BeadStack(:,:,ii)));
end

% 3. Setup PSFStruct
PSFStruct = smi_psf.PointSpreadFunction.createPSFStruct();
PSFStruct.Lambda = 0.69;
PSFStruct.NA = 1.35;
PSFStruct.N = 1.4;
PSFStruct.PixelSize = 0.108;
PSFStruct.Z = Z_um;
PSFStruct.SZ = size(BeadStack, 1);

% 4. Phase retrieval
[PSFStruct, PSF] = smi_psf.PointSpreadFunction.phaseRetrieval(...
    PSFStruct, BeadStack, 22, 81);

% 5. Validate
figure;
for ii = 1:size(PSF, 3)
    subplot(1,2,1); imagesc(BeadStack(:,:,ii)); title('Data');
    subplot(1,2,2); imagesc(PSF(:,:,ii)); title('Model');
    sgtitle(sprintf('Z = %.2f μm', Z_um(ii)));
    pause(0.2);
end

% 6. Save calibration
save('calibrated_psf.mat', 'PSFStruct', 'PSF');
```

### Comparing PSF Designs

Evaluating different 3D PSF strategies:

```matlab
% Common parameters
Lambda = 0.69;
NA = 1.4;
N = 1.518;
PixelSize = 0.108;
Z = -0.8:0.05:0.8;
Photons = 1000;
Bg = 20;

% Design 1: Simple astigmatism
PSF1 = smi_psf.PointSpreadFunction.createPSFStruct();
PSF1.Lambda = Lambda; PSF1.NA = NA; PSF1.N = N;
PSF1.PixelSize = PixelSize; PSF1.Z = Z;
PSF1.ZC_Phase = zeros(15, 1);
PSF1.ZC_Phase(5) = 2.0;
[~, PSF1] = smi_psf.PointSpreadFunction.scalarPSFZernike(PSF1);
[CRLB1, ~] = smi_psf.PointSpreadFunction.crlbPSFPupil(PSF1, Photons, Bg, 0);

% Design 2: Tetrapod PSF
PSF2 = smi_psf.PointSpreadFunction.createPSFStruct();
PSF2.Lambda = Lambda; PSF2.NA = NA; PSF2.N = N;
PSF2.PixelSize = PixelSize; PSF2.Z = Z;
PSF2.ZC_Phase = zeros(15, 1);
PSF2.ZC_Phase(8) = 1.5;  % Trefoil
PSF2.ZC_Phase(9) = 1.5;
[~, PSF2] = smi_psf.PointSpreadFunction.scalarPSFZernike(PSF2);
[CRLB2, ~] = smi_psf.PointSpreadFunction.crlbPSFPupil(PSF2, Photons, Bg, 0);

% Design 3: Prasad zones
PSF3 = smi_psf.PointSpreadFunction.createPSFStruct();
PSF3.Lambda = Lambda; PSF3.NA = NA; PSF3.N = N;
PSF3.PixelSize = PixelSize; PSF3.Z = Z;
[~, PSF3] = smi_psf.PointSpreadFunction.scalarPSFPrasadZone(PSF3, 4);
[CRLB3, ~] = smi_psf.PointSpreadFunction.crlbPSFPupil(PSF3, Photons, Bg, 0);

% Compare
figure;
subplot(2,2,1);
plot(Z, sqrt(CRLB1(:,3))*1000, 'LineWidth', 2); hold on;
plot(Z, sqrt(CRLB2(:,3))*1000, 'LineWidth', 2);
plot(Z, sqrt(CRLB3(:,3))*1000, 'LineWidth', 2);
xlabel('Z (μm)'); ylabel('Z Precision (nm)');
legend('Astigmatism', 'Tetrapod', 'Prasad');
title('Axial Precision');
grid on;

subplot(2,2,2);
plot(Z, sqrt(mean(CRLB1(:,1:2),2))*1000, 'LineWidth', 2); hold on;
plot(Z, sqrt(mean(CRLB2(:,1:2),2))*1000, 'LineWidth', 2);
plot(Z, sqrt(mean(CRLB3(:,1:2),2))*1000, 'LineWidth', 2);
xlabel('Z (μm)'); ylabel('Lateral Precision (nm)');
title('Lateral Precision');
grid on;

% Summary statistics
fprintf('\nMean Z precision over range:\n');
fprintf('  Astigmatism: %.1f nm\n', mean(sqrt(CRLB1(:,3)))*1000);
fprintf('  Tetrapod:    %.1f nm\n', mean(sqrt(CRLB2(:,3)))*1000);
fprintf('  Prasad:      %.1f nm\n', mean(sqrt(CRLB3(:,3)))*1000);
```

---

## Performance and Requirements

### Hardware Requirements

**GPU Required:**
- NVIDIA GPU with compute capability ≥5.0
- Minimum 2 GB GPU memory
- CUDA toolkit compatible with MATLAB version

**Check GPU:**
```matlab
g = gpuDevice;
fprintf('GPU: %s\n', g.Name);
fprintf('Compute capability: %s\n', g.ComputeCapability);
fprintf('Memory: %.1f GB\n', g.TotalMemory / 1e9);
```

### Memory Considerations

PSF computations use GPU memory intensively:

```matlab
% For large PSF stacks, reduce oversampling
PSFStruct.OSZ = 128;  % Instead of default 256

% Or process Z planes in batches
NZ = length(PSFStruct.Z);
BatchSize = 10;
PSF = zeros(PSFStruct.SZ, PSFStruct.SZ, NZ);

for ii = 1:BatchSize:NZ
    EndIdx = min(ii + BatchSize - 1, NZ);
    PSFTemp = PSFStruct;
    PSFTemp.Z = PSFStruct.Z(ii:EndIdx);
    [PSFBatch, ~] = smi_psf.PointSpreadFunction.scalarPSF(PSFTemp);
    PSF(:,:,ii:EndIdx) = PSFBatch;

    % Clear GPU
    reset(gpuDevice);
end
```

### Computational Time

Typical computation times (NVIDIA RTX 3080):

- `scalarPSF`: ~0.1 seconds (21 Z planes, 64x64)
- `scalarPSFZernike`: ~0.2 seconds (15 coefficients)
- `phaseRetrieval`: ~30-60 seconds (1000 iterations)
- `crlbPSFPupil`: ~2 seconds (21 Z planes)
- `optimPSFZernike`: ~10-30 minutes (15 coefficients)

---

## Troubleshooting

### GPU Out of Memory

**Problem:** "Out of memory on device" during PSF generation

**Solutions:**
```matlab
% 1. Reduce oversampling size
PSFStruct.OSZ = 128;  % Default is 256

% 2. Reduce PSF size
PSFStruct.SZ = 32;  % Default is 64

% 3. Process fewer Z planes at once
% Split Z range into batches (see Performance section)

% 4. Clear GPU between operations
reset(gpuDevice);
```

### Phase Retrieval Convergence Issues

**Problem:** Phase retrieval produces poor fits or diverges

**Solutions:**
```matlab
% 1. Check data quality
Data = Data - min(Data(:));  % Remove offset
Data = Data / max(Data(:));  % Normalize

% Check for hot pixels
MaxVal = max(Data(:));
MeanVal = mean(Data(:));
if MaxVal > 10 * MeanVal
    warning('Hot pixels detected - consider median filtering');
end

% 2. Adjust Zernike smoothing
MaxZCMag = 15;    % More aggressive smoothing
MaxZCPhase = 36;  % Less smoothing for complex aberrations

% 3. Verify Z positions match data
fprintf('Data has %d planes\n', size(Data, 3));
fprintf('PSFStruct.Z has %d values\n', length(PSFStruct.Z));
% Must be equal!

% 4. Check optical parameters
fprintf('NA/N ratio: %.3f (should be < 1)\n', PSFStruct.NA / PSFStruct.N);
```

### Incorrect PSF Appearance

**Problem:** Generated PSF doesn't look right

**Solutions:**
```matlab
% 1. Verify optical parameters
fprintf('Wavelength: %.3f μm\n', PSFStruct.Lambda);
fprintf('NA: %.2f\n', PSFStruct.NA);
fprintf('n: %.3f\n', PSFStruct.N);
fprintf('Pixel size: %.3f μm\n', PSFStruct.PixelSize);

% 2. Check sampling
PSFSize_um = PSFStruct.SZ * PSFStruct.PixelSize;
AiryRadius_um = 0.61 * PSFStruct.Lambda / PSFStruct.NA;
fprintf('PSF physical size: %.2f μm\n', PSFSize_um);
fprintf('Airy disk radius: %.3f μm\n', AiryRadius_um);
fprintf('Airy disk in pixels: %.1f\n', AiryRadius_um / PSFStruct.PixelSize);
% Should have ~10+ pixels across Airy disk

% 3. Visualize pupil
if isfield(PSFStruct, 'Pupil')
    figure;
    subplot(1,2,1);
    imagesc(PSFStruct.Pupil(:,:,1)); axis equal tight;
    title('Pupil Magnitude'); colorbar;

    subplot(1,2,2);
    imagesc(PSFStruct.Pupil(:,:,2)); axis equal tight;
    title('Pupil Phase'); colorbar;
end
```

### CRLB Values Seem Wrong

**Problem:** Precision estimates don't match expectations

**Solutions:**
```matlab
% 1. Verify photon and background levels
fprintf('Photons: %d\n', Photons);
fprintf('Background per pixel: %.1f\n', Bg);

% Typical values:
% - Bright emitters: 1000-5000 photons
% - Background: 10-50 photons/pixel

% 2. Check CRLB units
% CRLB is in variance (squared units)
% Must take sqrt for standard deviation
fprintf('X precision: %.1f nm\n', sqrt(CRLB(1,2)) * 1000);  % Correct
fprintf('NOT: %.1f nm\n', CRLB(1,2) * 1000);  % Wrong!

% 3. Verify PSF is normalized
PSFIntegral = sum(PSF(:));
fprintf('PSF integral: %.3f (should be ~1)\n', PSFIntegral);
```

---

## See Also

### Core Concepts
- [Architecture Overview](../core-concepts/architecture.md)
- [SMF Structure Guide](../core-concepts/smf-structure.md)

### API References
- [+smi_core Namespace](./smi-core.md) - Core localization algorithms
- [+smi_sim Namespace](./smi-sim.md) - Simulation tools

### Workflows
- [SMLM Analysis Workflow](../workflows/smlm-analysis.md)

### How-To Guides
- [How to Localize Molecules](../how-to/localize-molecules.md)
- [How to Use GPU](../how-to/use-gpu.md)

### External Resources
- [Zernike Polynomials (Wikipedia)](https://en.wikipedia.org/wiki/Zernike_polynomials)
- [Point Spread Function (Wikipedia)](https://en.wikipedia.org/wiki/Point_spread_function)

---

## References

### Key Citations

1. **Phase Retrieval:**
   - Hanser, B. M., et al. "Phase retrieval for high-numerical-aperture optical systems." Optics Letters 28.10 (2003): 801-803.

2. **Astigmatic 3D Imaging:**
   - Huang, B., et al. "Three-dimensional super-resolution imaging by stochastic optical reconstruction microscopy." Science 319.5864 (2008): 810-813.

3. **CRLB Theory:**
   - Ober, R. J., Ram, S., & Ward, E. S. "Localization accuracy in single-molecule microscopy." Biophysical Journal 86.2 (2004): 1185-1200.

4. **Engineered PSFs:**
   - Prasad, S. "Rotating point spread function via pupil-phase engineering." Optics Letters 38.4 (2013): 585-587.

---

## Summary

The +smi_psf namespace provides comprehensive tools for PSF modeling and analysis in smite:

**PSF Generation:** Create theoretical PSF models using scalar diffraction theory with arbitrary pupil functions defined via Zernike polynomials or direct pupil specification.

**Experimental Calibration:** Retrieve pupil magnitude and phase from experimental bead measurements using iterative phase retrieval with Zernike regularization.

**Precision Analysis:** Compute Cramer-Rao lower bounds to predict theoretical localization precision for given imaging conditions and PSF designs.

**3D Imaging:** Design and optimize phase masks for 3D super-resolution, including astigmatism, tetrapod, and custom engineered PSFs.

**Zernike Utilities:** Convert between indexing conventions, generate polynomial images, and decompose arbitrary phase patterns into Zernike basis.

These classes are essential for advanced 3D localization, PSF calibration, and optimizing microscope performance for single molecule imaging. Understanding +smi_psf enables custom PSF design, precise 3D localization, and quantitative analysis of imaging system performance.
