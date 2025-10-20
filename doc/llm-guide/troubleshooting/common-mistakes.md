---
title: "Common Mistakes and How to Fix Them"
category: "troubleshooting"
level: "beginner"
tags: ["troubleshooting", "errors", "mistakes", "debugging", "coordinate-system", "units", "camera", "smf", "smd"]
prerequisites: ["../getting-started/quickstart.md", "../core-concepts/architecture.md"]
related: ["../core-concepts/coordinate-system.md", "../core-concepts/smf-structure.md", "../workflows/frame-connection.md", "../workflows/spt-tracking.md"]
summary: "Comprehensive guide to frequently made errors in smite, why they happen, and how to fix them"
estimated_time: "25 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# Common Mistakes and How to Fix Them

## Purpose

This guide documents the most frequent errors users make when working with smite, explains why they occur, and provides clear solutions. Understanding these common pitfalls will save you hours of debugging and help you develop correct analysis habits from the start. Whether you're mixing up coordinates, using wrong units, or misunderstanding data structures, this guide will get you back on track.

## Prerequisites

- Basic understanding of [smite architecture](../core-concepts/architecture.md)
- Completion of [Quick Start](../getting-started/quickstart.md)
- Familiarity with MATLAB basics

## Overview

Common mistakes in smite fall into several categories:

- **Coordinate system errors**: Swapping X/Y, wrong indexing, forgetting pixel centers
- **Unit confusion**: Mixing pixels with micrometers or nanometers
- **Camera parameter mistakes**: Wrong gain, offset, or pixel size values
- **Setup issues**: Forgetting to call setupSMITE, incorrect paths
- **SMF configuration errors**: Misunderstanding parameter effects
- **Data structure misuse**: Incorrect SMD/TR field access
- **Workflow confusion**: Mixing frame connection with SPT tracking

This guide covers each category with examples of what goes wrong and how to fix it.

## Coordinate System Mistakes

### Mistake 1: Swapping X and Y Coordinates

**What happens:**
Your localizations appear transposed or rotated compared to the original image.

**Why it's wrong:**
MATLAB uses row-column indexing: `image(row, column)` = `image(Y, X)`. Many users reverse this, using `image(X, Y)`, which swaps the coordinates.

**Example of wrong code:**
```matlab
% WRONG: X and Y are swapped
for i = 1:length(SMD.X)
    pixel_value = image(SMD.X(i), SMD.Y(i));  % X is first index!
end
```

**Why this fails:**
- `SMD.X` represents horizontal position (column index)
- `SMD.Y` represents vertical position (row index)
- Array indexing is `image(row, column)`, so `Y` must come first

**Correct solution:**
```matlab
% CORRECT: Y (row) first, X (column) second
for i = 1:length(SMD.X)
    % For integer indexing
    row = round(SMD.Y(i));
    col = round(SMD.X(i));
    pixel_value = image(row, col);

    % OR for sub-pixel interpolation
    pixel_value = interp2(image, SMD.X(i), SMD.Y(i));
end
```

**Verification:**
```matlab
% Check if your coordinates make sense
fprintf('Image size: %d rows × %d columns\n', size(image, 1), size(image, 2));
fprintf('X range: %.1f to %.1f (should match columns)\n', min(SMD.X), max(SMD.X));
fprintf('Y range: %.1f to %.1f (should match rows)\n', min(SMD.Y), max(SMD.Y));
```

### Mistake 2: Forgetting Pixel-Centered Coordinates

**What happens:**
Visualizations are off by half a pixel, ROIs don't align correctly, or images appear shifted.

**Why it's wrong:**
In MATLAB (and smite), coordinate (1,1) refers to the **center** of the top-left pixel, not its corner. Pixel corners are at half-integer positions (0.5, 1.5, 2.5, etc.).

**Example of wrong code:**
```matlab
% WRONG: Assumes pixels go from 1 to N
figure;
imagesc(image);
hold on;
plot(SMD.X, SMD.Y, 'r.', 'MarkerSize', 10);
% Localizations appear shifted!
```

**Why this fails:**
`imagesc(image)` sets axis limits to [1, N], but pixels actually extend from [0.5, N+0.5]. Your localizations at pixel centers appear offset.

**Correct solution:**
```matlab
% CORRECT: Account for pixel-centered coordinates
[rows, cols] = size(image);
x_extent = [0.5, cols + 0.5];
y_extent = [0.5, rows + 0.5];

figure;
imagesc(x_extent, y_extent, image);
hold on;
plot(SMD.X, SMD.Y, 'r.', 'MarkerSize', 10);
axis equal tight;
xlabel('X (pixels)'); ylabel('Y (pixels)');
```

**Key concept:**
```
Pixel edges:  0.5    1.5    2.5    3.5    4.5    5.5
               |      |      |      |      |      |
Pixel centers:   1.0    2.0    3.0    4.0    5.0
```

### Mistake 3: Rounding Sub-Pixel Positions

**What happens:**
You lose the precision advantage of super-resolution microscopy.

**Why it's wrong:**
Localization positions are real numbers with sub-pixel precision (e.g., 10.347 pixels), not integers. Rounding destroys this precision and defeats the purpose of SMLM.

**Example of wrong code:**
```matlab
% WRONG: Loses sub-pixel precision
X_int = round(SMD.X);
Y_int = round(SMD.Y);
SR_image(Y_int, X_int) = 1;  % Only uses integer positions!
```

**Why this fails:**
- Typical localization precision: 0.05-0.2 pixels (5-20 nm for 100 nm pixels)
- `round()` quantizes to nearest pixel (100 nm for 100 nm pixels)
- You've just thrown away 5-20× precision improvement

**Correct solution:**
```matlab
% CORRECT: Use sub-pixel rendering
SR_zoom = 20;  % 20× magnification
SR_size = [size(image, 1), size(image, 2)] * SR_zoom;
SR_image = zeros(SR_size);

% Map to SR grid (keeps sub-pixel precision)
X_SR = SMD.X * SR_zoom;
Y_SR = SMD.Y * SR_zoom;

% Render with Gaussian
sigma_SR = 0.5 * SR_zoom;  % Sub-pixel rendering
for i = 1:length(SMD.X)
    SR_image = SR_image + smi_vis.gaussianImage(...
        [Y_SR(i), X_SR(i)], sigma_SR, SR_size);
end
```

**Or use smite's built-in visualization:**
```matlab
% Even better: Use GenerateImages
GI = smi_vis.GenerateImages(SMD, SMF);
SR_image = GI.genSRImage();  % Handles sub-pixel rendering correctly
```

## Unit Confusion

### Mistake 4: Mixing Pixels and Physical Units

**What happens:**
Your ROI selections are wrong, distance measurements are meaningless, or thresholds don't work.

**Why it's wrong:**
SMD stores positions in **pixels**, but scientific measurements need **micrometers** or **nanometers**. Mixing units leads to incorrect calculations.

**Example of wrong code:**
```matlab
% WRONG: Mixing units
ROI_nm = [100, 100, 500, 500];  % Intended as nanometers
in_ROI = (SMD.X >= ROI_nm(1)) & (SMD.X <= ROI_nm(3));  % But SMD.X is in pixels!
```

**Why this fails:**
If pixels are 100 nm each, `SMD.X = 1` means 100 nm from origin, not 1 nm. Your ROI from 100-500 nm is actually selecting 100-500 pixels (10,000-50,000 nm)!

**Correct solution:**
```matlab
% CORRECT: Convert units consistently
ROI_nm = [100, 100, 500, 500];  % [X_min, Y_min, X_max, Y_max] in nm
pixel_size_nm = SMD.PixelSize * 1000;  % Convert μm to nm

% Convert ROI to pixels
ROI_pixels = ROI_nm / pixel_size_nm;

% Now select in same units
in_ROI = (SMD.X >= ROI_pixels(1)) & (SMD.X <= ROI_pixels(3)) & ...
         (SMD.Y >= ROI_pixels(2)) & (SMD.Y <= ROI_pixels(4));

fprintf('ROI: [%.1f, %.1f] to [%.1f, %.1f] pixels\n', ROI_pixels);
fprintf('Selected %d localizations\n', sum(in_ROI));
```

**Unit conversion reference:**
```matlab
% Pixels → Physical
X_um = SMD.X * SMD.PixelSize;              % micrometers
X_nm = SMD.X * SMD.PixelSize * 1000;       % nanometers

% Physical → Pixels
X_pixels = X_um / SMD.PixelSize;           % from micrometers
X_pixels = X_nm / (SMD.PixelSize * 1000);  % from nanometers
```

### Mistake 5: Wrong PixelSize Units

**What happens:**
All your distance measurements are off by factors of 1000 or more.

**Why it's wrong:**
`SMF.Data.PixelSize` expects **micrometers**, but users often enter nanometers or millimeters.

**Example of wrong code:**
```matlab
% WRONG: Using nanometers instead of micrometers
SMF.Data.PixelSize = 100;  % Intended as 100 nm, but interpreted as 100 μm!
```

**Why this fails:**
smite interprets 100 as 100 μm = 100,000 nm. Your 100 nm pixels are now 1000× too large, making all distances wrong.

**Correct solution:**
```matlab
% CORRECT: Always use micrometers for PixelSize
SMF.Data.PixelSize = 0.1;    % 100 nm = 0.1 μm ✓
SMF.Data.PixelSize = 0.108;  % 108 nm = 0.108 μm ✓
SMF.Data.PixelSize = 0.065;  % 65 nm = 0.065 μm ✓
```

**Verification:**
```matlab
% Check if PixelSize is reasonable
if SMF.Data.PixelSize < 0.05 || SMF.Data.PixelSize > 0.5
    warning('PixelSize %.3f μm seems unusual. Check units!', SMF.Data.PixelSize);
    fprintf('Common values: 0.1 μm (100 nm), 0.108 μm (108 nm), 0.065 μm (65 nm)\n');
end

% Verify against known structures
fprintf('With PixelSize = %.3f μm:\n', SMF.Data.PixelSize);
fprintf('  100 pixels = %.1f nm\n', 100 * SMF.Data.PixelSize * 1000);
fprintf('  50 nm = %.1f pixels\n', 50 / (SMF.Data.PixelSize * 1000));
```

## Camera Parameter Mistakes

### Mistake 6: Using Wrong Camera Gain

**What happens:**
Photon counts are incorrect, thresholds don't work, localization precision is wrong.

**Why it's wrong:**
The gain converts camera output (ADU) to physical photons. Wrong gain means wrong photon counts for everything.

**Example of wrong code:**
```matlab
% WRONG: Using default gain without calibration
SMF.Data.CameraGain = 1;  % Default, but your camera has gain ≠ 1
SMF.Data.CameraType = 'EMCCD';
```

**Why this fails:**
- Typical EMCCD gain: 2-5 ADU/photon (not 1)
- Typical sCMOS gain: 0.3-0.6 ADU/photon (not 1)
- With gain = 1, photon counts are wrong by factors of 2-10×

**Correct solution:**
```matlab
% CORRECT: Measure or obtain actual gain from calibration

% For EMCCD (uniform gain)
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;  % Measured from calibration
SMF.Data.CameraOffset = 100;  % Measured background

% For sCMOS (pixel-wise gain)
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = '/path/to/calibration.mat';
% calibration.mat contains CameraGain, CameraOffset, CameraVar as images
```

**How to get the right gain:**
```matlab
% Method 1: Use calibration file (best)
% See doc/llm-guide/how-to/calibrate-camera.md

% Method 2: Quick estimate (if you must)
% Take dark frames and illuminated frames
dark_mean = mean(dark_frames(:));  % Offset
illum_mean = mean(illuminated_frames(:));

% Expected photons per pixel (from light source)
expected_photons = 100;
gain_estimate = (illum_mean - dark_mean) / expected_photons;
fprintf('Estimated gain: %.2f ADU/photon\n', gain_estimate);
```

### Mistake 7: Forgetting Camera Offset

**What happens:**
Dark images aren't properly subtracted, causing bias in photon estimates and wrong background levels.

**Why it's wrong:**
Cameras add a baseline offset to prevent negative values. If not accounted for, this offset is counted as signal.

**Example of wrong code:**
```matlab
% WRONG: Zero offset when camera actually has ~100 ADU baseline
SMF.Data.CameraOffset = 0;  % Default, but incorrect
```

**Why this fails:**
With 100 ADU offset, every pixel has 100 ADU that's not signal. When gain = 2, this looks like 50 extra photons everywhere!

**Correct solution:**
```matlab
% CORRECT: Measure offset from dark frames
% Take 100+ frames with no light
dark_frames = ... % Your dark acquisition
SMF.Data.CameraOffset = mean(dark_frames(:));

fprintf('Measured camera offset: %.1f ADU\n', SMF.Data.CameraOffset);

% Typical values:
% EMCCD: 100-500 ADU
% sCMOS: 100-200 ADU per pixel (varies by pixel)
```

### Mistake 8: EMCCD vs sCMOS Configuration

**What happens:**
Fitting fails, precision estimates are wrong, or performance is terrible.

**Why it's wrong:**
EMCCD and sCMOS have fundamentally different noise models and require different configuration.

**Example of wrong code:**
```matlab
% WRONG: sCMOS camera configured as EMCCD
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 0.45;  % sCMOS-like gain but wrong type
```

**Why this fails:**
- EMCCD: Uses Poisson + Gaussian noise model, uniform gain
- sCMOS: Uses pixel-specific gain/variance maps
- Wrong type → wrong noise model → wrong fitting → wrong results

**Correct solution:**
```matlab
% CORRECT: sCMOS with calibration
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = '/path/to/scmos_calibration.mat';
% Calibration file contains:
%   CameraGain:   2D array of gain per pixel
%   CameraOffset: 2D array of offset per pixel
%   CameraVar:    2D array of variance per pixel

% CORRECT: EMCCD with scalar parameters
SMF.Data.CameraType = 'EMCCD';
SMF.Data.CameraGain = 2.5;      % Single value for whole sensor
SMF.Data.CameraOffset = 100;    % Single value
SMF.Data.CameraNoise = 8.5;     % Read noise in electrons
```

**How to tell which camera you have:**
- **EMCCD**: Photometrics Evolve, Andor iXon, Hamamatsu ImagEM
- **sCMOS**: Hamamatsu Orca Flash, PCO Edge, Andor Zyla/Sona

## Setup and Path Issues

### Mistake 9: Forgetting to Call setupSMITE

**What happens:**
Commands fail with "Undefined function or variable" errors for smite classes.

**Why it's wrong:**
MATLAB doesn't know where smite is unless you add it to the path and run `setupSMITE`.

**Example error:**
```matlab
>> SMF = smi_core.SingleMoleculeFitting()
Undefined function or variable 'smi_core.SingleMoleculeFitting'.
```

**Why this fails:**
The `+smi_core` folder isn't on MATLAB's path, so MATLAB can't find the namespace.

**Correct solution:**
```matlab
% One-time setup: Add to startup.m
% Location: ~/Documents/MATLAB/startup.m (create if doesn't exist)

% Add these lines to startup.m:
addpath '~/Documents/MATLAB/smite/MATLAB'  % Or your smite location
setupSMITE

% Then restart MATLAB

% Verify it worked:
>> which setupSMITE
~/Documents/MATLAB/smite/MATLAB/setupSMITE.m

>> SMF = smi_core.SingleMoleculeFitting()
SMF =
  struct with fields:
    Data: [1×1 struct]
    ...
```

**Quick fix for current session:**
```matlab
% If you forgot to set up startup.m
addpath '~/Documents/MATLAB/smite/MATLAB'
setupSMITE

% Now smite commands will work (until you close MATLAB)
```

### Mistake 10: Wrong File Paths

**What happens:**
Data doesn't load, results don't save, or mysterious "file not found" errors appear.

**Why it's wrong:**
Incorrect path separators (/ vs \), relative paths that change with directory, or typos in file names.

**Example of wrong code:**
```matlab
% WRONG: Relative path (breaks if you change directory)
SMF.Data.FileDir = '../data';
SMF.Data.FileName = {'experiment1.h5'};
```

**Why this fails:**
If you `cd` to a different directory or run from a script elsewhere, `../data` points to the wrong location.

**Correct solution:**
```matlab
% CORRECT: Use absolute paths
SMF.Data.FileDir = 'C:\Users\YourName\Documents\data';  % Windows
% OR
SMF.Data.FileDir = '/home/username/data';  % Linux/Mac

SMF.Data.FileName = {'experiment1.h5'};

% Verify before running
if ~exist(fullfile(SMF.Data.FileDir, SMF.Data.FileName{1}), 'file')
    error('File not found: %s', fullfile(SMF.Data.FileDir, SMF.Data.FileName{1}));
end
```

**Platform-independent paths:**
```matlab
% BEST: Use fullfile() for cross-platform compatibility
base_dir = 'C:\Users\YourName\Documents';  % or '/home/username'
SMF.Data.FileDir = fullfile(base_dir, 'data', 'experiment1');
SMF.Data.FileName = {'Cell1.h5', 'Cell2.h5'};

% fullfile() handles / vs \ automatically
```

## SMF Configuration Errors

### Mistake 11: Conflicting Parameter Settings

**What happens:**
Unexpected behavior, no localizations found, or incorrect results.

**Why it's wrong:**
Some SMF parameters interact or override each other. Setting conflicting values causes confusion.

**Example of wrong code:**
```matlab
% WRONG: AutoThreshold enabled but threshold manually set
SMF.Thresholding.AutoThreshLogL = true;  % Auto-calculate threshold
SMF.Thresholding.MinLogL = 50;          % Manual threshold (ignored!)
```

**Why this fails:**
When `AutoThreshLogL = true`, smite calculates the threshold automatically and ignores your manual `MinLogL` setting. Your carefully chosen value of 50 never gets used.

**Correct solution:**
```matlab
% Option 1: Use auto-threshold (recommended)
SMF.Thresholding.AutoThreshLogL = true;
% Don't set MinLogL - it will be calculated

% Option 2: Use manual threshold
SMF.Thresholding.AutoThreshLogL = false;
SMF.Thresholding.MinLogL = 50;  % Now this is used

% Check what's actually being used
if SMF.Thresholding.AutoThreshLogL
    fprintf('Using auto-calculated threshold\n');
else
    fprintf('Using manual threshold: %.1f\n', SMF.Thresholding.MinLogL);
end
```

**Other conflicting parameters:**
```matlab
% Frame connection vs tracking
SMF.FrameConnection.On = true;   % For SMLM (blinking)
SMF.Tracking.On = true;          % For SPT (motion)
% WARNING: Don't enable both! Choose based on experiment type

% ROI specification
SMF.Data.DataROI = [10, 10, 100, 100];  % Specific ROI
% But if you use ROITools, it may override this
```

### Mistake 12: BoxSize Too Small or Too Large

**What happens:**
Missing localizations (too small) or slow performance with bad fits (too large).

**Why it's wrong:**
`BoxSize` defines the fitting window around each molecule. Too small truncates the PSF; too large includes multiple molecules and slows fitting.

**Example of wrong code:**
```matlab
% WRONG: BoxSize too small for PSF
SMF.Fitting.PSFSigma = 1.5;  % PSF width
SMF.BoxFinding.BoxSize = 3;  % Only 3×3 pixels!
```

**Why this fails:**
A Gaussian PSF with σ=1.5 pixels extends ±3σ = ±4.5 pixels. A 3×3 box captures only the peak, missing most of the PSF. Fits will be biased and imprecise.

**Correct solution:**
```matlab
% CORRECT: BoxSize should capture full PSF
SMF.Fitting.PSFSigma = 1.5;  % PSF width

% Rule of thumb: BoxSize ≥ 6*PSFSigma + 1 (make it odd)
BoxSize_min = 6 * SMF.Fitting.PSFSigma + 1;  % = 10
SMF.BoxFinding.BoxSize = 2*floor(BoxSize_min/2) + 1;  % Make odd: 11

% Typical values:
% Dense data: 7-9 pixels (faster, okay if well-separated)
% Sparse data: 11-15 pixels (better fits, no overlap issues)
```

**Check your setting:**
```matlab
% Verify BoxSize is appropriate
PSF_extent = 3 * SMF.Fitting.PSFSigma;  % 99% of PSF within ±3σ
fprintf('PSF extends ±%.1f pixels\n', PSF_extent);
fprintf('BoxSize %d captures ±%.1f pixels\n', ...
    SMF.BoxFinding.BoxSize, SMF.BoxFinding.BoxSize/2);

if SMF.BoxFinding.BoxSize < 6*SMF.Fitting.PSFSigma
    warning('BoxSize may be too small! Consider increasing.');
end
```

## Data Structure Misunderstandings

### Mistake 13: Modifying SMD Fields Incorrectly

**What happens:**
Results become inconsistent, analysis breaks, or data is corrupted.

**Why it's wrong:**
SMD fields are interdependent. Changing one without updating others breaks the structure.

**Example of wrong code:**
```matlab
% WRONG: Remove some localizations without updating other fields
SMD.X = SMD.X(1:1000);  % Keep first 1000
% But SMD.Y, SMD.Photons, SMD.FrameNum, etc. still have all data!
```

**Why this fails:**
Now `SMD.X` has 1000 elements but `SMD.Y` has (say) 5000. Trying to plot `plot(SMD.X, SMD.Y)` fails because lengths don't match.

**Correct solution:**
```matlab
% CORRECT: Use logical indexing to filter consistently
keep = 1:1000;  % Or any logical condition

% Apply to ALL vector fields simultaneously
SMD.X = SMD.X(keep);
SMD.Y = SMD.Y(keep);
SMD.Photons = SMD.Photons(keep);
SMD.Bg = SMD.Bg(keep);
SMD.X_SE = SMD.X_SE(keep);
SMD.Y_SE = SMD.Y_SE(keep);
SMD.FrameNum = SMD.FrameNum(keep);
SMD.DatasetNum = SMD.DatasetNum(keep);
% ... and all other vector fields

% Better: Use smi_core.SelectSMD utility (if available)
% Or create a helper function:
```

**Helper function for safe filtering:**
```matlab
function SMD = filterSMD(SMD, keep)
    % Apply logical index to all SMD vector fields
    fields = fieldnames(SMD);
    for i = 1:length(fields)
        field = fields{i};
        if isvector(SMD.(field)) && length(SMD.(field)) == length(keep)
            SMD.(field) = SMD.(field)(keep);
        end
    end
end

% Usage:
keep = (SMD.Photons > 500) & (SMD.X_SE < 0.2);
SMD = filterSMD(SMD, keep);
```

### Mistake 14: Confusing SMD and TR Structures

**What happens:**
Code fails when you try to access fields that don't exist or are organized differently.

**Why it's wrong:**
SMD is a single structure with all localizations. TR is an array of SMD structures, one per trajectory.

**Example of wrong code:**
```matlab
% WRONG: Treating TR like SMD
TR = smi_core.TrackingResults(...);
plot(TR.X, TR.Y);  % ERROR: TR doesn't have X and Y fields!
```

**Why this fails:**
`TR` is an array: `TR(1)`, `TR(2)`, etc. are individual trajectories (each is an SMD). There's no single `TR.X` field.

**Correct solution:**
```matlab
% CORRECT: TR is an array of trajectories
% Each TR(i) is an SMD for one trajectory

% Plot all trajectories
figure; hold on;
for i = 1:length(TR)
    plot(TR(i).X, TR(i).Y, '-o');
end
xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('%d Trajectories', length(TR)));

% Access specific trajectory
traj_3 = TR(3);  % Third trajectory (this is an SMD)
plot(traj_3.X, traj_3.Y, 'r-o', 'LineWidth', 2);

% Get all positions (concatenate across trajectories)
all_X = [];
all_Y = [];
for i = 1:length(TR)
    all_X = [all_X; TR(i).X];
    all_Y = [all_Y; TR(i).Y];
end
```

**Key differences:**
```matlab
% SMD: Single structure
% Fields: X, Y, Photons, ... (all vectors)
fprintf('SMD has %d localizations\n', length(SMD.X));

% TR: Array of SMD structures
% TR(i): One trajectory (is an SMD)
fprintf('TR has %d trajectories\n', length(TR));
fprintf('Trajectory 1 has %d localizations\n', length(TR(1).X));
```

## Workflow Confusion

### Mistake 15: Using Frame Connection for SPT Data

**What happens:**
Moving particles get incorrectly linked, trajectories are broken, motion is lost.

**Why it's wrong:**
Frame connection assumes stationary emitters (SMLM). It links nearby localizations across frames, which is wrong for moving particles.

**Example of wrong code:**
```matlab
% WRONG: Frame connection enabled for tracking data
SMF.FrameConnection.On = true;     % For stationary emitters
% But your data has moving particles!

% Run SMLM pipeline
SMLMobj = smi.SMLM(imageStack, SMF);
SMLMobj.run();  % Wrong: treats motion as noise
```

**Why this fails:**
A particle at (10, 10) in frame 1 moves to (12, 15) in frame 2. Frame connection sees these as different emitters (too far apart) or incorrectly links them, destroying the trajectory.

**Correct solution:**
```matlab
% CORRECT: Use SPT workflow for moving particles
SMF.FrameConnection.On = false;  % Disable for SPT
SMF.Tracking.On = true;           % Enable tracking instead

% Configure tracking
SMF.Tracking.D = 0.1;            % Expected diffusion (pixels²/frame)
SMF.Tracking.MaxDistFF = 5;      # Max distance frame-to-frame
SMF.Tracking.MaxFrameGap = 5;    # Allow gaps (blinking)

% Run SPT pipeline
SPTobj = smi.SPT(imageStack, SMF);
[SMD, ~, TR] = SPTobj.run();      # Returns trajectories in TR

% Analyze motion
MSD = smi_stat.computeMSD(TR);
```

**When to use which:**
```matlab
% Use Frame Connection (FrameConnection.On = true):
% - DNA-PAINT, dSTORM, PALM (stationary emitters)
% - Goal: Improve precision by combining blinks
% - Output: Combined SMD (one localization per emitter)

% Use SPT Tracking (Tracking.On = true):
% - Live-cell tracking, single particle tracking
% - Goal: Follow particle trajectories, measure dynamics
% - Output: TR array (trajectories over time)
```

### Mistake 16: Wrong MaxFrameGap for Blinking

**What happens:**
Frame connection misses blinks (too small) or incorrectly merges different emitters (too large).

**Why it's wrong:**
`MaxFrameGap` in frame connection limits how many frames a molecule can be dark before it's considered a new emitter.

**Example of wrong code:**
```matlab
% WRONG: MaxFrameGap = 1 for highly blinking dye
SMF.FrameConnection.MaxFrameGap = 1;  # Only allows 1 frame gap
% But your dye blinks off for 5-10 frames regularly!
```

**Why this fails:**
Your molecule blinks off for 8 frames. Frame connection sees the reappearance as a new emitter (gap too large). One molecule becomes many in your results.

**Correct solution:**
```matlab
% CORRECT: Set MaxFrameGap based on dye kinetics
% For DNA-PAINT: typically 1-3 frames (fast on/off)
SMF.FrameConnection.MaxFrameGap = 3;

% For dSTORM: typically 5-20 frames (longer blinks)
SMF.FrameConnection.MaxFrameGap = 10;

% For PALM: varies by photoactivation
SMF.FrameConnection.MaxFrameGap = 5;

% Analyze your data to decide:
% Look at LocalizeData results before frame connection
% Check how long molecules stay off
```

**Estimating from data:**
```matlab
% After initial localization (before frame connection)
SMD_raw = LD.genLocalizations();

% Find same molecule in consecutive frames (within 1 pixel)
same_mol = [];
for i = 1:length(SMD_raw.X)-1
    dist = sqrt((SMD_raw.X(i) - SMD_raw.X(i+1))^2 + ...
                (SMD_raw.Y(i) - SMD_raw.Y(i+1))^2);
    if dist < 1.0
        gap = SMD_raw.FrameNum(i+1) - SMD_raw.FrameNum(i);
        same_mol = [same_mol; gap];
    end
end

% Analyze gaps
fprintf('Frame gap distribution:\n');
fprintf('  Mean: %.1f frames\n', mean(same_mol));
fprintf('  Median: %.1f frames\n', median(same_mol));
fprintf('  95th percentile: %.1f frames\n', prctile(same_mol, 95));

% Set MaxFrameGap to ~95th percentile
SMF.FrameConnection.MaxFrameGap = ceil(prctile(same_mol, 95));
```

## Key Takeaways

1. **Coordinate system**: MATLAB uses `(row, column)` = `(Y, X)` indexing. Always put Y first when accessing images.

2. **Pixel centers**: Coordinate (1,1) is at the **center** of the top-left pixel. Image extent is [0.5, size+0.5].

3. **Sub-pixel precision**: Never round localization positions. Use sub-pixel rendering or interpolation.

4. **Units matter**: SMD positions are in **pixels**. PixelSize is in **micrometers**. Always convert consistently.

5. **PixelSize units**: Always enter in micrometers (0.1 for 100 nm, not 100).

6. **Camera calibration**: Use measured gain and offset, not defaults. EMCCD vs sCMOS requires different setup.

7. **Path setup**: Call `setupSMITE` in startup.m. Use absolute paths for data files.

8. **SMF parameters**: Understand interactions. AutoThreshold overrides manual values. BoxSize must fit PSF.

9. **SMD consistency**: When filtering SMD, update **all** vector fields together. Don't modify fields independently.

10. **SMD vs TR**: SMD is a single structure with vectors. TR is an **array** of SMD structures (one per trajectory).

11. **Frame connection vs tracking**: Frame connection for stationary emitters (SMLM). Tracking for moving particles (SPT). Never both.

12. **MaxFrameGap**: Set based on actual blinking behavior. Analyze your data to find the right value.

## Debugging Checklist

When something goes wrong, check these in order:

**Setup issues:**
- [ ] Did you call `setupSMITE`? Check with `which setupSMITE`
- [ ] Are paths correct? Verify file exists with `exist(fullfile(...))`
- [ ] Is smite version current? Check with `setupSMITE` output

**Coordinate issues:**
- [ ] Are you using `image(Y, X)` not `image(X, Y)`?
- [ ] Do plot axes account for pixel centers (0.5 offset)?
- [ ] Are you preserving sub-pixel precision (no rounding)?

**Unit issues:**
- [ ] Is PixelSize in micrometers? (0.1 for 100 nm, not 100)
- [ ] Are you converting units consistently throughout?
- [ ] Do distances make physical sense? (nm, μm scale)

**Camera issues:**
- [ ] Is CameraType correct (EMCCD vs SCMOS)?
- [ ] Do you have measured Gain and Offset (not defaults)?
- [ ] For sCMOS, do you have a calibration file?

**SMF configuration:**
- [ ] Is BoxSize appropriate for PSFSigma? (≥ 6×PSFSigma)
- [ ] Are thresholding parameters reasonable? (not too strict)
- [ ] Are you using the right workflow (Frame Connection vs Tracking)?

**Data structure issues:**
- [ ] When filtering SMD, did you update all fields?
- [ ] Are you treating TR as an array, not a single structure?
- [ ] Are field lengths consistent after modifications?

## See Also

- [Coordinate System Guide](../core-concepts/coordinate-system.md) - Detailed coordinate system explanation
- [SMF Structure Reference](../core-concepts/smf-structure.md) - Complete SMF field documentation
- [SMD Structure Reference](../core-concepts/smd-structure.md) - Complete SMD field documentation
- [Frame Connection Workflow](../workflows/frame-connection.md) - When and how to use frame connection
- [SPT Tracking Workflow](../workflows/spt-tracking.md) - Single particle tracking details
- [Camera Calibration Guide](../how-to/calibrate-camera.md) - How to calibrate your camera properly
- [Installation Guide](../getting-started/installation.md) - Setup and path configuration

## Next Steps

After avoiding these mistakes:
1. **Verify your setup**: Run through the debugging checklist
2. **Calibrate properly**: Follow [Camera Calibration](../how-to/calibrate-camera.md)
3. **Understand coordinates**: Read [Coordinate System](../core-concepts/coordinate-system.md) in depth
4. **Choose right workflow**: Study [Frame Connection](../workflows/frame-connection.md) vs [SPT Tracking](../workflows/spt-tracking.md)
5. **Build good habits**: Always verify units, check parameter interactions, test on small datasets first
