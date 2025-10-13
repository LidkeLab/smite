---
title: "File Formats Reference"
category: "reference"
level: "intermediate"
tags: ["file-formats", "hdf5", "mat", "calibration", "registration", "results", "data-io"]
prerequisites: ["../core-concepts/smf-structure.md", "../core-concepts/smd-structure.md"]
related: ["../how-to/load-data.md", "../getting-started/first-analysis.md", "dependencies.md"]
summary: "Complete reference for all file formats used in smite including HDF5, MAT, calibration, registration, and results files"
estimated_time: "20-25 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# File Formats Reference

## Purpose

smite works with multiple file formats for different purposes: raw data acquisition, camera calibration, channel registration, and analysis results. Understanding these formats is essential for:

- Loading and organizing raw microscopy data
- Configuring camera parameters correctly
- Setting up multi-channel experiments
- Saving and sharing analysis results
- Troubleshooting data loading issues
- Integrating smite with other analysis tools

This document provides complete specifications for all file formats used in smite workflows.

## Prerequisites

- Understanding of [SMF structure](../core-concepts/smf-structure.md)
- Understanding of [SMD structure](../core-concepts/smd-structure.md)
- Basic familiarity with MATLAB data types

## File Format Overview

smite uses five main categories of files:

| Format | Purpose | Extension | Required |
|--------|---------|-----------|----------|
| **Raw Data** | Microscopy image stacks | `.h5`, `.mat` | Yes |
| **Calibration** | Camera gain/offset/noise maps | `.mat` | Recommended |
| **Channel Registration** | Multi-channel alignment | `.mat` | For multi-channel |
| **Results** | Localization results | `.mat` | Output |
| **Parameters** | Analysis configuration | `.mat` | Optional |

---

## Raw Data Formats

Raw microscopy data contains the image sequences acquired from your camera. smite supports two formats.

### HDF5 Format (.h5)

**Recommended format** for new acquisitions due to efficient storage and standardized structure.

#### File Structure

HDF5 files use a hierarchical group structure:

```
file.h5
├── /Channel01/
│   └── /Zposition001/
│       ├── /Data0001/          # First dataset
│       │   └── Data0001        # 3D array [Frames × Y × X]
│       ├── /Data0002/          # Second dataset
│       │   └── Data0002
│       └── ...
├── /Metadata/                  # SEQv2 format only
│   └── Attributes
└── /Calibration/               # SEQv2 format only
    └── Calibration data
```

#### Hierarchical Organization

**Level 1: Channel**
- Group name: `/Channel01`, `/Channel02`, etc.
- Contains different imaging channels (wavelengths)
- Most experiments use `/Channel01` only

**Level 2: Z-Position**
- Group name: `/Zposition001`, `/Zposition002`, etc.
- Contains different focal planes
- Most 2D experiments use `/Zposition001` only

**Level 3: Datasets**
- Group names: `/Data0001`, `/Data0002`, etc.
- Each group contains one dataset (acquisition sequence)
- Multiple datasets = multiple acquisitions from same sample

**Level 4: Data Arrays**
- Dataset name matches group name: `Data0001`
- 3D array: `[NFrames × NYPixels × NXPixels]`
- Data type: typically `uint16` (camera ADU values)
- Order: Time (frames) × Y (rows) × X (columns)

#### HDF5 File Versions

smite recognizes three HDF5 structure versions:

**SEQv0** (Legacy)
```
/Channel01/Zposition001/
    Datasets directly under Zposition001
```

**SEQv1** (Standard)
```
/Channel01/Zposition001/Data0001/Data0001
    Each dataset in its own group
```

**SEQv2** (Extended)
```
/Channel01/Zposition001/Data0001/Data0001
/Metadata/
/Calibration/
    Includes metadata and calibration groups
```

#### Group Attributes

Each group can have HDF5 attributes containing metadata:

```matlab
% Example attributes for /Channel01/Zposition001/Data0001/
Attributes:
    FrameRate: 100.0           % Hz
    ExposureTime: 0.01         % seconds
    LaserPower: 50.0           % mW
    Wavelength: 647            % nm
    PixelSize: 0.1             % micrometers
```

#### Loading HDF5 Files

**Using smite's LoadData:**
```matlab
% Method 1: Via SMF structure
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = 'experiment.h5';

LD = smi_core.LoadData();
[~, Data, SMF] = LD.loadRawData(SMF, 1);  % Load dataset 1
```

**Using MATLAB built-in functions:**
```matlab
% View file structure
Info = h5info('experiment.h5');

% Count datasets
NDatasets = numel(Info.Groups.Groups.Datasets);

% Load specific dataset
Data = h5read('experiment.h5', '/Channel01/Zposition001/Data0001');
% Data is [NFrames × NYPixels × NXPixels] array
```

**Reading with smite's readH5File:**
```matlab
% Extract entire file structure
H5Structure = smi_core.LoadData.readH5File('experiment.h5');

% Extract specific group
H5Group = smi_core.LoadData.readH5File('experiment.h5', 'Data0001');

% Access data
ImageStack = H5Structure.Data.Data0001;
```

#### HDF5 Best Practices

**Advantages:**
- Efficient storage with compression
- Standardized scientific data format
- Self-documenting with attributes
- Platform-independent
- Handles large datasets well
- Can store multiple datasets in one file

**Naming Conventions:**
- Use descriptive filenames: `CellType_Condition_Date.h5`
- Keep dataset numbering sequential
- Store related acquisitions in same file

**Performance Tips:**
- HDF5 allows chunked reading (partial loading)
- Use appropriate compression (balance size vs. speed)
- For very large files, consider loading frames in batches

**Common Issues:**
- **File corruption**: Use HDF5 tools (`h5ls`, `h5dump`) to verify
- **Version mismatch**: Check if file is SEQv0/v1/v2
- **Missing datasets**: Use `h5info()` to list available datasets

---

### MAT Format (.mat)

MATLAB's native format, convenient but less efficient for large datasets.

#### File Structure

MAT files contain MATLAB variables saved with `save()`:

```matlab
% A typical smite .mat data file contains:
sequence    % [NFrames × NYPixels × NXPixels] array
            % Variable name is configurable via SMF.Data.DataVariable
```

#### Data Variable

The image stack variable name defaults to `'sequence'` but is configurable:

```matlab
% Set custom variable name in SMF
SMF.Data.DataVariable = 'imageStack';

% Or for multiple files with different names
SMF.Data.DataVariable = {'sequence', 'frames', 'data'};
```

#### Loading MAT Files

**Using smite's LoadData:**
```matlab
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = '/path/to/data';
SMF.Data.FileName = 'experiment.mat';
SMF.Data.DataVariable = 'sequence';  % Variable name in .mat file

LD = smi_core.LoadData();
[~, Data, SMF] = LD.loadRawData(SMF, 1);
```

**Using MATLAB built-in functions:**
```matlab
% Load entire file
load('experiment.mat');  % Loads 'sequence' variable

% Load specific variable
Data = load('experiment.mat', 'sequence');
ImageStack = Data.sequence;

% Check what's in the file
FileContents = who('-file', 'experiment.mat');
```

#### Multiple MAT Files

For multiple acquisitions, use cell arrays:

```matlab
% Multiple files in SMF
SMF.Data.FileName = {'data1.mat', 'data2.mat', 'data3.mat'};
SMF.Data.FileDir = '/path/to/data';

% All files use same variable name
SMF.Data.DataVariable = 'sequence';

% Load specific file (1, 2, or 3)
[~, Data, SMF] = LD.loadRawData(SMF, 2);  % Loads data2.mat
```

#### MAT File Array Format

Image data should be a 3D numeric array:

```matlab
% Correct format
sequence = zeros(NFrames, NYPixels, NXPixels, 'uint16');
sequence(frame, y, x) = value;

% Dimensions
size(sequence)
% [2000, 256, 256]  = 2000 frames of 256×256 images

% Data type
class(sequence)
% 'uint16' (camera ADU values, typically 0-65535)
```

#### MAT Format Best Practices

**Advantages:**
- Native MATLAB format (no conversion needed)
- Simple structure (just variables)
- Easy to inspect and modify
- Good for small datasets

**Disadvantages:**
- Less efficient storage than HDF5
- Larger file sizes
- One dataset per file typically
- No built-in compression (use v7.3 format)

**When to Use:**
- Small datasets (< 1 GB)
- Quick tests and simulations
- Legacy compatibility
- Simple workflows

**File Size Management:**
```matlab
% Use compressed format for large files
save('data.mat', 'sequence', '-v7.3');  % Uses HDF5 internally

% Check file size before saving
DataSize = whos('sequence');
FileSizeMB = DataSize.bytes / 1024^2;
```

---

## Camera Calibration Files

Calibration files contain pixel-wise camera properties essential for accurate localization. **Highly recommended for sCMOS cameras**, optional for EMCCD.

### Purpose

Camera calibration accounts for:
- **Gain**: Electron-to-ADU conversion per pixel
- **Offset**: Baseline signal per pixel
- **Variance**: Read noise per pixel

Without calibration, localization precision and photon counts will be inaccurate.

### File Structure

Calibration files are `.mat` files containing a `Params` structure:

```matlab
>> load('GainCalibration.mat')
>> Params

Params =
  struct with fields:
           MeanLevel: [256×256×20 single]  % Mean at each light level
            VarLevel: [256×256×20 single]  % Variance at each light level
      DarkImagesMean: [256×256×5 single]   % Mean dark frames
       DarkImagesVar: [256×256×5 single]   % Variance dark frames
 DarkImagesExpTime: [0.01 0.13 0.26 0.38 0.50]  % Exposure times
DarkImagesSequenceLength: 1000
              CCDVar: [256×256 single]      % Measured variance
             CCDGain: [256×256 single]      % Deprecated, use Gain
           CCDOffset: [256×256 single]      % Measured offset
           CameraObj: [1×1 struct]          % Camera settings
      LampPowerRange: [1×20 single]         % Light levels used
                Gain: [256×256 single]      % *** PRIMARY: Gain map ***
```

### Required Fields

The three essential fields extracted by smite:

```matlab
Params.Gain         % [NYPixels × NXPixels single] - ADU per electron
Params.CCDOffset    % [NYPixels × NXPixels single] - Baseline ADU
Params.CCDVar       % [NYPixels × NXPixels single] - Read noise variance
```

These must match the ROI size of your data or be the full chip size.

### Generating Calibration Files

Complete calibration procedure using MIC (MATLAB Instrument Control):

```matlab
% 1. Connect to camera
CameraSCMOS = MIC_HamamatsuCamera();
CameraSCMOS.ReturnType = 'matlab';
CameraSCMOS.DefectCorrection = 1;  % No correction
CameraSCMOS.ScanMode = 1;          % Slow scan
CameraSCMOS.ExpTime_Sequence = 0.01;
CameraSCMOS.SequenceLength = 1000;
CameraSCMOS.ROI = [897, 1152, 897, 1152];  % [XStart, XEnd, YStart, YEnd]

% 2. Connect to light source
Lamp660 = MIC_ThorlabsLED('Dev1', 'ao0');

% 3. Define light levels (target ~600 max ADU)
LampPowerRange = linspace(0, 3.6, 20);

% 4. Collect data at each light level
MeanLevel = [];
VarLevel = [];
CameraSCMOS.AcquisitionType = 'sequence';
CameraSCMOS.setup_acquisition();

for ii = 1:numel(LampPowerRange)
    fprintf('Light level %i of %i\n', ii, numel(LampPowerRange));
    Lamp660.setPower(LampPowerRange(ii));
    pause(1);  % Allow stabilization
    CameraSCMOS.start_sequence();
    MeanLevel = cat(3, MeanLevel, mean(CameraSCMOS.Data, 3));
    VarLevel = cat(3, VarLevel, var(single(CameraSCMOS.Data), [], 3));
end
Lamp660.setPower(0);

% 5. Perform least-squares fit (variance vs. mean)
Beta = NaN(size(MeanLevel, 1), size(MeanLevel, 2), 2);
for ii = 1:size(MeanLevel, 1)
    for jj = 1:size(MeanLevel, 2)
        % Beta(:,:,1) = offset, Beta(:,:,2) = gain
        Beta(ii, jj, 1:2) = smi_stat.leastSquaresFit(...
            squeeze(MeanLevel(ii, jj, :)), ...
            squeeze(VarLevel(ii, jj, :)), ...
            1 ./ squeeze(VarLevel(ii, jj, :)));
    end
end

% 6. Store results
Params.MeanLevel = single(MeanLevel);
Params.VarLevel = single(VarLevel);
Params.CCDVar = single(VarLevel(:,:,1));     % Variance at zero light
Params.Gain = single(Beta(:,:,2));           % Slope of variance vs mean
Params.CCDOffset = single(MeanLevel(:,:,1)); % Mean at zero light
Params.CameraObj.ROI = CameraSCMOS.ROI;
Params.LampPowerRange = single(LampPowerRange);

% 7. Save calibration file
SaveDir = '/path/to/calibrations';
FileName = fullfile(SaveDir, ...
    ['GainCalibration-', smi_helpers.genTimeString()]);
save(FileName, 'Params', '-v7.3');
```

### Using Calibration Files

**Set path in SMF:**
```matlab
SMF.Data.CameraType = 'SCMOS';
SMF.Data.CalibrationFilePath = '/path/to/GainCalibration.mat';
```

**Automatic loading:**
When `CalibrationFilePath` is set, smite automatically:
1. Loads `Params.Gain`, `Params.CCDOffset`, `Params.CCDVar`
2. Extracts ROI-matched region if needed
3. Applies calibration during localization

**Manual inspection:**
```matlab
load('/path/to/GainCalibration.mat', 'Params');

% Visualize gain map
figure; imagesc(Params.Gain);
colorbar; title('Camera Gain (ADU/e-)');

% Check statistics
fprintf('Mean Gain: %.2f ADU/e-\n', mean(Params.Gain(:)));
fprintf('Gain Std: %.2f ADU/e-\n', std(Params.Gain(:)));
fprintf('Mean Offset: %.1f ADU\n', mean(Params.CCDOffset(:)));
```

### Calibration Best Practices

**When to Calibrate:**
- Before each new experiment series
- After camera maintenance
- If temperature changes significantly
- Every 3-6 months for stability

**Quality Checks:**
```matlab
% Gain should be relatively uniform
histogram(Params.Gain(:));
% Typical sCMOS: 0.5-5 ADU/e-, narrow distribution

% Offset should be stable
histogram(Params.CCDOffset(:));
% Should be tight distribution around baseline

% Check for bad pixels
BadPixels = (Params.Gain < 0.1) | (Params.Gain > 10);
fprintf('Bad pixels: %d (%.2f%%)\n', sum(BadPixels(:)), ...
    100*sum(BadPixels(:))/numel(BadPixels));
```

**Storage and Organization:**
```matlab
% Naming convention
% GainCalibration-CameraModel-ROI-Date.mat
% Examples:
'GainCalibration-OrcaFlash-256x256-20250110.mat'
'GainCalibration-HamamatsuC13440-FullChip-20250110.mat'
```

---

## Channel Registration Files

For multi-channel imaging, registration files correct spatial misalignment between color channels.

### Purpose

Different wavelength channels often have:
- **Translation**: XY offset between channels
- **Rotation**: Angular misalignment
- **Scaling**: Different magnifications
- **Distortion**: Non-linear warping

Registration files store transformation parameters to align channels.

### File Structure

Registration files are `.mat` files containing transformation information:

```matlab
>> load('ChannelRegistration.mat')

% Typical contents (format varies by registration method)
RegistrationStruct =
  struct with fields:
    TransformType: 'affine'           % or 'polynomial', 'projective'
         Transform: [3×3 double]      % Transformation matrix
      FixedChannel: 1                 % Reference channel
     MovingChannel: 2                 % Channel to transform
          PixelSize: 0.1               % Micrometers
     ControlPoints: [N×4 double]      % Optional: [x1 y1 x2 y2] pairs
```

### Common Registration Formats

**Affine Transformation:**
```matlab
% 3×3 transformation matrix
T = [cos(θ)  -sin(θ)  tx;
     sin(θ)   cos(θ)  ty;
     0        0       1];

% Apply transformation
[x_reg, y_reg] = transformPointsForward(T, x, y);
```

**Polynomial Transformation:**
```matlab
% Higher-order correction for distortion
PolyParams = struct();
PolyParams.XCoeffs = [a0, a1, a2, ...];  % X = f(x,y)
PolyParams.YCoeffs = [b0, b1, b2, ...];  % Y = g(x,y)
```

### Using Registration Files

**Set path in SMF:**
```matlab
SMF.Data.RegistrationFilePath = '/path/to/ChannelRegistration.mat';
```

**Manual application:**
```matlab
% Load registration
load('ChannelRegistration.mat', 'RegistrationStruct');

% Transform coordinates from channel 2 to channel 1
SMD_Ch2_Registered = SMD_Ch2;  % Copy structure
if strcmp(RegistrationStruct.TransformType, 'affine')
    T = affine2d(RegistrationStruct.Transform);
    [X_reg, Y_reg] = transformPointsForward(T, SMD_Ch2.X, SMD_Ch2.Y);
    SMD_Ch2_Registered.X = X_reg;
    SMD_Ch2_Registered.Y = Y_reg;
end
```

### Creating Registration Files

**Using fiducial markers (beads):**
```matlab
% 1. Acquire images of fluorescent beads in both channels
% 2. Localize beads in each channel
% 3. Match bead pairs between channels
% 4. Compute transformation

% Example using MATLAB's fitgeotrans
ControlPoints_Ch1 = [x1, y1];  % Bead positions in channel 1
ControlPoints_Ch2 = [x2, y2];  % Same beads in channel 2

% Fit transformation (affine recommended)
tform = fitgeotrans(ControlPoints_Ch2, ControlPoints_Ch1, 'affine');

% Save registration
RegistrationStruct.TransformType = 'affine';
RegistrationStruct.Transform = tform.T;
RegistrationStruct.FixedChannel = 1;
RegistrationStruct.MovingChannel = 2;
RegistrationStruct.ControlPoints = [ControlPoints_Ch2, ControlPoints_Ch1];
save('ChannelRegistration.mat', 'RegistrationStruct');
```

### Registration Best Practices

**When to Register:**
- Before multi-channel experiments
- After optical alignment changes
- Monthly for stability checks

**Quality Assessment:**
```matlab
% Compute registration error
[X_pred, Y_pred] = transformPointsForward(tform, ControlPoints_Ch2(:,1), ...
    ControlPoints_Ch2(:,2));
Errors = sqrt((X_pred - ControlPoints_Ch1(:,1)).^2 + ...
              (Y_pred - ControlPoints_Ch1(:,2)).^2);

fprintf('Mean registration error: %.2f pixels\n', mean(Errors));
fprintf('RMS registration error: %.2f pixels\n', rms(Errors));
% Target: < 0.5 pixels for good registration
```

---

## Results Files

Analysis results are saved as `.mat` files containing SMD structures and metadata.

### Default Results Location

```matlab
% Set in SMF.Data.ResultsDir
SMF.Data.FileDir = '/path/to/data';
SMF.Data.ResultsDir = fullfile(SMF.Data.FileDir, 'Results');
% Results saved to: /path/to/data/Results/
```

### Results File Contents

A typical results file contains:

```matlab
>> load('Data0001_Results.mat')

% Variables in file:
SMD             % Single Molecule Data structure (main results)
SMF             % Single Molecule Fitting structure (parameters used)
SMR             % Single Molecule Results (additional analysis)
FrameConnect    % Frame connection results (if applicable)
```

### SMD Structure (Main Results)

Complete localization results:

```matlab
SMD =
  struct with fields:
          X: [N×1 single]         % X positions (pixels)
          Y: [N×1 single]         % Y positions (pixels)
          Z: [N×1 single]         % Z positions (micrometers, if 3D)
    Photons: [N×1 single]         % Total photons per localization
         Bg: [N×1 single]         % Background (photons/pixel)
        X_SE: [N×1 single]         % X standard error (pixels)
        Y_SE: [N×1 single]         % Y standard error (pixels)
        Z_SE: [N×1 single]         % Z standard error (micrometers)
  Photons_SE: [N×1 single]         % Photon uncertainty
       LogL: [N×1 single]         % Log-likelihood of fit
    FrameNum: [N×1 int32]         % Frame number of localization
  DatasetNum: [N×1 int16]         % Dataset number
      ConnectID: [N×1 int32]      % Frame-connection cluster ID
    NFrameConns: [N×1 int16]      % Number of frame connections
      ThreshFlag: [N×1 logical]   % Passed thresholding?
```

### Reading Results

**Load and inspect:**
```matlab
% Load results file
load('/path/to/Results/Data0001_Results.mat', 'SMD', 'SMF');

% How many localizations?
NLocs = length(SMD.X);
fprintf('Found %d localizations\n', NLocs);

% Basic statistics
fprintf('Mean photons: %.0f\n', mean(SMD.Photons));
fprintf('Mean precision: %.2f pixels\n', mean(SMD.X_SE));

% Plot results
figure;
plot(SMD.X, SMD.Y, '.', 'MarkerSize', 1);
axis equal; xlabel('X (pixels)'); ylabel('Y (pixels)');
title(sprintf('%d Localizations', NLocs));
```

**Filter results:**
```matlab
% Select high-quality localizations
GoodLocs = (SMD.Photons > 500) & ...
           (SMD.X_SE < 0.15) & ...
           (SMD.Y_SE < 0.15);

% Create filtered SMD
SMD_Filtered = smi_core.SelectLocalizations.selectLocalizations(SMD, GoodLocs);

fprintf('Kept %d / %d localizations\n', ...
    sum(GoodLocs), length(GoodLocs));
```

### Results File Naming

Default naming convention:

```
Data0001_Results.mat            % Basic results
Data0001_AnalysisID_Results.mat % With AnalysisID
Data0001_AllFrames_Results.mat  % Example with custom ID
```

Controlled by SMF parameters:
```matlab
SMF.Data.AnalysisID = 'HighPhotons';  % Optional identifier
% Results saved as: Data0001_HighPhotons_Results.mat
```

### Batch Results Organization

For batch processing with multiple datasets:

```
ExperimentFolder/
├── Cell01/
│   ├── Label01/
│   │   ├── Data0001.h5
│   │   └── Results/
│   │       ├── Data0001_Results.mat
│   │       └── Data0001_SR.mat         % SR image
│   └── Label02/
│       └── ...
├── Cell02/
│   └── ...
└── CombinedResults/                     % Optional
    └── AllCells_Combined.mat
```

### Exporting Results

**To CSV:**
```matlab
% Create table from SMD
T = table(SMD.X, SMD.Y, SMD.Photons, SMD.FrameNum, ...
    'VariableNames', {'X_pixels', 'Y_pixels', 'Photons', 'Frame'});

% Write to CSV
writetable(T, 'localizations.csv');
```

**To other formats:**
```matlab
% Export coordinates only
Coords = [SMD.X, SMD.Y, SMD.Photons];
save('coordinates.txt', 'Coords', '-ascii');

% Export to struct for JSON
ExportData.x = SMD.X;
ExportData.y = SMD.Y;
ExportData.photons = SMD.Photons;
ExportData.frame = SMD.FrameNum;
jsonStr = jsonencode(ExportData);
fid = fopen('results.json', 'w');
fprintf(fid, '%s', jsonStr);
fclose(fid);
```

---

## Parameter Files

SMF structures can be saved and loaded for reproducibility.

### Saving SMF Parameters

```matlab
% After configuring SMF
SMF = smi_core.SingleMoleculeFitting();
% ... set parameters ...

% Save for reuse
save('MyAnalysisParams.mat', 'SMF');
```

### Loading SMF Parameters

```matlab
% Load saved parameters
load('MyAnalysisParams.mat', 'SMF');

% Update data-specific fields
SMF.Data.FileName = 'NewData.h5';
SMF.Data.FileDir = '/new/path';

% Run analysis with loaded parameters
SMLMObj = smi.SMLM(SMF);
```

### Parameter File Best Practices

**Naming conventions:**
```
SMF_Standard2D.mat              % Standard 2D SMLM
SMF_Tracking_FastDiffusion.mat  % SPT for fast particles
SMF_3D_Astigmatism.mat          % 3D with astigmatism
```

**Version control:**
```matlab
% Include metadata
SMF.Data.AnalysisID = 'StandardAnalysis_v2.1';
SMF.Metadata.Description = 'Standard SMLM analysis for fixed cells';
SMF.Metadata.Created = datetime('now');
save('SMF_StandardAnalysis_v2.1.mat', 'SMF');
```

---

## Data ROI Specification

Control which region and frames to analyze using `SMF.Data.DataROI`.

### ROI Format

```matlab
% Full ROI specification (6 elements)
SMF.Data.DataROI = [YStart, XStart, YSize, XSize, ZStart, FrameStart];

% Common examples:
SMF.Data.DataROI = [1, 1, 256, 256, 1, 1];      % Full 256×256 image, all frames
SMF.Data.DataROI = [50, 50, 100, 100, 1, 1];    % 100×100 crop starting at (50,50)
SMF.Data.DataROI = [1, 1, 256, 256, 1, 100];    % Skip first 99 frames
```

### Automatic ROI

If `DataROI` is empty or invalid, smite sets it automatically:

```matlab
SMF.Data.DataROI = [];  % Auto-detect from data
% Automatically set to: [1, 1, NYPixels, NXPixels, 1, 1]
```

### Partial Loading for Large Files

```matlab
% Analyze only center 128×128 region
ImageSize = [256, 256];
CropSize = 128;
Start = (ImageSize - CropSize) / 2 + 1;
SMF.Data.DataROI = [Start(1), Start(2), CropSize, CropSize, 1, 1];

% Or analyze only frames 100-500
SMF.Data.DataROI = [1, 1, 256, 256, 1, 100];
% Then limit frames in analysis separately
```

---

## Coordinate Systems and Units

Understanding coordinate conventions is crucial for multi-software workflows.

### Pixel Coordinates

smite follows MATLAB image coordinate conventions:

```matlab
% Image indexing: Image(row, column) = Image(y, x)
% Coordinate (1, 1) = center of top-left pixel
% Coordinate (2, 1) = one pixel down, same column

% Example: 256×256 image
%   X range: 1 to 256 (left to right)
%   Y range: 1 to 256 (top to bottom)
%   Origin: top-left corner at (1, 1)
```

### Physical Units

Convert between pixels and physical units:

```matlab
% SMF.Data.PixelSize is in micrometers
PixelSize = 0.1;  % μm/pixel (100 nm)

% Pixel to micrometers
X_um = SMD.X * PixelSize;
Y_um = SMD.Y * PixelSize;

% Micrometers to pixels
X_pixels = X_um / PixelSize;
```

### Field of View

Calculate imaging area:

```matlab
% Image size in pixels
[NYPixels, NXPixels, NFrames] = size(Data);

% Physical field of view
FOV_X_um = NXPixels * SMF.Data.PixelSize;
FOV_Y_um = NYPixels * SMF.Data.PixelSize;

fprintf('Field of view: %.1f × %.1f μm\n', FOV_X_um, FOV_Y_um);
% Example: "Field of view: 25.6 × 25.6 μm"
```

---

## Troubleshooting File Format Issues

### HDF5 Files

**Problem: "No such file or directory"**
```matlab
% Check file exists
exist('/path/to/file.h5', 'file')  % Should return 2

% Check absolute vs relative paths
SMF.Data.FileDir = pwd;  % Use absolute paths
```

**Problem: "Dataset not found"**
```matlab
% List available datasets
Info = h5info('file.h5');
GroupNames = {Info.Groups.Groups.Groups.Name}

% Check version
[~, Version] = smi_core.LoadData.seqH5Data('file.h5');
fprintf('HDF5 version: %s\n', Version);
```

**Problem: "Unable to read data"**
```matlab
% Verify HDF5 integrity
h5disp('file.h5');  % Display entire structure

% Try reading with h5read directly
TestData = h5read('file.h5', '/Channel01/Zposition001/Data0001');
```

### MAT Files

**Problem: "Variable not found"**
```matlab
% Check variable names in file
who('-file', 'data.mat')

% Load and inspect
FileContents = load('data.mat');
fieldnames(FileContents)

% Set correct variable name
SMF.Data.DataVariable = 'ActualVariableName';
```

**Problem: "Out of memory"**
```matlab
% Check file size
FileInfo = dir('data.mat');
FileSizeGB = FileInfo.bytes / 1024^3;
fprintf('File size: %.2f GB\n', FileSizeGB);

% Use matfile() for partial loading
m = matfile('data.mat');
Data = m.sequence(1:100, :, :);  % Load first 100 frames only
```

### Calibration Files

**Problem: "ROI size mismatch"**
```matlab
% Check calibration ROI
load('Calibration.mat', 'Params');
CalibSize = size(Params.Gain);

% Check data ROI
DataSize = size(ImageStack, [1,2]);

% They must match!
if ~isequal(CalibSize, DataSize)
    warning('Size mismatch: Calib %dx%d vs Data %dx%d', ...
        CalibSize(1), CalibSize(2), DataSize(1), DataSize(2));
end
```

**Problem: "Missing required fields"**
```matlab
% Verify required fields exist
RequiredFields = {'Gain', 'CCDOffset', 'CCDVar'};
MissingFields = setdiff(RequiredFields, fieldnames(Params));
if ~isempty(MissingFields)
    error('Missing calibration fields: %s', strjoin(MissingFields, ', '));
end
```

---

## Quick Reference Tables

### File Extensions

| Extension | Type | Required | MATLAB Function |
|-----------|------|----------|-----------------|
| `.h5` | HDF5 raw data | Yes (data) | `h5read()`, `h5info()` |
| `.mat` | MAT raw data | Yes (data) | `load()`, `matfile()` |
| `.mat` | Calibration | Recommended | `load()` |
| `.mat` | Registration | For multi-channel | `load()` |
| `.mat` | Results | Output | `load()` |

### SMF Data Fields

| Field | Type | Example | Purpose |
|-------|------|---------|---------|
| `FileName` | cell | `{'data.h5'}` | Data filename(s) |
| `FileDir` | char | `'/path/to/data'` | Data directory |
| `ResultsDir` | char | `'/path/to/data/Results'` | Output directory |
| `FileType` | char | `'.h5'` | File extension |
| `DataVariable` | char | `'sequence'` | MAT variable name |
| `CalibrationFilePath` | char | `'/path/to/cal.mat'` | Camera calibration |
| `RegistrationFilePath` | char | `'/path/to/reg.mat'` | Channel registration |
| `DataROI` | 1×6 | `[1,1,256,256,1,1]` | Region to analyze |
| `PixelSize` | double | `0.1` | Micrometers per pixel |
| `FrameRate` | double | `100` | Frames per second |

### Data Array Dimensions

| Format | Dimensions | Indexing | Example |
|--------|------------|----------|---------|
| HDF5 | `[Frames, Y, X]` | `Data(frame, row, col)` | `[2000, 256, 256]` |
| MAT | `[Frames, Y, X]` | `Data(frame, row, col)` | `[2000, 256, 256]` |
| SMD.X | `[N, 1]` | `SMD.X(localization)` | `[50000, 1]` |

---

## Summary

smite uses well-defined file formats for each stage of analysis:

1. **Raw Data** (`.h5` or `.mat`): Image stacks from microscope
2. **Calibration** (`.mat`): Camera gain/offset/noise maps
3. **Registration** (`.mat`): Multi-channel alignment transformations
4. **Results** (`.mat`): SMD localization data and SMF parameters
5. **Parameters** (`.mat`): Reusable SMF configurations

**Key Takeaways:**

- HDF5 (`.h5`) is recommended for new data acquisition
- MAT format works well for smaller datasets
- Calibration files are essential for sCMOS cameras
- Registration enables multi-channel analysis
- Results files contain complete SMD + SMF for reproducibility
- All file paths should be absolute for reliability

**Next Steps:**

- [How to Load Data](../how-to/load-data.md) - Practical loading examples
- [SMF Structure](../core-concepts/smf-structure.md) - Complete parameter reference
- [SMD Structure](../core-concepts/smd-structure.md) - Results data format
- [First Analysis](../getting-started/first-analysis.md) - See formats in action

---

*This reference covers all standard file formats in smite. For custom formats or specific integration needs, consult the source code in `MATLAB/+smi_core/@LoadData/` or open a GitHub issue.*
