---
title: "API Reference: +smi_helpers Namespace"
category: "api-reference"
level: "beginner"
tags: ["api", "smi_helpers", "utilities", "filters", "helpers", "array-tools"]
prerequisites: ["../core-concepts/smd-structure.md"]
related: ["../workflows/smlm-analysis.md", "smi-core.md"]
summary: "Complete API reference for the +smi_helpers namespace covering utility functions, filtering tools, and helper classes for data manipulation"
estimated_time: "20 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# API Reference: +smi_helpers Namespace

## Purpose

The +smi_helpers namespace provides utility functions and helper classes that support smite analyses. These tools handle common tasks like filtering SMD structures, manipulating arrays, subdividing images and data, working with file paths, and managing structures. While not required for basic analyses, these utilities simplify advanced workflows and custom processing pipelines.

## Prerequisites

- Basic understanding of [SMD structure](../core-concepts/smd-structure.md)
- Familiarity with MATLAB structures and arrays
- Knowledge of image and data processing concepts

## Overview

The +smi_helpers namespace provides:

**Filtering Tools:**
- `Filters` - Quality filtering for BaGoL and SMLM analyses

**Data Subdivision:**
- `subdivideImage` - Split images into spatial sub-regions
- `subdivideSMD` - Split localizations into spatial sub-regions

**Array and Data Manipulation:**
- `arrayMUX` - Generic array multiplexer
- `compressToRange` - Compress integers to continuous range
- `padStruct` - Merge structures with defaults
- `removeBorder` - Remove image borders

**File and Path Utilities:**
- `selectFiles` - GUI file selection
- `getFileNames` - List files matching pattern
- `getDirectoryNames` - List directories matching pattern
- `filenameWindows2Unix` - Cross-platform path conversion

**Analysis Helpers:**
- `triangle_threshold` - Automatic threshold detection
- `findStartEndInds` - Find event boundaries in boolean arrays
- `genTimeString` - Generate timestamp strings

**System Utilities:**
- `versionSMITE` - Report smite version
- `requiredToolboxes` - List required MATLAB toolboxes
- `mkSMITETmpDir` - Create temporary directories for testing

## Class and Function Reference

### Filtering Tools

#### Filters

**Purpose:** Quality filtering of SMD structures, especially for BaGoL analyses.

**Key Concept:** BaGoL (Bayesian Grouping of Localizations) requires carefully filtered data. The Filters class provides a standard pipeline: remove negative coordinates, filter by intensity, inflate standard errors, remove short frame connections, and apply nearest-neighbor filtering.

**Class Definition:**
```matlab
classdef Filters < handle
```

**Static Methods:**

`filterNonNeg(SMD, Verbose)` - Remove localizations with negative coordinates
```matlab
SMD = smi_helpers.Filters.filterNonNeg(SMD, true);
% Removes any localizations where X or Y < 0
```

`filterIntensity(SMD, Verbose, MeanMultiplier)` - Filter by photon intensity
```matlab
% Remove localizations with photons > 3 * mean
SMD = smi_helpers.Filters.filterIntensity(SMD, true, 3);
% Removes bright outliers that may be artifacts
```

`inflateSE(SMD, Verbose, SEAdjust)` - Inflate standard errors
```matlab
% Multiply all standard errors by 1.5
SMD = smi_helpers.Filters.inflateSE(SMD, true, 1.5);
% Conservative SE estimates for robust clustering
```

`filterFC(SMD, Verbose, nFC)` - Filter by frame connection count
```matlab
% Keep only localizations representing ≥2 frame connections
SMD = smi_helpers.Filters.filterFC(SMD, true, 2);
% Removes single-frame blinks
```

`filterNN(SMD, Verbose, n_NN, MedianMultiplier)` - Nearest neighbor density filter
```matlab
% Require 5 neighbors within 3 * median precision
SMD = smi_helpers.Filters.filterNN(SMD, true, 5, 3);
% Removes isolated localizations (do NOT use on dSTORM!)
```

**Standard BaGoL Filter Pipeline:**

```matlab
% Load localized data
load('Results.mat', 'SMD');

fprintf('Starting with %d localizations\n', length(SMD.X));

% Step 1: Remove negative coordinates
SMD = smi_helpers.Filters.filterNonNeg(SMD, true);

% Step 2: Remove intensity outliers
SMD = smi_helpers.Filters.filterIntensity(SMD, true, 3);

% Step 3: Inflate standard errors (conservative)
SMD = smi_helpers.Filters.inflateSE(SMD, true, 1.5);

% Step 4: Require at least 2 frame connections
SMD = smi_helpers.Filters.filterFC(SMD, true, 2);

% Step 5: Nearest neighbor filter (not for dSTORM!)
% Only for PAINT or photoactivation-based techniques
SMD = smi_helpers.Filters.filterNN(SMD, true, 5, 3);

fprintf('After filtering: %d localizations\n', length(SMD.X));

% Now ready for BaGoL
```

**Usage Notes:**

- **filterNN warning:** Do NOT use nearest-neighbor filtering on dSTORM data where emitters can be densely packed. This filter assumes sparse labeling.
- Verbose flag controls console output (true = show progress)
- All filters operate by calling `smi_core.SingleMoleculeData.isolateSubSMD` internally
- Filters can be applied independently or in sequence

**Selective Filtering:**

```matlab
% More aggressive intensity filter
SMD = smi_helpers.Filters.filterIntensity(SMD, false, 2);
% Removes everything above 2x mean (stricter)

% Require more neighbors for dense structures
SMD = smi_helpers.Filters.filterNN(SMD, true, 10, 2);
% Requires 10 neighbors within 2x median precision

% Conservative frame connection requirement
SMD = smi_helpers.Filters.filterFC(SMD, true, 5);
% Only keep emitters seen in ≥5 frames
```

**See Also:**
- BaGoL workflow documentation
- [SMD Structure](../core-concepts/smd-structure.md)

---

### Data Subdivision Tools

#### subdivideImage

**Purpose:** Divides images into spatial sub-regions for parallel processing or local analysis.

**Key Concept:** Large images can be processed more efficiently by dividing into smaller tiles. Edge tiles are automatically sized to fit remaining space.

**Function Signature:**
```matlab
[DividedImages, ImageROIs] = smi_helpers.subdivideImage(Image, SubROISize)
```

**Inputs:**
- `Image`: Image to subdivide (M×N or M×N×F array)
- `SubROISize`: Nominal sub-ROI size in pixels ([Height, Width])

**Outputs:**
- `DividedImages`: Cell array of sub-images
- `ImageROIs`: ROI coordinates (NROIs×4: [YStart, XStart, YEnd, XEnd])

**Usage Examples:**

Basic image subdivision:
```matlab
% Load image
Image = imread('large_image.tif');
fprintf('Original size: %d × %d pixels\n', size(Image));

% Divide into 256×256 pixel tiles
SubROISize = [256, 256];
[Tiles, ROIs] = smi_helpers.subdivideImage(Image, SubROISize);

fprintf('Created %d tiles\n', length(Tiles));

% Process each tile
for ii = 1:length(Tiles)
    % Apply processing to Tiles{ii}
    ProcessedTile = myProcessing(Tiles{ii});
end
```

Processing image stack:
```matlab
% 3D image stack (Y × X × Frames)
DataStack = rand(512, 512, 1000);

% Subdivide spatially (preserves all frames in each tile)
[SubStacks, ROIs] = smi_helpers.subdivideImage(DataStack, [128, 128]);

% Each SubStacks{ii} is 128×128×1000 (or smaller for edge tiles)
fprintf('Divided into %d spatial regions\n', length(SubStacks));

% Parallel processing
parfor ii = 1:length(SubStacks)
    % Localize each region independently
    SMDRegion(ii) = processRegion(SubStacks{ii}, ROIs(ii, :));
end
```

Handling edge tiles:
```matlab
% 500×500 image divided into 128×128 tiles
Image = rand(500, 500);
[Tiles, ROIs] = smi_helpers.subdivideImage(Image, [128, 128]);

% Check tile sizes
for ii = 1:length(Tiles)
    TileSize = size(Tiles{ii});
    fprintf('Tile %d: %d × %d pixels\n', ii, TileSize(1), TileSize(2));
end

% Tiles 1-12 are 128×128
% Edge tiles are smaller (e.g., 116×128, 128×116, 116×116)
```

**See Also:**
- `subdivideSMD` - Subdivide localization data
- Parallel processing with parfor

---

#### subdivideSMD

**Purpose:** Divides SMD localizations into spatial sub-regions matching image subdivision.

**Key Concept:** After subdividing an image for processing, subdivide the SMD results to maintain spatial correspondence. Localizations near tile boundaries may appear in multiple sub-SMDs.

**Function Signature:**
```matlab
[SMDSub, SMDROIs] = smi_helpers.subdivideSMD(SMD, SubROISize)
```

**Inputs:**
- `SMD`: Single Molecule Data structure
- `SubROISize`: Nominal sub-ROI size ([Height, Width] in pixels)

**Outputs:**
- `SMDSub`: Array of SMD structures (one per sub-region)
- `SMDROIs`: ROI coordinates (NROIs×4: [YStart, XStart, YEnd, XEnd])

**Usage Examples:**

Basic SMD subdivision:
```matlab
% Load results
load('Results.mat', 'SMD');
fprintf('Total localizations: %d\n', length(SMD.X));

% Divide into 100×100 pixel regions
[SMDSub, ROIs] = smi_helpers.subdivideSMD(SMD, [100, 100]);

fprintf('Created %d sub-regions\n', length(SMDSub));

% Analyze each region
for ii = 1:length(SMDSub)
    fprintf('Region %d: %d localizations\n', ii, length(SMDSub(ii).X));
end
```

Parallel region analysis:
```matlab
% Subdivide data
[SMDSub, ROIs] = smi_helpers.subdivideSMD(SMD, [50, 50]);

% Process regions in parallel
Density = zeros(length(SMDSub), 1);
parfor ii = 1:length(SMDSub)
    % Compute localization density for each region
    AreaUM2 = (50 * SMD.PixelSize)^2;  % Area in µm²
    Density(ii) = length(SMDSub(ii).X) / AreaUM2;
end

% Visualize density map
NRows = ceil(sqrt(length(SMDSub)));
NCols = ceil(length(SMDSub) / NRows);
DensityMap = reshape(Density, NCols, NRows)';

figure;
imagesc(DensityMap);
colorbar;
title('Localization Density (localizations/µm²)');
```

Matched image and SMD subdivision:
```matlab
% Subdivide image and localizations with same ROI size
SubSize = [64, 64];

% Image subdivision
[ImageTiles, ImageROIs] = smi_helpers.subdivideImage(Image, SubSize);

% SMD subdivision (same ROIs)
[SMDSub, SMDROIs] = smi_helpers.subdivideSMD(SMD, SubSize);

% Process corresponding tiles
for ii = 1:length(ImageTiles)
    % Overlay localizations on corresponding image tile
    figure;
    imagesc(ImageTiles{ii});
    hold on;
    plot(SMDSub(ii).X - SMDROIs(ii, 2) + 1, ...
         SMDSub(ii).Y - SMDROIs(ii, 1) + 1, 'r.', 'MarkerSize', 5);
    title(sprintf('Region %d', ii));
end
```

**Important Notes:**

- Localizations near tile boundaries may appear in multiple sub-SMDs (by design)
- ROI coordinates use MATLAB convention: (1,1) is top-left pixel center
- Default SubROISize uses full SMD extent ([SMD.YSize, SMD.XSize])

**See Also:**
- `subdivideImage` - Subdivide images
- [SMD Structure](../core-concepts/smd-structure.md)

---

### Array and Data Manipulation

#### arrayMUX

**Purpose:** Generic array multiplexer for selecting outputs based on integer index.

**Key Concept:** Simplifies conditional output selection in a compact, readable way. Commonly used for verbosity-dependent function arguments.

**Function Signature:**
```matlab
Output = smi_helpers.arrayMUX(OutputOptions, Select)
```

**Inputs:**
- `OutputOptions`: Array or cell array of options
- `Select`: Zero-based index (0 = first option, 1 = second, etc.)

**Output:**
- `Output`: Selected element from OutputOptions

**Usage Examples:**

Verbosity-based display control:
```matlab
% Set display option based on verbosity level
Verbose = 2;  % 0=none, 1=iter, 2=final

DisplayOption = smi_helpers.arrayMUX({'none', 'iter', 'final'}, Verbose);

% Use in optimization
Options = optimset('fmincon', 'Display', DisplayOption);
% Verbose=0 → 'none', Verbose=1 → 'iter', Verbose=2 → 'final'
```

Selecting plot markers:
```matlab
% Different markers for different conditions
ConditionType = 1;  % 0, 1, or 2

Marker = smi_helpers.arrayMUX({'o', 's', '^'}, ConditionType);
Color = smi_helpers.arrayMUX({'r', 'g', 'b'}, ConditionType);

plot(X, Y, [Color, Marker]);
```

Numeric arrays:
```matlab
% Select threshold based on data quality
Quality = 2;  % 0=low, 1=medium, 2=high

Threshold = smi_helpers.arrayMUX([50, 100, 200], Quality);
% Quality=0 → 50, Quality=1 → 100, Quality=2 → 200
```

**Convention:** Uses zero-based indexing to match multiplexer hardware conventions.

---

#### compressToRange

**Purpose:** Compresses an integer array to a continuous range with no missing values.

**Key Concept:** Useful for converting sparse or non-sequential IDs (like ConnectID) into dense 1:N arrays for indexing or grouping.

**Function Signature:**
```matlab
Range = smi_helpers.compressToRange(IntegerArray)
```

**Input:**
- `IntegerArray`: Array of integers (may have gaps)

**Output:**
- `Range`: Compressed array with values 1:NUnique (preserves original order)

**Usage Examples:**

Basic compression:
```matlab
% Sparse IDs with gaps
IDs = [2; 4; 6; 6; 7; 10];

% Compress to continuous range
Compressed = smi_helpers.compressToRange(IDs);
% Result: [1; 2; 3; 3; 4; 5]

% Same relative ordering, but no gaps
```

Converting ConnectID for indexing:
```matlab
% Frame-connected data with non-sequential IDs
load('Results.mat', 'SMD');

% ConnectID might be [5, 5, 7, 7, 7, 12, 12, ...]
% Compress for use as array indices
CompactID = smi_helpers.compressToRange(SMD.ConnectID);
% Now: [1, 1, 2, 2, 2, 3, 3, ...]

% Count localizations per emitter
NEmitters = max(CompactID);
LocsPerEmitter = histcounts(CompactID, NEmitters);

% Access data by compact ID
for ii = 1:NEmitters
    EmitterLocs = find(CompactID == ii);
    % Process localizations for emitter ii
end
```

Preserves order:
```matlab
% Non-sorted input
IDs = [1; 5; 8; 8; 11; 3];

Compressed = smi_helpers.compressToRange(IDs);
% Result: [1; 3; 4; 4; 5; 2]

% Order preserved: 1 → 1, 5 → 3, 8 → 4, 11 → 5, 3 → 2
```

**See Also:**
- Frame connection and ConnectID
- Grouping and clustering operations

---

#### padStruct

**Purpose:** Merges input structure with default structure, filling missing fields.

**Key Concept:** Ensures structure has all required fields without overwriting existing values. Essential for backward compatibility and optional parameters.

**Function Signature:**
```matlab
Struct = smi_helpers.padStruct(Struct, DefaultStruct)
```

**Inputs:**
- `Struct`: Input structure (may be incomplete)
- `DefaultStruct`: Structure with default field values

**Output:**
- `Struct`: Input padded with missing fields from DefaultStruct

**Usage Examples:**

Parameter structure with defaults:
```matlab
% Define defaults
Default.MaxIter = 100;
Default.Tolerance = 1e-6;
Default.Verbose = true;
Default.OutputDir = pwd();

% User provides partial configuration
UserParams.MaxIter = 50;
UserParams.Verbose = false;
% Missing: Tolerance, OutputDir

% Merge with defaults
Params = smi_helpers.padStruct(UserParams, Default);

% Result has all fields:
% Params.MaxIter = 50 (from UserParams)
% Params.Tolerance = 1e-6 (from Default)
% Params.Verbose = false (from UserParams)
% Params.OutputDir = pwd() (from Default)
```

Backward compatibility:
```matlab
% Old SMF structure loaded from file
load('old_analysis.mat', 'SMF');

% Pad with current defaults
SMFDefault = smi_core.SingleMoleculeFitting();
SMF = smi_helpers.padStruct(SMF, SMFDefault);

% SMF now has all current fields with old values preserved
```

Nested structures:
```matlab
% Note: padStruct is shallow (does not recurse into sub-structures)
% For SMF-like nested structures, use smi_core.SingleMoleculeFitting.padSMF
```

**See Also:**
- `smi_core.SingleMoleculeFitting.padSMF` - SMF-specific padding

---

#### removeBorder

**Purpose:** Removes border pixels from images.

**Key Concept:** Edge artifacts in images can affect analysis. This function cleanly crops borders while handling arbitrary dimensions.

**Function Signature:**
```matlab
Image = smi_helpers.removeBorder(Image, Border, Direction)
```

**Inputs:**
- `Image`: Input image (1-4D array)
- `Border`: Border width in pixels (scalar or vector per dimension)
- `Direction`: 'both', 'pre', or 'post' (which edges to remove)

**Output:**
- `Image`: Cropped image

**Usage Examples:**

Simple border removal:
```matlab
% Remove 10-pixel border from all edges
Image = rand(100, 100);
ImageCropped = smi_helpers.removeBorder(Image, 10);
% Result: 80×80 image

% Same border on all dimensions
```

Asymmetric borders:
```matlab
% Different border per dimension
Image = rand(200, 300, 50);  % Y × X × Frames

% Remove [Y_border, X_border, Frame_border]
Border = [20, 30, 5];
ImageCropped = smi_helpers.removeBorder(Image, Border);
% Result: 160×240×40
```

Directional cropping:
```matlab
% Remove border only from start (top/left)
ImageCropped = smi_helpers.removeBorder(Image, 10, 'pre');

% Remove border only from end (bottom/right)
ImageCropped = smi_helpers.removeBorder(Image, 10, 'post');

% Remove from both (default)
ImageCropped = smi_helpers.removeBorder(Image, 10, 'both');
```

Removing edge artifacts:
```matlab
% Load image with bright edges
load('image_with_edges.mat', 'Image');

% Remove 5-pixel border to eliminate edge effects
ImageClean = smi_helpers.removeBorder(Image, 5);

% Continue with analysis on clean image
```

**See Also:**
- Image preprocessing
- ROI extraction

---

### File and Path Utilities

#### selectFiles

**Purpose:** Interactive GUI for selecting multiple files.

**Key Concept:** Provides user-friendly multi-file selection with pattern filtering.

**Function Signature:**
```matlab
[pathname, files] = smi_helpers.selectFiles(startdatadir, txt, pattern)
```

**Inputs:**
- `startdatadir`: Starting directory for file browser
- `txt`: (Optional) Custom message to display in dialog
- `pattern`: (Optional) File pattern filter (default: '*.mat')

**Outputs:**
- `pathname`: Directory path where files reside
- `files`: Cell array of selected filenames

**Usage Examples:**

Basic file selection:
```matlab
% Select .h5 files from data directory
[Path, Files] = smi_helpers.selectFiles('/data', ...
    'Select datasets to process', '*.h5');

% Process selected files
for ii = 1:length(Files)
    FilePath = fullfile(Path, Files{ii});
    fprintf('Processing: %s\n', Files{ii});
    % Load and process
end
```

Multiple file types:
```matlab
% Select any image file
[Path, Files] = smi_helpers.selectFiles(pwd(), ...
    'Select images', '*.*');

% Check file types
for ii = 1:length(Files)
    [~, ~, Ext] = fileparts(Files{ii});
    fprintf('File %d has extension: %s\n', ii, Ext);
end
```

Batch processing:
```matlab
% Let user select results files to combine
[Path, Files] = smi_helpers.selectFiles('C:\Results', ...
    'Select results to merge', '*Results.mat');

% Merge all SMD structures
SMDCombined = [];
for ii = 1:length(Files)
    load(fullfile(Path, Files{ii}), 'SMD');
    SMDCombined = smi_core.SingleMoleculeData.catSMD(SMDCombined, SMD);
end

fprintf('Combined %d datasets with %d localizations\n', ...
    length(Files), length(SMDCombined.X));
```

**See Also:**
- `getFileNames` - List files without GUI
- `getDirectoryNames` - List directories

---

#### getFileNames

**Purpose:** Programmatically lists files matching a pattern.

**Key Concept:** Non-interactive file listing with wildcard support. Faster than selectFiles for scripted workflows.

**Function Signature:**
```matlab
FileNames = smi_helpers.getFileNames(FileDir, NameString)
```

**Inputs:**
- `FileDir`: Directory to search (default: pwd())
- `NameString`: Pattern to match, supports wildcards (default: '*')

**Output:**
- `FileNames`: Cell array of matching filenames (not full paths)

**Usage Examples:**

List all .mat files:
```matlab
% Find all .mat files in directory
Files = smi_helpers.getFileNames('/data/results', '*.mat');

fprintf('Found %d .mat files\n', length(Files));
for ii = 1:length(Files)
    fprintf('  %s\n', Files{ii});
end
```

Pattern matching:
```matlab
% Find files with specific naming pattern
Files = smi_helpers.getFileNames('/data', 'Cell*_Results.h5');

% Process matching files
for ii = 1:length(Files)
    FilePath = fullfile('/data', Files{ii});
    % Load and process
end
```

Batch file processing:
```matlab
% Find all datasets
DataDir = 'C:\Experiments\2025-01-10';
Files = smi_helpers.getFileNames(DataDir, '*.h5');

% Setup for batch processing
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = DataDir;

for ii = 1:length(Files)
    fprintf('Processing %d/%d: %s\n', ii, length(Files), Files{ii});
    SMF.Data.FileName = Files(ii);

    % Run analysis
    % ... (localization code)
end
```

**See Also:**
- `getDirectoryNames` - List directories
- `selectFiles` - GUI file selection

---

#### getDirectoryNames

**Purpose:** Lists subdirectories matching a pattern.

**Function Signature:**
```matlab
DirNames = smi_helpers.getDirectoryNames(ParentDir, NameString)
```

**Usage Example:**
```matlab
% Find all "Cell*" directories
Dirs = smi_helpers.getDirectoryNames('/data/experiment', 'Cell*');

% Process each cell directory
for ii = 1:length(Dirs)
    CellDir = fullfile('/data/experiment', Dirs{ii});
    % Process cell data
end
```

---

### Analysis Helper Functions

#### triangle_threshold

**Purpose:** Automatic threshold detection using the triangle method.

**Key Concept:** Finds optimal threshold for bimodal distributions by fitting a triangle to the histogram. Robust and assumption-free.

**Function Signature:**
```matlab
Cutoff = smi_helpers.triangle_threshold(Data, Percentile, ShowFig)
```

**Inputs:**
- `Data`: Data array (any size, flattened internally)
- `Percentile`: Fraction of extreme values to ignore (default: 1e-4)
- `ShowFig`: Display diagnostic figure (default: false)

**Output:**
- `Cutoff`: Threshold value

**Usage Examples:**

Automatic thresholding for image segmentation:
```matlab
% Load image
Image = imread('fluorescence.tif');

% Find optimal threshold
Threshold = smi_helpers.triangle_threshold(Image, 1e-4, true);
% ShowFig=true displays histogram with threshold line

% Apply threshold
BinaryMask = Image > Threshold;

figure;
subplot(1,2,1); imagesc(Image); title('Original');
subplot(1,2,2); imagesc(BinaryMask); title('Thresholded');
```

Filtering localizations:
```matlab
% Load localization data
load('Results.mat', 'SMD');

% Find photon threshold
PhotonThreshold = smi_helpers.triangle_threshold(SMD.Photons, 1e-3, true);

% Keep only bright localizations
BrightLocs = SMD.Photons > PhotonThreshold;
SMDBright = smi_core.SingleMoleculeData.isolateSubSMD(SMD, find(BrightLocs));

fprintf('Threshold: %.1f photons\n', PhotonThreshold);
fprintf('Kept %d / %d localizations\n', ...
    length(SMDBright.X), length(SMD.X));
```

Quality control:
```matlab
% Find outliers in localization precision
Precision = sqrt(SMD.X_SE.^2 + SMD.Y_SE.^2);

% Triangle threshold for precision cutoff
MaxPrecision = smi_helpers.triangle_threshold(Precision, 1e-4, false);

% Filter
GoodLocs = Precision < MaxPrecision;
SMDGood = smi_core.SingleMoleculeData.isolateSubSMD(SMD, find(GoodLocs));
```

**Algorithm:** Based on Zack GW, Rogers WE, Latt SA (1977), "Automatic measurement of sister chromatid exchange frequency", J. Histochem. Cytochem. 25(7):741–53.

---

#### findStartEndInds

**Purpose:** Finds start and end indices of continuous true sequences in boolean arrays.

**Key Concept:** Useful for identifying events, regions, or continuous states in time series or spatial data.

**Function Signature:**
```matlab
[StartInds, EndInds] = smi_helpers.findStartEndInds(BoolArray)
```

**Input:**
- `BoolArray`: 1D logical array

**Outputs:**
- `StartInds`: Indices where true sequences begin
- `EndInds`: Indices where true sequences end

**Usage Examples:**

Finding events:
```matlab
% Boolean array marking events
Events = [0; 1; 1; 1; 0; 0; 1; 0; 1; 1];

[Starts, Ends] = smi_helpers.findStartEndInds(Events);
% Starts = [2; 7; 9]
% Ends = [4; 7; 10]

% Event 1: indices 2-4 (3 frames)
% Event 2: index 7 (1 frame)
% Event 3: indices 9-10 (2 frames)
```

Identifying trajectory gaps:
```matlab
% Trajectory with missing frames
Trajectory = [1; 1; 1; 0; 0; 1; 1; 0; 1; 1; 1];
% 1 = molecule present, 0 = gap

[Starts, Ends] = smi_helpers.findStartEndInds(Trajectory);

% Analyze segments
for ii = 1:length(Starts)
    SegmentLength = Ends(ii) - Starts(ii) + 1;
    fprintf('Segment %d: frames %d-%d (%d frames)\n', ...
        ii, Starts(ii), Ends(ii), SegmentLength);
end
```

Spatial regions:
```matlab
% Find connected regions above threshold
Image = rand(100, 100);
Threshold = 0.7;

% For each row, find high-intensity regions
for row = 1:size(Image, 1)
    BoolRow = Image(row, :) > Threshold;
    [Starts, Ends] = smi_helpers.findStartEndInds(BoolRow);

    % Mark regions
    for ii = 1:length(Starts)
        fprintf('Row %d: high region at columns %d-%d\n', ...
            row, Starts(ii), Ends(ii));
    end
end
```

**See Also:**
- Event detection
- Trajectory segmentation

---

#### genTimeString

**Purpose:** Generates formatted timestamp strings.

**Key Concept:** Creates consistent, filesystem-safe timestamp strings for file naming and logging.

**Function Signature:**
```matlab
TimeString = smi_helpers.genTimeString(Delimiter, MinFieldWidth)
```

**Inputs:**
- `Delimiter`: Character between time fields (default: '_')
- `MinFieldWidth`: Minimum digits per field (default: 2)

**Output:**
- `TimeString`: Formatted timestamp (e.g., '2025_01_10_14_23_05')

**Usage Examples:**

File naming with timestamps:
```matlab
% Generate unique filename
TimeStamp = smi_helpers.genTimeString('_', 2);
FileName = sprintf('Results_%s.mat', TimeStamp);
% e.g., 'Results_2025_01_10_14_23_05.mat'

save(FileName, 'SMD', 'SMF');
```

Different delimiters:
```matlab
% Use dash delimiter
TimeStamp = smi_helpers.genTimeString('-', 2);
% Result: '2025-01-10-14-23-05'

% Use no delimiter (compact)
TimeStamp = smi_helpers.genTimeString('', 2);
% Result: '20250110142305'
```

Logging:
```matlab
% Add timestamp to log messages
TimeStamp = smi_helpers.genTimeString(':', 2);
fprintf('[%s] Starting analysis...\n', TimeStamp);
% [2025:01:10:14:23:05] Starting analysis...
```

**See Also:**
- File naming conventions
- Logging utilities

---

### System Utilities

#### versionSMITE

**Purpose:** Reports current smite version.

**Usage:**
```matlab
Version = smi_helpers.versionSMITE();
fprintf('smite version: %s\n', Version);
```

---

#### requiredToolboxes

**Purpose:** Lists MATLAB toolboxes required by smite components.

**Usage:**
```matlab
smi_helpers.requiredToolboxes();
% Prints required toolboxes for each namespace
```

**Common Requirements:**
- Image Processing Toolbox - Core image operations
- Statistics and Machine Learning Toolbox - Statistical methods
- Parallel Computing Toolbox - GPU acceleration, parfor
- Curve Fitting Toolbox - FRC and fitting operations (optional)
- Optimization Toolbox - Advanced fitting (optional)

---

#### mkSMITETmpDir

**Purpose:** Creates temporary directories for smite testing and examples.

**Usage:**
```matlab
TmpDir = smi_helpers.mkSMITETmpDir('TestName');
% Creates: tempdir/smite/TestName/
```

---

## Common Usage Patterns

### Complete Filtering Pipeline for BaGoL

```matlab
% Load localized data
load('LocalizedData.mat', 'SMD');
fprintf('Starting: %d localizations\n', length(SMD.X));

% Apply standard BaGoL filter sequence
SMD = smi_helpers.Filters.filterNonNeg(SMD, false);
SMD = smi_helpers.Filters.filterIntensity(SMD, false, 3);
SMD = smi_helpers.Filters.inflateSE(SMD, false, 1.5);
SMD = smi_helpers.Filters.filterFC(SMD, false, 2);

% Only for PAINT/photoactivation (not dSTORM!)
SMD = smi_helpers.Filters.filterNN(SMD, false, 5, 3);

fprintf('After filtering: %d localizations\n', length(SMD.X));

% Ready for BaGoL analysis
```

### Parallel Processing with Subdivision

```matlab
% Load data and image
load('Data.mat', 'Data', 'SMD');

% Subdivide for parallel processing
TileSize = [64, 64];
[ImageTiles, ImageROIs] = smi_helpers.subdivideImage(Data, TileSize);
[SMDTiles, SMDROIs] = smi_helpers.subdivideSMD(SMD, TileSize);

% Process in parallel
Results = cell(length(ImageTiles), 1);
parfor ii = 1:length(ImageTiles)
    % Each worker processes one tile
    Results{ii} = processTile(ImageTiles{ii}, SMDTiles(ii));
end

% Combine results
FinalResult = combineResults(Results);
```

### Batch File Processing

```matlab
% Find all datasets
DataDir = 'C:\Experiments\Batch1';
Files = smi_helpers.getFileNames(DataDir, 'Cell*.h5');

fprintf('Found %d datasets to process\n', length(Files));

% Setup parameters
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.FileDir = DataDir;
SMF.Data.PixelSize = 0.108;
SMF.Data.FrameRate = 100;

% Process each file
for ii = 1:length(Files)
    fprintf('[%d/%d] Processing %s\n', ii, length(Files), Files{ii});

    SMF.Data.FileName = Files(ii);

    % Load and localize
    % ... (analysis code)

    % Save with timestamp
    TimeStamp = smi_helpers.genTimeString('_', 2);
    ResultFile = sprintf('%s_Results_%s.mat', Files{ii}(1:end-3), TimeStamp);
    save(fullfile(DataDir, 'Results', ResultFile), 'SMD', 'SMF');
end
```

### Adaptive Thresholding Workflow

```matlab
% Load data
load('Results.mat', 'SMD');

% Automatically determine photon threshold
PhotonCutoff = smi_helpers.triangle_threshold(SMD.Photons, 1e-4, true);

% Automatically determine precision threshold
Precision = sqrt(SMD.X_SE.^2 + SMD.Y_SE.^2);
PrecisionCutoff = smi_helpers.triangle_threshold(Precision, 1e-4, true);

% Apply thresholds
GoodLocs = (SMD.Photons > PhotonCutoff) & (Precision < PrecisionCutoff);
SMDFiltered = smi_core.SingleMoleculeData.isolateSubSMD(SMD, find(GoodLocs));

fprintf('Adaptive thresholds:\n');
fprintf('  Photons > %.1f\n', PhotonCutoff);
fprintf('  Precision < %.3f pixels\n', PrecisionCutoff);
fprintf('  Kept %d / %d localizations (%.1f%%)\n', ...
    length(SMDFiltered.X), length(SMD.X), ...
    100 * length(SMDFiltered.X) / length(SMD.X));
```

---

## Performance Considerations

### Parallel Processing with Subdivision

When processing large datasets, subdivision enables parallel processing:

```matlab
% Determine optimal tile size based on available workers
NWorkers = feature('numcores');
TotalPixels = SMD.XSize * SMD.YSize;
PixelsPerTile = TotalPixels / (NWorkers * 2);  % 2 tiles per worker
TileSize = round(sqrt(PixelsPerTile));

% Subdivide
[SMDTiles, ~] = smi_helpers.subdivideSMD(SMD, [TileSize, TileSize]);

% Process with parfor
parfor ii = 1:length(SMDTiles)
    % Process tile
end
```

### Memory Efficiency

For large arrays, avoid unnecessary copies:

```matlab
% Efficient: operate in place
SMD = smi_helpers.Filters.filterNonNeg(SMD, false);
SMD = smi_helpers.Filters.filterIntensity(SMD, false, 3);

% Less efficient: creates intermediate variables
SMD1 = smi_helpers.Filters.filterNonNeg(SMD, false);
SMD2 = smi_helpers.Filters.filterIntensity(SMD1, false, 3);
```

---

## Troubleshooting

### Filters Removing Too Many Localizations

**Problem:** Filter pipeline removes most localizations

**Solutions:**
```matlab
% Check each filter step individually
fprintf('Start: %d\n', length(SMD.X));

SMD = smi_helpers.Filters.filterNonNeg(SMD, false);
fprintf('After filterNonNeg: %d\n', length(SMD.X));

SMD = smi_helpers.Filters.filterIntensity(SMD, false, 3);
fprintf('After filterIntensity: %d\n', length(SMD.X));

% Identify problematic filter and adjust parameters
```

### Subdivision Creates Empty Tiles

**Problem:** Some tiles have no localizations

**Solution:** This is expected for sparse data. Filter empty tiles:
```matlab
[SMDTiles, ROIs] = smi_helpers.subdivideSMD(SMD, [100, 100]);

% Find non-empty tiles
NonEmpty = arrayfun(@(s) length(s.X) > 0, SMDTiles);

% Process only non-empty
for ii = find(NonEmpty)'
    % Process SMDTiles(ii)
end
```

---

## See Also

### Core Concepts
- [SMD Structure](../core-concepts/smd-structure.md)
- [Architecture Overview](../core-concepts/architecture.md)

### Workflows
- [SMLM Analysis](../workflows/smlm-analysis.md)
- BaGoL workflow documentation

### Related APIs
- [+smi_core Namespace](smi-core.md) - Core processing
- ROITools documentation

---

## Summary

The +smi_helpers namespace provides essential utility functions that support smite analyses:

**Filtering:** `Filters` class for quality filtering, especially for BaGoL preparation.

**Subdivision:** `subdivideImage` and `subdivideSMD` for spatial tiling and parallel processing.

**Array Tools:** `arrayMUX`, `compressToRange`, `padStruct`, `removeBorder` for data manipulation.

**File Tools:** `selectFiles`, `getFileNames`, `getDirectoryNames` for file management.

**Analysis Helpers:** `triangle_threshold`, `findStartEndInds`, `genTimeString` for common analysis tasks.

These utilities simplify workflows, improve code readability, and enable advanced processing patterns. While not required for basic analyses, mastering these tools enables efficient custom pipelines and batch processing.
