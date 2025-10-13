---
title: "Testing Reference"
category: "reference"
level: "intermediate"
tags: ["testing", "unit-tests", "verification", "quality-assurance", "debugging"]
prerequisites: ["../getting-started/installation.md", "../core-concepts/architecture.md"]
related: ["../troubleshooting/installation-issues.md", "../troubleshooting/compilation-errors.md", "../getting-started/quickstart.md"]
summary: "Complete reference for running tests in smite, from full test suites to individual component tests, with guidance on interpreting results"
estimated_time: "20 minutes"
last_updated: "2025-10-11"
status: "complete"
---

# Testing Reference

## Purpose

This reference provides comprehensive guidance for testing smite functionality, including:
- Running the complete test suite with `run_tests`
- Executing individual class unit tests
- Key tests for verifying installation
- Understanding test outputs and results
- Interpreting test failures and debugging
- Working with expected results and test data

Testing is essential for verifying that smite is installed correctly, that modifications haven't broken functionality, and that GPU acceleration is working as expected.

## Prerequisites

- smite successfully installed (see [Installation Guide](../getting-started/installation.md))
- Understanding of [smite architecture](../core-concepts/architecture.md)
- MATLAB R2021a or later
- Required MATLAB toolboxes installed
- (Optional) GPU configured for GPU-accelerated tests

## Overview

smite includes an extensive test suite covering all major components:
- **Main workflows:** SMLM, SPT, BaGoL
- **Core processing:** LocalizeData, FindROI, FrameConnection, DriftCorrection
- **Clustering:** DBSCAN, Voronoi, H-SET, PairCorrelation
- **Statistical methods:** HMM, ChangeDetection, DiffusionEstimator
- **Simulation:** SimSMLM, SimSPT, GaussBlobs
- **Visualization:** GenerateImages, GenerateMovies
- **PSF tools:** PointSpreadFunction, Zernike

Tests generate their own data when needed, ensuring reproducibility across installations. Results are saved to temporary directories for inspection.

## Running the Complete Test Suite

### The run_tests Script

The `run_tests` function executes all major unit tests in sequence, providing comprehensive verification of smite functionality.

**Location:** `smite/MATLAB/run_tests.m`

**Basic Usage:**

```matlab
run_tests
```

**What it Does:**

1. Tests all major namespaces (+smi, +smi_core, +smi_cluster, +smi_psf, +smi_sim, +smi_stat, +smi_vis)
2. Executes unit tests for each class in sequence
3. Uses try-catch blocks to continue testing even if individual tests fail
4. Prints progress messages for each test
5. Saves results to `tempdir/smite/unitTest/name_of_test/`

**Expected Runtime:** 5-15 minutes depending on hardware (GPU vs CPU)

### Complete Test List

The full test suite includes these tests (in execution order):

**Main Workflows:**
```matlab
smi.SMLM.unitTest()
smi.SPT.unitTestFFGC()
```

**Clustering (+smi_cluster):**
```matlab
smi_cluster.Clustering.unitTest()
smi_cluster.PairCorrelation.unitTest()
smi_cluster.StatisticsClustering.unitTest()
```

**Core Processing (+smi_core):**
```matlab
smi_core.ChannelRegistration.unitTest()
smi_core.DataToPhotons.unitTest()
smi_core.DriftCorrection.unitTest()
smi_core.FrameConnection.unitTest()
smi_core.LocalizeData.unitTest()
smi_core.Threshold.unitTest()
```

**PSF Tools (+smi_psf):**
```matlab
smi_psf.PointSpreadFunction.unitTest()
smi_psf.Zernike.unitTest()
```

**Simulation (+smi_sim):**
```matlab
smi_sim.GaussBlobs.unitTest()
smi_sim.SimSMLM.unitTest()
```

**Statistical Methods (+smi_stat):**
```matlab
smi_stat.ChangeDetection.unitTest()
smi_stat.DiffusionEstimator.unitTest()
smi_stat.HMM.unitTest()
```

**Visualization (+smi_vis):**
```matlab
smi_vis.GenerateImages.blobColorOverlay_unitTest()
smi_vis.GenerateImages.circleImage_unitTest()
smi_vis.GenerateImages.colorImage_unitTest()
smi_vis.GenerateImages.driftImage_unitTest()
smi_vis.GenerateImages.gaussianImage_unitTest()
smi_vis.GenerateImages.histogramImage_unitTest()
```

### Understanding Test Output

**Typical Output:**

```matlab
>> run_tests

smi.SMLM.unitTest
Running smi.SMLM.unitTest.
Testing all smi.SMLM functionality.
Testing loading of various datatypes.
Simulating data.
Saving data as mat files.
Loading and analyzing data saved as mat files.
   (Only doing box finding and fitting.)
Loading and analyzing data saved as mat file done.

Saving data as h5 files.
Loading and analyzing data saved as h5 files.
   (Only doing box finding and fitting.)
Loading and analyzing data saved as h5 file done.

Simulating realistic 2D SMLM data
Saving realistic SMSR data.
Testing 2D analysis.
-> ResultsDir = C:\Users\...\Temp\smite\unitTest\SMLM
Files produced in C:\Users\...\Temp\smite\unitTest\SMLM:
Done.

smi.SPT.unitTestFFGC (frame-to-frame and gap closing processes)
[progress output...]
```

**Success Indicators:**
- Tests complete without throwing unhandled errors
- "Done" messages appear
- Files are created in output directories
- No MATLAB error messages in red

**Failure Indicators:**
- Error messages in red
- Tests abort before completion
- Missing output files
- NaN or Inf values in results

## Running Individual Unit Tests

### Why Run Individual Tests

Individual tests are useful for:
- **Targeted verification:** Test specific functionality after changes
- **Debugging:** Isolate problems to specific components
- **Development:** Quick iteration when modifying code
- **Installation verification:** Test critical components without full suite

### Syntax

All unit test methods follow the same pattern:

```matlab
Success = ClassName.unitTest()
```

**Returns:** Boolean array `Success` indicating which sub-tests passed (1 = pass, 0 = fail)

### Core Component Tests

**LocalizeData (Localization Pipeline):**

```matlab
Success = smi_core.LocalizeData.unitTest()
```

Tests vital localization functionality:
- Simulates data with known emitter positions
- Runs complete localization pipeline
- Verifies localizations match input
- Checks thresholding consistency
- Saves diagnostic plots to `tempdir/smite/unitTest/LocalizeData/`

**Expected Output:**
```matlab
Success =
     1
```

**Runtime:** ~30 seconds

---

**FrameConnection (Temporal Clustering):**

```matlab
Success = smi_core.FrameConnection.unitTest()
```

Tests frame connection algorithm:
- Creates synthetic blinking data
- Tests LAP-FC algorithm
- Verifies temporal clustering
- Checks precision improvement

**Runtime:** ~1 minute

---

**DriftCorrection (Stage Drift):**

```matlab
Success = smi_core.DriftCorrection.unitTest()
```

Tests drift correction algorithms:
- Simulates data with known drift
- Tests KNN drift correction (intra and inter-dataset)
- Verifies drift recovery accuracy
- Tests brightfield registration (if available)

**Runtime:** ~2 minutes

---

**DataToPhotons (Camera Calibration):**

```matlab
Success = smi_core.DataToPhotons.unitTest()
```

Tests conversion from ADU to photons:
- Tests EMCCD calibration (scalar gain/offset)
- Tests sCMOS calibration (per-pixel maps)
- Verifies correct photon counting
- Tests read noise propagation

**Runtime:** ~10 seconds

---

**Threshold (Quality Filtering):**

```matlab
Success = smi_core.Threshold.unitTest()
```

Tests thresholding functionality:
- Tests threshold flag setting
- Verifies threshold application
- Tests multiple threshold combinations
- Checks rejection visualization

**Runtime:** ~15 seconds

---

**ChannelRegistration (Multi-Channel Alignment):**

```matlab
Success = smi_core.ChannelRegistration.unitTest()
```

Tests channel registration:
- Creates synthetic multi-channel data
- Tests transform computation
- Verifies coordinate transformation
- Checks registration accuracy

**Runtime:** ~1 minute

### Main Workflow Tests

**SMLM (Complete SMLM Pipeline):**

```matlab
Success = smi.SMLM.unitTest()
```

Comprehensive SMLM workflow test:
- Tests loading .mat and .h5 files
- Tests complete analysis pipeline
- Tests box finding, fitting, thresholding
- Tests frame connection and drift correction
- Generates realistic test data

**Sub-tests:**
- `Success(1)`: .mat file loading and analysis
- `Success(2)`: .h5 file loading and analysis
- `Success(3)`: Full pipeline with corrections

**Expected Output:**
```matlab
Success =
     1     1     1
```

**Runtime:** ~3-5 minutes

**Output Location:** `tempdir/smite/unitTest/SMLM/`

**Generated Files:**
- `SMLM_testData.mat` - Results file
- Various diagnostic plots

---

**SPT (Single Particle Tracking):**

```matlab
Success = smi.SPT.unitTestFFGC()
```

Tests frame-to-frame and gap-closing tracking:
- Simulates particle trajectories
- Tests frame-to-frame connection
- Tests gap closing over frame gaps
- Converts to trajectory format (TR)
- Generates tracking movie

**Expected Output:**
```matlab
Success =
     1
```

**Runtime:** ~2 minutes (GUI opens with movie player)

**Note:** This test opens a GUI showing the tracking results. Close the GUI to complete the test.

### Clustering Tests

**Clustering (Algorithms Suite):**

```matlab
Success = smi_cluster.Clustering.unitTest()
```

Tests clustering algorithms:
- DBSCAN
- Voronoi tessellation
- H-SET hierarchical clustering
- Nearest neighbor analysis
- Tests on 2D and 3D data

**Runtime:** ~2 minutes

**Output:** Generates numerous plots in `tempdir/smite/unitTest/Clustering/`

---

**PairCorrelation (Spatial Statistics):**

```matlab
Success = smi_cluster.PairCorrelation.unitTest()
```

Tests pair correlation analysis:
- Tests 2D and 3D correlation functions
- Tests Hopkins statistic
- Tests Ripley's K function
- Verifies statistical significance testing

**Runtime:** ~1 minute

---

**StatisticsClustering (Statistical Clustering):**

```matlab
Success = smi_cluster.StatisticsClustering.unitTest()
```

Tests statistical clustering methods:
- Tests various clustering algorithms
- Tests cluster validation metrics
- Tests parameter optimization

**Runtime:** ~1 minute

### Simulation Tests

**SimSMLM (SMLM Data Simulation):**

```matlab
Success = smi_sim.SimSMLM.unitTest()
```

Tests SMLM simulation:
- Creates realistic SMLM datasets
- Tests various blinking models
- Tests photobleaching
- Tests multi-color simulation

**Runtime:** ~30 seconds

---

**GaussBlobs (Gaussian Blob Generation):**

```matlab
Success = smi_sim.GaussBlobs.unitTest()
```

Tests Gaussian blob image generation:
- Tests random blob generation
- Tests specified blob positions
- Tests Poisson noise addition
- Verifies photon counting

**Runtime:** ~15 seconds

### PSF Tests

**PointSpreadFunction (PSF Modeling):**

```matlab
Success = smi_psf.PointSpreadFunction.unitTest()
```

Tests PSF calculation methods:
- Runs multiple sub-tests:
  - `crlbPSFPupil_unitTest` - CRLB calculation
  - `optimPSFZernike_unitTest` - Zernike optimization
  - `psfROIStack_unitTest` - PSF ROI generation
  - `scalarPSFPrasadZone_unitTest` - Scalar PSF calculation
  - `zernikeImage_unitTest` - Zernike image generation

**Runtime:** ~2-3 minutes

**Note:** Some sub-tests require additional toolboxes.

---

**Zernike (Zernike Polynomials):**

```matlab
Success = smi_psf.Zernike.unitTest()
```

Tests Zernike polynomial calculations:
- Tests polynomial generation
- Tests orthogonality
- Tests coefficient fitting

**Runtime:** ~30 seconds

### Statistical Method Tests

**HMM (Hidden Markov Model):**

```matlab
Success = smi_stat.HMM.unitTest()
```

Tests HMM state detection:
- Tests Viterbi algorithm
- Tests Baum-Welch training
- Tests state sequence recovery

**Runtime:** ~30 seconds

---

**ChangeDetection (Change Point Detection):**

```matlab
Success = smi_stat.ChangeDetection.unitTest()
```

Tests change point detection:
- Tests various detection algorithms
- Tests on simulated state changes
- Verifies change point accuracy

**Runtime:** ~30 seconds

---

**DiffusionEstimator (Diffusion Analysis):**

```matlab
Success = smi_stat.DiffusionEstimator.unitTest()
```

Tests diffusion coefficient estimation:
- Tests MSD calculation
- Tests various diffusion models
- Tests state identification

**Runtime:** ~30 seconds

### Visualization Tests

Visualization tests verify image and overlay generation:

```matlab
Success = smi_vis.GenerateImages.blobColorOverlay_unitTest()
Success = smi_vis.GenerateImages.circleImage_unitTest()
Success = smi_vis.GenerateImages.colorImage_unitTest()
Success = smi_vis.GenerateImages.driftImage_unitTest()
Success = smi_vis.GenerateImages.gaussianImage_unitTest()
Success = smi_vis.GenerateImages.histogramImage_unitTest()
```

**Runtime:** ~10 seconds each

**Output:** Generates test images in respective test directories

## Key Tests for Installation Verification

After installing smite, run these critical tests to verify proper installation:

### 1. LocalizeData Test (GPU Verification)

**Why:** Tests GPU-accelerated fitting (GaussMLE) and complete localization pipeline.

```matlab
Success = smi_core.LocalizeData.unitTest()
```

**What it Verifies:**
- GPU is detected and working
- CUDA files are compiled correctly
- Localization pipeline functions end-to-end
- Thresholding works properly

**If this passes:** GPU acceleration is working correctly.

**If this fails:** Check GPU installation, CUDA compilation, or recompile PTX files.

### 2. SMLM Test (Complete Workflow)

**Why:** Tests entire SMLM workflow from file loading to final results.

```matlab
Success = smi.SMLM.unitTest()
```

**What it Verifies:**
- File I/O (.mat and .h5)
- Data conversion (ADU to photons)
- Localization
- Thresholding
- Frame connection
- Drift correction
- Results saving

**If this passes:** smite is fully functional for SMLM analysis.

### 3. SPT Test (Tracking Workflow)

**Why:** Tests tracking algorithms and trajectory analysis.

```matlab
Success = smi.SPT.unitTestFFGC()
```

**What it Verifies:**
- Frame-to-frame connection
- Gap closing
- Trajectory formation
- Cost matrix calculation
- LAP solver

**If this passes:** Tracking functionality is working.

### Minimal Verification (Quick Check)

For a quick sanity check after installation:

```matlab
% Quick GPU test
gpuDevice

% Quick localization test
SMF = smi_core.SingleMoleculeFitting();
SMF.Data.CameraType = 'EMCCD';
B = smi_sim.GaussBlobs.genRandomBlobImage();
LD = smi_core.LocalizeData(poissrnd(B), SMF);
[SMD] = LD.genLocalizations();
fprintf('Found %d localizations\n', length(SMD.X));
```

**Expected:** Should find ~250 localizations without errors.

## Test Output Locations

### Temporary Directory Structure

Tests save results to platform-specific temporary directories:

**Windows:**
```
C:\Users\YourName\AppData\Local\Temp\smite\unitTest\
```

**Linux/Mac:**
```
/tmp/smite/unitTest/
```

**Finding Test Results:**

```matlab
% Get base temp directory
tempDir = tempdir

% Navigate to smite test results
testDir = fullfile(tempdir, 'smite', 'unitTest')

% List all test result directories
dir(testDir)
```

### Per-Test Directories

Each unit test creates its own subdirectory:

```
tempdir/smite/unitTest/
├── SMLM/                      % smi.SMLM.unitTest results
├── LocalizeData/              % smi_core.LocalizeData.unitTest results
├── Clustering/                % smi_cluster.Clustering.unitTest results
├── FrameConnection/           % smi_core.FrameConnection.unitTest results
├── DriftCorrection/           % smi_core.DriftCorrection.unitTest results
└── [other test directories]
```

### Typical Test Output Files

**Data Files:**
- `*.mat` - MATLAB data files (SMD, SMF, results)
- `*.h5` - HDF5 data files (raw datasets)

**Diagnostic Plots:**
- `*.png` - Saved figure images
- `*.fig` - MATLAB figure files
- `LD1.png` - LocalizeData diagnostic plot

**Movies (if generated):**
- `*.avi` - Tracking or localization movies

### Expected Results Directory

smite includes a directory with expected results for comparison:

**Location:** `smite/MATLAB/ExpectedResults/`

**Contents:**
- Reference output from unit tests
- Expected values for validation
- Note: Very large files excluded to reduce repository size

**Usage:**

```matlab
% Compare your results to expected
yourResults = load(fullfile(tempdir, 'smite', 'unitTest', 'LocalizeData', 'results.mat'));
expectedResults = load(fullfile(smite_path, 'MATLAB', 'ExpectedResults', 'unitTest', 'LocalizeData', 'results.mat'));

% Compare fields
fprintf('Your localizations: %d\n', length(yourResults.SMD.X));
fprintf('Expected localizations: %d\n', length(expectedResults.SMD.X));
```

**Note:** Due to stochastic nature of tests, exact numerical agreement is not expected. Tests should be "close" but not identical.

## Interpreting Test Results

### Understanding Success Arrays

Many tests return boolean arrays indicating sub-test results:

```matlab
Success = smi.SMLM.unitTest()
% Success = [1 1 1] means all 3 sub-tests passed
% Success = [1 0 1] means sub-test 2 failed
```

**Checking Results:**

```matlab
Success = smi.SMLM.unitTest();

if all(Success)
    fprintf('All tests passed!\n');
else
    fprintf('Failed tests: %s\n', mat2str(find(~Success)));
end
```

### Common Test Failure Patterns

**All Tests Fail:**
- **Cause:** Fundamental installation issue
- **Check:** MATLAB path, GPU, required toolboxes
- **Action:** Review [Installation Issues](../troubleshooting/installation-issues.md)

**GPU-Related Tests Fail:**
- **Cause:** GPU not detected or CUDA compilation issue
- **Check:** `gpuDevice` output
- **Action:** Recompile CUDA files or check [GPU Problems](../troubleshooting/gpu-problems.md)

**Specific Namespace Fails:**
- **Cause:** Missing toolbox or dependency
- **Check:** `ver` to list installed toolboxes
- **Action:** Install required toolbox

**Intermittent Failures:**
- **Cause:** Stochastic test behavior
- **Action:** Re-run test. If persists, investigate further.

### Test Verbosity

Many tests have adjustable verbosity:

```matlab
% Verbose levels:
% 0 - Silent (no output)
% 1 - Basic progress (default)
% 2 - Detailed progress
% 3 - Full diagnostics with saved figures

% Example with LocalizeData
SMF = smi_core.SingleMoleculeFitting();
Data = smi_sim.GaussBlobs.genRandomBlobImage();

% Silent
LD = smi_core.LocalizeData(Data, SMF, 0);

% Detailed
LD = smi_core.LocalizeData(Data, SMF, 2);

% Full diagnostics
LD = smi_core.LocalizeData(Data, SMF, 3);
```

**Verbosity 3 Benefits:**
- Saves diagnostic images
- Generates progress plots
- Creates detailed logs
- Useful for debugging failures

## Stochastic Nature of Tests

### Why Tests Vary

Many smite tests involve stochastic processes:
- Random number generation for synthetic data
- Poisson noise in simulated images
- Random initial conditions
- Monte Carlo simulations

**Consequence:** Running the same test multiple times produces slightly different numerical results.

### Seeding for Reproducibility

Some tests seed the random number generator:

```matlab
% From LocalizeData.unitTest
rng(1234)  % Set seed for reproducibility
```

**When tests are seeded:** Results should be identical across runs.

**When tests are not seeded:** Small variations are expected and normal.

### Acceptable Variation

**Expected variations:**
- Localization positions: ±0.1 pixels
- Photon counts: ±5%
- Number of localizations: ±2%
- Fit quality metrics: ±10%

**Unacceptable variations:**
- Order-of-magnitude differences
- Consistently NaN or Inf results
- Zero localizations when many expected
- Complete test failure

### Comparing Results

When comparing test results:

```matlab
% Load two test runs
run1 = load('test_run1.mat');
run2 = load('test_run2.mat');

% Compare number of localizations
n1 = length(run1.SMD.X);
n2 = length(run2.SMD.X);
fprintf('Run 1: %d localizations\n', n1);
fprintf('Run 2: %d localizations\n', n2);
fprintf('Difference: %.1f%%\n', 100 * abs(n1 - n2) / n1);

% Compare photon counts
photon_diff = abs(median(run1.SMD.Photons) - median(run2.SMD.Photons));
fprintf('Median photon difference: %.1f (%.1f%%)\n', ...
    photon_diff, 100 * photon_diff / median(run1.SMD.Photons));
```

**Rule of Thumb:** If differences are <5%, tests are consistent.

## Debugging Test Failures

### Step-by-Step Debugging Process

**1. Identify Which Test Failed:**

```matlab
Success = smi.SMLM.unitTest()
% Success = [1 0 1] - second sub-test failed
```

**2. Run Test with High Verbosity:**

```matlab
% Modify test file temporarily to increase verbosity
% Or run test components manually with verbose output
```

**3. Check Test Output Directory:**

```matlab
testDir = fullfile(tempdir, 'smite', 'unitTest', 'SMLM')
dir(testDir)
% Look for error logs, incomplete files, or diagnostic images
```

**4. Examine Error Messages:**

Look for:
- GPU errors: "Invalid PTX file", "Out of GPU memory"
- File errors: "File not found", "Permission denied"
- Numerical errors: "NaN or Inf values detected"
- Fitting errors: "Maximum iterations reached"

**5. Test Dependencies:**

```matlab
% Test GPU
gpuDevice

% Test required toolboxes
ver

% Test file permissions
testFile = fullfile(tempdir, 'test.mat');
save(testFile, 'testvar', 1);
delete(testFile);
```

### Common Failure Modes

**GPU Not Found:**

```
Error: GPU device not found or not supported
```

**Solution:**
```matlab
% Check GPU
gpuDevice

% Verify compute capability
g = gpuDevice;
fprintf('Compute capability: %s\n', g.ComputeCapability);
% Must be >= 5.0
```

**PTX Compilation Error:**

```
Error: Invalid PTX file
```

**Solution:**
```matlab
% Recompile CUDA files
cd(fullfile(smite_path, 'MATLAB', 'source', 'cuda'))
cuda_Make
```

**Out of Memory:**

```
Error: Out of GPU memory
```

**Solution:**
```matlab
% Reset GPU
reset(gpuDevice)

% Or reduce data size in test
% Edit test file to use smaller datasets
```

**Missing Toolbox:**

```
Error: Undefined function or variable 'fitdist'
```

**Solution:** Install Statistics and Machine Learning Toolbox

**File Permission Error:**

```
Error: Permission denied
```

**Solution:**
```matlab
% Check temp directory permissions
tempDir = tempdir;
[status, msg] = fileattrib(tempDir);
fprintf('Writable: %d\n', msg.UserWrite);

% Use alternative temp directory
setenv('TMPDIR', 'C:\alternate_temp');
```

### Manual Test Reproduction

To debug a test failure, run test code manually:

```matlab
% Example: Manually run LocalizeData test steps

% 1. Setup
SaveDir = fullfile(tempdir, 'smite', 'unitTest', 'LocalizeData');
if ~exist(SaveDir, 'dir')
    mkdir(SaveDir);
end

% 2. Generate test data (from unitTest source)
rng(1234)
FrameSizeFull = 256;
NFrames = 10;
SMD = smi_core.SingleMoleculeData.createSMD();
SMD.X = repmat(128 + 64*[0; 1; 1; -1; -1], [NFrames, 1]);
SMD.Y = repmat(128 + 64*[0; 1; -1; 1; -1], [NFrames, 1]);
SMD.Photons = 1e3 * ones(5*NFrames, 1);
SMD.PSFSigma = 1.3;
SMD.FrameNum = repelem((1:NFrames).', 5);
SMD.Bg = zeros(5*NFrames, 1);
SMD.NFrames = NFrames;
SMD.YSize = FrameSizeFull;
SMD.XSize = FrameSizeFull;
[~, ScaledData] = smi_sim.GaussBlobs.gaussBlobImage(SMD);

% 3. Generate SMF
SMF = smi_core.SingleMoleculeFitting;
SMF.BoxFinding.BoxSize = 10;
SMF.Fitting.PSFSigma = SMD.PSFSigma;

% 4. Run localization with high verbosity
LD = smi_core.LocalizeData(ScaledData, SMF, 3);
[SMDout, SMDPreThresh] = LD.genLocalizations();

% 5. Check results
fprintf('Input emitters: %d\n', length(SMD.X));
fprintf('Output localizations: %d\n', length(SMDout.X));
fprintf('Pre-threshold localizations: %d\n', length(SMDPreThresh.X));

% 6. Verify success
NEmitters = 5;
Success = all(SMDout.FrameNum==repelem((1:NFrames).', NEmitters)) ...
    && (numel(SMDout.X)==NEmitters*NFrames);
fprintf('Test success: %d\n', Success);
```

## Advanced Testing Topics

### Writing Custom Tests

You can write custom tests following smite patterns:

```matlab
function Success = myCustomTest()
% myCustomTest - Custom test for my modifications
%
% OUTPUTS:
%   Success: Boolean array of test results

Success = false(1, 3);  % Three sub-tests

% Test 1: Basic functionality
try
    % Your test code here
    result1 = myFunction(testInput);
    Success(1) = true;
catch ME
    fprintf('Test 1 failed: %s\n', ME.message);
end

% Test 2: Edge cases
try
    % Test edge cases
    result2 = myFunction(edgeCaseInput);
    Success(2) = true;
catch ME
    fprintf('Test 2 failed: %s\n', ME.message);
end

% Test 3: Integration
try
    % Test integration with other components
    result3 = integratedWorkflow();
    Success(3) = true;
catch ME
    fprintf('Test 3 failed: %s\n', ME.message);
end

fprintf('Tests passed: %d / %d\n', sum(Success), length(Success));
end
```

### Performance Testing

Test execution speed:

```matlab
% Time individual operations
tic
Success = smi_core.LocalizeData.unitTest();
elapsed = toc;
fprintf('LocalizeData test completed in %.1f seconds\n', elapsed);

% Compare GPU vs CPU
gpuDevice(1)  % Use GPU
tic
[SMD1] = performLocalization(Data, SMF);
gpu_time = toc;

% Disable GPU (not directly possible, but compare with CPU-only methods)
fprintf('GPU time: %.2f seconds\n', gpu_time);
```

### Regression Testing

After making modifications, ensure existing functionality still works:

```matlab
% Before modifications: Run tests and save results
run_tests
copyfile(fullfile(tempdir, 'smite', 'unitTest'), ...
         fullfile(tempdir, 'smite', 'unitTest_baseline'));

% After modifications: Run tests and compare
run_tests

% Compare results
baseline = load(fullfile(tempdir, 'smite', 'unitTest_baseline', 'SMLM', 'results.mat'));
current = load(fullfile(tempdir, 'smite', 'unitTest', 'SMLM', 'results.mat'));

% Detailed comparison
compare_results(baseline, current);
```

### Continuous Integration

For automated testing in development:

```matlab
% Script: run_tests_ci.m
% Automated test runner for CI/CD

diary('test_log.txt')
diary on

fprintf('Starting smite test suite...\n');
fprintf('Date: %s\n', datetime('now'));

try
    run_tests
    fprintf('\nAll tests completed successfully!\n');
    exit(0)
catch ME
    fprintf('\nTest suite failed:\n');
    fprintf('Error: %s\n', ME.message);
    fprintf('Location: %s\n', ME.stack(1).file);
    diary off
    exit(1)
end
```

## Summary

The smite test suite provides comprehensive verification of functionality:

**Complete Test Suite:** `run_tests` executes all unit tests (~5-15 minutes)

**Individual Tests:** Each class has a `.unitTest()` method for targeted testing

**Key Installation Tests:**
- `smi_core.LocalizeData.unitTest()` - Verifies GPU and localization
- `smi.SMLM.unitTest()` - Verifies complete SMLM workflow
- `smi.SPT.unitTestFFGC()` - Verifies tracking functionality

**Test Results:** Saved to `tempdir/smite/unitTest/[test_name]/`

**Stochastic Nature:** Small variations between runs are normal and expected

**Debugging:** Use high verbosity, check output directories, and run test components manually

Regular testing ensures smite functions correctly across installations, platforms, and modifications. When reporting issues, include test results and error messages to help diagnose problems.

## See Also

### Getting Started
- [Installation Guide](../getting-started/installation.md) - Install smite properly
- [Quick Start](../getting-started/quickstart.md) - Get running quickly
- [First Analysis](../getting-started/first-analysis.md) - Complete analysis workflow

### Troubleshooting
- [Installation Issues](../troubleshooting/installation-issues.md) - Resolve installation problems
- [Compilation Errors](../troubleshooting/compilation-errors.md) - Fix mex/CUDA issues
- [GPU Problems](../troubleshooting/gpu-problems.md) - GPU-specific troubleshooting

### Core Concepts
- [Architecture](../core-concepts/architecture.md) - Understanding smite structure
- [SMF Structure](../core-concepts/smf-structure.md) - Analysis parameters
- [SMD Structure](../core-concepts/smd-structure.md) - Results format

### API Reference
- [+smi_core Namespace](../api-reference/smi-core.md) - Core class reference
- [+smi Namespace](../api-reference/smi-namespace.md) - Main workflow classes
