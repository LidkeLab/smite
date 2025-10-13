---
title: "Compilation Reference"
category: "reference"
level: "advanced"
tags: ["compilation", "mex", "cuda", "build", "development", "nvcc", "compiler"]
prerequisites: ["../getting-started/installation.md", "../how-to/use-gpu.md"]
related: ["../troubleshooting/compilation-errors.md", "../troubleshooting/gpu-problems.md"]
summary: "Complete reference for compiling mex and CUDA files in smite"
estimated_time: "30-45 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Compilation Reference

## Purpose

This reference provides complete documentation for compiling smite's mex and CUDA files from source. smite includes precompiled binaries for all supported platforms (Linux, MacOS, Windows), so compilation is typically unnecessary. This guide is intended for developers modifying source code, users with platform-specific issues preventing use of precompiled binaries, or users needing to target specific GPU architectures.

## When Compilation is Required

**Important:** Compilation is only needed in specific circumstances:

1. **Source code modification**: You changed `.cpp` or `.cu` files in `MATLAB/source/`
2. **Precompiled files fail**: Platform-specific issues with included mex or PTX files
3. **GPU architecture targeting**: Your GPU architecture requires different compile flags
4. **Platform mismatch**: ARM vs x86_64 architecture differences
5. **Development work**: Contributing to smite development

**Before compiling:** Test whether smite works with precompiled files. Most users never need to compile.

## Overview of Compilation Process

### Mex Files

**Mex (MATLAB Executable)** files are compiled C++ code that MATLAB can call directly. smite uses mex for performance-critical operations:

- `smi_c_FrameConnection.cpp`: Frame-to-frame particle linking
- `c_HistRecon.cpp`: Histogram reconstruction for super-resolution images
- `c_HistImTime.cpp`: Time-series histogram generation
- `c_GenColorChannels.cpp`: Multi-color channel image generation
- `c_lap.cpp`: Linear assignment problem solver for tracking

**Output location:** `MATLAB/mex/`
**Output extensions:** `.mexa64` (Linux), `.mexmaci64` (MacOS), `.mexw64` (Windows)

### CUDA Files

**CUDA files** are GPU kernels compiled to PTX (Parallel Thread Execution) intermediate representation. smite uses CUDA for:

- `smi_cuda_gaussMLEv2.cu`: Maximum likelihood estimation for Gaussian fitting
- `smi_cuda_FindROI.cu`: Region of interest detection
- `smi_cuda_gaussBlobROIStack.cu`: Gaussian blob generation on ROI stacks
- `smi_cuda_PSFSample3DBlob.cu`: 3D point spread function sampling

**Output location:** `MATLAB/ptx/`
**Output extension:** `.ptx` (platform-independent intermediate code)

**PTX portability:** PTX files are designed to be architecture-independent. Compiling for compute capability 5.0 typically works on all newer GPUs due to forward compatibility, though newer GPUs may benefit from architecture-specific compilation.

## Mex Compilation

### Requirements

**All platforms:**
- MATLAB R2021a or later
- C/C++ compiler compatible with your MATLAB version
- Configured mex environment (`mex -setup` must succeed)

**Platform-specific compilers:**

**Linux:**
- GCC (GNU Compiler Collection)
- Typically pre-installed or available via package manager
- Supported versions depend on MATLAB release (see [compiler compatibility](https://www.mathworks.com/support/requirements/supported-compilers.html))

**MacOS:**
- Xcode Command Line Tools
- Clang/LLVM compiler
- Free download from Apple

**Windows:**
- Visual Studio (Community edition is free)
- OR MinGW-w64 C/C++ Compiler (MATLAB Add-On, mex-only, not for CUDA)
- Visual Studio required if also compiling CUDA files

### Compiler Setup

#### Linux Compiler Setup

**Ubuntu/Debian:**
```bash
# Install build tools
sudo apt-get update
sudo apt-get install build-essential

# Verify installation
gcc --version
g++ --version
```

**RHEL/CentOS/Fedora:**
```bash
# Install development tools
sudo yum groupinstall "Development Tools"

# Verify
gcc --version
g++ --version
```

**Configure mex:**
```matlab
% In MATLAB
mex -setup C++
% Expected output: "MEX configured to use 'g++' for C++ language compilation"
```

#### MacOS Compiler Setup

**Install Xcode Command Line Tools:**
```bash
# Terminal
xcode-select --install
```

This opens installation dialog. Accept and complete installation.

**Verify installation:**
```bash
gcc --version
clang --version
# Should show Apple clang version
```

**Configure mex:**
```matlab
% In MATLAB
mex -setup C++
% Expected: "MEX configured to use 'Xcode with Clang' for C++ language compilation"
```

**Apple Silicon (M1/M2/M3) considerations:**

MATLAB on Apple Silicon must run under Rosetta 2 for x86_64 compatibility:
```bash
# Install Rosetta 2 if not present
softwareupdate --install-rosetta

# Launch MATLAB under Rosetta
arch -x86_64 /Applications/MATLAB_R2021a.app/bin/matlab

# Or set MATLAB.app to always use Rosetta:
# Finder > Applications > Right-click MATLAB > Get Info > "Open using Rosetta"
```

#### Windows Compiler Setup

**Option 1: MinGW-w64 (mex only, easiest)**

Suitable if only compiling mex files (not CUDA):

```matlab
% In MATLAB
% Home tab > Add-Ons > Get Add-Ons
% Search: "MinGW-w64 Compiler"
% Click "Install"

% After installation
mex -setup C++
% Expected: "MEX configured to use 'MinGW64 Compiler (C++)' for C++ language compilation"
```

**Option 2: Visual Studio (required for CUDA)**

Required if compiling CUDA files:

1. Download Visual Studio Community (free): https://visualstudio.microsoft.com/
2. Run installer
3. Select "Desktop development with C++"
4. Ensure these components are checked:
   - MSVC v143 - VS 2022 C++ x64/x86 build tools (or appropriate version)
   - Windows 10 SDK (or Windows 11 SDK)
   - C++ CMake tools (optional but useful)
5. Complete installation (requires ~7 GB disk space)
6. Restart computer

**Configure mex:**
```matlab
% In MATLAB
mex -setup C++
% Expected: "MEX configured to use 'Microsoft Visual C++ 2022' for C++ language compilation"
```

### Compiler Version Compatibility

MATLAB versions have specific supported compiler versions. Using unsupported versions may produce warnings but often works.

**Check compatibility:**
- https://www.mathworks.com/support/requirements/supported-compilers.html
- Select your MATLAB version to see supported compilers

**Common compatibility:**
- MATLAB R2021a: GCC 8.x-9.x, Visual Studio 2017-2019, Xcode 11.x
- MATLAB R2022a: GCC 9.x-10.x, Visual Studio 2019-2022, Xcode 12.x-13.x
- MATLAB R2023a: GCC 10.x-11.x, Visual Studio 2019-2022, Xcode 13.x-14.x

**Handling version warnings:**

If you see warnings like:
```
Warning: You are using gcc version '11.0.0'. The version of gcc is not supported.
```

**Action:** Try compilation anyway. Warnings often don't prevent successful compilation. If compilation succeeds and mex files work, ignore the warning.

### Compilation Procedure

**Step 1: Navigate to source directory**

```matlab
% From MATLAB
cd ~/Documents/MATLAB/smite/MATLAB/source/c
% Or Windows
cd C:\Users\username\Documents\MATLAB\smite\MATLAB\source\c
```

**Step 2: Run mex_Make script**

```matlab
mex_Make
```

**What mex_Make does:**

The script compiles each C++ source file and outputs to `MATLAB/mex/`:

```matlab
% Simplified view of mex_Make.m

% Set up paths
baseFile = which('smi_core.FrameConnection');
[smite_CorePath] = fileparts(fileparts(baseFile));
[basePath] = fileparts(smite_CorePath);
sourcePath = fullfile(basePath,'source','c');
mexFilePath = fullfile(basePath,'mex');

% Compile each mex file
mex(fullfile(sourcePath,'smi_c_FrameConnection.cpp'), '-outdir', mexFilePath);
mex(fullfile(sourcePath,'c_HistRecon.cpp'), '-outdir', mexFilePath);
mex(fullfile(sourcePath,'c_HistImTime.cpp'), '-outdir', mexFilePath);
mex(fullfile(sourcePath,'c_GenColorChannels.cpp'), '-outdir', mexFilePath);
mex(fullfile(sourcePath,'c_lap.cpp'), '-outdir', mexFilePath);
```

**Expected output:**
```
Building with 'g++'.
MEX completed successfully.
Building with 'g++'.
MEX completed successfully.
[... etc for each file ...]
```

**Step 3: Verify compilation**

```matlab
% Check output directory
mex_dir = fullfile(smite_root, 'MATLAB', 'mex');
dir(fullfile(mex_dir, 'smi_c_FrameConnection.*'))
% Should show .mexa64 (Linux), .mexmaci64 (Mac), or .mexw64 (Windows)

% Test a mex function
smi_core.FrameConnection.unitTest()
```

### Manual Compilation

If `mex_Make` fails, compile individual files manually:

```matlab
% Navigate to source directory
cd ~/Documents/MATLAB/smite/MATLAB/source/c

% Compile single file
mex smi_c_FrameConnection.cpp -outdir ../../mex

% With explicit include paths (if needed)
mex -I"fullpath/to/matlab/extern/include" smi_c_FrameConnection.cpp -outdir ../../mex

% With optimization flags
mex -O smi_c_FrameConnection.cpp -outdir ../../mex

% With debugging symbols
mex -g smi_c_FrameConnection.cpp -outdir ../../mex
```

### Platform-Specific Mex Notes

**Linux:**
- Standard compilation usually works without modification
- Ensure write permissions on `MATLAB/mex/` directory
- GCC version warnings are common but usually harmless

**MacOS:**
- First run may require security approval: System Preferences > Security & Privacy > Allow
- Or use terminal: `xattr -d com.apple.quarantine ~/Documents/MATLAB/smite/MATLAB/mex/*.mexmaci64`
- Apple Silicon: Must run MATLAB under Rosetta 2

**Windows:**
- Long path names can cause issues; use short paths or enable long path support
- MinGW sufficient for mex; Visual Studio required only for CUDA
- Precompiled files usually work; compilation less frequently needed on Windows

## CUDA Compilation

### Requirements

**Hardware:**
- NVIDIA GPU (AMD, Intel, Apple GPUs not supported by MATLAB)
- GPU compute capability 5.0 or higher
- Sufficient GPU memory (2+ GB minimum)

**Software:**

**All platforms:**
- MATLAB R2021a or later with Parallel Computing Toolbox
- NVIDIA CUDA Toolkit (version compatible with MATLAB)
- nvcc compiler (included with CUDA Toolkit)

**Windows only:**
- Visual Studio with C++ compiler (cl.exe)
- Both nvcc AND cl.exe must be on system PATH

**Linux/Mac:**
- nvcc compiler on PATH
- Compatible with system GCC/Clang

### CUDA Toolkit Installation

**Check CUDA Toolkit compatibility:**

MATLAB versions support specific CUDA Toolkit versions. See: https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html

**Compatibility matrix:**
- MATLAB R2021a: CUDA 11.0 - 11.2
- MATLAB R2021b: CUDA 11.0 - 11.4
- MATLAB R2022a: CUDA 11.2 - 11.6
- MATLAB R2022b: CUDA 11.2 - 11.8
- MATLAB R2023a: CUDA 11.2 - 12.0
- MATLAB R2023b: CUDA 11.2 - 12.2

**Important:** Use CUDA version within your MATLAB's supported range.

#### Linux CUDA Installation

**Download CUDA Toolkit:**

Visit: https://developer.nvidia.com/cuda-downloads

Select:
- Linux
- Your distribution (Ubuntu, RHEL, etc.)
- Architecture (x86_64)
- Distribution version
- Installer type (runfile recommended for control)

**Example installation (Ubuntu):**
```bash
# Download installer (check website for current version)
wget https://developer.download.nvidia.com/compute/cuda/12.2.0/local_installers/cuda_12.2.0_535.54.03_linux.run

# Run installer
sudo sh cuda_12.2.0_535.54.03_linux.run

# Follow prompts:
# - Accept EULA
# - Deselect "Driver" if you already have recent NVIDIA driver
# - Select "CUDA Toolkit"
# - Accept default installation location: /usr/local/cuda-12.2
```

**Add to PATH:**

Edit `~/.bashrc` (or `~/.zshrc` for zsh):
```bash
# Add CUDA to PATH
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
```

Apply changes:
```bash
source ~/.bashrc
```

**Verify installation:**
```bash
nvcc --version
# Should show CUDA compilation tools version
```

#### MacOS CUDA Installation

**Important:** CUDA support on MacOS is limited and deprecated. Only older Intel Macs with NVIDIA GPUs support CUDA. Modern Macs (Apple Silicon M1/M2/M3) do not have NVIDIA GPUs and cannot use CUDA.

For Intel Macs with NVIDIA GPUs (rare):

**Download older CUDA Toolkit:**

Visit: https://developer.nvidia.com/cuda-toolkit-archive

Select CUDA 10.1 or 10.2 (later versions not available for Mac).

**Install from DMG:**
1. Open downloaded .dmg file
2. Run installer package
3. Follow installation prompts
4. Default location: `/Developer/NVIDIA/CUDA-10.1`

**Add to PATH:**

Edit `~/.bash_profile` or `~/.zshrc`:
```bash
export PATH=/Developer/NVIDIA/CUDA-10.1/bin:$PATH
export DYLD_LIBRARY_PATH=/Developer/NVIDIA/CUDA-10.1/lib:$DYLD_LIBRARY_PATH
```

Apply:
```bash
source ~/.bash_profile
```

**Verify:**
```bash
nvcc --version
```

#### Windows CUDA Installation

**Download CUDA Toolkit:**

Visit: https://developer.nvidia.com/cuda-downloads

Select:
- Windows
- x86_64
- Windows version (10 or 11)
- exe (local) recommended

**Install:**
1. Run installer executable
2. Choose "Custom" installation
3. Select components:
   - CUDA Toolkit (required)
   - CUDA Documentation (optional)
   - CUDA Samples (optional)
   - Driver components (only if updating driver)
4. Accept default installation path: `C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2`
5. Complete installation
6. Restart computer

**Verify nvcc installation:**
```cmd
:: Open Command Prompt
nvcc --version
```

If command not found, add to PATH manually:

**Add CUDA to System PATH:**
1. Control Panel > System and Security > System
2. Advanced system settings
3. Environment Variables
4. Under "System variables", select "Path", click Edit
5. Click New, add: `C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\bin`
6. Click OK, OK, OK
7. Restart Command Prompt and MATLAB

**Verify cl.exe (Visual Studio compiler):**

Windows requires both nvcc AND cl.exe for CUDA compilation:

```cmd
cl
:: Should show Microsoft C/C++ Optimizing Compiler version info
```

If not found:

**Add Visual Studio to PATH:**
1. Locate cl.exe in Visual Studio installation:
   ```
   C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.35.32215\bin\Hostx64\x64\cl.exe
   ```
   (Version number varies; search for `cl.exe`)

2. Add to System PATH (same process as CUDA above)

**Alternative: Use Visual Studio Command Prompt**

Start MATLAB from "x64 Native Tools Command Prompt for VS 2022":
- Start Menu > Visual Studio 2022 > x64 Native Tools Command Prompt
- From that prompt: `matlab`
- Both nvcc and cl.exe will be on PATH automatically

### CUDA Compilation Procedure

**Step 1: Configure PATH in cuda_Make.m**

Before first compilation, edit `cuda_Make.m` to set correct paths.

**Edit:** `~/Documents/MATLAB/smite/MATLAB/source/cuda/cuda_Make.m`

**Linux/Mac:**

Lines 30-34:
```matlab
else % Linux/MacOS
   % Prepend the system path with the NVIDIA GPU Computing Toolkit binaries
   % directory for compiling with nvcc.  Update with the correct path for
   % newer versions of the toolkit.  The toolkit is compatible with gcc.
   setenv('PATH', ['/usr/local/cuda/bin:' getenv('PATH')]);
end
```

**Update path** if your CUDA installation is not `/usr/local/cuda`:
```matlab
setenv('PATH', ['/usr/local/cuda-12.2/bin:' getenv('PATH')]);
```

**Windows:**

Lines 14-29:
```matlab
if ispc % Windows
   % Prepend the system path with the NVIDIA GPU Computing Toolkit binaries
   % directory for compiling with nvcc.  Prepend (rather than append) in case
   % multiple versions of the toolkit are installed (which automatically are
   % added to the system path on installation).  Update with the correct path
   % for newer versions of the toolkit.
   setenv('PATH', ['C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\bin;' getenv('PATH')]);

   % Prepend the system path with the Visual Studio binaries directory for
   % compiling with cl.  See comment above about prepend versus append.
   % Update with the correct path for newer versions of Visual Studio.
   % [VS2019]
   setenv('PATH', ['C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.29.30133\bin\Hostx64\x64;' getenv('PATH')]);
```

**Update paths** to match your installations:
- CUDA path: Update version number (`v12.2` → your version)
- Visual Studio path: Update year (2019 → 2022) and MSVC version number

**To find your MSVC version:**
```cmd
dir "C:\Program Files (x86)\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC"
:: Shows version number directories, e.g., 14.35.32215
```

**Step 2: Navigate to CUDA source directory**

```matlab
% From MATLAB
cd ~/Documents/MATLAB/smite/MATLAB/source/cuda
% Or Windows
cd C:\Users\username\Documents\MATLAB\smite\MATLAB\source\cuda
```

**Step 3: Run cuda_Make script**

```matlab
cuda_Make
```

**What cuda_Make does:**

For each CUDA source file:
1. Compiles `.cu` file to `.ptx` using nvcc
2. Copies `.ptx` to `MATLAB/ptx/`
3. Copies `.cu` source to `MATLAB/ptx/` (for reference)
4. Copies `.m` wrapper files to `MATLAB/ptx/` (if present)

**Compilation command used:**
```matlab
nvcc_cmd = 'nvcc -arch=sm_50 -ptx %s -o %s\n';
```

Flags:
- `-arch=sm_50`: Target compute capability 5.0 (Maxwell architecture and newer)
- `-ptx`: Output PTX intermediate representation (platform-independent)

**Expected output:**
```
Compiling smi_cuda_gaussMLEv2 ...
s = 0
r = [compilation output from nvcc]

Compiling smi_cuda_FindROI ...
s = 0
r = [compilation output]

Compiling smi_cuda_gaussBlobROIStack ...
s = 0
r = [compilation output]

Compiling smi_cuda_PSFSample3DBlob ...
s = 0
r = [compilation output]
```

`s = 0` indicates successful compilation (exit status 0).

**Step 4: Verify compilation**

```matlab
% Check PTX files created
ptx_dir = fullfile(smite_root, 'MATLAB', 'ptx');
dir(fullfile(ptx_dir, '*.ptx'))
% Should show:
% smi_cuda_gaussMLEv2.ptx
% smi_cuda_FindROI.ptx
% smi_cuda_gaussBlobROIStack.ptx
% smi_cuda_PSFSample3DBlob.ptx

% Test GPU functionality
gpu = gpuDevice;
fprintf('GPU: %s (compute capability %.1f)\n', gpu.Name, gpu.ComputeCapability);

% Run localization unit test
smi_core.LocalizeData.unitTest()
```

### Targeting Specific GPU Architectures

The default compilation targets compute capability 5.0 (`-arch=sm_50`), which provides forward compatibility with newer GPUs. However, targeting your specific GPU architecture can improve performance.

**Check your GPU compute capability:**
```matlab
gpu = gpuDevice;
fprintf('Compute capability: %.1f\n', gpu.ComputeCapability);
```

**Edit cuda_Make.m** line 38 to target your architecture:

**Single architecture (smaller file, faster compile):**
```matlab
% For Maxwell (compute capability 5.0)
nvcc_cmd = 'nvcc -arch=sm_50 -ptx %s -o %s\n';

% For Pascal (6.0)
nvcc_cmd = 'nvcc -arch=sm_60 -ptx %s -o %s\n';

% For Volta (7.0)
nvcc_cmd = 'nvcc -arch=sm_70 -ptx %s -o %s\n';

% For Turing (7.5)
nvcc_cmd = 'nvcc -arch=sm_75 -ptx %s -o %s\n';

% For Ampere (8.0)
nvcc_cmd = 'nvcc -arch=sm_80 -ptx %s -o %s\n';

% For Hopper (9.0)
nvcc_cmd = 'nvcc -arch=sm_90 -ptx %s -o %s\n';
```

**Multiple architectures (larger file, more compatible):**
```matlab
% Support both compute capability 5.0 and 8.0
nvcc_cmd = 'nvcc -gencode arch=compute_50,code=sm_50 -gencode arch=compute_80,code=sm_80 -ptx %s -o %s\n';
```

**Architecture reference:**

| Compute Capability | Architecture | Example GPUs |
|--------------------|--------------|--------------|
| 5.0 | Maxwell | GTX 750 Ti, GTX 980, Quadro M4000 |
| 6.0 | Pascal | GTX 1070/1080, Tesla P100 |
| 6.1 | Pascal | GTX 1050/1060, Titan X |
| 7.0 | Volta | Tesla V100, Titan V |
| 7.5 | Turing | RTX 2060/2070/2080, Quadro RTX 5000 |
| 8.0 | Ampere | RTX 3090, A100 |
| 8.6 | Ampere | RTX 3060/3070/3080 |
| 8.9 | Ada Lovelace | RTX 4090 |
| 9.0 | Hopper | H100 |

Check your GPU's compute capability: https://developer.nvidia.com/cuda-gpus

**Forward compatibility note:** PTX compiled for older compute capability works on newer GPUs (e.g., sm_50 works on sm_80 GPU), but not vice versa. Targeting newer architectures may enable additional GPU features.

### Manual CUDA Compilation

Compile individual CUDA files manually:

```bash
# From terminal (Linux/Mac)
cd ~/Documents/MATLAB/smite/MATLAB/source/cuda

# Compile single file
nvcc -arch=sm_50 -ptx smi_cuda_gaussMLEv2/smi_cuda_gaussMLEv2.cu -o ../../ptx/smi_cuda_gaussMLEv2.ptx

# Copy source to ptx directory
cp smi_cuda_gaussMLEv2/smi_cuda_gaussMLEv2.cu ../../ptx/
```

```cmd
:: From Command Prompt (Windows)
cd C:\Users\username\Documents\MATLAB\smite\MATLAB\source\cuda

:: Compile
nvcc -arch=sm_50 -ptx smi_cuda_gaussMLEv2\smi_cuda_gaussMLEv2.cu -o ..\..\ptx\smi_cuda_gaussMLEv2.ptx

:: Copy source
copy smi_cuda_gaussMLEv2\smi_cuda_gaussMLEv2.cu ..\..\ptx\
```

## Testing Compiled Files

### Test Mex Files

**Basic test - check file exists and loads:**
```matlab
% Navigate to mex directory
mex_dir = fullfile(smite_root, 'MATLAB', 'mex');

% Check smi_c_FrameConnection
mex_file = dir(fullfile(mex_dir, 'smi_c_FrameConnection.*'));
if ~isempty(mex_file)
    fprintf('Found: %s\n', mex_file.name);
else
    error('smi_c_FrameConnection mex file not found');
end

% Try to load (this succeeds if mex file is valid)
which smi_c_FrameConnection
% Should show path to mex file
```

**Functional test - run unit test:**
```matlab
% Test frame connection (uses smi_c_FrameConnection mex)
smi_core.FrameConnection.unitTest()
% Should complete without errors

% Test histogram reconstruction (uses c_HistRecon mex)
% Requires generating data first
smi_core.LoadData.unitTest()
```

**Test each mex file individually:**
```matlab
% c_lap (Linear Assignment Problem solver, used in tracking)
smi.SPT.unitTestFFGC()  % Tests tracking, which uses c_lap

% Other mex files tested through their respective classes
run_tests  % Comprehensive test suite
```

### Test CUDA Files

**Check GPU availability:**
```matlab
% Test GPU device
try
    gpu = gpuDevice;
    fprintf('GPU detected: %s\n', gpu.Name);
    fprintf('Compute capability: %.1f\n', gpu.ComputeCapability);
    fprintf('Driver version: %s\n', gpu.DriverVersion);
    fprintf('Toolkit version: %s\n', gpu.ToolkitVersion);
    fprintf('Available memory: %.2f GB\n', gpu.AvailableMemory / 1e9);
catch ME
    error('GPU not available: %s', ME.message);
end
```

**Test PTX files load:**
```matlab
% Check PTX files exist
ptx_dir = fullfile(smite_root, 'MATLAB', 'ptx');
ptx_files = dir(fullfile(ptx_dir, '*.ptx'));
fprintf('Found %d PTX files:\n', length(ptx_files));
for i = 1:length(ptx_files)
    fprintf('  %s\n', ptx_files(i).name);
end

% Test PTX loads (happens automatically when calling GPU functions)
% If PTX file is invalid, MATLAB will error when trying to use it
```

**Functional test - run GPU localization:**
```matlab
% Quick GPU test - localize synthetic data
smi_core.LocalizeData.unitTest()
% This test:
% 1. Generates synthetic data
% 2. Uses GPU for box finding (smi_cuda_FindROI)
% 3. Uses GPU for Gaussian fitting (smi_cuda_gaussMLEv2)
% 4. Validates results
% If successful, GPU compilation is working correctly

% Full SMLM pipeline test
smi.SMLM.unitTest()
% Tests complete SMLM workflow including all GPU operations
```

**Test specific GPU functions:**
```matlab
% Test gaussMLEv2 directly
% Create test data
testData = 100 * rand(11, 11, 100, 'single') + 10;  % 100 ROIs, 11x11 pixels
testData = gpuArray(testData);  % Move to GPU

% This will use smi_cuda_gaussMLEv2.ptx if it loads
% (Actual call requires proper setup through LocalizeData class)

% If errors occur, check:
% - PTX file exists
% - GPU has sufficient memory
% - Compute capability compatibility
```

### Comprehensive Test Suite

**Run full test suite:**
```matlab
% Complete unit test suite (tests all functionality including mex and CUDA)
run_tests

% Saves results to: tempdir/smite/unitTest/
% Check console output for PASSED/FAILED status
```

**Key tests that validate compilation:**
- `smi.SMLM.unitTest` - Tests SMLM pipeline (GPU localization, mex histogram)
- `smi.SPT.unitTestFFGC` - Tests tracking (c_lap mex for frame connection)
- `smi_core.LocalizeData.unitTest` - Tests GPU fitting algorithms
- `smi_core.FrameConnection.unitTest` - Tests smi_c_FrameConnection mex

## Troubleshooting

### Common Issues

**Mex compilation fails:**
- Check compiler is installed: `mex -setup C++`
- Verify compiler version compatibility with MATLAB
- Check write permissions on `MATLAB/mex/` directory
- See [Compilation Errors troubleshooting guide](../troubleshooting/compilation-errors.md)

**nvcc not found:**
- CUDA Toolkit not installed or not on PATH
- Edit `cuda_Make.m` to add CUDA bin directory to PATH
- Verify: `system('nvcc --version')` from MATLAB

**cl.exe not found (Windows):**
- Visual Studio not installed or cl.exe not on PATH
- Edit `cuda_Make.m` to add Visual Studio bin directory to PATH
- Or launch MATLAB from Visual Studio Command Prompt

**Compute capability mismatch:**
- Edit `cuda_Make.m` line 38 to target your GPU architecture
- Check GPU compute capability: `gpuDevice.ComputeCapability`
- Use `-arch=sm_XX` matching your GPU (e.g., `-arch=sm_80` for Ampere)

**CUDA version incompatible with MATLAB:**
- Check supported CUDA versions for your MATLAB release
- Install compatible CUDA Toolkit version
- See: https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html

**Compilation succeeds but runtime errors:**
- Test GPU availability: `gpuDevice`
- Check GPU memory: `gpuDevice.AvailableMemory`
- Verify PTX files in correct location: `dir('~/Documents/MATLAB/smite/MATLAB/ptx/*.ptx')`
- May need to target specific GPU architecture

### Getting Detailed Help

For compilation issues, collect diagnostic information:

```matlab
% MATLAB and platform info
version
computer
computer('arch')

% Compiler info
mex -setup C++

% GPU info (for CUDA issues)
try
    gpu = gpuDevice;
    fprintf('GPU: %s\n', gpu.Name);
    fprintf('Compute capability: %.1f\n', gpu.ComputeCapability);
    fprintf('Driver: %s\n', gpu.DriverVersion);
    fprintf('Toolkit: %s\n', gpu.ToolkitVersion);
catch ME
    fprintf('GPU error: %s\n', ME.message);
end

% PATH
getenv('PATH')
```

Report issues with this information at: https://github.com/LidkeLab/smite/issues

## See Also

- [Installation Guide](../getting-started/installation.md) - Initial smite setup
- [GPU Usage](../how-to/use-gpu.md) - GPU requirements and configuration
- [Compilation Errors](../troubleshooting/compilation-errors.md) - Detailed troubleshooting
- [GPU Problems](../troubleshooting/gpu-problems.md) - Runtime GPU issues
- [MATLAB Compiler Support](https://www.mathworks.com/support/requirements/supported-compilers.html)
- [NVIDIA CUDA Installation](https://docs.nvidia.com/cuda/)

## Summary

**Key points:**

1. **Compilation is rarely needed** - Precompiled files work for most users
2. **Mex compilation** requires C++ compiler configured via `mex -setup`
3. **CUDA compilation** requires NVIDIA CUDA Toolkit and (on Windows) Visual Studio
4. **Run `mex_Make`** in `source/c/` for mex files
5. **Run `cuda_Make`** in `source/cuda/` for CUDA files (after editing paths)
6. **Test compilation** with unit tests: `run_tests` or specific tests like `smi.SMLM.unitTest`
7. **PTX forward compatibility** - Compiling for older architecture works on newer GPUs
8. **Platform differences** - Linux easiest, Windows requires more setup, MacOS limited CUDA support

**When to compile:**
- Modifying source code
- Precompiled files don't work on your system
- Targeting specific GPU architecture for performance
- Platform or architecture mismatch

**When NOT to compile:**
- Fresh installation (try precompiled files first)
- smite working without errors
- No development work planned

Compilation enables development and customization but is not required for standard smite usage.
