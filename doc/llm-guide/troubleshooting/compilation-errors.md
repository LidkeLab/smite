---
title: "Troubleshooting Compilation Errors"
category: "troubleshooting"
level: "advanced"
tags: ["compilation", "mex", "cuda", "nvcc", "errors", "troubleshooting"]
prerequisites: ["../getting-started/installation.md", "../how-to/use-gpu.md"]
related: ["gpu-problems.md", "installation-issues.md", "../reference/compilation.md"]
summary: "Comprehensive troubleshooting guide for mex and CUDA compilation errors in smite"
estimated_time: "20-30 minutes"
last_updated: "2025-01-10"
status: "complete"
---

# Troubleshooting Compilation Errors

## Purpose

This guide helps diagnose and resolve compilation errors when building smite's mex and CUDA files. Compilation is only needed when precompiled binaries fail or when modifying source code. This troubleshooting guide covers compiler setup issues, missing dependencies, platform-specific errors, and linking problems.

## Prerequisites

- Understanding of [installation process](../getting-started/installation.md)
- Basic command line knowledge
- Administrator/sudo access for installing compilers
- Familiarity with [GPU requirements](../how-to/use-gpu.md) for CUDA compilation

## When Compilation is Needed

**Important:** smite includes precompiled mex and CUDA files for Linux, MacOS, and Windows. Compilation is only required if:

1. **Precompiled files don't work** on your system
2. **You modified source code** in `source/c/` or `source/cuda/`
3. **Your GPU architecture is newer** than supported by precompiled PTX files
4. **Platform mismatch** (e.g., ARM vs x86_64)

**Before compiling:** Try using precompiled files first. If smite works without errors, compilation is unnecessary.

## Mex Compilation Errors

### Error: "mex: command not found" or "mex not recognized"

**Symptoms:**
```matlab
>> mex -setup
Unrecognized function or variable 'mex'.
```

Or from command line:
```bash
$ mex
mex: command not found
```

**Diagnosis:**

MATLAB is not properly installed, or mex command is not accessible. This indicates a fundamental MATLAB installation issue.

**Solution:**

```matlab
% 1. Check MATLAB is installed and running
version

% 2. Check if mex function exists
which mex -all

% 3. Check MATLAB bin directory is on system PATH
getenv('PATH')

% 4. Verify mex binary exists
% Linux/Mac: /usr/local/MATLAB/R2021a/bin/mex
% Windows: C:\Program Files\MATLAB\R2021a\bin\mex.bat
```

**If mex is not found:**

- **Linux/Mac:** Add MATLAB bin directory to PATH
  ```bash
  export PATH="/usr/local/MATLAB/R2021a/bin:$PATH"
  ```

- **Windows:** Add MATLAB bin directory via System Environment Variables
  - Control Panel > System > Advanced > Environment Variables
  - Add to PATH: `C:\Program Files\MATLAB\R2021a\bin`

- **Restart MATLAB** after modifying PATH

**Still not working:** Reinstall MATLAB ensuring mex utilities are selected during installation.

### Error: "No supported compiler was found"

**Symptoms:**
```matlab
>> mex -setup
Error using mex
No supported compiler was found. For options, visit
https://www.mathworks.com/support/requirements/supported-compilers.html
```

**Diagnosis:**

No C/C++ compiler is installed or configured for mex. This is the most common issue preventing mex compilation.

**Solution by Platform:**

#### Windows

**Option 1: MinGW-w64 (Recommended for mex only)**

```matlab
% Install MinGW-w64 Add-On via MATLAB
% Home > Add-Ons > Get Add-Ons
% Search for "MinGW-w64 Compiler"
% Install

% Then setup mex
mex -setup C++
```

**Option 2: Visual Studio (Required for CUDA compilation)**

1. Download Visual Studio Community (free): https://visualstudio.microsoft.com/
2. During installation, select "Desktop development with C++"
3. Install C++ compiler components
4. Restart MATLAB
5. Configure mex:
   ```matlab
   mex -setup C++
   ```

**Verify installation:**
```matlab
mex -setup C++
% Should show: "MEX configured to use 'Microsoft Visual C++ 2022' for C++ language compilation"
```

#### Linux

**Install GCC:**

```bash
# Ubuntu/Debian
sudo apt-get update
sudo apt-get install build-essential

# RHEL/CentOS/Fedora
sudo yum groupinstall "Development Tools"

# Verify installation
gcc --version
g++ --version
```

**Configure mex:**
```matlab
mex -setup C++
% Should show: "MEX configured to use 'g++' for C++ language compilation"
```

#### MacOS

**Install Xcode Command Line Tools:**

```bash
xcode-select --install
```

This opens a dialog to install developer tools. Accept and wait for installation.

**Verify installation:**
```bash
gcc --version
clang --version
```

**Configure mex:**
```matlab
mex -setup C++
% Should show: "MEX configured to use 'Xcode with Clang' for C++ language compilation"
```

### Error: Compiler version mismatch

**Symptoms:**
```matlab
>> mex_Make
Warning: You are using gcc version '11.0.0'. The version of gcc is not supported.
The version currently supported with MEX is '9.3.0'.
For a list of currently supported compilers see:
https://www.mathworks.com/support/compilers/current_release
```

**Diagnosis:**

Your compiler version is newer or older than officially supported by your MATLAB version. This is usually a warning, not a fatal error.

**Solution:**

**Option 1: Ignore warning and try compilation**

Often works despite warning:
```matlab
cd ~/Documents/MATLAB/smite/MATLAB/source/c
mex_Make  % Try anyway
```

If compilation succeeds, the mex files work despite the version mismatch warning.

**Option 2: Install supported compiler version**

```bash
# Linux: Install specific gcc version
sudo apt-get install gcc-9 g++-9

# Set as default for current session
export CC=gcc-9
export CXX=g++-9

# Then compile in MATLAB
```

**Option 3: Update MATLAB**

Newer MATLAB versions support newer compilers. Check compatibility: https://www.mathworks.com/support/requirements/supported-compilers.html

### Error: Missing headers during mex compilation

**Symptoms:**
```matlab
>> mex_Make
Error using mex
...smi_c_FrameConnection.cpp:8:10: fatal error: 'matrix.h' file not found
#include "matrix.h"
         ^~~~~~~~~~
```

**Diagnosis:**

MATLAB include directories are not accessible. This indicates mex configuration issue or incomplete MATLAB installation.

**Solution:**

```matlab
% 1. Verify MATLAB root
matlabroot

% 2. Check if extern/include exists
dir(fullfile(matlabroot, 'extern', 'include'))

% 3. Manually specify include path
mex_cmd = sprintf('mex -I"%s" smi_c_FrameConnection.cpp -outdir ../mex', ...
    fullfile(matlabroot, 'extern', 'include'));
eval(mex_cmd);
```

**If extern/include is missing:**

MATLAB installation is incomplete. Reinstall MATLAB with development libraries.

**Persistent issues:**

Edit `mex_Make.m` to explicitly add include paths:

```matlab
includePath = fullfile(matlabroot, 'extern', 'include');
mexFlags = sprintf('-I"%s"', includePath);
mex([mexFlags ' ' fullfile(sourcePath,'smi_c_FrameConnection.cpp')], ...
    '-outdir', mexFilePath);
```

### Error: Linking errors during mex compilation

**Symptoms:**
```matlab
>> mex_Make
Undefined symbols for architecture x86_64:
  "_mxCreateNumericArray", referenced from:
      frameConnect() in smi_c_FrameConnection.o
```

**Diagnosis:**

Linker cannot find MATLAB libraries. Common on Mac with multiple MATLAB versions or incorrect architecture.

**Solution:**

```matlab
% Check MATLAB architecture
computer('arch')
% Should match your system (maci64, glnxa64, win64)

% Verify library directory exists
libdir = fullfile(matlabroot, 'bin', computer('arch'))
dir(libdir)

% Manually specify library path
mex_cmd = sprintf('mex -L"%s" -lmx -lmex smi_c_FrameConnection.cpp -outdir ../mex', libdir);
eval(mex_cmd);
```

**MacOS specific - multiple MATLAB versions:**

```bash
# Check which MATLAB is active
which matlab

# Use specific MATLAB version
/Applications/MATLAB_R2021a.app/bin/matlab
```

**Architecture mismatch (Mac M1/M2/M3):**

MATLAB on Apple Silicon must run under Rosetta 2:
```bash
# Check if running under Rosetta
pgrep -lf MATLAB

# If native ARM: MATLAB needs to run under Rosetta
# Get Rosetta version of MATLAB or run with arch -x86_64
arch -x86_64 /Applications/MATLAB_R2021a.app/bin/matlab
```

## CUDA Compilation Errors

### Error: "nvcc: command not found" or "nvcc not recognized"

**Symptoms:**

From MATLAB:
```matlab
>> cuda_Make
...
[status, result] = system('nvcc -arch=sm_50 -ptx ...')
status = 1
result = 'nvcc' is not recognized as an internal or external command
```

From command line:
```bash
$ nvcc --version
nvcc: command not found
```

**Diagnosis:**

NVIDIA CUDA Toolkit is not installed, or nvcc compiler is not on system PATH.

**Solution:**

#### Step 1: Check CUDA Toolkit Installation

**Linux:**
```bash
# Check if CUDA is installed
ls /usr/local/cuda

# Check nvcc
ls /usr/local/cuda/bin/nvcc

# Check version
/usr/local/cuda/bin/nvcc --version
```

**Windows:**
```cmd
:: Check if CUDA is installed
dir "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA"

:: Check nvcc
dir "C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\bin\nvcc.exe"

:: Check version
"C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\bin\nvcc.exe" --version
```

**MacOS:**
```bash
# Check CUDA installation
ls /Developer/NVIDIA/CUDA-*/bin/nvcc
nvcc --version
```

#### Step 2: Install CUDA Toolkit (if not installed)

**Linux:**

Download from: https://developer.nvidia.com/cuda-downloads

```bash
# Ubuntu example (check website for your distribution)
wget https://developer.download.nvidia.com/compute/cuda/12.2.0/local_installers/cuda_12.2.0_535.54.03_linux.run
sudo sh cuda_12.2.0_535.54.03_linux.run

# Follow installer prompts
# Select: CUDA Toolkit
# Do NOT install driver if already have recent NVIDIA driver
```

**Windows:**

Download from: https://developer.nvidia.com/cuda-downloads

1. Run installer
2. Choose "Custom installation"
3. Select "CUDA Toolkit"
4. Complete installation
5. Restart computer

**MacOS:**

CUDA support on MacOS is deprecated. NVIDIA GPUs are not available in modern Macs (Apple Silicon). Only Intel Macs with NVIDIA GPUs can use CUDA.

Download from: https://developer.nvidia.com/cuda-downloads (older versions only)

#### Step 3: Add nvcc to PATH

**Linux:**

Edit `~/.bashrc` or `~/.zshrc`:
```bash
export PATH=/usr/local/cuda/bin:$PATH
export LD_LIBRARY_PATH=/usr/local/cuda/lib64:$LD_LIBRARY_PATH
```

Apply changes:
```bash
source ~/.bashrc
```

**Windows:**

Add to System PATH via Environment Variables:
1. Control Panel > System > Advanced System Settings
2. Environment Variables
3. System variables > PATH > Edit
4. Add: `C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\bin`
5. OK, OK, OK
6. Restart MATLAB

**MacOS:**

Edit `~/.bash_profile` or `~/.zshrc`:
```bash
export PATH=/Developer/NVIDIA/CUDA-10.1/bin:$PATH
export DYLD_LIBRARY_PATH=/Developer/NVIDIA/CUDA-10.1/lib:$DYLD_LIBRARY_PATH
```

Apply:
```bash
source ~/.bash_profile
```

#### Step 4: Verify PATH in MATLAB

```matlab
% Check PATH includes CUDA
getenv('PATH')

% Try nvcc from MATLAB
[status, result] = system('nvcc --version')
```

If still not found, manually set PATH in `cuda_Make.m`:

**Linux/Mac:**
```matlab
setenv('PATH', ['/usr/local/cuda/bin:' getenv('PATH')]);
```

**Windows:**
```matlab
setenv('PATH', ['C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\bin;' getenv('PATH')]);
```

### Error: "cl.exe not found" (Windows only)

**Symptoms:**
```matlab
>> cuda_Make
...
nvcc fatal   : Cannot find compiler 'cl.exe' in PATH
```

**Diagnosis:**

Windows requires Visual Studio C++ compiler (`cl.exe`) in addition to nvcc. CUDA uses Visual Studio's compiler as the host compiler.

**Solution:**

#### Step 1: Install Visual Studio

Download Visual Studio Community: https://visualstudio.microsoft.com/

During installation:
1. Select "Desktop development with C++"
2. Ensure "MSVC C++ build tools" is checked
3. Ensure "Windows SDK" is checked
4. Complete installation

#### Step 2: Locate cl.exe

```cmd
:: Common locations
dir "C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\*\bin\Hostx64\x64\cl.exe"
dir "C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\*\bin\Hostx64\x64\cl.exe"
```

#### Step 3: Add cl.exe to PATH

**Option 1: Use Visual Studio Command Prompt**

Start MATLAB from "x64 Native Tools Command Prompt for VS 2022":
1. Start Menu > Visual Studio 2022 > x64 Native Tools Command Prompt
2. From that prompt: `matlab`
3. Run `cuda_Make`

**Option 2: Modify cuda_Make.m**

Edit `smite/MATLAB/source/cuda/cuda_Make.m`, line ~29:

```matlab
% Update path to match your Visual Studio installation
setenv('PATH', ['C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.35.32215\bin\Hostx64\x64;' getenv('PATH')]);
```

Replace version number (14.35.32215) with your actual MSVC version from Step 2.

**Option 3: Add to System PATH**

Add cl.exe location to System PATH (like Step 3 for nvcc).

#### Step 4: Verify

```matlab
% Check PATH includes cl.exe
getenv('PATH')

% Test cl.exe
[status, result] = system('cl')
% Should show Microsoft C++ compiler version info
```

### Error: Compute capability mismatch

**Symptoms:**
```matlab
>> cuda_Make
nvcc warning : The 'compute_50' target architecture is deprecated and may be removed in a future release

% Or at runtime:
>> gpuDevice
Error: CUDA_ERROR_NO_BINARY_FOR_GPU: no kernel image available for execution on device
```

**Diagnosis:**

The PTX file was compiled for a different GPU architecture than your GPU supports. The compilation target (`-arch=sm_50`) may be incompatible with your GPU or CUDA version.

**Solution:**

#### Step 1: Check your GPU compute capability

```matlab
gpu = gpuDevice;
fprintf('Your GPU compute capability: %.1f\n', gpu.ComputeCapability);
```

#### Step 2: Edit cuda_Make.m to target your GPU

Open `smite/MATLAB/source/cuda/cuda_Make.m`, line ~38:

```matlab
% Original (targets compute capability 5.0+)
nvcc_cmd = 'nvcc -arch=sm_50 -ptx %s -o %s\n';

% For newer GPUs (Ampere, compute capability 8.0+):
nvcc_cmd = 'nvcc -arch=sm_80 -ptx %s -o %s\n';

% For even newer (Hopper, compute capability 9.0+):
nvcc_cmd = 'nvcc -arch=sm_90 -ptx %s -o %s\n';

% For multiple architectures (larger file but more compatible):
nvcc_cmd = 'nvcc -gencode arch=compute_50,code=sm_50 -gencode arch=compute_80,code=sm_80 -ptx %s -o %s\n';
```

**Architecture codes:**
- `sm_50`: Maxwell (GTX 750, GTX 980, Quadro M4000)
- `sm_60`: Pascal (GTX 1080, Tesla P100)
- `sm_70`: Volta (Tesla V100)
- `sm_75`: Turing (RTX 2080, Quadro RTX 5000)
- `sm_80`: Ampere (RTX 3090, A100)
- `sm_90`: Hopper (H100)

#### Step 3: Recompile

```matlab
cd ~/Documents/MATLAB/smite/MATLAB/source/cuda
cuda_Make
```

**Note:** PTX is designed to be forward-compatible. Compiling for `sm_50` usually works on newer GPUs, but not vice versa.

### Error: CUDA Toolkit version incompatible with MATLAB

**Symptoms:**
```matlab
>> cuda_Make
...compilation succeeds...

% But at runtime:
>> gpuDevice
Error using gpuDevice
The CUDA driver version is not compatible with the installed CUDA Toolkit version.
```

**Diagnosis:**

MATLAB's Parallel Computing Toolbox has specific CUDA version requirements. Your installed CUDA Toolkit may be too new or too old.

**Solution:**

#### Step 1: Check MATLAB's supported CUDA version

```matlab
% Check MATLAB version
version

% Get GPU device info
gpuDevice
% Shows ToolkitVersion supported by MATLAB
```

See: https://www.mathworks.com/help/parallel-computing/gpu-support-by-release.html

**MATLAB version compatibility:**
- R2021a: CUDA 11.0 - 11.2
- R2021b: CUDA 11.0 - 11.4
- R2022a: CUDA 11.2 - 11.6
- R2022b: CUDA 11.2 - 11.8
- R2023a: CUDA 11.2 - 12.0
- R2023b: CUDA 11.2 - 12.2

#### Step 2: Options

**Option 1: Install compatible CUDA version**

Download specific version from: https://developer.nvidia.com/cuda-toolkit-archive

**Option 2: Update MATLAB**

Download newer MATLAB that supports your CUDA version.

**Option 3: Use precompiled PTX files**

If you're compiling unnecessarily, use the precompiled PTX files included with smite. They work across CUDA versions.

### Error: Out of memory during compilation

**Symptoms:**
```matlab
>> cuda_Make
nvcc fatal   : Memory allocation error during compilation
```

Or:
```matlab
system error: Out of memory
```

**Diagnosis:**

GPU ran out of memory during compilation, or system RAM is insufficient. This is rare but can occur on systems with limited resources.

**Solution:**

```matlab
% 1. Close other applications using GPU
% Check GPU memory usage
gpu = gpuDevice;
fprintf('Available GPU memory: %.2f GB\n', gpu.AvailableMemory / 1e9);

% 2. Reset GPU
reset(gpu);

% 3. Close other MATLAB sessions

% 4. Try compilation again
cd ~/Documents/MATLAB/smite/MATLAB/source/cuda
cuda_Make
```

**If persistent:**

Compile one file at a time by commenting out sections in `cuda_Make.m`:

```matlab
% Compile only gaussMLEv2 first
%% smi_cuda_gaussMLEv2
cuda_dir = 'smi_cuda_gaussMLEv2';
fprintf('Compiling %s ...\n', cuda_dir);
% ... rest of section ...

% Comment out other sections temporarily
% %% smi_cuda_FindROI
% % ...
```

## Platform-Specific Issues

### Windows: Long path issues

**Symptoms:**
```matlab
>> mex_Make
Error: The system cannot find the path specified
```

**Solution:**

```matlab
% Use short paths
cd C:\Users\klidke\DOCUME~1\MATLAB\smite\MATLAB\source\c
mex_Make
```

Or enable long paths in Windows 10+:
1. Registry Editor: `HKEY_LOCAL_MACHINE\SYSTEM\CurrentControlSet\Control\FileSystem`
2. Set `LongPathsEnabled` to 1
3. Restart

### Linux: Permission denied

**Symptoms:**
```bash
$ cd ~/Documents/MATLAB/smite/MATLAB/source/c
$ matlab -r "mex_Make; quit"
Error: Permission denied writing to ../mex/smi_c_FrameConnection.mexa64
```

**Solution:**

```bash
# Check ownership of mex directory
ls -ld ~/Documents/MATLAB/smite/MATLAB/mex

# Fix permissions
chmod -R u+w ~/Documents/MATLAB/smite/MATLAB/mex

# Or run with sudo (not recommended)
sudo matlab -r "mex_Make; quit"
```

### MacOS: Code signing issues

**Symptoms:**

On first run of compiled mex file:
```
"smi_c_FrameConnection.mexmaci64" cannot be opened because the developer cannot be verified
```

**Solution:**

```bash
# Allow unsigned mex files
xattr -d com.apple.quarantine ~/Documents/MATLAB/smite/MATLAB/mex/*.mexmaci64

# Or allow in System Preferences:
# Security & Privacy > General > Allow apps downloaded from: Anywhere
```

### MacOS: Rosetta 2 required (Apple Silicon)

**Symptoms:**
```
Bad CPU type in executable
```

**Diagnosis:**

MATLAB on Apple Silicon (M1/M2/M3) must run under Rosetta 2 for x86_64 compatibility.

**Solution:**

```bash
# Install Rosetta 2 if not already installed
softwareupdate --install-rosetta

# Launch MATLAB under Rosetta
arch -x86_64 /Applications/MATLAB_R2021a.app/bin/matlab

# Or set MATLAB to always open with Rosetta:
# Finder > Applications > Right-click MATLAB > Get Info
# Check "Open using Rosetta"
```

**Note:** GPU operations are not available on Apple Silicon (no NVIDIA GPU). Compilation will succeed but GPU functions will fail at runtime.

## Verification After Compilation

### Verify Mex Files

```matlab
% Check mex files were created
mex_dir = fullfile(smite_root, 'MATLAB', 'mex');
dir(fullfile(mex_dir, 'smi_c_FrameConnection.*'))
dir(fullfile(mex_dir, 'c_HistRecon.*'))
dir(fullfile(mex_dir, 'c_lap.*'))

% Test a mex function
% (Requires data, see unit tests)
smi_core.FrameConnection.unitTest()
```

### Verify CUDA Files

```matlab
% Check PTX files were created
ptx_dir = fullfile(smite_root, 'MATLAB', 'ptx');
dir(fullfile(ptx_dir, 'smi_cuda_gaussMLEv2.ptx'))
dir(fullfile(ptx_dir, 'smi_cuda_FindROI.ptx'))

% Test GPU functions
gpu = gpuDevice;
fprintf('GPU: %s (compute capability %.1f)\n', gpu.Name, gpu.ComputeCapability);

% Quick GPU localization test
smi_core.LocalizeData.unitTest()
```

### Full Test Suite

```matlab
% Run all unit tests
run_tests

% Or test specific components
smi.SMLM.unitTest()
smi.SPT.unitTestFFGC()
```

## Getting Help

If compilation issues persist after trying these solutions:

### Collect Diagnostic Information

```matlab
% MATLAB version
version

% Platform info
computer
computer('arch')

% Compiler info
mex -setup C++

% GPU info (if applicable)
try
    gpu = gpuDevice;
    fprintf('GPU: %s\n', gpu.Name);
    fprintf('Compute capability: %.1f\n', gpu.ComputeCapability);
    fprintf('Driver: %s\n', gpu.DriverVersion);
    fprintf('Toolkit: %s\n', gpu.ToolkitVersion);
catch
    fprintf('No GPU available\n');
end

% PATH
getenv('PATH')

% Compiler versions (run in terminal)
% Linux/Mac: gcc --version, nvcc --version
% Windows: cl, nvcc --version
```

### Report Issue

GitHub Issues: https://github.com/LidkeLab/smite/issues

Include:
1. Operating system and version
2. MATLAB version
3. Compiler type and version
4. GPU model and CUDA version (for CUDA issues)
5. Complete error message
6. Output of diagnostic commands above

## See Also

- [Installation Guide](../getting-started/installation.md) - Initial setup
- [GPU Usage Guide](../how-to/use-gpu.md) - GPU requirements and configuration
- [GPU Problems](gpu-problems.md) - Runtime GPU issues
- [Installation Issues](installation-issues.md) - General installation troubleshooting
- [doc/mex+CUDA.md](../../mex+CUDA.md) - Compilation reference documentation
- [MATLAB Compiler Support](https://www.mathworks.com/support/requirements/supported-compilers.html)
- [NVIDIA CUDA Documentation](https://docs.nvidia.com/cuda/)

## Summary

**Most common issues:**

1. **Mex compiler not found** → Install GCC/Visual Studio/Xcode
2. **nvcc not found** → Install CUDA Toolkit and add to PATH
3. **cl.exe not found (Windows)** → Install Visual Studio, add to PATH
4. **Compute capability mismatch** → Edit `cuda_Make.m` to target your GPU architecture
5. **CUDA version incompatible** → Install CUDA version compatible with your MATLAB

**Remember:** Compilation is usually unnecessary. Try precompiled files first. Only compile if errors occur or you modify source code.
