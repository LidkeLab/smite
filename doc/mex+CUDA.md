# Compiling ***smite*** mex and CUDA files

## Prerequisites

To compile ***smite*** mex files, mex must have been set up in MATLAB with
an installed C/C++ compiler (type `mex -setup` to MATLAB to see if this is
the case).
To compile CUDA files, the NVIDIA GPU Computing Toolkit needs to be installed,
with the `nvcc` compiler located on the file PATH.  (For Windows, `cl.exe`
from a Visual Studio installation [freely available] also needs to be on
the file PATH.)
See
[Linux](https://docs.nvidia.com/cuda/cuda-installation-guide-linux/) and
[Windows](https://docs.nvidia.com/cuda/cuda-installation-guide-microsoft-windows/)
documentation for installing NVIDIA CUDA.
Note that due to forward compatibility rules, if the GPU architecture is
newer than what is
[supported by the MATLAB version](https://www.mathworks.com/help/releases/R2021b/parallel-computing/gpu-support-by-release.html),
CUDA compilations may fail.

## mex

To compile the mex files in `~/Documents/MATLAB/smite/MATLAB/source/c`
(for example, if they have been modified), run `mex_Make` in the
source directory:
```
   cd ~/Documents/MATLAB/smite/MATLAB/source/c
   mex_Make
```
Here, mex must have been set up (check in MATLAB with `mex -setup`)
to handle C/C++ integration.  The compiled files will be placed in the
`smite/MATLAB/mex` directory with extension mexa64 (Linux), mexmaci64
(MacOS) or mexw64 (Windows), assuming a 64-bit machine.  NOTE:
compilation should only be needed if source files change or the current
mex files in the mex directory fail to work properly.

## CUDA

Simularly, to compile the CUDA files in
`~/Documents/MATLAB/smite/MATLAB/source/cuda`, run `cuda_Make` in the
source directory:
```
   cd ~/Documents/MATLAB/smite/MATLAB/source/cuda
   cuda_Make
```
To generate CUDA files, the `nvcc` compiler (part of the NVIDIA GPU
Computing Toolkit installation), and in addition for Windows, the `cl.exe`
compiler from Visual Studio, must be on the user's `PATH` (see the
`setenv` lines at the top of the
[cuda_Make.m](../MATLAB/source/cuda/cuda_Make.m) file).  The user may need
to modify the line(s) here to `cuda_Make` of the form:
```
   setenv('PATH', [getenv('PATH') ':nvcc_PATH']);
```
(MacOS/Linux), or
```
   setenv('PATH', [getenv('PATH') ';nvcc_PATH']);
   setenv('PATH', [getenv('PATH') ';cl_PATH']);
```
(Windows)
where `nvcc_PATH` should be replaced by the actual path to where the
`nvcc` binary resides and cl_PATH should be replaced by the actual path
to where the Visual Studio `cl.exe` compiler resides.  The compiled
files will be placed in the `smite/MATLAB/ptx` directory.

NOTES: ptx is is designed to be independent of hardware architecture.
Compilation should only be needed if source files change or the current
ptx files in the ptx directory fail to work properly.  Also, MATLAB
does not support machines (like Macs) with non-NVIDIA GPUs when using
the Parallel Computing Toolbox (which ***smite*** uses).
