# Compiling ***smite*** mex and CUDA files

## mex

To compile the mex files in `~/Documents/MATLAB/smite/MATLAB/source/c`
(for example, if they have been modified), run `mex_Make` in the
source directory:
```
   cd ~/Documents/MATLAB/smite/MATLAB/source/c
   mex_Make
```
Here, mex must have been initially set up (via `mex -setup C++`)
to handle C++ integration.  The compiled files will be placed in the
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
To generate CUDA files, the `nvcc` compiler (part of the NVIDIA CUDA
installation) must be on the user's `PATH` (see the `setenv` lines at the
top of the `cuda_Make.m` file).  The user may need to add a line here to
`cuda_Make` of the form:
```
   setenv('PATH', [getenv('PATH') ':nvcc_PATH']);
```
(MacOS/Linux), or
```
   setenv('PATH', [getenv('PATH') ';nvcc_PATH']);
```
(Windows)
where `nvcc_PATH` should be replaced by the actual path to where the
`nvcc` binary resides.  The compiled files will be placed in the
`smite/MATLAB/ptx` directory.

NOTES: ptx is is designed to be independent of hardware architecture.
Compilation should only be needed if source files change or the current
mex files in the ptx directory fail to work properly.  Also, MATLAB
does not support machines (like Macs) with non-NVIDIA GPUs when using
the Parallel Computing Toolbox.
