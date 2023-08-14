# Compiling ***smite** mex and CUDA files

To compile the mex files in `~/Documents/MATLAB/smite/MATLAB/source/c`
(for example, if they have been modified), use `mex_Make` in the
source directory:
```
   cd ~/Documents/MATLAB/smite/MATLAB/source/c
   mex_Make
```
Here, mex must have been initially set up (via `mex -setup C++`) to handle
C++ integration.

Simularly, to compile the CUDA files in
`~/Documents/MATLAB/smite/MATLAB/source/cuda`, use `cuda_Make` in th3
source directory:
```
   cd ~/Documents/MATLAB/smite/MATLAB/source/cuda
   cuda_Make
```
To generate CUDA files, the `nvcc` compiler must be on the user's PATH
(see the lines at the top of the `cuda_Make.m` file).  The user may need
to add a line here to `cuda_Make` of the form:
```
   setenv('PATH', [getenv('PATH') ':nvcc_PATH']);
```
(MacOS/Linux)
```
   setenv('PATH', [getenv('PATH') ';nvcc_PATH']);
```
(Windows)
where `nvcc_PATH` is the actual path to where the `nvcc` binary resides.

The compiled files will be placed in the `smite/MATLAB/mex` or
`smite/MATLAB/ptx` directories, respectively.
