%% This acts as a Makefile for the ptx files.

% This only needs to be run once.  ptx and cu files will be saved to the
% smite/MATLAB/ptx directory.
%
% IMPORTANT:
%    cuda_Make MUST be run while in the smite/MATLAB/source/cuda directory.
%
% REQUIREMENTS:
%    Need to have a CUDA toolkit and VS2013 installed (Windows).

clc

if ispc % Windows
   % Prepend the system path with the NVIDIA GPU Computing Toolkit binaries
   % directory for compiling with nvcc.  Prepend (rather than append) in case
   % multiple versions of the toolkit are installed (which automatically are
   % added to the system path on installation).  Update with the correct path
   % for newer versions of the toolkit.
%  setenv('PATH', ['C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v10.1\bin;' getenv('PATH')]);
   setenv('PATH', ['C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v12.2\bin;' getenv('PATH')]);

   % Prepend the system path with the Visual Studio binaries directory for
   % compiling with cl.  See comment above about prepend versus append.
   % Update with the correct path for newer versions of Visual Studio.
   % [VS2013]
%  setenv('PATH', ['C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin;' getenv('PATH')]);
   % [VS2019]
   setenv('PATH', [';C:\Program Files (x86)\Microsoft Visual Studio\2019\Community\VC\Tools\MSVC\14.29.30133\bin\Hostx64\x64;' getenv('PATH')]);
else % Linux/MacOS
   % Prepend the system path with the NVIDIA GPU Computing Toolkit binaries
   % directory for compiling with nvcc.  Update with the correct path for
   % newer versions of the toolkit.  The toolkit is compatible with gcc.
   setenv('PATH', ['/usr/local/cuda-10.1/bin:' getenv('PATH')]);
end

% -arch=sm_50 says build for GPUs with computing capability at least 5.0 .
nvcc_cmd = 'nvcc -arch=sm_50 -ptx %s -o %s\n';

%% smi_cuda_gaussMLEv2
clc
cuda_dir = 'smi_cuda_gaussMLEv2';
fprintf('Compiling %s ...\n', cuda_dir);
addpath(cuda_dir);

[s, r] = system(sprintf(nvcc_cmd, ...
                        fullfile(cuda_dir, [cuda_dir, '.cu']), ...
                        fullfile('..', '..', 'ptx', [cuda_dir, '.ptx'])))
copyfile(fullfile(cuda_dir, [cuda_dir, '.cu']), fullfile('..', '..', 'ptx'));

%% smi_cuda_FindROI

cuda_dir = 'smi_cuda_FindROI';
fprintf('Compiling %s ...\n', cuda_dir);
addpath(cuda_dir);

[s, r] = system(sprintf(nvcc_cmd, ...
                        fullfile(cuda_dir, [cuda_dir, '.cu']), ...
                        fullfile('..', '..', 'ptx', [cuda_dir, '.ptx'])))
copyfile(fullfile(cuda_dir, [cuda_dir, '.cu']), fullfile('..', '..', 'ptx'));

%% smi_cuda_gaussBlobROIStack

cuda_dir = 'smi_cuda_gaussBlobROIStack';
fprintf('Compiling %s ...\n', cuda_dir);
addpath(cuda_dir);

[s, r] = system(sprintf(nvcc_cmd, ...
                        fullfile(cuda_dir, [cuda_dir, '.cu']), ...
                        fullfile('..', '..', 'ptx', [cuda_dir, '.ptx'])))
copyfile(fullfile(cuda_dir, [cuda_dir, '.cu']), fullfile('..', '..', 'ptx'));
copyfile(fullfile(cuda_dir, [cuda_dir, '.m']),  fullfile('..', '..', 'ptx'));

%% smi_cuda_PSFSample3DBlob

cuda_dir = 'smi_cuda_PSFSample3DBlob';
fprintf('Compiling %s ...\n', cuda_dir);
addpath(cuda_dir);

[s, r] = system(sprintf(nvcc_cmd, ...
                        fullfile(cuda_dir, [cuda_dir, '.cu']), ...
                        fullfile('..', '..', 'ptx', [cuda_dir, '.ptx'])))
                    
copyfile(fullfile(cuda_dir, [cuda_dir, '.cu']), fullfile('..', '..', 'ptx'));
