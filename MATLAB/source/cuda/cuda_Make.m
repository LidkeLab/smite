%% This acts as a Makefile for the ptx files.

% This only needs to be run once.  ptx and cu files will be saved to the
% SMA/ptx directory.
%
% IMPORTANT:
%    cuda_Make MUST be run while in the SMA/source/cuda directory.
%
% REQUIREMENTS:
%    Need to have a CUDA toolkit and VS2013 installed (Windows).

clc

if ispc
   % Adding system path for nvcc to compile with nvcc
   setenv('PATH', [getenv('PATH') ';C:\Program Files\NVIDIA GPU Computing Toolkit\CUDA\v11.5\bin']);
   % Adding system path for VS2013 to compile with cl
   % setenv('PATH', [getenv('PATH') ';C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin']);
   setenv('PATH', [getenv('PATH') ';C:\Program Files\Microsoft Visual Studio\2022\Community\VC\Tools\MSVC\14.30.30705\bin\Hostx64\x64'])     
else % Linux/MacOS
   % Adding system path for nvcc to compile with nvcc
   setenv('PATH', [getenv('PATH') ':/usr/local/cuda-10.1/bin']);
end

%% smi_cuda_gaussMLEv2
clc
cuda_dir = 'smi_cuda_gaussMLEv2';
fprintf('Compiling %s ...\n', cuda_dir);
addpath(cuda_dir);

%note -allow-unsupported-compiler is used to allow 2022 version.  DANGER!
[s, r] = system(sprintf('nvcc -allow-unsupported-compiler -ptx %s -o %s\n',...
                        fullfile(cuda_dir, [cuda_dir, '.cu']), ...
                        fullfile('..', '..', 'ptx', [cuda_dir, '.ptx'])))
copyfile(fullfile(cuda_dir, [cuda_dir, '.cu']), fullfile('..', '..', 'ptx'));

%% smi_cuda_FindROI

cuda_dir = 'smi_cuda_FindROI';
fprintf('Compiling %s ...\n', cuda_dir);
addpath(cuda_dir);

[s, r] = system(sprintf('nvcc -allow-unsupported-compiler -ptx %s -o %s\n',                ...
                        fullfile(cuda_dir, [cuda_dir, '.cu']), ...
                        fullfile('..', '..', 'ptx', [cuda_dir, '.ptx'])))
copyfile(fullfile(cuda_dir, [cuda_dir, '.cu']), fullfile('..', '..', 'ptx'));

%% smi_cuda_gaussBlobROIStack

cuda_dir = 'smi_cuda_gaussBlobROIStack';
fprintf('Compiling %s ...\n', cuda_dir);
addpath(cuda_dir);

[s, r] = system(sprintf('nvcc -allow-unsupported-compiler -ptx %s -o %s\n',                ...
                        fullfile(cuda_dir, [cuda_dir, '.cu']), ...
                        fullfile('..', '..', 'ptx', [cuda_dir, '.ptx'])));
copyfile(fullfile(cuda_dir, [cuda_dir, '.cu']), fullfile('..', '..', 'ptx'));
copyfile(fullfile(cuda_dir, [cuda_dir, '.m']),  fullfile('..', '..', 'ptx'));

%% smi_cuda_PSFSample3DBlob

cuda_dir = 'smi_cuda_PSFSample3DBlob';
fprintf('Compiling %s ...\n', cuda_dir);
addpath(cuda_dir);

[s, r] = system(sprintf('nvcc -allow-unsupported-compiler -ptx %s -o %s\n',                ...
                        fullfile(cuda_dir, [cuda_dir, '.cu']), ...
                        fullfile('..', '..', 'ptx', [cuda_dir, '.ptx'])))
                    
copyfile(fullfile(cuda_dir, [cuda_dir, '.cu']), fullfile('..', '..', 'ptx'));
