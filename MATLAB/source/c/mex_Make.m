%% This acts as a make file for mex functions. 

% This should be run to generate mex files for each architecture.
% mex_Make should be run while in the smite/MATLAB/source/c directory.

%% get path

baseFile = which('smi_core.FrameConnection');
[smite_CorePath] = fileparts(fileparts(baseFile));
[basePath] = fileparts(smite_CorePath);
sourcePath = fullfile(basePath,'source','c');
mexFilePath = fullfile(basePath,'mex');

%% c_FrameConnect

mex(fullfile(sourcePath,'smi_c_FrameConnection.cpp'), '-outdir', mexFilePath);
 
%% c_HistRecon

mex(fullfile(sourcePath,'c_HistRecon.cpp'), '-outdir', mexFilePath);
%copyfile(fullfile(sourcePath,'c_HistRecon.m'),mexFilePath);

%% c_HistImTime

mex(fullfile(sourcePath,'c_HistImTime.cpp'), '-outdir', mexFilePath);
%copyfile(fullfile(sourcePath,'c_HistImTime.m'),mexFilePath);

%% c_GenColorChannels

mex(fullfile(sourcePath,'c_GenColorChannels.cpp'), '-outdir', mexFilePath);
%copyfile(fullfile(sourcePath,'c_GenColorChannels.m'),mexFilePath);

%% c_lap

mex(fullfile(sourcePath,'c_lap.cpp'), '-outdir', mexFilePath);
