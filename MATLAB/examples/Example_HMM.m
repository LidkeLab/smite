% This script demonstrates the basic usage of smi_stat.HMM.

%% Reload results from smi.SPT tracking and search for dimer candidates.
% NOTE: 'FileDir' is a directory containing both channel 1 and channel 2
%        results files.  findDimerCandidatesFromFiles() will search this
%        directory for files matching FilePatterns{1} and treat those files
%        as channel 1 results.  It will then do the same for
%        FilePatterns{2} (channel 2 results) and then attempt to match the
%        channel 1 and channel 2 files based on their file names (e.g.,
%        data01_Channel1_Results.mat would be paired to
%        data01_Channel2_Results.mat).
FileDir = 'C:\Users\David\Documents\MATLAB\spt_demos\HMM_demo';
FilePatterns = {'*Channel1_Results.mat'; '*Channel2_Results.mat'};
MaxDimerSeparation = 2; % pixels
MaxSeparation = 5; % pixels
[TRArray, FileList] = smi_stat.HMM.findDimerCandidatesFromFiles(...
    FileDir, FilePatterns, MaxDimerSeparation, MaxSeparation);

%% Prepare the HMM class and run the analysis.
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.FrameRate = 20; % fps, specific to the loaded TRs above
SMF.Data.PixelSize = 0.1; % micrometers, specific to the loaded TRs above
HMM = smi_stat.HMM(TRArray, SMF);
HMM.MaxSeparation = MaxSeparation;
HMM.DimerSeparation = 0.5; % pixels, TRUE physical separation between dimers
HMM.DiffusionCoefficient = 0.0615; % px^2 / frame, specific to the loaded TRs above
HMM.RegistrationError = 0; % pixels, specific to the loaded TRs above
HMM.SaveDir = fullfile(FileDir, 'HMM_Results');
HMM.GeneratePlots = true;
HMM.UnitFlag = true;
HMM.performFullAnalysis()
