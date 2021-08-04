% This script demonstrates the basic usage of smi_stat.HMM.

%% Reload some results from smi.SPT tracking.
% Select the desired tracking results.
FileListChannel1 = uipickfiles('Prompt', 'Pick the channel 1 files');
FileListChannel2 = uipickfiles('Prompt', 'Pick the channel 2 files');

% Match the files based on their time stamps.
% NOTE: The channel 1 and channel 2 tracking results differ only by the
%       filename tags 'Channel1' and 'Channel2', so we can pair the
%       selected files by just removing those strings.
IgnoreText = {'Channel1', 'Channel2'};
[PairedChannel1, PairedChannel2] = ...
    smi_helpers.pairText(FileListChannel1, FileListChannel2, IgnoreText);

%% Isolate dimer candidate events from the TR structures.
MaxDimerSeparation = 2; % pixels
MaxSeparation = 5; % pixels
TRArray = struct([]);
for ff = 1:numel(PairedChannel1)
    % Load the channel 1 TR structure.
    load(PairedChannel1{ff}, 'TR')
    TR1 = TR;
    
    % Load the channel 2 TR structure.
    load(PairedChannel2{ff}, 'TR')
    TR2 = TR;
    
    % Find dimerization event candidates between TR1 and TR2 trajectories
    % and add those candidates to the concatenated 'TRArray'.
    TRArray = [TRArray; smi_stat.HMM.findDimerCandidates(TR1, TR2, ...
        MaxDimerSeparation, MaxSeparation)];
end

%% Prepare the HMM class and run the analysis.
SMF = smi_core.SingleMoleculeFitting;
SMF.Data.FrameRate = 20; % fps, specific to the loaded TRs above
SMF.Data.PixelSize = 0.1; % micrometers, specific to the loaded TRs above
HMM = smi_stat.HMM(TRArray, SMF);
HMM.MaxSeparation = MaxSeparation;
HMM.DimerSeparation = 0.5; % pixels, TRUE physical separation between dimers
HMM.DiffusionCoefficient = 0.0615; % px^2 / frame, specific to the loaded TRs above
HMM.RegistrationError = 0; % pixels, specific to the loaded TRs above
HMM.SaveDir = 'C:\Users\David\Documents\MATLAB\spt_demos\HMM_demo\smite_test';
HMM.GeneratePlots = true;
HMM.UnitFlag = true;
HMM.performFullAnalysis()