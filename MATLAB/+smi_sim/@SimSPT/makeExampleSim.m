function [Files, SimParams, DataParams] = makeExampleSim(...
    SimParams, DataParams, SaveDir)
%makeExampleSim prepares a basic two channel SPT simulation.
% This method prepares a basic example simulation of two-channel SPT data.
%
% INPUTS:
%   SimParams: Structure of simulation parameters (see
%              SimSPT.defineDefaultParams()).
%   DataParams: Structure of additional data related parameters, such as 
%               number of simulated datasets and noise characteristics, not
%               covered by 'SimParams'.
%               NDatasets: Number of simulated datasets. (Default = 10)
%               Background: Background photon counts. (Default = 0)
%   SaveDir: Base directory in which simulated data will be saved. 
%            (Default = pwd())
%
% OUTPUTS:
%   Files: Cell array of the filenames saved in 'SaveDir'.  Data are saved
%          as side-by-side images, e.g., for a 128x128 ROI size in each
%          channel, saved data will be 128x256 pixels, with channel 1 in
%          columns 1:128 and channel 2 in columns 129:256.

% Created by:
%   David J. Schodt (Lidke Lab, 2023)


% Set defaults.
if (~exist('SimParams', 'var') || isempty(SimParams))
    SimParams = smi_sim.SimSPT.defineDefaultParams();
else
    SimParams = smi_helpers.padStruct(SimParams, ...
        smi_sim.SimSPT.defineDefaultParams());
end
if (~exist('SaveDir', 'var') || isempty(SaveDir))
    SaveDir = pwd();
end
if ~isfolder(SaveDir)
    mkdir(SaveDir)
end
DefaultDataParams = struct();
DefaultDataParams.NDatasets = 10;
DefaultDataParams.Background = 0;
if (~exist('DataParams', 'var') || isempty(DataParams))
    DataParams = DefaultDataParams;
else
    DataParams = smi_helpers.padStruct(DataParams, DefaultDataParams);
end

% Loop over datasets and prepare the NDatasets simulations.
Files = cell(DataParams.NDatasets, 1);
for nn = 1:DataParams.NDatasets
    % Simulate trajectories.
    SPTSim = smi_sim.SimSPT;
    SPTSim.SimParams = SimParams;
    SPTSim.createSimulation()

    % Split the simulated trajectories randomly into two channels.
    NTraj = numel(SPTSim.TR);
    Inds = (1:NTraj);
    Channel1Inds = randperm(NTraj, ceil(NTraj/2));
    Channel2Inds = setdiff(Inds, Channel1Inds);
    SMD1 = smi_core.SingleMoleculeData.isolateSubSMD(...
        SPTSim.SMD, ismember(SPTSim.SMD.ConnectID, Channel1Inds));
    SMD2 = smi_core.SingleMoleculeData.isolateSubSMD(...
        SPTSim.SMD, ismember(SPTSim.SMD.ConnectID, Channel2Inds));

    % Simulate raw data.
    SMF = smi_core.SingleMoleculeFitting;
    SMF.Data.DataROI = [1, 1, SPTSim.SMD.YSize, SPTSim.SMD.XSize];
    SMF.Fitting.PSFSigma = 1.3;
    [~, RawData1] = smi_sim.GaussBlobs.gaussBlobImage(...
        SMD1, SMF, DefaultDataParams.Background);
    [~, RawData2] = smi_sim.GaussBlobs.gaussBlobImage(...
        SMD2, SMF, DefaultDataParams.Background);
    sequence = [RawData1, RawData2];
    pause(1.1) % pause before next iteration so time stamp is incremented
    Files{nn} = sprintf('Data_%03i_%s.mat', nn, smi_helpers.genTimeString());

    % Save the simulated data.
    save(fullfile(SaveDir, Files{nn}), 'sequence')
end


end
