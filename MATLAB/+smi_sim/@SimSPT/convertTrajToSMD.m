function [SMD] = convertTrajToSMD(TrajStruct, SimParams)
%convertTrajToSMD converts a TrajStruct to an SMD.
% This method reorganizes some vital fields in the input 'TrajStruct' and
% stores them into a Single Molecule Data structure 'SMD'.
%
% INPUTS:
%   TrajStruct: Structure containing trajectory data (see
%               smi_sim.SimSPT.simTrajectories())
%   SimParams: Structure of simulation parameters (see
%              smi_sim.SimSPT.defineDefaultParams()).
%
% OUTPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Pad the input 'SimParams'.
if (~exist('SimParams', 'var') || isempty(SimParams))
    SimParams = smi_sim.SimSPT.defineDefaultParams();
else
    DefaultParams = smi_sim.SimSPT.defineDefaultParams();
    SimParams = smi_helpers.padStruct(SimParams, DefaultParams);
end

% Prepare the SMD.
SMD = smi_core.SingleMoleculeData.createSMD();
SMD.NFrames = SimParams.NFrames; 
SMD.XSize = SimParams.FrameSize(2);
SMD.YSize = SimParams.FrameSize(1);
DataSize = size(TrajStruct.Photons);
NObservations = sum(TrajStruct.IsOn, 2);
SMD.ConnectID = repelem((1:DataSize(1)).', NObservations, 1);
ValidTrajInd = find(TrajStruct.IsOn.');
SMD.FrameNum = 1 + mod(ValidTrajInd-1, DataSize(2));
SMD.DatasetNum = ones(size(SMD.FrameNum));
SMD.PSFSigma = SimParams.PSFSigma * ones(size(SMD.FrameNum));
X = TrajStruct.Trajectories(:, :, 1).';
Y = TrajStruct.Trajectories(:, :, 2).';
SMD.X = X(ValidTrajInd);
SMD.Y = Y(ValidTrajInd);
X_SE = TrajStruct.Trajectories_SE(:, :, 1).';
Y_SE = TrajStruct.Trajectories_SE(:, :, 2).';
SMD.X_SE = X_SE(ValidTrajInd);
SMD.Y_SE = Y_SE(ValidTrajInd);
Photons = TrajStruct.Photons.';
SMD.Photons = Photons(ValidTrajInd);
Photons_SE = TrajStruct.Photons_SE.';
SMD.Photons_SE = Photons_SE(ValidTrajInd);
Bg = TrajStruct.Bg.';
SMD.Bg = Bg(ValidTrajInd);
Bg_SE = TrajStruct.Bg_SE.';
SMD.Bg_SE = Bg_SE(ValidTrajInd);


end