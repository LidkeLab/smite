function [SMD] = convertTrajToSMD(TrajStruct)
%convertTrajToSMD converts a TrajStruct to an SMD.
% This method reorganizes some vital fields in the input 'TrajStruct' and
% stores them into a Single Molecule Data structure 'SMD'.
%
% INPUTS:
%   TrajStruct: Structure containing trajectory data (see
%               smi_sim.SimSPT.simTrajectories())
%
% OUTPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Prepare the SMD.
SMD = smi_core.SingleMoleculeData.createSMD();
DataSize = size(TrajStruct.Photons);
NObservations = sum(TrajStruct.IsOn, 2);
SMD.ConnectID = repelem((1:DataSize(1)).', NObservations);
ValidTrajInd = find(TrajStruct.IsOn.');
SMD.FrameNum = mod(ValidTrajInd, DataSize(2));
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