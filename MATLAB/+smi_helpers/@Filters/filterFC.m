function SMD = filterFC(SMD, nFC)
%filterFC: Filter out localizations representing nFC or fewer frame connections.
%
% INPUTS:
%    SMD   Single Molecule Data structure
%    nFC   minimum number of frame connected localizations representing a
%          single localizaation allowed [DEFAULT = 1]
%
% OUTPUT:
%    SMD   modified Single Molecule Data structure

% Created by
%    David J. Schodt and Michael J. Wester (5/24/2022)

if ~exist('nFC', 'var')
   nFC = 1;
end

% Perform frame connection
SMF = smi_core.SingleMoleculeFitting();
SMF.FrameConnection.Method = 'Hypothesis test';
SMF.FrameConnection.LoS = 0.01;
SMF.FrameConnection.MaxSeparation = 1; % pixels
SMF.FrameConnection.MaxFrameGap = 5; % frames
FC = smi_core.FrameConnection(SMD, SMF);
FC.performFrameConnection();
SMD = FC.SMDCombined;

% Remove localizations representing fewer than nFC frame-connected
% localizations.
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, SMD.NCombined > nFC);

end
