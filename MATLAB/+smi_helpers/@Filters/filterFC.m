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
N = numel(SMD.FrameNum);
%SMD.DatasetNum = ones(size(SMD.FrameNum), 'uint32');
if numel(SMD.PValue) ~= N
   SMD.PValue = zeros(size(SMD.FrameNum));
end
if numel(SMD.ThreshFlag) ~= N
   SMD.ThreshFlag = zeros(size(SMD.FrameNum));
end
if numel(SMD.XBoxCorner) ~= N || numel(SMD.YBoxCorner) ~= N
   SMD.XBoxCorner = zeros(size(SMD.FrameNum));
   SMD.YBoxCorner = zeros(size(SMD.FrameNum));
end
%SMD.Photons = zeros(size(SMD.FrameNum));
%SMD.Bg = zeros(size(SMD.FrameNum));
%SMD.LogLikelihood = zeros(size(SMD.FrameNum));
%SMD.NDims = 2;
%SMD.NDatasets = 10;
%SMD.NFrames = max(SMD.FrameNum);
%SMD.XSize = 2 ^ nextpow2(max([SMD.X; SMD.Y]));
%SMD.YSize = SMD.XSize;

SMF = smi_core.SingleMoleculeFitting();
% Method is 'Hypothesis test' for the BaGoL paper calculations only, otherwise
% use 'LAP-FC' for new results.
%SMF.FrameConnection.Method = 'Hypothesis test';
SMF.FrameConnection.Method = 'LAP-FC';
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
