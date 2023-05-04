function [Success] = unitTest()
%unitTest tests vital functionality of the LocalizeData class.
% This method performs various tests to ensure that vital functionality of
% the smi_core.LocalizeData class is working as intended.
%
% NOTE: Failure of this unit test may in fact be caused by failure of some
%       other class methods OUTSIDE of the LocalizeData class.  A more
%       direct/better test of this class could not be conceptualized at the
%       time of writing.
%
% INPUTS:
%
% OUTPUTS:
%   Success: An array of boolean flags to indicate success of various tests
%            performed.
%            Success(1): genLocalizations - SMD generated successfully.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)

SaveDir = fullfile(tempdir, 'smite', 'unitTest', 'LocalizeData');
if ~isfolder(SaveDir)
   mkdir(fullfile(tempdir, 'smite'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest'));
   mkdir(fullfile(tempdir, 'smite', 'unitTest', 'LocalizeData'));
end

% Seed the random number generator so that we always get the same results.
rng(1234)

% Generate some simulated data in units of photons.
FrameSizeFull = 256; % don't change this! other numbers assume = 256
NFrames = 10;
SMD = smi_core.SingleMoleculeData.createSMD();
SMD.X = repmat(128 + 64*[0; 1; 1; -1; -1], [NFrames, 1]);
SMD.Y = repmat(128 + 64*[0; 1; -1; 1; -1], [NFrames, 1]);
SMD.Photons = 1e3 * ones(5*NFrames, 1);
SMD.PSFSigma = 1.3;
SMD.FrameNum = repelem((1:NFrames).', 5);
SMD.Bg = zeros(5*NFrames, 1);
SMD.NFrames = NFrames;
SMD.YSize = FrameSizeFull;
SMD.XSize = FrameSizeFull;
[~, ScaledData] = smi_sim.GaussBlobs.gaussBlobImage(SMD);

% Generate an SMF structure.
SMF = smi_core.SingleMoleculeFitting;
SMF.BoxFinding.BoxSize = 10;
SMF.Fitting.PSFSigma = SMD.PSFSigma;

% Attempt to generate localizations from the simulated data.
LD = smi_core.LocalizeData(ScaledData, SMF, 3);
[SMDout, SMDPreThresh] = LD.genLocalizations();
saveas(gcf, fullfile(SaveDir, 'LD1.png'));

% Check that SMD and SMDPreThresh make sense and are consistent with each
% other (this isn't meant to be an exact check of individual fields).
NEmitters = 5;
Success(1) = ...
    (all(SMDout.FrameNum==repelem((1:NFrames).', NEmitters)) ...
    && (numel(SMDout.X)==NEmitters*NFrames) ...
    && (numel(SMDout.X)==sum(SMDPreThresh.ThreshFlag==0)));


end
