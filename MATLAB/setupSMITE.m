function setupSMITE()
%Run this function in startup.m to setup required paths for smite.
% If the smite folder is located in userpath, then use the following: 
% 
% MATLAB 2017a and later:
%   run(fullfile(userpath,'smite','setupSMITE'))
% MATLAB 2016b and ealier:
%   run(fullfile(userpath(1:end-1),'smite','setupSMITE'))

SMITEPath=fileparts(which('setupSMITE'));

addpath(SMITEPath)
addpath(fullfile(SMITEPath,'mex'))
addpath(fullfile(SMITEPath,'ptx'))

ExtSoftPath = fileparts(SMITEPath);
run(fullfile(ExtSoftPath, 'ExternalSoftware', 'setupExternalSoftware'));
fprintf('SMITE version: %s\n', smi_helpers.versionSMITE());

end
