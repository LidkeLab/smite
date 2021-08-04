function setupExternalSoftware()
%Run function in startup.m to setup required paths for smite external software.
% If the smite folder is located in userpath, then use the following: 
% 
% MATLAB 2017a and later:
%   run(fullfile(userpath, 'smite', 'ExternalSoftWare', ...
%                'setupExternalSoftware'))
%
% MATLAB 2016b and ealier:
%   run(fullfile(userpath(1:end-1), 'smite', 'ExternalSoftWare', ...
%                'setupExternalSoftware'))

ExternalSoftwarePath = fileparts(which('setupExternalSoftware'));

%addpath(ExternalSoftwarePath)
addpath(fullfile(ExternalSoftwarePath, 'FRCresolution_software', ...
                 'matlabdistribution', 'FRCresolutionfunctions'));
addpath(fullfile(ExternalSoftwarePath, 'PlotSpread'));
addpath(fullfile(ExternalSoftwarePath, 'uipickfiles'));

end
