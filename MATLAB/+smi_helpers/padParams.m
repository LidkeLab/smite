function [Params] = padParams(Params, DefaultParams)
%padParams pads the input 'Params' with defaults in 'DefaultParams'.
% This method merges the two structures 'Params' and 'DefaultParams', with
% values in 'Params' taking precedent unless not available. 
%
% INPUTS:
%   Params: Structure array which is to be padded.
%   DefaultParams: Structure array of "default" parameters that will be
%                  added to 'Params' unless already present in 'Params', in
%                  which case the "default" values are ignored.
%
% OUTPUTS:
%   Params: Input 'Params' padded with values from 'DefaultParams'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Loop through the fields of the default parameters and add them to
% 'Params' as needed.
DefaultParameterNames = fieldnames(DefaultParams);
InputParameterNames = fieldnames(Params);
for pp = 1:numel(DefaultParameterNames)
    if ~any(ismember(DefaultParameterNames{pp}, InputParameterNames))
        % The field DefaultParameterNames{pp} is not present in the
        % Params structure and so the default must be added.
        Params.(DefaultParameterNames{pp}) = ...
            DefaultParams.(DefaultParameterNames{pp});
    end
end


end