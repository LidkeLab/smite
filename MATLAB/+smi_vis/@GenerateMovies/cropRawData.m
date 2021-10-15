function [RawData] = cropRawData(RawData, Params)
%cropRawData crops the input data to the region specified in Params.
% This method takes the input 'RawData' and crops it to the size specified
% by the parameters in 'Params'.
%
% INPUTS:
%   RawData: Raw data to be cropped. (YSizexXSizexNFrames or
%            YSizexXSizex3xNFrames numeric array)
%   Params: Structure of movie parameters (see
%           smi_vis.GenerateMovies.prepDefaults())
%
% OUTPUTS:
%   RawData: Cropped input 'RawData'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Crop the data.
if ~(isempty(Params.XPixels) && isempty(Params.YPixels) ...
        && isempty(Params.ZFrames))
    DataSize = size(RawData, 1:2);
    Params.XPixels = [max(1, round(Params.XPixels(1))), ...
        min(DataSize(2), round(Params.XPixels(2)))];
    Params.YPixels = [max(1, round(Params.YPixels(1))), ...
        min(DataSize(1), round(Params.YPixels(2)))];
    RawData = RawData(Params.YPixels(1):Params.YPixels(2), ...
        Params.XPixels(1):Params.XPixels(2), :, :);
end


end