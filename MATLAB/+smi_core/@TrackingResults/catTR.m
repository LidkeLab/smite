function [TR] = catTR(TR1, TR2, CheckDims)
%catTR concatenates two TR structures.
% This method pads each of the inputs 'TR1' and 'TR2' to ensure they share
% the same fields and then concatenates them into a larger structure 'TR'.
%
% INPUTS:
%   TR1: Tracking Results structure.
%   TR2: Tracking Results structure.
%   CheckDims: Flag to indicate manual checking of dimensions should be
%              done (e.g., try to force inputs to be columns). 
%              (Default = true)
%
% OUTPUTS:
%   TR: Concatenation of 'TR1' and 'TR2'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('CheckDims', 'var') || isempty(CheckDims))
    CheckDims = true;
end

% If one of the arrays is empty, we'll just return the non-empty array.
if isempty(TR1)
    TR = TR2;
    return
elseif isempty(TR2)
    TR = TR1;
    return
end

% Pad 'TR1' and 'TR2' to ensure they share the same fields.
TR1 = smi_core.TrackingResults.padTR(TR1, TR2);
TR2 = smi_core.TrackingResults.padTR(TR2, TR1);

% Ensure both 'TR1' and 'TR2' are column structs.
if CheckDims
    if isrow(TR1)
        TR1 = TR1.';
    end
    if isrow(TR2)
        TR2 = TR2.';
    end
end

% Concatenate the two TR structures.
TR = [TR1; TR2];


end