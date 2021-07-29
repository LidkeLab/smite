function [TR] = catTR(TR1, TR2)
%catTR concatenates two TR structures.
% This method pads each of the inputs 'TR1' and 'TR2' to ensure they share
% the same fields and then concatenates them into a larger structure 'TR'.
%
% INPUTS:
%   TR1: Tracking Results structure.
%   TR2: Tracking Results structure.
%
% OUTPUTS:
%   TR: Concatenation of 'TR1' and 'TR2'.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Pad 'TR1' and 'TR2' to ensure they share the same fields.
TR1 = smi_core.TrackingResults.padTR(TR1, TR2);
TR2 = smi_core.TrackingResults.padTR(TR2, TR1);

% Ensure both 'TR1' and 'TR2' are column structs.
if isrow(TR1)
    TR1 = TR1.';
end
if isrow(TR2)
    TR2 = TR2.';
end

% Concatenate the two TR structures.
TR = [TR1; TR2];


end