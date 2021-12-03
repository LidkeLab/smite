function [Shift, IntShift, CorrData, CorrParams, ShiftParams] = ...
    findOffsetIter(RefStack, MovingStack, NIterMax, Tolerance, ...
    CorrParams, ShiftParams)
%findOffsetIter estimates the sub-pixel shift between images.
% This method esimates the shift between 'RefStack' and 'MovingStack' by
% iteratively fitting a polynomial to a scaled cross-correlation between
% the stacks, shifting the stacks via FFT, and then re-computing the shift.
%
% INPUTS:
%   RefStack: Reference stack. (YSizexXSizexNImages)
%   MovingStack: Stack that has moved w.r.t. RefStack (YSizexXSizexNImages)
%   NIterMax: Maximum number of iterations. (Default = 10)
%   Tolerance: Tolerance of the shifts allowing early stopping before
%              NIterMax. That is, we stop before NIterMax when the newest
%              estimated shift is less than 'Tolerance'. 
%              (3x1 float)(Default = [0; 0; 0])
%   CorrParams: Structure of parameters passed to 
%               smi_stat.findStackOffset()
%   ShiftParams: Structure of parameters passed to smi_stat.shiftImage()
%
% OUTPUTS:
%   Shift: Shift between the image stacks. ([Y; X; Z])
%   IntShift: Integer shift between the image stacks. ([Y; X; Z])
%
% CITATION: 
%   Cross-correlation shift finding method: 
%       https://doi.org/10.1101/2021.03.26.437196
%   Iterative idea motivated by https://doi.org/10.1117/12.603304 , noting
%       that the xcorr bias should approach zero as the shift approaches
%       zero.

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Set defaults/validate inputs.
if (~exist('NIterMax', 'var') || isempty(NIterMax))
    NIterMax = 10;
end
if (~exist('Tolerance', 'var') || isempty(Tolerance))
    Tolerance = [0; 0; 0];
end
if (~exist('CorrParams', 'var') || isempty(CorrParams))
    CorrParams = struct([]);
end
if (~exist('ShiftParams', 'var') || isempty(ShiftParams))
    ShiftParams = struct([]);
end
StackSize = size(RefStack);
Tolerance = padarray(Tolerance, max(0, sum(StackSize>1)-numel(Tolerance)), ...
    'post');

% Iteratively estimate the shift.
[Shift, IntShift, CorrData, CorrParams] = ...
    smi_stat.findStackOffset(RefStack, MovingStack, CorrParams);
NewShift = Shift;
ii = 1;
while (all(abs(NewShift)>Tolerance) && (ii<NIterMax))  
    % Shift the image stack.
    ii = ii + 1;
    MovingStack =smi_stat.shiftImage(MovingStack, NewShift, ShiftParams);
    
    % Compute the shift.
    [NewShift, NewIntShift, CorrData, CorrParams] = ...
        smi_stat.findStackOffset(RefStack, MovingStack, CorrParams);
    Shift = Shift + NewShift;
    IntShift = IntShift + NewIntShift;
end


end