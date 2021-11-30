function [Shift] = findOffsetIter(RefStack, MovingStack, ...
    NIterMax, Tolerance, UseGPU)
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
%              NIterMax. (3x1 float)(Default = [0; 0; 0])
%   UseGPU: Flag indicating GPU should be used. (Default = false)
%
% OUTPUTS:
%   Shift: Shift between the image stacks. ([Y; X; Z])
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
if (~exist('UseGPU', 'var') || isempty(UseGPU))
    UseGPU = false;
end
StackSize = size(RefStack);
Tolerance = padarray(Tolerance, max(0, sum(StackSize>1)-numel(Tolerance)), ...
    'post');

% Iteratively estimate the shift.
PreviousShift = [inf; inf; inf];
Shift = [0; 0; 0];
ii = 1;
while (all(abs(Shift-PreviousShift)>Tolerance) && (ii<NIterMax))
    % Compute the shift.
    ii = ii + 1;
    [~, NewShift] = smi_stat.findStackOffset(RefStack, MovingStack, ...
        [], [], [], false, UseGPU);
    Shift = Shift + NewShift;
    
    % Shift the image stack.
    MovingStack = smi_stat.shiftImage(MovingStack, NewShift, UseGPU);
end


end