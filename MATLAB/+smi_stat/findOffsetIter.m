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
%   NIterMax: Maximum number of iterations. (Default = 5)
%   Tolerance: Tolerance of the shifts allowing early stopping before
%              NIterMax. That is, we stop before NIterMax when the newest
%              estimated shift is less than 'Tolerance'.
%              (3x1 float)(Default = [0; 0; 0])
%   CorrParams: Structure of parameters passed to smi_stat.findOffset().
%   ShiftParams: Structure of parameters passed to smi_stat.shiftImage()
%
% OUTPUTS:
%   Shift: Shift between the image stacks. ([Y; X; Z])
%   IntShift: Integer shift between the image stacks. ([Y; X; Z])
%
% CITATION:
%   Cross-correlation shift finding method:
%       Wester, M.J., Schodt, D.J., Mazloom-Farsibaf, H. et al. Robust,
%       fiducial-free drift correction for super-resolution imaging.
%       Sci Rep 11, 23672 (2021).
%       https://doi.org/10.1038/s41598-021-02850-7
%   Iterative idea motivated by https://doi.org/10.1117/12.603304 , noting
%       that the xcorr bias should approach zero as the shift approaches
%       zero.

% Created by:
%   David J. Schodt (Lidke Lab 2021)


% Set defaults/validate inputs.
if (~exist('NIterMax', 'var') || isempty(NIterMax))
    NIterMax = 5;
end
if (~exist('Tolerance', 'var') || isempty(Tolerance))
    Tolerance = [0; 0; 0];
end
StackSize = size(RefStack, 1:3);
if (~exist('CorrParams', 'var') || isempty(CorrParams) ...
        || ~isfield(CorrParams, 'FTSize') || isempty(CorrParams.FTSize))
    CorrParams.FTSize = 2 .^ nextpow2(StackSize);
end
if (~exist('ShiftParams', 'var') || isempty(ShiftParams))
    ShiftParams = struct([]);
end
Tolerance = padarray(Tolerance, max(0, sum(StackSize>1)-numel(Tolerance)), ...
    'post');

% Prepare a Nyquist mask for findStackOffset().
% NOTE: The extra scaling of FNyquist (which is 0.5) accounts for the size
%       difference between the Fourier transform and the input images.
if (~isfield(CorrParams, 'FTMask') || isempty(CorrParams.FTMask))
    FNyquist = 0.5 * (StackSize/CorrParams.FTSize);
    CorrParams.FTMask = smi_stat.frequencyMask(CorrParams.FTSize, FNyquist);
end

% Iteratively estimate the shift.
[Shift, IntShift, CorrData, CorrParams] = ...
    smi_stat.findOffset(RefStack, MovingStack, CorrParams);
NewShift = Shift;
ii = 1;
while (any(abs(NewShift)>Tolerance) && (ii<NIterMax))
    % Shift the image stack.
    ii = ii + 1;
    MovingStack =smi_stat.shiftImage(MovingStack, NewShift, ShiftParams);

    % Remove the border behind the shift direction (this border now
    % contains junk data that we don't want influencing the
    % cross-correlation).
    StackHalfWidth = floor(size(MovingStack, 1:3) / 2);
    for nn = 1:3
        BorderDirection = ...
            smi_helpers.arrayMUX({'pre', 'post'}, NewShift(nn) < 0);
        Border = [0, 0, 0];
        Border(nn) = ceil(abs(NewShift(nn)));
        if ((StackHalfWidth(nn)-Border(nn)) >= CorrParams.FitOffset(nn))
            % We should only remove the border if there is enough data to
            % still allow for a fit.
            MovingStack = smi_helpers.removeBorder(MovingStack, Border, ...
                BorderDirection);
            RefStack = smi_helpers.removeBorder(RefStack, Border, ...
                BorderDirection);
            CorrParams.BinaryMask = ...
                smi_helpers.removeBorder(CorrParams.BinaryMask, Border, ...
                BorderDirection);
        end
    end

    % Compute the shift.
    [NewShift, NewIntShift, CorrData, CorrParams] = ...
        smi_stat.findOffset(RefStack, MovingStack, CorrParams);
    Shift = Shift + NewShift;
    IntShift = IntShift + NewIntShift;
end


end