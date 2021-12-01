function [Shift, IntShift, CorrData, Params] = ...
    findStackOffset(Stack1, Stack2, Params)
%findStackOffset estimates a sub-pixel offset between two stacks of images.
% findStackoffset() will estimate the offset between two 3D stacks of
% images.  This method computes an integer pixel offset between the two
% stacks via a cross-correlation and then fits 2nd order polynomials to the 
% resulting cross-correlation.  An estimate of a sub-pixel offset is then 
% produced by determining the location of the peaks of the three (x, y, z)
% 2nd order polynomial fits.
%
% NOTE: The convention used here for the offset is based on indices as
%       follows: If Stack is a 3D stack of images, and
%       Stack1 = Stack(m:n, m:n, m:n) 
%       Stack2 = Stack((m:n)+x, (m:n)+y, (m:n)+z)
%       then PixelOffset = findStackOffset(Stack1, Stack2) == [x; y; z]
% NOTE: All inputs besides Stack1 and Stack2 are optional and can be
%       replaced by [] (an empty array).
% NOTE: Stack1 and Stack2 must be the same size in all 3 dimensions (x, y,
%       and z)
%
% INPUTS:
%   Stack1: The stack to which Stack2 is compared to, i.e. Stack1 is the
%           reference stack. (MxNxO) 
%   Stack2: The stack for which the offset relative to Stack1 is to be 
%           determined. (MxNxO)
%   Params: Structure of parameters.
%           MaxOffset: Maximum offset between Stack1 and Stack2 to be 
%                      considered in the calculation of 'Shift' and
%                      'IntShift'.
%                      (3x1)(Default = ceil(size(Stack1)/2)) 
%           FitOffset: Maximum offset from the peak of the 
%                      cross-correlation curve for which data will be fit
%                      to determine 'Shift'. 
%                      (3x1)(Default = [2; 2; 2]) 
%           BinaryMask: Mask to multiply the stacks with before computing
%                       to cross-correlation. 
%                       (MxNxO)(Default = ones(M, N, O)) 
%           PlotFlag: Specifies whether or not plot(s) will be generated.
%           UseGPU: Flag indicating GPU should be used.
%                   (boolean)(Default = logical(gpuDeviceCount()))
%           SuppressWarnings: Flag indicating we should suppress all
%                             warnings. (Default = false)
%
% OUTPUTS:
%   Shift: The sub-pixel offset of Stack2 relative to Stack1, approximated 
%          based on a 2nd order polynomial fit(s) to a scaled 
%          cross-correlation. (3x1)
%   IntShift: The integer pixel offset of Stack2 relative to Stack1,
%             determined based on the location of the peak of the xcorr 
%             coefficient field between the two stacks. (3x1)
%   CorrData: Structure array containing the scaled cross-correlation 
%             corresponding to MaxOffset as well as the fitting results 
%             which were used to determine 'Shift'.
%   Params: Input 'Params' padded with defaults/with values modified based
%           on the data.
% 
% REQUIRES:
%   MATLAB Parallel Computing Toolbox (if setting UseGPU = 1)
%   Supported NVIDIA CUDA GPU (if setting UseGPU = 1)
%
% CITATION:

% Created by:
%   David J. Schodt (LidkeLab 2018)


% Set default parameter values if needed.
DefaultParams.MaxOffset = ceil(size(Stack1) / 2);
DefaultParams.FitOffset = [2; 2; 2];
DefaultParams.BinaryMask = ones(size(Stack1));
DefaultParams.PlotFlag = false;
DefaultParams.UseGPU = logical(gpuDeviceCount());
DefaultParams.SuppressWarnings = false;
if (~exist('Params', 'var') || isempty(Params))
    Params = DefaultParams;
else
    Params = smi_helpers.padStruct(Params, DefaultParams);
end

% If requested, turn of warnings.
if Params.SuppressWarnings
    WarningState = warning('off');
end

% Ensure MaxOffset and FitOffset are column vectors.
if isrow(Params.MaxOffset)
    Params.MaxOffset = Params.MaxOffset.';
end
if isrow(Params.FitOffset)
    Params.FitOffset = Params.FitOffset.';
end

% Check if the stacks are actually stacks (i.e. multiple images in each
% stack).  If they are not, copy the given image to match the stack size of
% the other stack.
NImagesStack1 = size(Stack1, 3);
NImagesStack2 = size(Stack2, 3);
if ((NImagesStack1==1) || (NImagesStack2==1))
    % One of the two stacks is just a single image, copy that image to
    % match the size of the other stack.
    if (NImagesStack1 < NImagesStack2)
        % Stack1 is a single image: create copies of this image to ensure
        % Stack1 and Stack2 are the same size.
        Stack1 = repmat(Stack1, [1, 1, NImagesStack2]);
    elseif (NImagesStack1 > NImagesStack2)
        % Stack2 is a single image: create copies of this image to ensure
        % Stack1 and Stack2 are the same size.
        Stack2 = repmat(Stack2, [1, 1, NImagesStack1]);
    else
        % Both stacks contain only a single image.  Make these stacks of
        % two images so that we can still proceed with the shift finding in
        % the x and y dimensions.
        Stack1 = repmat(Stack1, [1, 1, 2]);
        Stack2 = repmat(Stack2, [1, 1, 2]);
    end
    
    % We no longer care about the z shift, so change the third element of
    % MaxOffset to 0 to reduce computation time.
    Params.MaxOffset(3) = 0;
end

% Convert the stacks to gpuArrays if needed.
if Params.UseGPU
    Stack1 = gpuArray(Stack1);
    Stack2 = gpuArray(Stack2);
    Params.BinaryMask = gpuArray(Params.BinaryMask);
end

% Ensure the stacks are floating point arrays.
% NOTE: If using the GPU, we should convert to single after sending the
%       stacks to the GPU with gpuArray() (it's faster this way).
Stack1 = single(Stack1);
Stack2 = single(Stack2);

% Determine dimensions relevant to the problem to improve code readability.
Stack1Size = size(Stack1).';
Stack2Size = size(Stack2).';
SizeOfFullXCorr = Stack1Size + Stack2Size - 1;

% Ensure that MaxOffset and FitOffset are valid, modifying their values if
% needed.
% NOTE: This is just ensuring that the MaxOffset corresponds to shifts
%       between the two stacks that still maintain some overlap.
IndicesToModify = find(Params.MaxOffset > floor(SizeOfFullXCorr/2)).';
for ii = IndicesToModify
    warning('MaxOffset(%i) = %g is too big and was reset to %i', ...
        ii, Params.MaxOffset(ii), floor(SizeOfFullXCorr(ii) / 2))
    Params.MaxOffset(ii) = floor(SizeOfFullXCorr(ii) / 2);
end
BadFitOffset = (Params.FitOffset > Params.MaxOffset);
Params.FitOffset(BadFitOffset) = Params.MaxOffset(BadFitOffset);

% Define the indices within a full cross-correlation (size SizeOfFullXCorr)
% that we wish to inspect.
CorrOffsetIndicesY = ...
    max(ceil(SizeOfFullXCorr(1)/2) - Params.MaxOffset(1), 1) ...
    : ceil(SizeOfFullXCorr(1)/2) + Params.MaxOffset(1);
CorrOffsetIndicesX = ...
    max(ceil(SizeOfFullXCorr(2)/2) - Params.MaxOffset(2), 1) ...
    : ceil(SizeOfFullXCorr(2)/2) + Params.MaxOffset(2);
CorrOffsetIndicesZ = ...
    max(ceil(SizeOfFullXCorr(3)/2) - Params.MaxOffset(3), 1) ...
    : ceil(SizeOfFullXCorr(3)/2) + Params.MaxOffset(3);

% Scale each image in each stack by intensity to reduce linear trends in 
% the cross-correlation.
for ii = 1:Stack1Size(3)
    Stack1(:, :, ii) = Stack1(:, :, ii) / sum(sum(Stack1(:, :, ii)));
end
for ii = 1:Stack2Size(3)
    Stack2(:, :, ii) = Stack2(:, :, ii) / sum(sum(Stack2(:, :, ii)));
end

% Whiten each image in the stack with respect to the entire stack, ignoring
% the parts which are covered by the BinaryMask when computing mean, std., 
% etc.
Stack1Masked = Stack1(logical(Params.BinaryMask));
Stack2Masked = Stack2(logical(Params.BinaryMask));
Stack1Whitened = (Stack1-mean(Stack1Masked)) ...
    / (std(Stack1Masked) * sqrt(numel(Stack1Masked)-1));
Stack2Whitened = (Stack2-mean(Stack2Masked)) ...
    / (std(Stack2Masked) * sqrt(numel(Stack2Masked)-1));

% Re-apply the binary mask to ensure the masked points cannot contribute to
% the cross-correlation.
Stack1Whitened = Params.BinaryMask .* Stack1Whitened;
Stack2Whitened = Params.BinaryMask .* Stack2Whitened;

% Compute the 3D FFT's of each stack, padding with zeros before computing.
% The padding size selected such that the result is approximately 
% equivalent to the brute forced cross-correlation.
% NOTE: Typically, we would pad to 2*size(Stack)-1, however using
%       2*size(Stack) will improve the performance of the FFT when the 
%       dimensions of Stack are powers of 2.
Stack1PaddedFFT = fftn(Stack1Whitened, 2 * size(Stack1Whitened));
Stack2PaddedFFT = fftn(Stack2Whitened, 2 * size(Stack2Whitened));

% Compute the 3D cross-correlation in the Fourier domain.
XCorr3D = ifftn(conj(Stack1PaddedFFT) .* Stack2PaddedFFT);

% Compute the binary cross-correlation for later use in scaling.
Stack1Binary = Params.BinaryMask .* ones(size(Stack1Whitened));
Stack2Binary = Params.BinaryMask .* ones(size(Stack2Whitened));
Stack1BinaryFFT = fftn(Stack1Binary, 2 * size(Stack1Whitened));
Stack2BinaryFFT = fftn(Stack2Binary, 2 * size(Stack2Whitened));
XCorr3DBinary = ifftn(conj(Stack1BinaryFFT) .* Stack2BinaryFFT);

% Scale the 3D cross-correlation by the cross-correlation of the
% zero-padded binary images (an attempt to reduce the bias to a [0, 0, 0]
% offset introduced by the zero-padded edge effects), scaling by
% max(XCorr3DBinary(:)) to re-convert to a correlation coefficient.
XCorr3D = (XCorr3D./XCorr3DBinary) * max(XCorr3DBinary(:));

% Shift the cross-correlation image such that an auto-correlation image 
% will have it's energy peak at the center of the 3D image.
XCorr3D = circshift(XCorr3D, size(Stack1Whitened) - 1);

% Isolate the central chunk of the cross-correlation.
XCorr3D = real(XCorr3D(CorrOffsetIndicesY, ...
    CorrOffsetIndicesX, ...
    CorrOffsetIndicesZ));

% Fetch the cross-correlation result from the GPU (if needed).
if Params.UseGPU
    XCorr3D = gather(XCorr3D);
end

% Stack the 3D xcorr cube into a 1D array.  
% NOTE: MATLAB stacks columns in each 2D array (dimensions 1 and 2) then 
%       stacks the resulting columns along the third dimension, e.g. when
%       Array(:, :, 1) = [1, 3; 2, 4] and Array(:, :, 2) = [5, 7; 6, 8], 
%       Array(:) = [1, 2, 3, 4, 5, 6, 7, 8].' 
StackedCorrCube = XCorr3D(:);

% Determine the integer offset between the two stacks.
[~, IndexOfMax] = max(StackedCorrCube);
[PeakRow, PeakColumn, PeakHeight] = ind2sub(size(XCorr3D), IndexOfMax);
RawOffsetIndices = [PeakRow; PeakColumn; PeakHeight];

% Compute the integer offset between the two stacks.  The additional minus
% sign was just a chosen convention for the output of this code (see note 
% at top of this code).
% NOTE: We subtract MaxOffset+1 because that is the location of the 
%       [0, 0, 0] offset (the center of the cross-correlation).
IntShift = -(RawOffsetIndices - Params.MaxOffset - 1);

% Fit a second order polynomial through a line varying with y at the peak
% of the cross-correlation in x, z, and use that polynomial to predict an
% offset.  If possible, center the fit around the integer peak of the
% cross-correlation.
YArray = (max(1, RawOffsetIndices(1)-Params.FitOffset(1)) ...
    : min(2*Params.MaxOffset(1)+1, RawOffsetIndices(1)+Params.FitOffset(1))).';
YData = XCorr3D(YArray, RawOffsetIndices(2), RawOffsetIndices(3));
X = [ones(numel(YArray), 1), YArray, YArray.^2];
Beta = ((X.'*X) \ X.') * YData;
RawOffsetFitY = -Beta(2) / (2 * Beta(3));
PolyFitFunctionY = @(R) Beta(1) + Beta(2)*R + Beta(3)*R.^2;

% Fit a second order polynomial through a line varying with x at the peak 
% of the cross-correlation in y, z.
XArray = (max(1, RawOffsetIndices(2)-Params.FitOffset(2)) ...
    : min(2*Params.MaxOffset(2)+1, RawOffsetIndices(2)+Params.FitOffset(2))).';
XData = ...
    XCorr3D(RawOffsetIndices(1), XArray, RawOffsetIndices(3)).';
X = [ones(numel(XArray), 1), XArray, XArray.^2];
Beta = ((X.'*X) \ X.') * XData;
RawOffsetFitX = -Beta(2) / (2 * Beta(3));
PolyFitFunctionX = @(R) Beta(1) + Beta(2)*R + Beta(3)*R.^2;

% Fit a second order polynomial through a line varying with z
% at the peak of the cross-correlation in x, y.
ZArray = (max(1, RawOffsetIndices(3)-Params.FitOffset(3)) ...
    : min(2*Params.MaxOffset(3)+1, RawOffsetIndices(3)+Params.FitOffset(3))).';
ZData = squeeze(...
    XCorr3D(RawOffsetIndices(1), RawOffsetIndices(2), ZArray));
X = [ones(numel(ZArray), 1), ZArray, ZArray.^2];
Beta = ((X.'*X) \ X.') * ZData;
RawOffsetFitZ = -Beta(2) / (2 * Beta(3));
PolyFitFunctionZ = @(R) Beta(1) + Beta(2)*R + Beta(3)*R.^2;

% Create arrays of the polynomial fits to use for visualization later on.
XArrayDense = linspace(XArray(1), XArray(end), size(Stack1, 1));
YArrayDense = linspace(YArray(1), YArray(end), size(Stack1, 1));
ZArrayDense = linspace(ZArray(1), ZArray(end), size(Stack1, 1));
XFitAtPeak = PolyFitFunctionX(XArrayDense);
YFitAtPeak = PolyFitFunctionY(YArrayDense);
ZFitAtPeak = PolyFitFunctionZ(ZArrayDense);

% Compute the predicted offset based on the polynomial fits.
RawOffsetFit = [RawOffsetFitY; RawOffsetFitX; RawOffsetFitZ];
        
% Determine the predicted offset between the stack.  The additional minus
% sign was just a chosen convention for the output of this code (see note 
% at top of this code).
% NOTE: We subtract MaxOffset+1 because that is the location of the 
%       [0, 0, 0] offset (the center of the cross-correlation).
Shift = -(RawOffsetFit - Params.MaxOffset - 1);

% Populate the CorrData struct with information that we might wish to use
% later.
CorrData.XCorr3D = XCorr3D;
CorrData.XFitAtPeak = XFitAtPeak;
CorrData.YFitAtPeak = YFitAtPeak;
CorrData.ZFitAtPeak = ZFitAtPeak;

% Display line sections through the integer location of the
% cross-correlation, overlain on the fit along those lines.
if Params.PlotFlag
    PlotFigure = findobj('Tag', 'CorrWindow');
    if isempty(PlotFigure)
        PlotFigure = figure('Tag', 'CorrWindow');
    end
    clf(PlotFigure);
    PlotAxes = subplot(3, 1, 1, 'Parent', PlotFigure);
    plot(PlotAxes, -Params.MaxOffset(1):Params.MaxOffset(1), ...
        XCorr3D(:, RawOffsetIndices(2), RawOffsetIndices(3)), 'x')
    hold(PlotAxes, 'on')
    plot(PlotAxes, YArrayDense-Params.MaxOffset(1)-1, YFitAtPeak)
    title(PlotAxes, 'Y Correlation')
    PlotAxes = subplot(3, 1, 2, 'Parent', PlotFigure);
    plot(PlotAxes, -Params.MaxOffset(2):Params.MaxOffset(2), ...
        XCorr3D(RawOffsetIndices(1), :, RawOffsetIndices(3)), 'x')
    hold(PlotAxes, 'on')
    plot(PlotAxes, XArrayDense-Params.MaxOffset(2)-1, XFitAtPeak)
    title(PlotAxes, 'X Correlation')
    ylabel(PlotAxes, 'Correlation Coefficient')
    PlotAxes = subplot(3, 1, 3, 'Parent', PlotFigure);
    plot(PlotAxes, -Params.MaxOffset(3):Params.MaxOffset(3), ...
        squeeze(XCorr3D(RawOffsetIndices(1), RawOffsetIndices(2), :)), 'x')
    hold(PlotAxes, 'on')
    plot(PlotAxes, ZArrayDense-Params.MaxOffset(3)-1, ZFitAtPeak)
    title(PlotAxes, 'Z Correlation')
    xlabel(PlotAxes, 'Pixel Offset')
end

% Restore the warning state.
if Params.SuppressWarnings
    warning(WarningState)
end


end