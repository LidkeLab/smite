function [PixelOffset, SubPixelOffset, CorrData, MaxOffset] = ...
    findStackOffset(Stack1, Stack2, MaxOffset, FitOffset, ...
    BinaryMask, PlotFlag, UseGPU)
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
%   Stack1: (mxnxo) The stack to which Stack2 is compared to, i.e. 
%           Stack1 is the reference stack.
%   Stack2: (mxnxo) The stack for which the offset relative to Stack1 
%           is to be determined.
%   MaxOffset: (3x1 or 1x3)(Default = [5; 5; 5]) Maximum offset between 
%           Stack1 and Stack2 to be considered in the calculation of
%           PixelOffset and SubPixelOffset.
%   FitOffset: (3x1 or 1x3)(Default = [2; 2; 2]) Maximum offset from the
%           peak of the cross-correlation curve for which data will be fit
%           to determine SubPixelOffset.
%   BinaryMask: (mxnxo)(Default = ones(m, n, o)) Mask to multiply the
%           stacks with before computing to cross-correlation.
%   PlotFlag: (boolean)(Default = true) Specifies whether or not the 1D 
%           line plots through the peak of the xcorr will be shown.  
%           PlotFlag = true will allow the plots to be shown, 
%           PlotFlag = false will not allow plots to be displayed.
%   UseGPU: (boolean)(Default = logical(gpuDeviceCount()))
%
% OUTPUTS:
%   PixelOffset: (3x1)(integer) The integer pixel offset of Stack2 relative
%           to Stack1, determined based on the location of the peak of the
%           xcorr coefficient field between the two stacks.
%   SubPixelOffset: (3x1)(float) The sub-pixel offset of Stack2 relative to
%           Stack1, approximated based on a 2nd order polynomial fit(s) to 
%           the cross-correlation.
%   CorrData: (struct) Structured array containing the scaled
%             cross-correlation corresponding to MaxOffset as well as the
%             fitting results which were used to determine SubPixelOffset.
%   MaxOffset: (3x1 or 1x3) Maximum offset between Stack1 and Stack2 
%           considered in the calculation of PixelOffset and 
%           SubPixelOffset.  This is returned because the user input value
%           of MaxOffset is truncated if too large. 
% 
% REQUIRES:
%   MATLAB Parallel Computing Toolbox (if setting UseGPU = 1)
%   Supported NVIDIA CUDA GPU (if setting UseGPU = 1)
%
% CITATION:

% Created by:
%   David J. Schodt (LidkeLab 2018)


% Set default parameter values if needed.
if (~exist('MaxOffset', 'var') || isempty(MaxOffset))
    MaxOffset = [5; 5; 5];
end
if (~exist('FitOffset', 'var') || isempty(FitOffset))
    FitOffset = [2; 2; 2];
end
if (~exist('PlotFlag', 'var') || isempty(PlotFlag))
    PlotFlag = true;
end
if (~exist('BinaryMask', 'var') || isempty(BinaryMask))
    BinaryMask = ones(size(Stack1));
    if ((size(Stack1, 3)==1) || (size(Stack2, 3)==1))
        % One or both of the stacks are just a single image, so we need to
        % change the size of the BinaryMask to account for that (single
        % images are converted to a 3D stack later on by copying the image
        % twice in the z dimension).
        BinaryMask = repmat(BinaryMask, [1, 1, 2]);
    end
end
if (~exist('UseGPU', 'var') || isempty(UseGPU))
    UseGPU = logical(gpuDeviceCount());
end

% Ensure MaxOffset is a column vector for consistency.
if isrow(MaxOffset)
    MaxOffset = MaxOffset.';
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
    MaxOffset(3) = 0;
end

% Convert the stacks to gpuArrays if needed.
if UseGPU
    Stack1 = gpuArray(Stack1);
    Stack2 = gpuArray(Stack2);
    BinaryMask = gpuArray(BinaryMask);
end

% Ensure the stacks are floating point arrays.
% NOTE: If using the GPU, we should convert to single after sending the
%       stacks to the GPU with gpuArray() (it's faster this way).
Stack1 = single(Stack1);
Stack2 = single(Stack2);

% Determine dimensions relevant to the problem to improve code readability.
Stack1Size = size(Stack1).';
Stack2Size = size(Stack2).';
SizeOfFullXCorr = Stack1Size + Stack2Size - 1; % size of a full xcorr stack

% Ensure that the MaxOffset input is valid, modifying it's values if
% needed.
% NOTE: This is just ensuring that the MaxOffset corresponds to shifts
%       between the two stacks that still maintain some overlap.
IndicesToModify = find(MaxOffset > floor(SizeOfFullXCorr/2)).';
for ii = IndicesToModify
    warning('MaxOffset(%i) = %g is too big and was reset to %i', ...
        ii, MaxOffset(ii), floor(SizeOfFullXCorr(ii) / 2))
    MaxOffset(ii) = floor(SizeOfFullXCorr(ii) / 2);
end

% Define the indices within a full cross-correlation (size SizeOfFullXCorr)
% that we wish to inspect.
CorrOffsetIndicesX = max(ceil(SizeOfFullXCorr(1)/2) - MaxOffset(1), 1) ...
    : ceil(SizeOfFullXCorr(1)/2) + MaxOffset(1);
CorrOffsetIndicesY = max(ceil(SizeOfFullXCorr(2)/2) - MaxOffset(2), 1) ...
    : ceil(SizeOfFullXCorr(2)/2) + MaxOffset(2);
CorrOffsetIndicesZ = max(ceil(SizeOfFullXCorr(3)/2) - MaxOffset(3), 1) ...
    : ceil(SizeOfFullXCorr(3)/2) + MaxOffset(3);

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
Stack1Masked = Stack1(logical(BinaryMask));
Stack2Masked = Stack2(logical(BinaryMask));
Stack1Whitened = (Stack1-mean(Stack1Masked)) ...
    / (std(Stack1Masked) * sqrt(numel(Stack1Masked)-1));
Stack2Whitened = (Stack2-mean(Stack2Masked)) ...
    / (std(Stack2Masked) * sqrt(numel(Stack2Masked)-1));

% Re-apply the binary mask to ensure the masked points cannot contribute to
% the cross-correlation.
Stack1Whitened = BinaryMask .* Stack1Whitened;
Stack2Whitened = BinaryMask .* Stack2Whitened;

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
Stack1Binary = BinaryMask .* ones(size(Stack1Whitened));
Stack2Binary = BinaryMask .* ones(size(Stack2Whitened));
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
XCorr3D = real(XCorr3D(CorrOffsetIndicesX, ...
    CorrOffsetIndicesY, ...
    CorrOffsetIndicesZ));

% Fetch the cross-correlation result from the GPU (if needed).
if UseGPU
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
PixelOffset = -(RawOffsetIndices - MaxOffset - 1);

% Fit a second order polynomial through a line varying with x at the peak
% of the cross-correlation in y, z, and use that polynomial to predict an
% offset.  If possible, center the fit around the integer peak of the
% cross-correlation.
XArray = (max(1, RawOffsetIndices(1)-FitOffset(1)) ...
    : min(2*MaxOffset(1)+1, RawOffsetIndices(1)+FitOffset(1))).';
XData = XCorr3D(XArray, RawOffsetIndices(2), RawOffsetIndices(3));
X = [ones(numel(XArray), 1), XArray, XArray.^2];
Beta = ((X.'*X) \ X.') * XData; % this is just least-squares fitting
RawOffsetFitX = -Beta(2) / (2 * Beta(3)); % x at peak of the 2D polynomial
PolyFitFunctionX = @(R) Beta(1) + Beta(2)*R + Beta(3)*R.^2;

% Fit a second order polynomial through a line varying with y
% at the peak of the cross-correlation in x, z.
YArray = (max(1, RawOffsetIndices(2)-FitOffset(2)) ...
    : min(2*MaxOffset(2)+1, RawOffsetIndices(2)+FitOffset(2))).';
YData = ...
    XCorr3D(RawOffsetIndices(1), YArray, RawOffsetIndices(3)).';
X = [ones(numel(YArray), 1), YArray, YArray.^2];
Beta = ((X.'*X) \ X.') * YData; % this is just least-squares fitting
RawOffsetFitY = -Beta(2) / (2 * Beta(3)); % y at peak of the 2D polynomial
PolyFitFunctionY = @(R) Beta(1) + Beta(2)*R + Beta(3)*R.^2;

% Fit a second order polynomial through a line varying with z
% at the peak of the cross-correlation in x, y.
ZArray = (max(1, RawOffsetIndices(3)-FitOffset(3)) ...
    : min(2*MaxOffset(3)+1, RawOffsetIndices(3)+FitOffset(3))).';
ZData = squeeze(...
    XCorr3D(RawOffsetIndices(1), RawOffsetIndices(2), ZArray));
X = [ones(numel(ZArray), 1), ZArray, ZArray.^2];
Beta = ((X.'*X) \ X.') * ZData; % this is just least-squares fitting
RawOffsetFitZ = -Beta(2) / (2 * Beta(3)); % z at peak of the 2D polynomial
PolyFitFunctionZ = @(R) Beta(1) + Beta(2)*R + Beta(3)*R.^2;

% Create arrays of the polynomial fits to use for visualization later on.
XArrayDense = linspace(XArray(1), XArray(end), size(Stack1, 1));
YArrayDense = linspace(YArray(1), YArray(end), size(Stack1, 1));
ZArrayDense = linspace(ZArray(1), ZArray(end), size(Stack1, 1));
XFitAtPeak = PolyFitFunctionX(XArrayDense);
YFitAtPeak = PolyFitFunctionY(YArrayDense);
ZFitAtPeak = PolyFitFunctionZ(ZArrayDense);

% Compute the predicted offset based on the polynomial fits.
RawOffsetFit = [RawOffsetFitX; RawOffsetFitY; RawOffsetFitZ];
        
% Determine the predicted offset between the stack.  The additional minus
% sign was just a chosen convention for the output of this code (see note 
% at top of this code).
% NOTE: We subtract MaxOffset+1 because that is the location of the 
%       [0, 0, 0] offset (the center of the cross-correlation).
SubPixelOffset = -(RawOffsetFit - MaxOffset - 1);

% Populate the CorrData struct with information that we might wish to use
% later.
CorrData.XCorr3D = XCorr3D;
CorrData.XFitAtPeak = XFitAtPeak;
CorrData.YFitAtPeak = YFitAtPeak;
CorrData.ZFitAtPeak = ZFitAtPeak;

% Display line sections through the integer location of the
% cross-correlation, overlain on the fit along those lines.
if PlotFlag
    PlotFigure = findobj('Tag', 'CorrWindow');
    if isempty(PlotFigure)
        PlotFigure = figure('Tag', 'CorrWindow');
    end
    clf(PlotFigure); % clear the figure window
    figure(PlotFigure); % ensure we plot into the correct figure
    subplot(3, 1, 1)
    plot(-MaxOffset(1):MaxOffset(1), ...
        XCorr3D(:, RawOffsetIndices(2), RawOffsetIndices(3)), 'x')
    hold('on')
    plot(XArrayDense-MaxOffset(1)-1, XFitAtPeak)
    title('X Correlation')
    subplot(3, 1, 2)
    plot(-MaxOffset(2):MaxOffset(2), ...
        XCorr3D(RawOffsetIndices(1), :, RawOffsetIndices(3)), 'x')
    hold('on')
    plot(YArrayDense-MaxOffset(2)-1, YFitAtPeak)
    title('Y Correlation')
    ylabel('Correlation Coefficient')
    subplot(3, 1, 3)
    plot(-MaxOffset(3):MaxOffset(3), ...
        squeeze(XCorr3D(RawOffsetIndices(1), RawOffsetIndices(2), :)), 'x')
    hold('on')
    plot(ZArrayDense-MaxOffset(3)-1, ZFitAtPeak)
    title('Z Correlation')
    xlabel('Pixel Offset')
end


end