function [Shift, IntShift, CorrData, Params] = ...
    findOffset(Stack1, Stack2, Params)
%findOffset estimates a sub-pixel offset between two stacks of images.
% findOffset() will estimate the offset between two 3D stacks of
% images.  This method computes an integer pixel offset between the two
% stacks via a cross-correlation and then fits 2nd order polynomials to the
% resulting cross-correlation.  An estimate of a sub-pixel offset is then
% produced by determining the location of the peaks of the three (y, x, z)
% 2nd order polynomial fits.
%
% NOTE: The convention used here for the offset is based on indices as
%       follows: If Stack is a 3D stack of images, and
%       Stack1 = Stack(m:n, m:n, m:n)
%       Stack2 = Stack((m:n)+y, (m:n)+x, (m:n)+z)
%       then PixelOffset = findStackOffset(Stack1, Stack2) == [y; x; z]
% NOTE: All inputs besides Stack1 and Stack2 are optional and can be
%       replaced by [] (an empty array).
% NOTE: Stack1 and Stack2 must be the same size in all 3 dimensions
%       (y, x, and z)
%
% INPUTS:
%   Stack1: The stack to which Stack2 is compared to, i.e. Stack1 is the
%           reference stack. (MxNxO)
%   Stack2: The stack for which the offset relative to Stack1 is to be
%           determined. (MxNxO)
%   Params: Structure of parameters.
%           MaxOffset: Maximum offset between Stack1 and Stack2 to be
%                      considered in the calculation of 'Shift' and
%                      'IntShift'.  This also gets applied as a maximum
%                      shift for the output 'Shift'.
%                      (1x3)(Default = ceil(size(Stack1)/2))
%           FitOffset: Maximum offset from the peak of the
%                      cross-correlation curve for which data will be fit
%                      to determine 'Shift'.
%                      (1x3)(Default = [2; 2; 2])
%           SymmetrizeFit: Flag indicating we should attempt to symmetrize
%                          the fit offset points (sometimes the peak is off
%                          by +-1 due to noise, and this will attempt to
%                          shift the fit points to reflect that).
%                          (Default = true)
%           BinaryMask: Mask to multiply the stacks with before computing
%                       to cross-correlation.
%                       (MxNxO)(Default = ones(M, N, O))
%           FTSize: Size of the Fourier transform used internally.  It's
%                   best to set this to be a power of 2 along each
%                   dimension as doing so improves FFT speed.
%                   (Default = 2 ^ nextpow2(Stack1Size))
%           FTMask: Mask to apply to the Fourier domain cross-correlation
%                   before inverting. (FTSize array)
%           PlotFlag: Specifies whether or not plot(s) will be generated.
%           UseGPU: Flag indicating GPU should be used.
%                   (boolean)(Default = false)
%           SuppressWarnings: Flag indicating we should suppress all
%                             warnings. (Default = false)
%
% OUTPUTS:
%   Shift: The sub-pixel offset of Stack2 relative to Stack1, approximated
%          based on a 2nd order polynomial fit(s) to a scaled
%          cross-correlation. (3x1)([y; x; z])
%   IntShift: The integer pixel offset of Stack2 relative to Stack1,
%             determined based on the location of the peak of the xcorr
%             coefficient field between the two stacks. (3x1)([y; x; z])
%   CorrData: Structure array containing the scaled cross-correlation
%             corresponding to MaxOffset as well as the fitting results
%             which were used to determine 'Shift'.
%   Params: Input 'Params' padded with defaults/with values modified based
%           on the data.
%
% REQUIRES:
%   MATLAB Statistics and Machine Learning Toolbox (if setting
%       SymmetrizeFit = true)
%   MATLAB Parallel Computing Toolbox (if setting UseGPU = true)
%   Supported NVIDIA CUDA GPU (if setting UseGPU = true)
%
% CITATION:
%   Wester, M.J., Schodt, D.J., Mazloom-Farsibaf, H. et al. Robust,
%   fiducial-free drift correction for super-resolution imaging.
%   Sci Rep 11, 23672 (2021). https://doi.org/10.1038/s41598-021-02850-7

% Created by:
%   David J. Schodt (LidkeLab 2018)


% Set default parameter values if needed.
Stack1Size = size(Stack1, 1:3);
Stack2Size = size(Stack2, 1:3);
DefaultParams.MaxOffset = ceil(Stack1Size / 2);
DefaultParams.FitOffset = [2, 2, 2];
DefaultParams.SymmetrizeFit = true;
DefaultParams.BinaryMask = ones(Stack1Size);
DefaultParams.FTSize = 2 .^ nextpow2(Stack1Size);
DefaultParams.FTMask = [];
DefaultParams.PlotFlag = false;
DefaultParams.UseGPU = logical(gpuDeviceCount());
DefaultParams.SuppressWarnings = false;
if (~exist('Params', 'var') || isempty(Params))
    Params = DefaultParams;
else
    Params = smi_helpers.padStruct(Params, DefaultParams);
end
if isempty(Params.FTMask)
    Params.FTMask = ones(Params.FTSize);
elseif ~all(size(Params.FTMask, 1:3) == Params.FTSize)
    if ~Params.SuppressWarnings
        warning(['Params.FTMask reset to ones: size not ', ...
            'compatible with Params.FTSize.'])
    end
    Params.FTMask = ones(Params.FTSize);
end

% If requested, turn of warnings.
if Params.SuppressWarnings
    WarningState = warning('off');
end

% Ensure MaxOffset and FitOffset are row vectors.
if iscolumn(Params.MaxOffset)
    Params.MaxOffset = Params.MaxOffset.';
end
if iscolumn(Params.FitOffset)
    Params.FitOffset = Params.FitOffset.';
end

% Check if the stacks are actually stacks (i.e. multiple images in each
% stack).  If they are not, copy the given image to match the stack size of
% the other stack.
if ((Stack1Size(3)==1) || (Stack2Size(3)==1))
    % One of the two stacks is just a single image, copy that image to
    % match the size of the other stack.
    if (Stack1Size(3) < Stack2Size(3))
        % Stack1 is a single image: create copies of this image to ensure
        % Stack1 and Stack2 are the same size.
        Stack1 = repmat(Stack1, [1, 1, Stack2Size(3)]);
    elseif (Stack1Size(3) > Stack2Size(3))
        % Stack2 is a single image: create copies of this image to ensure
        % Stack1 and Stack2 are the same size.
        Stack2 = repmat(Stack2, [1, 1, Stack1Size(3)]);
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
    Params.FitOffset(3) = 0;
end

% Convert the stacks to gpuArrays if needed.
if Params.UseGPU
    Stack1 = gpuArray(Stack1);
    Stack2 = gpuArray(Stack2);
    Params.BinaryMask = gpuArray(Params.BinaryMask);
    Params.FTMask = gpuArray(Params.FTMask);
end

% Ensure the stacks are floating point arrays.
% NOTE: If using the GPU, we should convert to single after sending the
%       stacks to the GPU with gpuArray() (it's faster this way).
Stack1 = single(Stack1);
Stack2 = single(Stack2);
Params.BinaryMask = single(Params.BinaryMask);
Params.FTMask = single(Params.FTMask);

% Ensure that MaxOffset is valid.
% NOTE: This is just ensuring that the MaxOffset corresponds to shifts
%       between the two stacks that still maintain some overlap.
MaxAllowedOffset = max(0, floor(Params.FTSize / 2) - 1);
IndicesToModify = find(Params.MaxOffset > MaxAllowedOffset);
for ii = IndicesToModify
    warning('MaxOffset(%i) = %g is too big and was reset to %i', ...
        ii, Params.MaxOffset(ii), MaxAllowedOffset(ii))
    Params.MaxOffset(ii) = MaxAllowedOffset(ii);
end

% Scale each image in each stack by intensity to reduce linear trends in
% the cross-correlation.
Stack1 = Stack1 ./ sum(Stack1, [1, 2]);
Stack2 = Stack2 ./ sum(Stack2, [1, 2]);

% Whiten each image in the stack with respect to the entire stack, ignoring
% the parts which are covered by the BinaryMask when computing mean, std.,
% etc.
Stack1Masked = Stack1 .* Params.BinaryMask;
Stack2Masked = Stack2 .* Params.BinaryMask;
Stack1Whitened = (Stack1-mean(Stack1Masked(:))) ...
    / (std(Stack1Masked(:)) * sqrt(numel(Stack1Masked)-1));
Stack2Whitened = (Stack2-mean(Stack2Masked(:))) ...
    / (std(Stack2Masked(:)) * sqrt(numel(Stack2Masked)-1));

% Compute the 3D FFT's of each stack, padding with zeros before computing.
% Also, ensure that the Params.BinaryMask is reapplied.
Stack1PaddedFFT = fftn(Params.BinaryMask .* Stack1Whitened, Params.FTSize);
Stack2PaddedFFT = fftn(Params.BinaryMask .* Stack2Whitened, Params.FTSize);

% Compute the 3D cross-correlation in the Fourier domain.
XCorr3D = ...
    real(ifftn(conj(Stack1PaddedFFT) .* Stack2PaddedFFT .* Params.FTMask));

% Compute the binary cross-correlation for later use in scaling.
BinaryStackFFT = fftn(Params.BinaryMask, Params.FTSize);
XCorr3DBinary = ...
    real(ifftn(conj(BinaryStackFFT) .* BinaryStackFFT .* Params.FTMask));

% Scale the 3D cross-correlation by the cross-correlation of the
% zero-padded binary images (an attempt to reduce the bias to a [0, 0, 0]
% offset introduced by the zero-padded edge effects), scaling by
% max(XCorr3DBinary(:)) to re-convert to a correlation coefficient.
XCorr3D = (XCorr3D./XCorr3DBinary) * max(XCorr3DBinary(:));

% Shift the cross-correlation image such that an auto-correlation image
% will have it's energy peak at the center of the 3D image.
XCorr3D = circshift(XCorr3D, ceil(Params.FTSize/2) - 1);

% Define the indices within the cross-correlation that we wish to inspect.
CorrCenter = floor(Params.FTSize / 2) + rem(Params.FTSize, 2);
CorrOffsetIndicesY = (-Params.MaxOffset(2):Params.MaxOffset(2)) ...
    + CorrCenter(2);
CorrOffsetIndicesX = (-Params.MaxOffset(1):Params.MaxOffset(1)) ...
    + CorrCenter(1);
CorrOffsetIndicesZ = (-Params.MaxOffset(3):Params.MaxOffset(3)) ...
    + CorrCenter(3);

% Isolate the central chunk of the cross-correlation.
XCorr3D = XCorr3D(CorrOffsetIndicesY, ...
    CorrOffsetIndicesX, ...
    CorrOffsetIndicesZ);

% Fetch the cross-correlation result from the GPU (if needed).
if Params.UseGPU
    XCorr3D = gather(XCorr3D);
end

% Determine the integer offset between the two stacks.
[~, IndexOfMax] = max(XCorr3D(:));
[PeakRow, PeakColumn, PeakHeight] = ind2sub(size(XCorr3D), IndexOfMax);
RawOffsetIndices = [PeakRow; PeakColumn; PeakHeight];

% Define some arrays w.r.t. the integer offset.
YArray = (max(1, RawOffsetIndices(1)-Params.FitOffset(1)) ...
    : min(size(XCorr3D, 1), RawOffsetIndices(1)+Params.FitOffset(1))).';
YData = XCorr3D(YArray, RawOffsetIndices(2), RawOffsetIndices(3));
XArray = (max(1, RawOffsetIndices(2)-Params.FitOffset(2)) ...
    : min(size(XCorr3D, 2), RawOffsetIndices(2)+Params.FitOffset(2))).';
XData = XCorr3D(RawOffsetIndices(1), XArray, RawOffsetIndices(3)).';
ZArray = (max(1, RawOffsetIndices(3)-Params.FitOffset(3)) ...
    : min(size(XCorr3D, 3), RawOffsetIndices(3)+Params.FitOffset(3))).';
ZData = squeeze(XCorr3D(RawOffsetIndices(1), RawOffsetIndices(2), ZArray));

% If requested, attempt to "symmetrize" our fit points (sometimes the peak
% xcorr is +-1 off from the visual symmetry center due to noise, which
% might negatively affect fit results).
if Params.SymmetrizeFit
    % Symmetrize Z arrays.
    % NOTE: Experimentally, this asymmetry issue is more prevalent for Z,
    %       so it's probably best to resymmetrize Z first.
    NeighborInds = min(numel(ZData), ...
        max(1, Params.FitOffset(3) + [0, 2]));
    [~, MaxNeighborInd] = max(ZData(NeighborInds));
    ProposedInd = RawOffsetIndices(3) ...
        + (MaxNeighborInd==2) - (MaxNeighborInd==1);
    ProposedInd = min(size(XCorr3D, 3), max(1, ProposedInd));
    ZArrayProposed = (max(1, ProposedInd-Params.FitOffset(3)) ...
        : min(size(XCorr3D, 3), ProposedInd+Params.FitOffset(3))).';
    ZDataProposed = squeeze(...
        XCorr3D(RawOffsetIndices(1), RawOffsetIndices(2), ZArrayProposed));
    if (abs(skewness(ZDataProposed)) < abs(skewness(ZData)))
        RawOffsetIndices(3) = ProposedInd;
        ZArray = ZArrayProposed;
        ZData = ZDataProposed;
        YArray = (max(1, RawOffsetIndices(1)-Params.FitOffset(1)) ...
            : min(size(XCorr3D, 1), RawOffsetIndices(1)+Params.FitOffset(1))).';
        YData = XCorr3D(YArray, RawOffsetIndices(2), RawOffsetIndices(3));
        XArray = (max(1, RawOffsetIndices(2)-Params.FitOffset(2)) ...
            : min(size(XCorr3D, 2), RawOffsetIndices(2)+Params.FitOffset(2))).';
        XData = XCorr3D(RawOffsetIndices(1), XArray, RawOffsetIndices(3)).';
    end

    % Symmetrize Y arrays.
    NeighborInds = min(numel(YData), ...
        max(1, Params.FitOffset(1) + [0, 2]));
    [~, MaxNeighborInd] = max(YData(NeighborInds));
    ProposedInd = RawOffsetIndices(1) ...
        + (MaxNeighborInd==2) - (MaxNeighborInd==1);
    ProposedInd = min(size(XCorr3D, 1), max(1, ProposedInd));
    YArrayProposed = (max(1, ProposedInd-Params.FitOffset(1)) ...
        : min(size(XCorr3D, 1), ProposedInd+Params.FitOffset(1))).';
    YDataProposed = ...
        XCorr3D(YArrayProposed, RawOffsetIndices(2), RawOffsetIndices(3));
    if (abs(skewness(YDataProposed)) < abs(skewness(YData)))
        RawOffsetIndices(1) = ProposedInd;
        YArray = YArrayProposed;
        YData = YDataProposed;
        XArray = (max(1, RawOffsetIndices(2)-Params.FitOffset(2)) ...
            : min(size(XCorr3D, 2), RawOffsetIndices(2)+Params.FitOffset(2))).';
        XData = XCorr3D(RawOffsetIndices(1), XArray, RawOffsetIndices(3)).';
        ZArray = (max(1, RawOffsetIndices(3)-Params.FitOffset(3)) ...
            : min(size(XCorr3D, 3), RawOffsetIndices(3)+Params.FitOffset(3))).';
        ZData = squeeze(XCorr3D(RawOffsetIndices(1), RawOffsetIndices(2), ...
            ZArray));
    end

    % Symmetrize X arrays.
    NeighborInds = min(numel(XData), ...
        max(1, Params.FitOffset(2) + [0, 2]));
    [~, MaxNeighborInd] = max(XData(NeighborInds));
    ProposedInd = RawOffsetIndices(2) ...
        + (MaxNeighborInd==2) - (MaxNeighborInd==1);
    ProposedInd = min(size(XCorr3D, 2), max(1, ProposedInd));
    XArrayProposed = (max(1, ProposedInd-Params.FitOffset(2)) ...
        : min(size(XCorr3D, 2), ProposedInd+Params.FitOffset(2))).';
    XDataProposed = ...
        XCorr3D(RawOffsetIndices(1), XArrayProposed, RawOffsetIndices(3)).';
    if (abs(skewness(XDataProposed)) < abs(skewness(XData)))
        RawOffsetIndices(2) = ProposedInd;
        XArray = XArrayProposed;
        XData = XDataProposed;
        YArray = (max(1, RawOffsetIndices(1)-Params.FitOffset(1)) ...
            : min(size(XCorr3D, 1), RawOffsetIndices(1)+Params.FitOffset(1))).';
        YData = XCorr3D(YArray, RawOffsetIndices(2), RawOffsetIndices(3));
        ZArray = (max(1, RawOffsetIndices(3)-Params.FitOffset(3)) ...
            : min(size(XCorr3D, 3), RawOffsetIndices(3)+Params.FitOffset(3))).';
        ZData = squeeze(XCorr3D(RawOffsetIndices(1), RawOffsetIndices(2), ...
            ZArray));
    end
end

% Compute the integer offset between the two stacks.
IntShift = Params.MaxOffset.' - RawOffsetIndices + 1;
IntShift(isnan(IntShift)) = 0;

% Fit a second order polynomial through a line varying with y at the peak
% of the cross-correlation in x, z, and use that polynomial to predict an
% offset.  If possible, center the fit around the integer peak of the
% cross-correlation.
X = [ones(numel(YArray), 1), YArray, YArray.^2];
Beta = ((X.'*X) \ X.') * YData;
RawOffsetFitY = -Beta(2) / (2 * Beta(3));
PolyFitFunctionY = @(R) Beta(1) + Beta(2)*R + Beta(3)*R.^2;

% Fit a second order polynomial through a line varying with x at the peak
% of the cross-correlation in y, z.
X = [ones(numel(XArray), 1), XArray, XArray.^2];
Beta = ((X.'*X) \ X.') * XData;
RawOffsetFitX = -Beta(2) / (2 * Beta(3));
PolyFitFunctionX = @(R) Beta(1) + Beta(2)*R + Beta(3)*R.^2;

% Fit a second order polynomial through a line varying with z
% at the peak of the cross-correlation in x, y.
X = [ones(numel(ZArray), 1), ZArray, ZArray.^2];
Beta = ((X.'*X) \ X.') * ZData;
RawOffsetFitZ = -Beta(2) / (2 * Beta(3));
PolyFitFunctionZ = @(R) Beta(1) + Beta(2)*R + Beta(3)*R.^2;

% Compute the predicted offset based on the polynomial fits.
RawOffsetFit = [RawOffsetFitY; RawOffsetFitX; RawOffsetFitZ];

% Determine the predicted offset between the stack.
Shift = Params.MaxOffset.' - RawOffsetFit + 1;
Shift(isnan(Shift)) = 0;
Shift = min(Params.MaxOffset.', max(-Params.MaxOffset.', Shift));

% Create arrays of the polynomial fits to use for visualization.
XArrayDense = linspace(XArray(1), XArray(end), 8*Params.FitOffset(2));
YArrayDense = linspace(YArray(1), YArray(end), 8*Params.FitOffset(1));
ZArrayDense = linspace(ZArray(1), ZArray(end), 8*Params.FitOffset(3));
XFitAtPeak = PolyFitFunctionX(XArrayDense);
YFitAtPeak = PolyFitFunctionY(YArrayDense);
ZFitAtPeak = PolyFitFunctionZ(ZArrayDense);

% Populate the CorrData struct.
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