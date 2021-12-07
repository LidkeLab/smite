function [Shift, IntShift, CorrData, Params] = ...
    findZOffset(Image, Stack, Params)
%findZOffset finds the offset between Image and Stack along Z.
% findZOffset() uses a cross-correlation fit by a polynomial to determine
% the z-shift between Image and Stack, where the central image in 'Stack' 
% defines the 0 coordinate.
%
% INPUTS:
%   Image: The image to be compared to 'Stack'. (MxNx1) 
%   Stack: The stack of images to be compared to Image. (MxNxNImages)
%   Params: Structure of parameters.
%           FitOffset: Maximum offset from the peak of the 
%                      cross-correlation curve for which data will be fit
%                      to determine 'Shift'. (Default = 2) 
%           BinaryMask: Mask to multiply the images before computing the
%                       shift. (Default = zeros(size(Image)))
%           PlotFlag: Specifies whether or not plot(s) will be generated.
%           UseGPU: Flag indicating GPU should be used.
%                   (boolean)(Default = false)
%           SuppressWarnings: Flag indicating we should suppress all
%                             warnings. (Default = false)
%
% OUTPUTS:
%   Shift: The approximate shift along z from 'Image' to the center image
%          in 'Stack'.
%   IntShift: The integer pixel offset of Image relative to Stack.
%   CorrData: Structure array containing the cross-correlation maxima
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
%   David J. Schodt (LidkeLab 2021)


% Set default parameter values if needed.
ImSize = size(Image);
DefaultParams.FitOffset = 2;
DefaultParams.BinaryMask = ones(ImSize);
DefaultParams.PlotFlag = false;
DefaultParams.UseGPU = false;
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

% Convert the stacks to gpuArrays if needed.
if Params.UseGPU
    Image = gpuArray(Image);
    Stack = gpuArray(Stack);
    Params.BinaryMask = gpuArray(Params.BinaryMask);
end

% Ensure the stacks are floating point arrays.
% NOTE: If using the GPU, we should convert to single after sending the
%       stacks to the GPU with gpuArray() (it's faster this way).
Image = single(Image);
Stack = single(Stack);
Params.BinaryMask = single(Params.BinaryMask);

% Ensure that FitOffset is valid, modifying if needed.
NImages = size(Stack, 3);
if (Params.FitOffset > floor(NImages/2))
    warning('FitOffset = %i is too big and was reset to %i', ...
        Params.FitOffset, floor(NImages/2))
    Params.FitOffset(ii) = floor(NImages / 2);
end

% Scale each image in the stack by intensity to reduce linear trends in 
% the cross-correlation.
Stack = Stack ./ sum(Stack, [1, 2]);

% Whiten the image and apply the binary mask.
ImageMasked = Image .* Params.BinaryMask;
StackMasked = Stack .* repmat(Params.BinaryMask, [1, 1, NImages]);
ImageWhitened = (ImageMasked-mean(ImageMasked(:))) ...
    / (std(ImageMasked(:)) * sqrt(numel(ImageMasked)-1));

% Loop through images in the stack and compute the maximum of the
% cross-correlation.
XCorrMax = zeros(NImages, 1);
for nn = 1:NImages
    % Whiten the current stack image.
    CurrentIm = StackMasked(:, :, nn);
    CurrentIm = (CurrentIm-mean(CurrentIm(:))) ...
        / (std(CurrentIm(:)) * sqrt(numel(CurrentIm)-1));

    % Compute the maximum of the cross-correlation.
    XCorrMax(nn) = sum(Params.BinaryMask .* ImageWhitened .* CurrentIm, 1:2);
end

% Compute the integer offset between the two stacks.
[~, IndexOfMax] = max(XCorrMax);
ZeroImageInd = ceil(NImages / 2);
IntShift = ZeroImageInd - IndexOfMax;

% Fit a second order polynomial through the xcorr maxima to determine a
% sub-pixel offset.
XArray = (max(1, IndexOfMax-Params.FitOffset) ...
    : min(NImages, IndexOfMax+Params.FitOffset)).';
X = [ones(numel(XArray), 1), XArray, XArray.^2];
Beta = ((X.'*X) \ X.') * XCorrMax(XArray);
Shift = ZeroImageInd + Beta(2) / (2*Beta(3));
PolyFitFunction = @(R) Beta(1) + Beta(2)*R + Beta(3)*R.^2;

% Create arrays of the polynomial fits to use for visualization.
XArrayDense = linspace(XArray(1), XArray(end), max(ImSize));
FitAtPeak = PolyFitFunction(XArrayDense);

% Populate the CorrData struct.
CorrData.XCorrMax = XCorrMax;
CorrData.FitAtPeak = FitAtPeak;

% Display line sections through the integer location of the
% cross-correlation, overlain on the fit along those lines.
if Params.PlotFlag
    PlotFigure = findobj('Tag', 'CorrWindow');
    if isempty(PlotFigure)
        PlotFigure = figure('Tag', 'CorrWindow');
    end
    clf(PlotFigure);
    PlotAxes = axes(PlotFigure);
    plot(PlotAxes, ZeroImageInd - XArray, XCorrMax(XArray), 'x')
    hold(PlotAxes, 'on')
    plot(PlotAxes, ZeroImageInd - XArrayDense, FitAtPeak)
    title(PlotAxes, 'Max. Correlation') 
end

% Restore the warning state.
if Params.SuppressWarnings
    warning(WarningState)
end


end