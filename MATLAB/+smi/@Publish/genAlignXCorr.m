function [XCorrStruct] = genAlignXCorr(AlignRegStruct, SaveDir)
%genAlignXCorr generates xcorr curve data for AlignRegStruct.
% This method will create various plot(s) and analysis data related to the
% brightfield registration process of the acquistion.  A plot will be
% created containing information about the cross-correlation process used
% to correct for sample drift during brightfield registration.  This plot
% will be saved in the SaveDir.
%
% INPUTS:
%   AlignRegStruct: Structured array containing the brightfield
%                   registration data structures for each dataset.
%   SaveDir: Directory associated with Dataset in which the analysis
%            results will be saved.
%
% OUTPUTS:
%   XCorrStruct: Structured array containing the cross-correlation
%                information computed in this method.

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Pre-allocate the XCorrStruct structured array.
XCorrStruct = struct();

% Determine the number of sequences contained in this dataset.
NSequences = numel(AlignRegStruct);

% Grab the reference stack from the first dataset in the AlignRegStruct
% structure (it doesn't have to be from the first dataset, they should all
% contain the same reference stack).
ReferenceStack = AlignRegStruct(1).Data.ReferenceStack;

% Loop through each dataset and grab the stack taken directly before the
% dataset was collected.
XCorrDirName = fullfile(SaveDir, 'XCorrPlots');
if ~isfolder(XCorrDirName)
    mkdir(XCorrDirName);
end
MaxCorr = zeros(numel(AlignRegStruct), 1); % pre-allocate
MaxCorrFit = zeros(numel(AlignRegStruct), 3); % pre-allocate
OffsetFitSuccess = zeros(NSequences, 1, 'logical'); % pre-allocate
MaxIterReached  = zeros(NSequences, 1, 'logical'); % pre-allocate
for ii = 1:NSequences
    % Determine if the last round of registration was succesful (to
    % emphasize on the plot which points might be suspicious/unreliable).
    OffsetFitSuccess(ii) = all(...
        AlignRegStruct(ii).Data.OffsetFitSuccess);
    MaxIterReached(ii) = AlignRegStruct(ii).Data.MaxIterReached;
    
    % Grab the current stack from the ii-th dataset.
    CurrentStack = AlignRegStruct(ii).Data.CurrentStack;
    
    % Grab some misc. parameters needed from the AlignRegStruct.
    % NOTE: these might change from one dataset to the next, so we need
    %       this inside the for loop.
    ZStackMaxDev = AlignRegStruct(ii).Attributes.ZStack_MaxDev;
    ZStackStep = AlignRegStruct(ii).Attributes.ZStack_Step;
    IsInitialRegistration = ...
        AlignRegStruct(ii).Attributes.IsInitialRegistration;
    
    % Define the indices of the ReferenceStack which correspond to the
    % CurrentStack (CurrentStack might be smaller than ReferenceStack, but
    % they share a center image along the z-direction of the stack).
    FocalInd = 1 + ZStackMaxDev/ZStackStep;
    StackSteps = ZStackMaxDev/ZStackStep;
    ZStackRefInds = (FocalInd-StackSteps):(FocalInd+StackSteps);
    
    % Isolate the portion of the ReferenceStack that we'll want to compare
    % our CurrentStack to.
    if IsInitialRegistration
        % If IsInitialRegistration flag was set, the CurrentStack is
        % already the same size as the ReferenceStack so we don't need to
        % isolate the central portion.
        ReferenceSubStack = ReferenceStack;
    else
        ReferenceSubStack = ReferenceStack(:, :, ZStackRefInds);
    end
    
    % Re-compute the scaled cross-correlations between the two stacks.
    CorrParams.PlotFlag = true;
    CorrParams.SuppressWarnings = true;
    NIter = 3;
    [~, ~, CorrData] = smi_stat.findOffsetIter(...
        ReferenceSubStack, CurrentStack, NIter, [], CorrParams);
    
    % Generate the x, y, z cross-correlation fitting plots and save them.
    FigHandle = findobj('Tag', 'CorrWindow');
    saveas(FigHandle, ...
        fullfile(XCorrDirName, sprintf('XCorrSequence%i', ii)), 'png');
    saveas(FigHandle, ...
        fullfile(XCorrDirName, sprintf('XCorrSequence%i', ii)), 'fig');
    close(FigHandle);
    
    % Determine the maximum value of the scaled cross-correlation between
    % the two stacks.
    MaxCorr(ii) = max(CorrData.XCorr3D(:));
    MaxCorrFit(ii, :) = [max(CorrData.XFitAtPeak), ...
        max(CorrData.YFitAtPeak), ...
        max(CorrData.ZFitAtPeak)];
end

% Save the above computed cross-correlation information in the
% XCorrStruct.
XCorrStruct.MaxCorr = MaxCorr;
XCorrStruct.MaxCorrFit = MaxCorrFit;
XCorrStruct.OffsetFitSuccess = OffsetFitSuccess;
XCorrStruct.MaxIterReached = MaxIterReached;

% Determine which sequences did not have a succesful registration and
% create an array of their max-correlations/fits.
% NOTE: I've added the NaN points to ensure something is always in the
%       plot, which will ensure we don't get legend() warnings (NaN's don't
%       actually show up in the plot, but are still given handles visible
%       to legend()).
SequenceArray = (1:NSequences).';
SequenceArrayUnsuccessfulFit = [NaN; SequenceArray(~OffsetFitSuccess)];
SequenceArrayMaxIterReached = [NaN; SequenceArray(MaxIterReached)];
MaxCorrUnsuccessfulFit = [NaN; MaxCorr(~OffsetFitSuccess)];
MaxCorrMaxIterReached = [NaN; MaxCorr(MaxIterReached)];

% Plot the maximum correlation coefficients found throughout the
% acquisition.
FigureHandle = figure();
PlotAxes = subplot(2, 1, 1, 'Parent', FigureHandle);
plot(PlotAxes, SequenceArray, MaxCorr, 'x')
hold(PlotAxes, 'on');
plot(PlotAxes, [1, NSequences], ones(2, 1) * mean(MaxCorr), ':')
plot(PlotAxes, SequenceArrayUnsuccessfulFit, MaxCorrUnsuccessfulFit, ...
    'o', 'MarkerSize', 10, 'LineWidth', 2)
plot(PlotAxes, SequenceArrayMaxIterReached, MaxCorrMaxIterReached, ...
    's', 'MarkerSize', 10, 'LineWidth', 2)
if (min(MaxCorr) ~= max(MaxCorr))
    % Change the YTicks displayed on the plot, unless the min and max are
    % the same value.
    PlotAxes.YTick = [min(MaxCorr), ...
        (max(MaxCorr)+min(MaxCorr)) / 2, ...
        max(MaxCorr)];
end
legend(PlotAxes, ...
    {'Corr. Coeff.', sprintf('Mean = %.4f', mean(MaxCorr)), ...
    'Failed Registration Fit', 'Max. Registration Iter. Reached'}, ...
    'Location', 'best')
xlabel(PlotAxes, 'Sequence Number')
ylabel(PlotAxes, 'Max. Corr.')

% Plot the maximum values of the X, Y, and Z line fits through the peak of
% the cross-correlation volume.
PlotAxes = subplot(2, 1, 2, 'Parent', FigureHandle);
plot(PlotAxes, SequenceArray, MaxCorrFit(:, 1), 'x')
hold(PlotAxes, 'on');
plot(PlotAxes, SequenceArray, MaxCorrFit(:, 2), 'x')
plot(PlotAxes, SequenceArray, MaxCorrFit(:, 3), 'x')
PlotAxes.YTick = [min(MaxCorrFit(:)), ...
    (max(MaxCorrFit(:))+min(MaxCorrFit(:))) / 2, ...
    max(MaxCorrFit(:))];
legend(PlotAxes, {'X Fit', 'Y Fit', 'Z Fit'}, 'Location', 'best');
xlabel(PlotAxes, 'Sequence Number')
ylabel(PlotAxes, 'Max. Corr. Fit')

% Save the plot in the SaveDir and then close the figure.
saveas(FigureHandle, fullfile(SaveDir, 'AlignRegXCorrMaxima.png'), 'png');
close(FigureHandle);


end