function [] = genOverlayResults(obj)
%genOverlayResults generates some two-color overlay results.
% This method will generate several results pertaining to two-color
% overlays of super-resolution data.


% Generate a list of all two color overlay images.
OverlayFileStruct = ...
    dir(fullfile(obj.SaveBaseDir, '*GaussianOverlay*'));
if isempty(OverlayFileStruct)
    % The overlays don't exist, so we can't produce results for the
    % overlays.
    return
end
OverlayFileNames = {OverlayFileStruct.name};
NOverlays = numel(OverlayFileNames);

% Loop through all overlay images and compute the shift between the color
% channels.
ImageShift = zeros(NOverlays, 2);
for ii = 1:NOverlays
    % Display a status message in the command line.
    if obj.Verbose
        fprintf(['Publish.performFullAnalysis(): ', ...
            'Computing shift for overlay image %i of %i\n'], ...
            ii, NOverlays);
    end
    
    % Load the image into the workspace.
    OverlayImage = ...
        imread(fullfile(obj.SaveBaseDir, OverlayFileNames{ii}));
    
    % Compute overlay stats. for this overlay image.
    % NOTE: OverlayImage should be a green-magenta overlay, but we can 
    %       ignore the third color channel because that information is 
    %       already contained in the first channel.
    if (obj.Verbose < 3)
        % Turn off some warnings that come from findStackOffset().
        warning('off')
    end
    [~, SubPixelShift] = MIC_Reg3DTrans.findStackOffset(...
        OverlayImage(:, :, 2), OverlayImage(:, :, 1), ...
        [size(OverlayImage, [1, 2])-1, 0], [], [], 0, 0);
    ImageShift(ii, :) = SubPixelShift(1:2);
    if (obj.Verbose < 3)
        warning('on')
    end
end

% Concatenate the max. correlation coefficients from the image in each 
% overlay image (keeping the labels separate) so that we have just one 
% array containing the max. correlation coefficients from all of the
% overlay images produced.  Repeat for the registration errors.
CellLabelFields = fieldnames(obj.ResultsStruct);
FirstLabelFields = ...
    CellLabelFields(contains(CellLabelFields, 'Label_01'));
SecondLabelFields = ...
    CellLabelFields(contains(CellLabelFields, 'Label_02'));
ConcatenatedMaxCorr = cell(numel(FirstLabelFields), 2);
ConcatenatedRegError = cell(numel(FirstLabelFields), 2);
for ii = 1:numel(FirstLabelFields)
    % Isolate the sub-structure for this cell/label pair.
    FirstLabelStructure = obj.ResultsStruct.(FirstLabelFields{ii});
    FirstLabelFieldName = fieldnames(FirstLabelStructure);
    SecondLabelStructure = obj.ResultsStruct.(SecondLabelFields{ii});
    SecondLabelFieldName = fieldnames(SecondLabelStructure);
    
    % Add the MaxCorr array for the ii-th cell/label pair to the 
    % concatenated array, assuming only one dataset existed for this
    % cell/label pair. Repeat for the RegError array.
    ConcatenatedMaxCorr{ii, 1} = ...
        FirstLabelStructure.(FirstLabelFieldName{1}).MaxCorr;
    ConcatenatedRegError{ii, 1} = ...
        FirstLabelStructure.(FirstLabelFieldName{1}).RegError;
    ConcatenatedMaxCorr{ii, 2} = ...
        SecondLabelStructure.(SecondLabelFieldName{1}).MaxCorr;
    ConcatenatedRegError{ii, 2} = ...
        SecondLabelStructure.(SecondLabelFieldName{1}).RegError;
end

% Generate the overlay plots across all cells.
SRPixelSize = obj.SMF.Data.PixelSize / obj.SRImageZoom;
obj.makeOverlayPlots(ImageShift, ...
    ConcatenatedRegError, ...
    ConcatenatedMaxCorr, ...
    SRPixelSize, obj.SMF.Data.PixelSize, ...
    obj.SaveBaseDir)


end