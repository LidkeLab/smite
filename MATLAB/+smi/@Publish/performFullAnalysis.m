function [] = performFullAnalysis(obj)
%performFullAnalysis is the main run method for the smi.Publish class.
% This method is the main run method for the smi.Publish class, meaning
% that it can be used to perform the standard analysis expected for use of
% this class.


% Define the results directory, which will be in the top level directory 
% obj.CoverslipDir for easy access.
if isempty(obj.SaveBaseDir)
    obj.SaveBaseDir = fullfile(obj.CoverslipDir, 'Results');
end

% Determine the names of the sub-directories of interest within
% obj.CoverslipDir.  These correspond to individual cells imaged during the
% experiment.
CellNames = smi_helpers.getDirectoryNames(obj.CoverslipDir, 'Cell*');
NCells = numel(CellNames);
if (obj.Verbose > 1)
    fprintf(['Publish.performFullAnalysis(): %i ', ...
        'cell directories found:\n'], NCells)
    for ii = 1:NCells
        fprintf('\t%s\n', CellNames{ii})
    end
elseif (obj.Verbose == 1)
    fprintf(['Publish.performFullAnalysis(): analyzing %i ', ...
        'cell directories...\n'], NCells)
end

% Loop through the cell directories and analyze the contents.
for ii = 1:NCells
    if (obj.Verbose > 1)
        fprintf(['Publish.performFullAnalysis(): analyzing ', ...
            'cell %i of %i...\n'], ii, NCells)
    end
    obj.processCell(CellNames{ii});
end

% Generate misc. stats. comparing all cells that were analyzed.
if obj.GenerateOverlayStats
    % Generate a list of all two color overlay images.
    OverlayFileStruct = dir(...
        fullfile(obj.SaveBaseDir, '*HistogramOverlay*'));
    if isempty(OverlayFileStruct)
        % The histogram overlays don't exist, so we can't
        % produce results for the overlays.
        return
    end
    OverlayFileNames = {OverlayFileStruct.name};
    NOverlays = numel(OverlayFileNames);
    
    % Loop through all overlay images and compute the shift
    % between the color channels using findshift().
    ImageShift = zeros(NOverlays, 2); % pre-allocate
    for ii = 1:NOverlays
        % Display a status message in the command line.
        fprintf('Computing shift for overlay image %i\n', ii);
        
        % Load the image into the workspace.
        OverlayImage = imread(...
            fullfile(obj.SaveBaseDir, OverlayFileNames{ii}));
        
        % Compute overlay stats. for this overlay image.
        % NOTE: OverlayImage should be a green-magenta overlay,
        %       but we can ignore the third color channel
        %       because that information is already contained
        %       in the first channel.
        ImageShift(ii, :) = findshift(...
            OverlayImage(:, :, 2), OverlayImage(:, :, 1), ...
            'iter').';
    end
    
    % Concatenate the max. correlation coefficients from the
    % image in each overlay image (keeping the labels separate)
    % so that we have just one array containing the max.
    % correlation coefficients from all of the overlay images
    % produced.  Repeat for the registration errors.
    CellLabelFields = fieldnames(obj.PublishedResultsStruct);
    FirstLabelFields = CellLabelFields(...
        contains(CellLabelFields, 'Label_01'));
    SecondLabelFields = CellLabelFields(...
        contains(CellLabelFields, 'Label_02'));
    ConcatenatedMaxCorr = cell(numel(FirstLabelFields), 2);
    ConcatenatedRegError = cell(numel(FirstLabelFields), 2);
    for ii = 1:numel(FirstLabelFields)
        % Isolate the sub-structure for this cell/label pair.
        FirstLabelStructure = obj.PublishedResultsStruct. ...
            (FirstLabelFields{ii});
        FirstLabelFieldName = fieldnames(...
            FirstLabelStructure);
        SecondLabelStructure = obj.PublishedResultsStruct. ...
            (SecondLabelFields{ii});
        SecondLabelFieldName = fieldnames(...
            SecondLabelStructure);
        
        % Add the MaxCorr array for the ii-th cell/label pair
        % to the concatenated array, assuming only one dataset
        % existed for this cell/label pair. Repeat for the
        % RegError array.
        ConcatenatedMaxCorr{ii, 1} = FirstLabelStructure. ...
            (FirstLabelFieldName{1}).MaxCorr;
        ConcatenatedRegError{ii, 1} = FirstLabelStructure. ...
            (FirstLabelFieldName{1}).RegError;
        ConcatenatedMaxCorr{ii, 2} = SecondLabelStructure. ...
            (SecondLabelFieldName{1}).MaxCorr;
        ConcatenatedRegError{ii, 2} = SecondLabelStructure. ...
            (SecondLabelFieldName{1}).RegError;
    end
    
    % Generate the overlay plots across all cells.
    obj.genOverlayPlots(ImageShift, ...
        ConcatenatedRegError, ...
        ConcatenatedMaxCorr, ...
        obj.SMF.SRPixelSize, obj.SMF.PixelSize, ...
        obj.SaveBaseDir)
end

% Indicate completion of the analysis/generation of results.
fprintf('Results have been published to %s \n', ...
    obj.CoverslipDir);


end