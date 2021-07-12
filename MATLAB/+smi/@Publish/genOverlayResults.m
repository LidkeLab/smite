function [] = genOverlayResults(obj)
%genOverlayResults generates some two-color overlay results.
% This method will generate several results pertaining to two-color
% overlays of super-resolution data.


% Generate a list of all two color overlay images.
OverlayFileStruct = dir(fullfile(obj.SaveBaseDir, '*GaussianOverlay*'));
if isempty(OverlayFileStruct)
    % The overlays don't exist, so we can't produce results for the
    % overlays.
    return
end
OverlayFileNames = {OverlayFileStruct.name}.';
NOverlays = numel(OverlayFileNames);

% Loop through all overlay images and compute the shift between the color
% channels.
ImageShift = zeros(NOverlays, 2);
ImageShiftLocal = cell(NOverlays, 4);
SubROIDivisor = 10;
MaxOffset = 20;
for ii = 1:NOverlays
    % Display a status message in the command line.
    if obj.Verbose
        fprintf(['Publish.performFullAnalysis(): ', ...
            'Computing shift for overlay image %i of %i...\n'], ...
            ii, NOverlays);
    end
    
    % Load the image into the workspace.
    OverlayImage = ...
        imread(fullfile(obj.SaveBaseDir, OverlayFileNames{ii}));
    
    % Compute the shift between labels for this image.
    % NOTE: OverlayImage should be a green-magenta overlay, but we can 
    %       ignore the third color channel because that information is 
    %       already contained in the first channel.
    if (obj.Verbose < 3)
        % Turn off some warnings that come from findStackOffset().
        warning('off')
    end
    PixelShift = MIC_Reg3DTrans.findStackOffset(...
        OverlayImage(:, :, 2), OverlayImage(:, :, 1), ...
        [MaxOffset, MaxOffset, 0], [], [], 0, 0);
    ImageShift(ii, :) = PixelShift(1:2);
    
    % Compute the local shift between labels for this image.
    SubROISize = ceil(size(OverlayImage, [1, 2]) / SubROIDivisor);
    [PixelOffsets, SubPixelOffsets, ImageROIs, ImageStats] = ...
        obj.estimateLocalShifts(...
        OverlayImage(:, :, 2), OverlayImage(:, :, 1), SubROISize, [], 1);
    ImageShiftLocal{ii, 1} = PixelOffsets;
    ImageShiftLocal{ii, 2} = SubPixelOffsets;
    ImageShiftLocal{ii, 3} = ImageROIs;
    ImageShiftLocal{ii, 4} = ImageStats;
    if (obj.Verbose < 3)
        warning('on')
    end
end

% Concatenate the max. correlation coefficients from the image in each 
% overlay image (keeping the labels separate) so that we have just one 
% array containing the max. correlation coefficients from all of the
% overlay images produced.  Repeat for the registration errors.
MakeShiftVsCorrPlots = ...
    all(isfield(obj.ResultsStruct, {'RegError', 'MaxCorr'}));
if MakeShiftVsCorrPlots
    ConcatenatedRegError = [{obj.ResultsStruct(:, 1).RegError}.',  ...
        {obj.ResultsStruct(:, 2).RegError}.'];
    ConcatenatedMaxCorr = [{obj.ResultsStruct(:, 1).MaxCorr}.', ...
        {obj.ResultsStruct(:, 2).MaxCorr}.'];
else
    ConcatenatedRegError = [];
    ConcatenatedMaxCorr = [];
end

% Create and save a structure which contains all of the data used to
% generate these plots.
OverlayInfoStruct.ImageShift = ImageShift;
OverlayInfoStruct.ImageShiftLocal = ImageShiftLocal;
OverlayInfoStruct.RegError = ConcatenatedRegError;
OverlayInfoStruct.MaxCorr = ConcatenatedMaxCorr;
save(fullfile(obj.SaveBaseDir, 'OverlayInfoStruct.mat'), ...
    'OverlayInfoStruct');

% % Estimate an affine transform from the coordinates of each label.
% ResultsStructDir = fullfile(obj.SaveBaseDir, 'ResultsStructs');
% Label1Results = dir(fullfile(ResultsStructDir, '*Label_01*'));
% Label1Paths = fullfile(ResultsStructDir, ...
%     {Label1Results(~[Label1Results.isdir]).name});
% NOverlays = numel(Label1Paths);
% AffineTransforms = cell(NOverlays, 1);
% for ii = 1:NOverlays
%     % Display a status message in the command line.
%     if obj.Verbose
%         fprintf(['Publish.performFullAnalysis(): ', ...
%             'Computing affine transform for overlay image %i of %i...\n'], ...
%             ii, NOverlays);
%     end
%     
%     % Load the results structs into the workspace.
%     load(Label1Paths{ii}, 'SMD')
%     SMD1 = SMD;
%     Label2Path = strrep(Label1Paths{ii}, 'Label_01', 'Label_02');
%     if exist(Label2Path, 'file')
%         load(Label2Path, 'SMD')
%         SMD2 = SMD;
%     else
%         continue
%     end
%     
%     % Compute an affine transform between the coordinates.
%     AffineTransforms{ii} = ...
%         smi_stat.findCoordAffine([SMD1.X, SMD1.Y], [SMD2.X, SMD2.Y], 50);
% end
% save(fullfile(obj.SaveBaseDir, 'AffineTransforms.mat'), ...
%     'AffineTransforms', 'Label1Paths')

% Generate the overlay plots across all cells.
if MakeShiftVsCorrPlots
    SRPixelSize = obj.SMF.Data.PixelSize / obj.SRImageZoom;
    obj.makeOverlayPlots(ImageShift, ...
        ConcatenatedRegError, ...
        ConcatenatedMaxCorr, ...
        SRPixelSize, obj.SMF.Data.PixelSize, ...
        obj.SaveBaseDir)
end


end