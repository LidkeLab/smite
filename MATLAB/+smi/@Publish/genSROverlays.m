function [] = genSROverlays(ResultsCellDir, SaveDir, AnalysisID, ...
    Mask, MaskName)
%genSROverlays generates various types of  SR overlay images.
% This method will generate different types of overlay images for
% super-resolution (SR) data.  A multicolor overlay of the Gaussian SR
% images will be generated and saved as a .png file with the identifier
% _GaussianOverlay_ placed in the filename.  A multicolor overlay of the
% histogram SR images will be generated and saved as a .png file with the
% identifier _HistogramOverlay_ in the filename.  A multicolor overlay
% containing circles with radii proportional to the average localization
% precision (the generalized mean between the X and Y precisions) for each
% localization will be generated and saved as a .png file with the
% identifier _CircleOverlay_ in the filename.
%
% INPUTS:
%   ResultsCellDir: Directory containing the sub-directories for each label
%   SaveDir: Directory in which the resulting .png overlay images will be
%            saved.
%   AnalysisID: Analysis ID used to help identify the correct files.
%               (Default = '')
%   Mask: Optional mask to be applied to the overlays. (Logical array, same
%         aspect ratio as SR images).
%   MaskName: Optional identifier for the mask (e.g., '20nm').
%             (Default = '')

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Set defaults/modify inputs.
if ~exist('AnalysisID', 'var')
    AnalysisID = '';
elseif ~isempty(AnalysisID)
    AnalysisID = ['_', AnalysisID];
end
if (~exist('Mask', 'var') || isempty(Mask))
    Mask = true;
end
if (~exist('MaskName', 'var') || isempty(MaskName))
    MaskName = '';
end

% Get the names of the directories containing the results for each label,
% throwing warnings if appropriate.
LabelDirNames = smi_helpers.getDirectoryNames(ResultsCellDir, 'Label*');
NLabels = numel(LabelDirNames);
if (NLabels > 4)
    warning('genSROverlay: 3 labels max')
    return
end
if (NLabels < 2)
    warning('genSROverlay: 2 labels min')
    return
end

% Load the Gaussian and histogram images for each label and store them in a
% single array.
GaussianImages = [];
CircleImages = [];
for ii = 1:NLabels
    % Create a list of sub-directories under the current label (there could
    % be multiple for a given label, e.g. an extra for a photobleaching
    % round of imaging).
    DatasetDirNames = smi_helpers.getDirectoryNames(...
        fullfile(ResultsCellDir, LabelDirNames{ii}), ...
        sprintf('Data*%s', AnalysisID));

    % If more than two datasets exists for this label, throw an error
    % (we can have two: one desired result, one photobleaching result)
    % since it's not clear which dataset to use.
    if (numel(DatasetDirNames) > 2)
        error('More than two datasets exist for %s', LabelDirNames{ii})
    end

    % Load our images into the appropriate arrays.
    DatasetDirNames = ...
        DatasetDirNames(contains(DatasetDirNames, AnalysisID));
    DatasetDirNames = ...
        DatasetDirNames(~contains(DatasetDirNames, 'bleach', ...
        'IgnoreCase', true));
    DatasetNames = erase(DatasetDirNames, AnalysisID);
    for jj = 1:numel(DatasetDirNames)
        % Create the appropriate filepaths and read in the images.
        FileDirectory = fullfile(ResultsCellDir, LabelDirNames{ii});
        FileNameCircle = sprintf('%s%s_CircleImage.png', ...
            DatasetNames{jj}, AnalysisID);
        FileNameGaussian = sprintf('%s%s_GaussImage.png', ...
            DatasetNames{jj}, AnalysisID);
        CircleImages = cat(3, CircleImages, sum(imread(fullfile(...
            FileDirectory, DatasetDirNames{jj}, FileNameCircle)), 3));
        GaussianImages = cat(3, GaussianImages, ...
            sum(imread(fullfile(FileDirectory, FileNameGaussian)), 3));
    end
end

% Generate our color overlay images (3 channel images).
[OverlayImageGaussian, ColorOrderTagGaussian] = ...
    smi_vis.GenerateImages.overlayNImages(GaussianImages);
[OverlayImageCircle, ColorOrderTagCircle] = ...
    smi_vis.GenerateImages.overlayNImages(CircleImages);

% Save the overlay images in the top level directory.
CellName = ResultsCellDir(regexp(ResultsCellDir, 'Cell*'):end);
CellNameClean = erase(CellName, '_');
OverlayImageGaussianName = sprintf('%s%s_GaussianOverlay_%s.png', ...
    CellNameClean, AnalysisID, ColorOrderTagGaussian);
imwrite(OverlayImageGaussian, fullfile(SaveDir, OverlayImageGaussianName));
OverlayImageCircleName = sprintf('%s%s_CircleOverlay_%s.png', ...
    CellNameClean, AnalysisID, ColorOrderTagCircle);
imwrite(OverlayImageCircle, fullfile(SaveDir, OverlayImageCircleName));

% Prepare and save masked overlays.
if ~all(Mask(:))
    Mask = imresize(Mask, size(OverlayImageGaussian, 1:2));
    MaskedGaussian = Mask .* OverlayImageGaussian;
    MaskedGaussianName = sprintf('%s%s%sMask_GaussianOverlay_%s.png', ...
        CellNameClean, AnalysisID, MaskName, ColorOrderTagGaussian);
    imwrite(MaskedGaussian, fullfile(SaveDir, MaskedGaussianName));
    Mask = imresize(Mask, size(OverlayImageCircle, 1:2));
    MaskedCircle = Mask .* OverlayImageCircle;
    MaskedCircleName = sprintf('%s%s%sMask_CircleOverlay_%s.png', ...
        CellNameClean, AnalysisID, MaskName, ColorOrderTagCircle);
    imwrite(MaskedCircle, fullfile(SaveDir, MaskedCircleName));
end


end