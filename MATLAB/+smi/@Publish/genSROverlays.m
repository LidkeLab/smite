function [] = genSROverlays(ResultsCellDir, SaveDir)
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

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


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
% NOTE: pre-allocation assumes 5120x5120 images.
GaussianImages = zeros(5120, 5120, NLabels);
HistogramImages = GaussianImages;
CircleImages = [];
for ii = 1:NLabels
    % Create a list of sub-directories under the current label (there could
    % be multiple for a given label, e.g. an extra for a photobleaching
    % round of imaging).
    DatasetDirNames = smi_helpers.getDirectoryNames(...
        fullfile(ResultsCellDir, LabelDirNames{ii}), 'Data*');
    
    % If more than two datasets exists for this label, throw an error
    % (we can have two: one desired result, one photobleaching result)
    % since it's not clear which dataset to use.
    if (numel(DatasetDirNames) > 2)
        error('More than two datasets exist for %s', LabelDirNames{ii})
    end
    
    % Load our images into the appropriate arrays.
    for jj = 1:numel(DatasetDirNames)
        % Ensure that we skip results from photobleaching rounds.
        if contains(DatasetDirNames{jj}, 'bleach', 'IgnoreCase', true)
            continue
        end
        
        % Create the appropriate filepaths and read in the images.
        FileDirectory = fullfile(ResultsCellDir, LabelDirNames{ii});
        FileNameCircle = sprintf('%s_CircleImage.png', ...
            DatasetDirNames{jj});
        FileNameGaussian = sprintf('%s_GaussImage.png', ...
            DatasetDirNames{jj});
        FileNameHistogram = sprintf('%s_HistImage.png', ...
            DatasetDirNames{jj});
        CircleImages = cat(3, CircleImages, sum(imread(fullfile(...
            FileDirectory, DatasetDirNames{jj}, FileNameCircle)), 3));
        GaussianImages(:, :, ii) = imread(...
            fullfile(FileDirectory, FileNameGaussian));
        HistogramImages(:, :, ii) = imread(...
            fullfile(FileDirectory, DatasetDirNames{jj}, FileNameHistogram));
    end
end

% Generate our color overlay images (3 channel images).
[OverlayImageGaussian, ColorOrderTagGaussian] = ...
    smi_vis.GenerateImages.overlayNImages(GaussianImages);
[OverlayImageHistogram, ColorOrderTagHistogram] = ...
    smi_vis.GenerateImages.overlayNImages(HistogramImages);
[OverlayImageCircle, ColorOrderTagCircle] = ...
    smi_vis.GenerateImages.overlayNImages(CircleImages);

% Save the overlay images in the top level directory.
CellName = ResultsCellDir(regexp(ResultsCellDir, 'Cell*'):end);
CellNameClean = erase(CellName, '_');
OverlayImageGaussianName = sprintf('%s_GaussianOverlay_%s.png', ...
    CellNameClean, ColorOrderTagGaussian);
imwrite(OverlayImageGaussian, ...
    fullfile(SaveDir, OverlayImageGaussianName));
OverlayImageHistogramName = sprintf('%s_HistogramOverlay_%s.png', ...
    CellNameClean, ColorOrderTagHistogram);
imwrite(OverlayImageHistogram, ...
    fullfile(SaveDir, OverlayImageHistogramName));
OverlayImageCircleName = sprintf('%s_CircleOverlay_%s.png', ...
    CellNameClean, ColorOrderTagCircle);
imwrite(OverlayImageCircle, ...
    fullfile(SaveDir, OverlayImageCircleName));


end