function [] = loadFiducials(obj)
%loadFiducials loads fiducial files and sets associated class properties.
% This method will load the fiducial files specified by the fields in
% obj.SMF.Data.  Once the files are loaded, this method will set various
% other class properties based on the fiducial images.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Load the first fiducial file and set obj.FiducialROI to a default if
% needed.
% NOTE: I'll be averaging over the "time" dimension of all fiducials used.
NFiles = numel(obj.SMF.Data.FileName);
if (obj.Verbose > 1)
    fprintf(['\tChannelRegistration.loadFiducials(): ', ...
        'Loading %i fiducial files from %s...\n'], ...
        NFiles, obj.SMF.Data.FileDir)
end
LoadData = smi_core.LoadData;
[~, TempImage] = ...
    LoadData.loadRawData(obj.SMF, 1, obj.SMF.Data.DataVariable);
TempImage = mean(TempImage, 3);
ImageSize = size(TempImage);
FullROI = [1, 1, ImageSize(1:2)];
if isempty(obj.SplitFormat)
    % When obj.SplitFormat is empty, the user must manually define
    % obj.FiducialROI!
    if isempty(obj.FiducialROI)
        error(['loadFiducials(): If you set obj.SplitFormat = [], ', ...
            'you must also manually define obj.FiducialROI'])
    end
else
    % We must define (or re-define) obj.FiducialROI based on
    % obj.SplitFormat and the number of files present.
    if (NFiles == 1)
        % Split the file into the specified ROIs.
        if (obj.SplitFormat == 1)
            % In this case, the user hasn't selected an appropriate
            % SplitFormat for the single file case, so we'll define a
            % (hopefully) useful default.
            if (ImageSize(1) > ImageSize(2))
                obj.SplitFormat = [1; 2];
            else
                obj.SplitFormat = [1, 2];
            end
            if (obj.Verbose > 0)
                warning(['loadFiducials(): obj.SplitFormat ', ...
                    'reset to [%i, %i]'], ...
                    obj.SplitFormat(1), obj.SplitFormat(2))
            end
        end
        obj.FiducialROI = obj.convertSplitFormatToROIs(...
            FullROI, obj.SplitFormat);
    elseif (NFiles > 1)
        % For multiple files, we'll set FiducialROI s.t. each image is used
        % in it's entirety.  If one of the images is too small, everything
        % will crash, but I'll assume that's the users fault!
        obj.SplitFormat = 1;
        obj.FiducialROI = FullROI;
    else
        error('loadFiducials(): no files defined in obj.SMF.Data.FileName')
    end
end

% Load the remaining fiducial images/split the fiducial into the specified
% ROIs for later use.
NROIs = size(obj.FiducialROI, 1);
if (NROIs == 1)
    % All fiducial images are stored in separate files, so we'll have to
    % load them one at a time.
    NFiducials = NFiles;
    FiducialImages = zeros([obj.FiducialROI(1, 3:4), NFiducials]);
    FiducialImages(:, :, 1) = TempImage(...
        obj.FiducialROI(1):obj.FiducialROI(3), ...
        obj.FiducialROI(2):obj.FiducialROI(4));
    for ii = 2:NFiducials
        [~, TempImage] = ...
            LoadData.loadRawData(obj.SMF, ii, obj.SMF.Data.DataVariable);
        FiducialImages(:, :, ii) = mean(...
            TempImage(obj.FiducialROI(1):obj.FiducialROI(3), ...
            obj.FiducialROI(2):obj.FiducialROI(4), :), ...
            3);
    end
else
    % There is only one fiducial image, which we'll need to split by ROI.
    if (obj.Verbose > 1)
        fprintf(['\tChannelRegistration.loadFiducials(): ', ...
            'Splitting fiducial image into channels.\n'])
    end
    NFiducials = NROIs;
    FiducialImages = zeros(...
        [obj.FiducialROI(1, 3:4)-obj.FiducialROI(1, 1:2)+1, NFiducials]);
    for ii = 1:NFiducials
        FiducialImages(:, :, ii) = TempImage(...
            obj.FiducialROI(ii, 1):obj.FiducialROI(ii, 3), ...
            obj.FiducialROI(ii, 2):obj.FiducialROI(ii, 4));
    end
end
obj.FiducialImages = FiducialImages;


end