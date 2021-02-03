function [RegistrationTransform] = findTransform(obj)
%findTransform finds a channel registration transform.
% This method will find a channel registration transform object that is
% intended to register coordinates from/features in the fiducial files
% specified by obj.SMF.Data.FileNames.
%
% OUTPUTS:
%   RegistrationTransform: A cell array of MATLAB tform objects (this
%                          will be the same cell array stored in
%                          obj.RegistrationTransform). (cell array)
%
% REQUIRES:
%   MATLAB 2013b or later
%   Image Processing Toolbox (to use fitgeotrans())

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Load the fiducial images (if needed). If the fiducials are already
% present, we should make sure obj.FiducialROI is set.
if isempty(obj.FiducialImages)
    obj.loadFiducials();
elseif isempty(obj.FiducialROI)
    obj.FiducialROI =  [1, 1, size(obj.FiducialImages, [1, 2])];
end
NFiducials = size(obj.FiducialImages, 3);

% Perform the gain and offset correction on each of the fiducial images
% (this probably doesn't matter for image transforms, but I'm going to do
% this up here anyways because it looks cleaner).  If
% obj.AutoscaleFiducials is set, we'll attempt to autoscale the fiducials
% instead of doing a proper gain/offset correction.
if obj.AutoscaleFiducials
    % Perform a full-scale histogram stretch on each image.
    ScaledData = obj.FiducialImages;
    for ii = 1:NFiducials
        CurrentImage = obj.FiducialImages(:, :, ii);
        ScaledData(:, :, ii) = ...
            (CurrentImage-min(CurrentImage(:))) ...
            ./ max(max(CurrentImage-min(CurrentImage(:))));
    end
    
    % Multiply by an extra factor so that the counts are somewhat realistic
    % for a camera (having the range [0, 1] messes up the fitting,
    % presumably because it seems like 0 photons but I'm not sure).
    ScaledData = ScaledData*100;
else
    % NOTE: We'll need to update this in the future to specify RawDataROI
    %       and CalibrationROI.
    [ScaledData] = smi_core.DataToPhotons.convertToPhotons(...
        obj.FiducialImages, ...
        obj.SMF.Data.CameraGain, obj.SMF.Data.CameraOffset, ...
        []);
end

% Proceed based on the setting of obj.TransformationBasis (which defines
% whether or not we need to find localizations from the fiducial data).
RegistrationTransform = cell(NFiducials, 1);
switch obj.TransformationBasis
    case 'coordinates'
        % Use smi_core.LocalizeData to find localizations in the fiducial
        % images.
        LocalizeData = smi_core.LocalizeData([], obj.SMF);
        FiducialSMD = cell(NFiducials, 1);
        FiducialSMDPreThresh = FiducialSMD;
        for ii = 1:NFiducials
            LocalizeData.ScaledData = ScaledData(:, :, ii);
            [FiducialSMD{ii}, FiducialSMDPreThresh{ii}] = ...
                LocalizeData.genLocalizations();
        end
        
        % Attempt to pair the sets of coordinates, treating the first set
        % as the reference.
        Coords1 = [FiducialSMD{1}.X, FiducialSMD{1}.Y];
        FinalCoordinates = cell(NFiducials, 1);
        FinalCoordinates{1} = repmat(Coords1, [1, 1, 2]);
        for ii = 2:NFiducials
            % Pair the coordinates.
            Coords2 = [FiducialSMD{ii}.X, FiducialSMD{ii}.Y];
            PairMap12 = obj.pairCoordinates(Coords1, Coords2, ...
                obj.SeparationThreshold);
            
            % Define PairedCoordinates{ii} as a 3D array, with the
            % third dimension corresponding to the fiducial number.
            IsPairedCoords1 = ~isnan(PairMap12);
            PairedCoordinates = cat(3, ...
                Coords1(IsPairedCoords1, :), ...
                Coords2(PairMap12(IsPairedCoords1), :));
            
            % If needed, let the user manually cull pairs of points (there
            % might be some non-sensical pairings/localizations).
            if obj.ManualCull
                [FinalCoordinates{ii}] = obj.performManualCull(...
                    obj.FiducialImages, PairedCoordinates);
            else
                FinalCoordinates{ii} = PairedCoordinates;
            end
            
            % Compute the transform.
            switch obj.TransformationType
                case 'lwm'
                    RegistrationTransform{ii} = fitgeotrans(...
                        FinalCoordinates{ii}(:, :, 2), ...
                        FinalCoordinates{ii}(:, :, 1), ...
                        obj.TransformationType, obj.NNeighborPoints);
                case 'polynomial'
                    RegistrationTransform{ii} = fitgeotrans(...
                        FinalCoordinates{ii}(:, :, 2), ...
                        FinalCoordinates{ii}(:, :, 1), ...
                        obj.TransformationType, obj.PolynomialDegree);
                otherwise
                    RegistrationTransform{ii} = fitgeotrans(...
                        FinalCoordinates{ii}(:, :, 2), ...
                        FinalCoordinates{ii}(:, :, 1), ...
                        obj.TransformationType);
            end
            obj.Coordinates = FinalCoordinates;
        end
    case 'images'
        % Find the transform directly from the images.
        [Optimizer, Metric] = imregconfig('multimodal');
        for ii = 2:NFiducials
            RegistrationTransform{ii} = imregtform(...
                ScaledData(:, :, ii), ScaledData(:, :, 1), ...
                obj.TransformationType, Optimizer, Metric);
        end
    otherwise
        fprintf(['Unrecognized transformation basis:\n', ...
            'obj.TransformationBasis = %s'], obj.TransformationBasis)
end
obj.RegistrationTransform = RegistrationTransform;


end