function [DensityIm] = computeDensityImage(SMD, SigmaSmooth, Mag, ...
    NEmitters, TimeSeries, FastGauss)
%computeDensityImage estimates the local density of observed emitters.
% This method prepares a smoothed Gaussian image of the localizations in
% SMD, normalized such that it represents an approximate density of
% emitters.  Specifically, for a single image in DensityIm, each pixel is
% in units of emitters/(pixel^2), so the number of emitters in, e.g.,
% DensityIm(:, :, n) is sum_{i,j} (1px^2) * DensityIm(i, j, n) 
%   = sum(sum(DensityIm(:, :, n))
%
% INPUTS:
%   SMD: Single Molecule Data structure.
%   SigmaSmooth: Standard deviation of the smoothing filter applied to the
%                Gaussian image to generate the density image, given in
%                units of SMD coordinates.
%                (pixels)(Default = 0)
%   Mag: Approximate magnification from SMD pixels to output image pixels.
%        (Default = 1)
%   NEmitters: Total number of emitters from which SMD was generated.
%              (NFramesx1 or 1x1)(Default based on setting of TimeSeries)
%   TimeSeries: Flag to indicate SMD.NFrames density images are desired
%               (i.e., one per frame) as opposed to single density image.
%               (Default = false)
%   FastGauss: Flag indicating an approximate Gaussian filter should be
%              used to improve speed for large stacks of data at the
%              expense of potential artifacts. (Default = false)
%
% OUTPUTS:
%   DensityIm: Gaussian image representing the approximate emitter
%              densities within each pixel. ...
%              (emitters / pixel^2)(SMD.YSize x SMD.XSize float).
%
% REQUIRES:
%   Image Processing Toolbox 2015a or sooner if FastGauss = false
%
% CITATION:

% Created by David J. Schodt (Lidke Lab, 2021)


% Set defaults.
if (~exist('SigmaSmooth', 'var') || isempty(SigmaSmooth))
    SigmaSmooth = 0;
end
if (~exist('Mag', 'var') || isempty(Mag))
    Mag = 20;
end
if (~exist('TimeSeries', 'var') || isempty(TimeSeries))
    TimeSeries = false;
end
if (~exist('NEmitters', 'var') || isempty(NEmitters))
    % For time series, we'll use the number of localizations per frame as
    % our number of emitters (e.g., the density in each frame is a
    % localization density, which for something like SPT ~represents the
    % emitter density).
    NEmitters = groupcounts(SMD.FrameNum, 1:(SMD.NFrames+1), ...
        'IncludeEmptyGroups', true, 'IncludedEdge', 'left');
    if ~TimeSeries
        % In this case, we'll take the number of emitters to be the maximum
        % number of localizations seen in any single frame.
        NEmitters = max(NEmitters);
    end
end
if (~exist('FastGauss', 'var') || isempty(FastGauss))
    FastGauss = false;
end

% Rescale the data in SMD based on 'Mag'.
% NOTE: gaussBlobImage() below defines the Gaussians in terms of
%       SMD.PSFSigma.  For the density image, we want to use the standard
%       errors for the Gaussian st. devs., thus we need to overwrite
%       PSFSigma.
ImSize = ceil(Mag * [SMD.YSize, SMD.XSize]);
SMD.XSize = ImSize(2);
SMD.YSize = ImSize(1);
SMD.X = Mag*(SMD.X - 0.5) + 0.5;
SMD.Y = Mag*(SMD.Y - 0.5) + 0.5;
SMD.PSFSigma = Mag * [SMD.Y_SE, SMD.X_SE];

% Generate the Gaussian images of the localizations in SMD.
if ~TimeSeries
    % If we set all of SMD.FrameNum to 1, it'll put every localization in
    % the same frame.
    SMD.FrameNum = ones(size(SMD.FrameNum));
    SMD.NFrames = 1;
end
GaussIm = smi_sim.GaussBlobs.gaussBlobImage(SMD);

% Apply a Gaussian filter to the GaussIm of localizations.
% NOTE: gaussInPlace() is an approximation to a Gaussian filter.  As a
%       result, some artifacts appear.
Sigma = SigmaSmooth * Mag;
if (Sigma <= 0)
    DensityIm = GaussIm;
elseif FastGauss
    DensityIm = gather(smi_core.FindROI.gaussInPlace(GaussIm, Sigma));
else
    DensityIm = imgaussfilt(GaussIm, Sigma);
end

% Rescale the smoothed Gaussian so that each frame represents an
% approximate emitter density.
if (numel(NEmitters) == size(DensityIm, 3))
    % In this case, we have NEmitters defined for every image of our
    % smoothed image, so we should normalize frame-to-frame.
    for ff = 1:numel(NEmitters)
        if any(any(DensityIm(:, :, ff)))
            DensityIm(:, :, ff) = DensityIm(:, :, ff) ...
                - min(min(DensityIm(:, :, ff)));
            DensityIm(:, :, ff) = DensityIm(:, :, ff) ...
                / sum(sum(DensityIm(:, :, ff)));
            DensityIm(:, :, ff) = ...
                NEmitters(ff) * DensityIm(:, :, ff);
        end
    end
else
    % If this happens, it's not clear what the user wants, so we'll issue a
    % warning and then apply the average of NEmitters to normalize the
    % entire stack at once.
    warning('Provided ''NEmitters'' does not match the number of frames!')
    if any(DensityIm(:))
        NEmitters = mean(NEmitters);
        DensityIm = DensityIm - min(DensityIm(:));
        DensityIm = DensityIm / sum(DensityIm(:));
        DensityIm = NEmitters * DensityIm;
    end
end


end