function [Density] = computeDensity(SMD)
%computeDensity estimates the per frame density of observed emitters.
% This method computes the number of observed localizations per pixel^2 
% present in each frame of data in the input SMD structure.  This is done
% by computing the area of the ROI in pixel^2, counting the number of
% localizations present in SMD for each frame, and then dividing the number
% of localizations in each frame by the area of the ROI.
%
% INPUTS:
%   SMD: Single Molecule Data structure.
%
% OUTPUTS:
%   Density: The density of observed emitters in each frame. 
%            (localizations / pixel^2)(NFrames x NDatasets numeric array)
%
% CITATION:

% Created by Derek Rinaldi (Diane Lidke Lab, 2018)
% Rewritten by David J. Schodt (Lidke Lab, 2020)


% Compute the area of the ROI in pixel^2.
ROIArea = SMD.XSize * SMD.YSize;

% Isolate some fields from the SMD/compute some new arrays as needed.
FrameNum = SMD.FrameNum;
UniqueFrames = sort(unique(FrameNum));
DatasetNum = SMD.DatasetNum;

% Loop through each data dataset present in the SMD structure and compute
% the per frame density of observed emitters.
Density = zeros(SMD.NFrames, SMD.NDatasets);
for dd = 1:SMD.NDatasets
    % Loop through each frame of this dataset.
    for ff = UniqueFrames.'
        % Identify all localizations in this frame and compute the density.
        Density(ff, dd) = sum((FrameNum==ff) & (DatasetNum==dd)) / ROIArea;
    end
end


end