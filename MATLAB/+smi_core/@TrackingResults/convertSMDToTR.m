function [TR] = convertSMDToTR(SMD, FileInfoStruct)
%convertSMDToTR converts an SMD into a TR structure.
% This method takes an SMD structure (see smi_core.SingleMoleculeData) and
% converts it into a TR structure.
% 
% INPUTS:
%   SMD: Single Molecule Data structure with a properly populated field
%        'ConnectID' (which at present, is converted to the field
%        'TrajectoryID' in the TR structure for the purpose of back
%        compatability/interpretability).
%   FileInfoStruct: A structure array with fields FileDir and FileName
% 
% OUTPUTS:
%   TR: Tracking Results.  The TR structure is a structure array containing
%       several fields present in the SMD structure but reorganized such
%       that each element of TR, e.g., TR(7), corresponds to a single
%       unique trajectory.  A TR structure should containg AT LEAST the
%       following fields (even if they are empty):
%           TrajectoryID: An integer identifier unique to each trajectory.
%           X: X position of the trajectory over time. (pixels)
%           Y: Y position of the trajectory over time. (pixels)
%           X_SE: The standard error of the X position estimate. (pixels)
%           Y_SE: The standard error of the Y position estimate. (pixels)
%           FrameNum: The frame in which the trajectory appears. (frames)
%           Photons: The integrated photons found in each localization of
%                    the trajectory. (photons)
%           Bg: The background found for each localization of the
%               trajectory. (photons)
%           IndSMD: An index back into the associated SMD structure, i.e., 
%                  TR(nn).X(ii) can be found in SMD.X(TR(nn).IndSMD(ii)).
%           FrameRate: The frame rate used to acquire this data. 
%                      (frames / second)
%           PixelSize: The pixel size of the sensor used to acquire this
%                      data BACK PROJECTED to the sample plane.
%                      (micrometers)
%           DataROI: ROI of original raw data in which the trajectories
%                    were found, formatted as [YStart, XStart, YEnd, XEnd].
%                    (pixels) 
%           XSize: The number of pixels along the X dimension of the ROI.
%           YSize: The number of pixels along the Y dimensino of the ROI.
%           FileDir: Directory containing the raw data corresponding to 
%                    this TR.
%           FileName: Name of the file corresponding to this TR.
% 
% REQUIRES:
%
% CITATION:

% Created by:
%   Hanieh Mazloom-Farsibaf (Lidke Lab, 2018)
%   Rewritten in smite, David J. Schodt (Lidke Lab, 2021)


% Set defaults as needed.
if (~exist('FileInfoStruct', 'var') || isempty(FileInfoStruct))
    FileInfoStruct.FileDir = '';
    FileInfoStruct.FileName = '';
end

% Create an empty TR structure, ending after this piece of code if no SMD
% was input (it might sometimes be useful to produce an empty TR
% structure).
TR = smi_core.TrackingResults.createTR();
[TR.FileDir] = deal(FileInfoStruct.FileDir);
[TR.FileName] = deal(FileInfoStruct.FileName);
if (~exist('SMD', 'var') || isempty(SMD) || isempty(fieldnames(SMD)))
    return
end

% Loop through each trajectory in SMD and place it in the output TR.
% NOTE: Some fields will be the same for all trajectories and are added
%       outside of the loop.
UniqueTrajIDs = unique(SMD.ConnectID);
for ii = numel(UniqueTrajIDs):-1:1
    % Create an index array corresponding only to the current trajectory.
    CurrentTrajIndices = find(SMD.ConnectID == UniqueTrajIDs(ii));

    % Sort the frame number array and then sort the CurrentTrajIndices in
    % the same manner.
    % NOTE: This isn't strictly necessary, but it's nice to have these
    %       fields ordered in time.
    [~, SortIndices] = sort(SMD.FrameNum(CurrentTrajIndices));
    CurrentTrajIndices = CurrentTrajIndices(SortIndices);
    
    % Populate the TR structure.
    TR(ii).TrajectoryID = SMD.ConnectID(CurrentTrajIndices(1));
    TR(ii).X = SMD.X(CurrentTrajIndices);
    TR(ii).Y = SMD.Y(CurrentTrajIndices);
    TR(ii).X_SE = SMD.X_SE(CurrentTrajIndices);
    TR(ii).Y_SE = SMD.Y_SE(CurrentTrajIndices);
    if (all(isfield(SMD, {'Z', 'Z_SE'})) ...
            && ~isempty(SMD.Z) && ~isempty(SMD.Z_SE))
        TR(ii).Z = SMD.Z(CurrentTrajIndices);
        TR(ii).Z_SE = SMD.Z_SE(CurrentTrajIndices);
    end
    TR(ii).FrameNum = SMD.FrameNum(CurrentTrajIndices);
    TR(ii).Photons = SMD.Photons(CurrentTrajIndices);
    TR(ii).Bg = SMD.Bg(CurrentTrajIndices);
    TR(ii).IndSMD = CurrentTrajIndices;
end
[TR.FrameRate] = deal(SMD.FrameRate);
[TR.PixelSize] = deal(SMD.PixelSize);
[TR.XSize] = deal(SMD.XSize);
[TR.YSize] = deal(SMD.YSize);
[TR.FileDir] = deal(FileInfoStruct.FileDir);
[TR.FileName] = deal(FileInfoStruct.FileName);


end