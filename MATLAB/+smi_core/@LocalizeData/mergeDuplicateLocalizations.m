function [SMD] = mergeDuplicateLocalizations(obj, SMD)
%mergeDuplicateLocalizations merges localizations at same position/frame.
% This method identifies and merges duplicate localizations that occur at
% the same position within the same frame. This can happen when FindROI
% detects multiple identical pixel maxima (common in noise-free simulation
% data). When duplicates are found, the localization with the higher
% LogLikelihood is kept.
%
% INPUTS:
%   SMD: Single Molecule Data structure with candidate localizations
%
% OUTPUTS:
%   SMD: Modified SMD structure with duplicates removed
%
% NOTE: This addresses Issue #21 where FindROI can generate multiple boxes
%       for the same emitter when pixels have identical maximum values.

% Created by:
%   Claude/Keith Lidke (Lidke Lab, 2025)


% Return early if no localizations to process
if isempty(SMD.X) || numel(SMD.X) < 2
    return;
end

% Distance threshold for considering duplicates (pixels)
% Two localizations are considered duplicates if they are within this
% distance in the same frame. 0.1 pixels is conservative - true duplicates
% from identical pixel maxima will be < 0.01 pixels apart.
dist_threshold = 0.1;

% Find duplicates frame-by-frame
unique_frames = unique(SMD.FrameNum);
keep_indices = true(size(SMD.X));
n_duplicates_found = 0;

for ff = 1:length(unique_frames)
    frame_mask = SMD.FrameNum == unique_frames(ff);
    frame_indices = find(frame_mask);

    if length(frame_indices) < 2
        continue;  % No duplicates possible with single localization
    end

    % Check pairwise distances within this frame
    for ii = 1:length(frame_indices)-1
        if ~keep_indices(frame_indices(ii))
            continue;  % Already marked for removal
        end

        for jj = ii+1:length(frame_indices)
            if ~keep_indices(frame_indices(jj))
                continue;  % Already marked for removal
            end

            % Calculate Euclidean distance
            dx = SMD.X(frame_indices(ii)) - SMD.X(frame_indices(jj));
            dy = SMD.Y(frame_indices(ii)) - SMD.Y(frame_indices(jj));
            dist = sqrt(dx^2 + dy^2);

            if dist < dist_threshold
                % Duplicates found - keep the one with higher LogLikelihood
                % (better fit quality)
                if SMD.LogLikelihood(frame_indices(ii)) >= ...
                        SMD.LogLikelihood(frame_indices(jj))
                    keep_indices(frame_indices(jj)) = false;
                    n_duplicates_found = n_duplicates_found + 1;
                else
                    keep_indices(frame_indices(ii)) = false;
                    n_duplicates_found = n_duplicates_found + 1;
                    break;  % Move to next ii since this one is removed
                end
            end
        end
    end
end

% Extract only the non-duplicate localizations
if n_duplicates_found > 0
    if obj.Verbose > 1
        fprintf(['\tLocalizeData.mergeDuplicateLocalizations(): ', ...
            'Removed %d duplicate localizations\n'], n_duplicates_found);
    end
    SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, keep_indices);
end


end
