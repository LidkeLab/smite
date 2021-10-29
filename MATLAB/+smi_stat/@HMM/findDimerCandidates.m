function [TRArray] = findDimerCandidates(TR1, TR2, ...
    MaxDimerSeparation, MaxSeparation, MinValidPoints, MinMedianPhotons, ...
    BorderPadding)
%findDimerCandidates finds dimer candidate pairs between TR1 and TR2.
% This method will search through all combinations of trajectories between
% TR1 and TR2 (one combination contains one trajectory from TR1, one from
% TR2) to find dimer candidates, i.e., trajectories that were observed
% within MaxDimerSeparation of each other at some point in time.
%
% INPUTS:
%   TR1: Tracking results structure (see smi_core.TrackingResults)
%   TR2: Tracking results structure (see smi_core.TrackingResults)
%   MaxDimerSeparation: Maximum value of the minimum trajectory separation
%                       for a pair to be considered a dimer candidate.
%                       (pixels)(Default = 1)
%   MaxSeparation: Maximum distance between candidate trajectories that
%                  will be saved in the output TRArray.  For example, if
%                  two trajectories were separated by 
%                  [10, 5, 3, 1, 1, 4, 12], and the MaxSeparation is 4, the
%                  output TRArray will contain the segments of these
%                  trajectories corresponding to the separations 
%                  [3, 1, 1, 4]. 
%                  (pixels)(Default = 10 * MaxDimerDistance).
%   MinValidPoints: Minimum number of points during which the candidates
%                   must be within MaxSeparation of each other.
%                   (frames)(Default = 0)
%   MinPhotons: Minimum value for the median photons a trajectory must have
%               to be considered a candidate. (photons)(Default = 0)
%   BorderPadding: Minimum number of pixels from the border the
%                  trajectories must be during their candidate duration.
%                  (pixels)(Default = 0)
%
% OUTPUTS:
%   TRArray: A structure array of TR structures, where the constituent TR
%            structures correspond to dimer candidate trajectories.  The 
%            first index corresponds to the "channel" of the trajectory and 
%            the second index corresponds to the pair number.  For example,
%            TRArray(j, 1) will contain information about a trajectory from 
%            TR1 that was observed within MaxDimerDistance of
%            TRArray(j, 2), a trajectory in TR2.
%            NOTE: The entirety of the two trajectories will be saved in
%                  the TRArray, and a new field DimerCandidateBool will be
%                  added to describe which datapoints correspond to the
%                  trajectories being close to each other (e.g., 
%                  TRArray(j, 1).FrameNum(TRArray(j, 1).DimerCandidateBool)
%                  will give the frames during which the trajectories were
%                  close).  Another field, ObservationStatus is defined as
%                  followed: ObservationStatus(1) is a boolean indicating
%                  whether or not the "off" to "on" transition
%                  (>MaxDimerSeparation to <=MaxDimerSeparation) was
%                  observed. ObservationStatus(2) is a boolean indicating
%                  whether or not the "on" to "off" transition was
%                  observed.

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Set default parameters if needed.
if (~exist('MaxDimerSeparation', 'var') || isempty(MaxDimerSeparation))
    MaxDimerSeparation = 1;
end
if (~exist('MaxSeparation', 'var') || isempty(MaxSeparation))
    MaxSeparation = 10 * MaxDimerSeparation;
end
if (~exist('MinValidPoints', 'var') || isempty(MinValidPoints))
    MinValidPoints = 0;
end
if (~exist('MinMedianPhotons', 'var') || isempty(MinMedianPhotons))
    MinMedianPhotons = 0;
end
if (~exist('BorderPadding', 'var') || isempty(BorderPadding))
    BorderPadding = 0;
end

% Make sure TR1 and TR2 share the same fields (this is a nice check which
% may allow inter-compability between old/new TR structures).
TR1 = smi_core.TrackingResults.padTR(TR1, TR2);
TR2 = smi_core.TrackingResults.padTR(TR2, TR1);

% Loop through TR1 and search for dimer candidates in TR2.
if BorderPadding
    % For certain simulations, it's easiest for us to not define these and
    % to just set BorderPadding = 0.
    XSize = TR1(1).XSize;
    YSize = TR1(1).YSize;
    UpperBorderX = XSize - BorderPadding;
    UpperBorderY = YSize - BorderPadding;
end
[TR1.DimerCandidateBool] = deal([]);
[TR1.ObservationStatus] = deal([]);
[TR1.Separation] = deal([]);
[TR1.AverageSE] = deal([]);
[TR2.DimerCandidateBool] = deal([]);
[TR2.ObservationStatus] = deal([]);
[TR2.Separation] = deal([]);
[TR2.AverageSE] = deal([]);
TRArray = struct([]);
for ii = 1:numel(TR1)
    % Isolate relevant fields from TR1 for the current channel 1 trajectory
    % and create other arrays as needed.
    XChannel1 = double(TR1(ii).X);
    YChannel1 = double(TR1(ii).Y);
    AverageSEChannel1 = mean(double([TR1(ii).X_SE, TR1(ii).Y_SE]), 2);
    FrameNumChannel1 = TR1(ii).FrameNum;
    NaNChannel1 = NaN(numel(XChannel1), 1);
    for jj = 1:numel(TR2)
        % Isolate relevant fields from TR2 for the current channel 2 
        % trajectory and create other arrays as needed.
        XChannel2 = double(TR2(jj).X);
        YChannel2 = double(TR2(jj).Y);
        AverageSEChannel2 = mean(double([TR2(jj).X_SE, TR2(jj).Y_SE]), 2);
        FrameNumChannel2 = TR2(jj).FrameNum;
        NaNChannel2 = NaN(numel(XChannel2), 1);
        
        % Determine which frames (if any) these two trajectories were
        % both observed in.
        OverlapBoolChannel1 = ismember(FrameNumChannel1, FrameNumChannel2);
        if sum(OverlapBoolChannel1) < MinValidPoints
            % If these trajectories didn't overlap in time for sufficiently
            % many frames, we can skip to the next iteration of the loop
            % now. 
            % NOTE: We'll still need to compare to MinValidPoints later
            %       on, but continuing here can speed up the code since we
            %       already know they don't overlap for very long.
            continue
        end
        OverlapBoolChannel2 = ismember(FrameNumChannel2, FrameNumChannel1);
        
        % Determine if these trajectories were ever sufficiently close to
        % be considered an interaction candidate.
        XChannel1Overlap = XChannel1(OverlapBoolChannel1);
        YChannel1Overlap = YChannel1(OverlapBoolChannel1);
        FrameNumChannel1Overlap = FrameNumChannel1(OverlapBoolChannel1);
        XChannel2Overlap = XChannel2(OverlapBoolChannel2);
        YChannel2Overlap = YChannel2(OverlapBoolChannel2);
        FrameNumChannel2Overlap = FrameNumChannel2(OverlapBoolChannel2);
        CenterToCenterSeparation = sqrt(...
            (XChannel1Overlap-XChannel2Overlap).^2 ...
            + (YChannel1Overlap-YChannel2Overlap).^2);
        WithinMaxSeparation = (CenterToCenterSeparation <= MaxSeparation);
        WithinMaxDimerDist = (CenterToCenterSeparation ...
            <= MaxDimerSeparation);
        if ~(any(WithinMaxSeparation) && any(WithinMaxDimerDist))
            % These trajectories were never sufficiently close to each
            % other so we can just skip the rest of the thresholding steps.
            continue
        end
                
        % Find the starting/ending indices of each interaction candidate 
        % for the current trajectory pair.
        % NOTE: It's easier to understand how these are found by
        %       testing them on an example.  I used
        %       WithinMaxSeparation = ...
        %           [1, 0, 1, 1, 0, 0, 1, 0, 1, 1, 0, 1].';
        StartIndices = find([WithinMaxSeparation(1); ...
            diff(WithinMaxSeparation)] == 1);
        EndIndices = find([0; diff(WithinMaxSeparation)] == -1) - 1;
        
        % If there is not an EndIndices element corresponding to every 
        % StartIndices element, that means there was an interaction that 
        % was not fully captured (WithinMaxSeparation(end) == 1).  We 
        % still want to keep that interaction so we should add the missing
        % EndIndices value.
        if (numel(StartIndices) ~= numel(EndIndices))
            EndIndices = [EndIndices; numel(WithinMaxSeparation)];
        end
                
        % Remove start/end indices that are the same, i.e. events with only
        % a single point captured within the MaxSeparation
        % NOTE: Just to be clear, it's okay to have a single point within
        %       MaxDimerSeparation, but if it's just one point within
        %       MaxSeparation then that data won't be very interpretable
        %       (if anything, this condition being met might mean that
        %       MaxSeparation is set too low!)
        SinglePointBool = (StartIndices == EndIndices);
        StartIndices = StartIndices(~SinglePointBool);
        EndIndices = EndIndices(~SinglePointBool);
        
        % Loop through each interaction that we've identified.
        NObservations = numel(FrameNumChannel1Overlap);
        for kk = 1:numel(StartIndices)
            % Define our indexing array for the current set of start/end 
            % indices.
            IndexArray = max(StartIndices(kk), 1) ...
                : min(EndIndices(kk), NObservations);
                        
            % Ensure that these trajectories were within MaxDimerDistance 
            % of each other during this portion of time.
            % NOTE: The thresholding on MaxDimerDistance above is just a
            %       preliminary threshold over the full trajectory, this 
            %       threshold is now thresholding over an isolated portion
            %       of the trajectory.
            if ~any(WithinMaxDimerDist(IndexArray))
                continue
            end
            
            % Map the start/end indices we've found back into the TR 
            % structures.
            % NOTE: This is a bit redudant for some TR fields, e.g.
            %       XChannel1Overlap(IndexArray) and
            %       XChannel1(DimerCandidateBoolChannel1) should give the
            %       same result.  I'm doing this so that we can easily 
            %       access other fields of the TR that we don't have an
            %       "overlap" version of.
            DimerCandidateBoolChannel1 = ismember(...
                FrameNumChannel1, ...
                FrameNumChannel1Overlap(IndexArray));
            DimerCandidateBoolChannel2 = ismember(...
                FrameNumChannel2, ...
                FrameNumChannel2Overlap(IndexArray));
            
            % Enforce the MinValidPoints threshold on this pair segment.
            if (sum(DimerCandidateBoolChannel1) < MinValidPoints)
                continue
            end
            
            % Enforce the border padding.
            CandidateXChannel1 = XChannel1(DimerCandidateBoolChannel1);
            CandidateYChannel1 = YChannel1(DimerCandidateBoolChannel1);
            CandidateXChannel2 = XChannel2(DimerCandidateBoolChannel2);
            CandidateYChannel2 = YChannel2(DimerCandidateBoolChannel2);
            if (BorderPadding ...
                    && (any(CandidateXChannel1 < BorderPadding) ...
                    || any(CandidateYChannel1 < BorderPadding) ...
                    || any(CandidateXChannel2 < BorderPadding) ...
                    || any(CandidateYChannel2 < BorderPadding) ...
                    || any(CandidateXChannel1 > UpperBorderX) ...
                    || any(CandidateYChannel1 > UpperBorderY) ...
                    || any(CandidateXChannel2 > UpperBorderX) ...
                    || any(CandidateYChannel2 > UpperBorderY)))
                % One of these trajectories got too close to the edge of
                % the field of view, so we'll discard these candidates.
                continue
            end
            
            % Enforce the min. photons.
            MedianPhotonsCh1 = ...
                median(TR1(ii).Photons(DimerCandidateBoolChannel1));
            MedianPhotonsCh2 = ...
                median(TR2(jj).Photons(DimerCandidateBoolChannel2));
            if ~((MedianPhotonsCh1>=MinMedianPhotons) ...
                    && (MedianPhotonsCh2>=MinMedianPhotons))
                continue
            end
            
            % Store the dimer candidate segments in the TRArray.
            TRArrayCurrent(1, 1) = TR1(ii);
            TRArrayCurrent(1).DimerCandidateBool = ...
                DimerCandidateBoolChannel1;
            SeparationChannel1 = NaNChannel1;
            SeparationChannel1(OverlapBoolChannel1) = ...
                CenterToCenterSeparation;
            CandidateSeparation = ...
                SeparationChannel1(DimerCandidateBoolChannel1);
            ObservationStatus = ...
                [(CandidateSeparation(1) > MaxDimerSeparation); 
                (CandidateSeparation(end) > MaxDimerSeparation)];
            TRArrayCurrent(1).ObservationStatus = ObservationStatus;
            TRArrayCurrent(1).Separation = SeparationChannel1;
            TRArrayCurrent(1).AverageSE = AverageSEChannel1;
            TRArrayCurrent(1, 2) = TR2(jj);
            TRArrayCurrent(2).DimerCandidateBool = ...
                DimerCandidateBoolChannel2;
            TRArrayCurrent(2).ObservationStatus = ObservationStatus;
            SeparationChannel2 = NaNChannel2;
            SeparationChannel2(OverlapBoolChannel2) = ...
                CenterToCenterSeparation;
            TRArrayCurrent(2).Separation = SeparationChannel2;
            TRArrayCurrent(2).AverageSE = AverageSEChannel2;
            TRArray = [TRArray; TRArrayCurrent];
            
            % Construct the ObservationStatus struct.
            CandidateSeparation = ...
                SeparationChannel1(DimerCandidateBoolChannel1);
            ObservationStatusCurrent.OffOnObserved = ...
                (CandidateSeparation(1) > MaxDimerSeparation);
            ObservationStatusCurrent.OnOffObserved = ...
                (CandidateSeparation(end) > MaxDimerSeparation);
        end
    end
end


end