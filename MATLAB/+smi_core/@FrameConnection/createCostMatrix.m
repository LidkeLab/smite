function [CostMatrix] = createCostMatrix(ClusterData, ...
    KOn, KOff, KBleach, PMiss, InitialDensity, MaxFrameGap, EndFrame, ...
    NonLinkMarker)
%createCostMatrix creates the cost matrix for frame connection
% This method creates a cost matrix whose elements represent the cost for
% connecting super-resolved localizations to one another as an alternative
% approach to frame connection.
%
% INPUTS:
%   ClusterData: Array of cluster data for a specific cluster, extracted
%                from the ClusterData cell array output from
%                organizeClusterData(). (NLocalizations x 6 numeric array)
%                [X, Y, X_SE, Y_SE, FrameNum, ConnectID]
%   KOn: Transition rate for emitters turning on from the dark state.
%        (1 / frame)
%   KOff: Transition rate for emitters reverting to the dark state.
%         (1 / frame)
%   KBleach: Transition rate for emitters bleaching.
%            (1 / frame)
%   PMiss: Probability of missing a localization (i.e., the probability
%          that an emitter is on but was missed during the localization
%          step).
%   InitialDensity: Initial emitter density associated with the current
%                   cluster. (emitters / pixel^2 at start of experiment)
%   MaxFrameGap: Maximum frame gap considered during the pre-clustering
%                step.
%   EndFrame: Last possible frame for the collected data (e.g., if the
%             experimental data consists of 2000 frames, EndFrame = 2000).
%   NonLinkMarker: A marker in the output CostMatrix that indicates we
%                  don't want to select that element in the linear
%                  assignment problem.
%                  (scalar, ~=inf, ~=nan, and typically < 0)
%
% OUTPUTS:
%   CostMatrix: The cost matrix whose elements represent the cost for
%               connecting one localization to another localization/cluster
%               of localizations.
%               (2*NLoc x 2*NLoc numeric array, where NLoc is the number
%               of localizations in the input SMD structure)
%
% CITATION:

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Extract some arrays from the ClusterData to make the code more readable.
X = ClusterData(:, 1);
Y = ClusterData(:, 2);
X_SE = ClusterData(:, 3);
Y_SE = ClusterData(:, 4);
FrameNum = ClusterData(:, 5);

% Define various parameters.
NLocalizations = numel(FrameNum);
CMSize = 2 * NLocalizations;

% Populate the upper-left connection block.
CostMatrix = NonLinkMarker * ones(CMSize);
for mm = 1:NLocalizations
    % For each localization, define the cost of linking to the other
    % localizations (divided by two, since a selection of one of these
    % costs leads to the same selection in the auxillary block).  In this
    % case, it's the negative log-likelihood, where the likelihood is given
    % as p(observed separation | localization error) ...
    %  * [p(blinking off w/o bleaching, staying off, then blinking back on)
    %   + p(missing last N frames localizations)p(not missing new loc.)]
    for nn = (mm+1):NLocalizations
        DeltaFrame = abs(FrameNum(mm) - FrameNum(nn));
        if (DeltaFrame == 0)
            continue
        end
        SigmaX = sqrt(X_SE(mm)^2 + X_SE(nn)^2);
        SigmaY = sqrt(Y_SE(mm)^2 + Y_SE(nn)^2);
        SeparationCost = log(2*pi*SigmaX*SigmaY) ...
            + (X(mm)-X(nn))^2 / (2*SigmaX^2) ...
            + (Y(mm)-Y(nn))^2 / (2*SigmaY^2);
        ObservationCost = -log((PMiss^(DeltaFrame-1)) * (1-PMiss));
        StillOnCost = (KOff + KBleach) * DeltaFrame;
        CostMatrix(mm, nn) = ...
            (SeparationCost + ObservationCost + StillOnCost) / 2;
        CostMatrix(nn+NLocalizations, mm+NLocalizations) = ...
            CostMatrix(mm, nn);
    end
end

% Populate the lower-left birth block and the upper-right death block.
StartFrame = 1;
DutyCycle = KOn / (KOn+KOff+KBleach);
FrameArray = (1:max(FrameNum)).';
RhoOn = InitialDensity * DutyCycle ...
    * exp(-KBleach*DutyCycle*(FrameArray-1));
RhoOff = RhoOn * KOff / KOn;
IndexArray = 1:NLocalizations;
for nn = 1:NLocalizations
    FrameGapPast = min(MaxFrameGap, FrameNum(nn) - StartFrame);
    FrameGapFuture = min(MaxFrameGap, EndFrame - FrameNum(nn));
    BirthCost = -log(1-PMiss) ...
        -log(RhoOff(FrameNum(nn))*(1-exp(-KOn))*exp(-FrameGapPast*KOn) ...
        + RhoOn(FrameNum(nn)-FrameGapPast)*(PMiss^FrameGapPast));
    DeathCost = -log((1-exp(-KOff)) ...
        + (1-exp(-KBleach)) ...
        + (PMiss^FrameGapFuture));
    CostMatrix(IndexArray+NLocalizations, nn) = BirthCost;
    CostMatrix(nn, IndexArray+NLocalizations) = DeathCost;
end

% Set infinite costs to the NonLinkMarker.
CostMatrix(isinf(CostMatrix) | isnan(CostMatrix)) = NonLinkMarker;


end