function [CostMatrix] = createCostMatrixFF(SMD, SMF, ...
    DiffusionConstants, FrameNumber, NonLinkMarker)
%createCostMatrixFF generates frame-to-frame connection cost matrix.
% This method creates the cost matrix for the frame-to-frame connection of
% localizations present in SMD.
%
% INPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%        containing the localizations that we wish to stitch into
%        trajectories.
%   SMF: Single Molecule Fitting structure defining many of the parameters
%        we'll just to populate cost matrix elements.
%       (see smi_core.SingleMoleculeFitting)
%   DiffusionConstants: Diffusion constants for each localization in SMD
%                       (column 1) and their SEs (column 2). If this is not
%                       provided, we'll use SMF.Tracking.D for all
%                       trajectories.
%                       (numel(SMD.FrameNum)x2 array)(px^2/frame)
%   FrameNumber: The frame containing the localizations for which we want
%                to construct a cost matrix (for connection to
%                FrameNumber+1).
%   NonLinkMarker: A marker in the output CostMatrix that indicates we
%                  don't want to select that element in the linear
%                  assignment problem.
%                  (scalar, ~=inf, ~=nan, and typically < 0)(Default = -1)
%
% OUTPUTS:
%   CostMatrix: The cost matrix whose elements represent the cost for
%               connecting a particle (whether it be a real localization
%               present in SMD, or a "ghost" particle that we don't see) in
%               frame FrameNumber to a particle in frame FrameNumber+1.
%               (m+n x m+n numeric array, where n is the number of
%               localizations in FrameNumber, m is the number of
%               localizations in FrameNumber+1)
%
% CITATION:

% Created by:
%   Will Kanagy (Lidke Lab, 2018)
%   Hanieh Mazloom-Farsibaf (Lidke Lab, 2019)
%   Revised and added comments, David J. Schodt (Lidke Lab, 2020)
%   Rewritten for smite + several new parameters, David J. Schodt 
%       (Lidke Lab, 2021)


% Set defaults as needed.
if (~exist('NonLinkMarker', 'var') || isempty(NonLinkMarker))
    NonLinkMarker = -1;
end
if (~exist('DiffusionConstants', 'var') || isempty(DiffusionConstants))
    DiffusionConstants = [SMF.Tracking.D, inf] ...
        .* ones(numel(SMD.FrameNum), 1);
end

% Determine which emitters in SMD were present in frames FrameNumber and
% FrameNumber+1.
EmitterIndicesFrame1 = find(SMD.FrameNum == FrameNumber);
EmitterIndicesFrame2 = find(SMD.FrameNum == (FrameNumber+1));
NFrame1 = numel(EmitterIndicesFrame1);
NFrame2 = numel(EmitterIndicesFrame2);

% Generate the upper left (the upper NxM) block of the cost matrix.
CMSize = NFrame2 + NFrame1;
CostMatrix = NonLinkMarker * ones(CMSize);
for ii = 1:NFrame1
    % Isolate the localizations in frame FrameNumber from SMD.
    CurrentIndicesFrame1 = EmitterIndicesFrame1(ii);
    XFrame1 = SMD.X(CurrentIndicesFrame1);
    X_SEFrame1 = SMD.X_SE(CurrentIndicesFrame1);
    YFrame1 = SMD.Y(CurrentIndicesFrame1);
    Y_SEFrame1 = SMD.Y_SE(CurrentIndicesFrame1);
    DFrame1 = DiffusionConstants(CurrentIndicesFrame1, :);
    for jj = 1:NFrame2
        % Isolate the localizations in frame FrameNumber+1 from SMD.
        CurrentIndicesFrame2 = EmitterIndicesFrame2(jj);
        XFrame2 = SMD.X(CurrentIndicesFrame2);
        X_SEFrame2 = SMD.X_SE(CurrentIndicesFrame2);
        YFrame2 = SMD.Y(CurrentIndicesFrame2);
        Y_SEFrame2 = SMD.Y_SE(CurrentIndicesFrame2);
        DFrame2 = DiffusionConstants(CurrentIndicesFrame2, :);
        
        % Define the standard deviations of the distribution of observed X
        % and Y separations between two localizations in consecutive
        % frames, assuming they came from the same emitter which
        % experienced Brownian motion with diffusion constant D.
        % NOTE: When using trajectory-wise D's, I'm not sure which approach
        %       is more valid: D==DFrame1 in the following equations, or
        %       D==mean(DFrame1+DFrame2)?  I think D==DFrame1 makes more
        %       sense, however that might be too tight of a restriction
        %       (e.g., if one of the localizations was previously estimated
        %       to be from a very slow trajectory, that might prevent us
        %       from connecting the two localizations).
        Sigma_X = sqrt(DFrame1(1) + DFrame2(1) ...
            + X_SEFrame1^2 + X_SEFrame2^2);
        Sigma_Y = sqrt(DFrame1(1) + DFrame2(1) ...
            + Y_SEFrame1^2 + Y_SEFrame2^2);
        
        % Compute the distance between the localizations and determine if
        % we need to proceed.
        XJump = XFrame1 - XFrame2;
        YJump = YFrame1 - YFrame2;
        Separation = sqrt(XJump^2 + YJump^2);
        ZScoreD = abs(DFrame1(1)-DFrame2(1)) ...
            / sqrt(DFrame1(2)^2 + DFrame2(2)^2);
        if ((Separation>SMF.Tracking.MaxDistFF) ...
                || (abs(XJump/Sigma_X)>SMF.Tracking.MaxZScoreDist) ...
                || (abs(YJump/Sigma_Y)>SMF.Tracking.MaxZScoreDist) ...
                || (ZScoreD>SMF.Tracking.MaxZScoreD))
            continue
        end
        
        % Define the negative log-likelihood of the observed X, Y from
        % FrameNumber and FrameNumber+1 having come from the Normal
        % distributions defined by Sigma_X and Sigma_Y (i.e., this is
        % -log(Normal(mean = 0, Sigma_X) * Normal(mean = 0, Sigma_Y))
        NegLikelihoodOfXandY = log(2*pi*Sigma_X*Sigma_Y) ...
            + (XFrame1-XFrame2)^2 / (2*Sigma_X^2) ...
            + (YFrame1-YFrame2)^2 / (2*Sigma_Y^2);
        
        % Compute the negative log-likelihood of these two localizations
        % being the same emitter (i.e., this is
        % -log(NormalX * NormalY * P(not turning off))
        % NOTE: There is an implicit DeltaT = 1 frame, i.e., cost of not
        %       turning off is K_off*DeltaT with DeltaT = 1 frame
        % NOTE: The additional 0.5 multiplying the link costs comes from
        %       our choice to set the "auxillary" (bottom right) block of
        %       the cost matrix equal to the transpose of the "link"
        %       block.  With that choice, we have argued (but not proven!)
        %       that any element selected in the link block during the
        %       linear assignment problem will be selected in the
        %       "auxillary" block as well (and vice versa), and thus the
        %       factor of 0.5 ensures that the total cost of elements
        %       selected in the LAP will be consistent with the physically
        %       reasonable costs we've defined otherwise.
        CostMatrix(ii, jj) = ...
            0.5 * (NegLikelihoodOfXandY+SMF.Tracking.K_off);
    end
end

% Define the costs of birth and death of particles between frames
% FrameNumber and FrameNumber+1 (the costs of introducing a new emitter
% appearing in FrameNumber+1/an emitter disappearing in FrameNumber+1,
% respectively).
if isscalar(SMF.Tracking.Rho_off)
    CostBirth = -log(SMF.Tracking.Rho_off * (1-exp(-SMF.Tracking.K_on)));
else
    % Rescale the data to match the size of our density image.
    Scale = size(SMF.Tracking.Rho_off) ./ [SMD.YSize, SMD.XSize];
    Y = min(size(SMF.Tracking.Rho_off, 1), ...
        max(1, ceil((SMD.Y(EmitterIndicesFrame2)-0.5)*Scale(1) + 0.5)));
    X = min(size(SMF.Tracking.Rho_off, 2), ...
        max(1, ceil((SMD.X(EmitterIndicesFrame2)-0.5)*Scale(2) + 0.5)));
    
    % Compute the birth costs
    CostBirth = inf(numel(Y));
    for ii = 1:numel(Y)
        CostBirth(:, ii) = -log(SMF.Tracking.Rho_off(Y(ii), X(ii)) ...
            * (1-exp(-SMF.Tracking.K_on)));
    end
end
CostDeath = -log(1-exp(-SMF.Tracking.K_off));

% Populate the remaining blocks of the cost matrix as appropriate.
% NOTE: The bottom right "auxillary" block is not physically meaningful (at
%       least not in an obvious way?).  Our choice to set it equal to the
%       transpose of the upper left block seems appropriate but we never
%       "proved" that this is always the best choice.
CostMatrix((NFrame1+1):CMSize, 1:NFrame2) = CostBirth;
CostMatrix(1:NFrame1, (NFrame2+1):CMSize) = CostDeath;
CostMatrix((NFrame1+1):CMSize, (NFrame2+1):CMSize) = ...
    CostMatrix(1:NFrame1, 1:NFrame2).';


end