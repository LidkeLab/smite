function [CostMatrix] = createCostMatrixGC(SMD, SMF, ...
    NonLinkMarker, CreateSparseMatrix)
%createCostMatrixGC creates the cost matrix for gap closing in SPT.
% This method creates a cost matrix whose elements represent the cost for
% connecting trajectory segments produced by the frame-to-frame step in the
% typical smi.SPT workflow.
%
% INPUTS:
%   SMD: Single Molecule Data structure (see smi_core.SingleMoleculeData)
%        containing the localizations that we wish to stitch into
%        trajectories.
%   SMF: Single Molecule Fitting structure defining many of the parameters
%        we'll just to populate cost matrix elements.
%       (see smi_core.SingleMoleculeFitting)
%   NonLinkMarker: A marker in the output CostMatrix that indicates we
%                  don't want to select that element in the linear
%                  assignment problem.
%                  (scalar, ~=inf, ~=nan, and typically < 0)(Default = -1)
%   CreateSparseMatrix: This boolean flag will determine whether or not we
%                       should define the output CM as a sparse matrix
%                       (and represent it in MATLAB as a sparse type) or
%                       as a "regular" matrix with a NonLinkMarker for
%                       uninteresting costs.
%                       (boolean flag, 0 or 1)(Default = 1)
%
% OUTPUTS:
%   CostMatrix: The cost matrix whose elements represent the cost for
%               connecting the end of one trajectory to the start of
%               another trajectory in a future frame.
%               (2*NTraj x 2*NTraj numeric, array, where NTraj is the
%               number of trajectories in the input SMD structure)
%
% CITATION:

% Created by
%   Elton Jhamba (Lidke Lab, 2018)
%   Hanieh Mazloom-Farsibaf, Michael Wester (Lidke Lab, 2019)
%   Rewritten for sparse matrices, David J. Schodt (Lidke Lab, Spring 2020)
%   Rewritten for smite, David J. Schodt (Lidke Lab, 2021)


% Set defaults as needed.
if (~exist('NonLinkMarker', 'var') || isempty(NonLinkMarker))
    NonLinkMarker = -1;
end
if (~exist('CreateSparseMatrix', 'var') || isempty(CreateSparseMatrix))
    CreateSparseMatrix = true;
end

% Extract some arrays from the SMD/SMF structure.  Doing this outside of
% for loops can speed things up for very large SMD structures (although for
% small SMD structures this has the potential to slow things down, but then
% we don't care about speed).
X = SMD.X;
Y = SMD.Y;
X_SE = SMD.X_SE;
Y_SE = SMD.Y_SE;
FrameNum = SMD.FrameNum;
ConnectID = SMD.ConnectID;
MaxFrameGap = SMF.Tracking.MaxFrameGap;
MaxDistGC = SMF.Tracking.MaxDistGC;
D = SMF.Tracking.D;
K_on = SMF.Tracking.K_on;
K_off = SMF.Tracking.K_off;
Rho_off = SMF.Tracking.Rho_off;

% Define various parameters that will be useful later on.
NTraj = max(ConnectID);
CMSize = 2 * NTraj * [1, 1];

% Determine the starting and ending indices within SMD of each trajectory
% (HMF noted that preparing these outside of the nested for loop below will
% speed up the code).
StartEndIndices = zeros(NTraj, 2);
for ii = 1:NTraj
    CurrentTrajBool = (ConnectID == ii);
    CurrentFrames = FrameNum(CurrentTrajBool);
    StartEndIndices(ii, :) = ...
        [find((FrameNum==min(CurrentFrames)) & CurrentTrajBool), ...
        find((FrameNum==max(CurrentFrames)) & CurrentTrajBool)];
end

% Initialize the cost matrix index array and element array (this is the
% best way to create a sparse matrix, but we can also use this for the
% non-sparse, "regular" matrix type).
CMIndices = [repelem((1:CMSize(1)).', CMSize(1)), ...
    repmat((1:CMSize(1)).', [CMSize(1), 1])];
CMElements = NonLinkMarker * ones(CMSize(1)^2, 1);

% Fill in the NTraj x NTraj upper-left block (the link block) and the
% auxillary (bottom-left) block.
for ee = 1:NTraj
    % Isolate the end of the ee-th trajectory.
    EndIndexCurrent = StartEndIndices(ee, 2);
    EndFrameCurrent = FrameNum(EndIndexCurrent);
    XEndCurrent = X(EndIndexCurrent);
    X_SEEndCurrent = X_SE(EndIndexCurrent);
    YEndCurrent = Y(EndIndexCurrent);
    Y_SEEndCurrent = Y_SE(EndIndexCurrent);
    
    % Loop through the trajectory beginnings and populate the cost matrix.
    for bb = 1:NTraj
        if (ee == bb)
            % No need to run the rest of the code, we don't care about this
            % cost.
            continue
        end
        
        % Isolate the index corresponding to the beginning of the bb-th
        % trajectory.
        % NOTE: We won't be granted any time savings by defining separate
        %       like XBeginCurrent here (as we had for XEndCurrent above)
        %       since each of those would only get used once.
        StartIndexCurrent = StartEndIndices(bb, 1);
        
        % Determine how many frames have elapsed between the end of
        % trajectory ee and the start of trajectory bb.
        DeltaFrame = FrameNum(StartIndexCurrent) - EndFrameCurrent;
        
        % Compute the distance between the end of trajectory ee and the
        % start of trajectory bb.
        DeltaXY = sqrt((XEndCurrent-X(StartIndexCurrent))^2 ...
            + (YEndCurrent-Y(StartIndexCurrent))^2);
        
        % Compute the cost corresponding to linking these two trajectories
        % (unless the MaxDist and MaxFrameGap aren't satisified, in which
        % case we don't care to calculate this cost).
        if ((DeltaFrame>0) ...
                && (DeltaFrame<=MaxFrameGap) ...
                && (DeltaXY<=MaxDistGC))
            % Define the standard deviations of the distribution of
            % observed X and Y separations between two localizations in
            % separated in time by DeltaFrame frames, assuming they came
            % from the same emitter which experienced Brownian motion with
            % diffusion constant D.
            Sigma_X = sqrt(2*D*DeltaFrame ...
                + X_SEEndCurrent^2 + X_SE(StartIndexCurrent)^2);
            Sigma_Y = sqrt(2*D*DeltaFrame ...
                + Y_SEEndCurrent^2 + Y_SE(StartIndexCurrent)^2);
            
            % Define the log-likelihood of the observed X, Y from
            % FrameNumber and FrameNumber+1 having come from the Normal
            % distributions defined by Sigma_X and Sigma_Y (i.e., this is
            % -log(Normal(mean = 0, Sigma_X) * Normal(mean = 0, Sigma_Y))
            NegLogLikelihoodOfXandY = log(2*pi*Sigma_X*Sigma_Y) ...
                + (XEndCurrent-X(StartIndexCurrent))^2 / (2*Sigma_X^2) ...
                + (YEndCurrent-Y(StartIndexCurrent))^2 / (2*Sigma_Y^2);
            
            % Compute the negative log-likelihood of these two trajectories
            % arising from the same emitter. To clarify, this is
            % -log(NormalX * NormalY
            %   * P(traj ee turning off)
            %   * P(traj ee not turning on for DeltaFrame-1 frames)
            %   * P(traj bb turning on)
            % NOTE: There are implicit DeltaT's multiplying all transition
            %       rates, where DeltaT = 1 frame.
            % NOTE: The additional 0.5 multiplying the gap-closing costs
            %       comes from our choice to set the "auxillary" (bottom
            %       right) block of the cost matrix equal to the transpose
            %       of the "gap-closing" block.  With that choice, we have
            %       argued (but not proven!) that any element selected in
            %       the gap-closing block during the linear assignment
            %       problem will be selected in the "auxillary" block as
            %       well (and vice versa), and thus the factor of 0.5
            %       ensures that the total cost of elements selected in the
            %       LAP will be consistent with the physically reasonable
            %       costs we've defined otherwise.
            % NOTE: The usage of sub2ind() has the rows and columns
            %       reversed from normal usage.  This is an artifact
            %       resulting from differing conventions used in MATLAB:
            %       usually, MATLAB is column oriented, but our usage of
            %       CMIndices here (and in sparse()) is row oriented.
            CMElements(sub2ind(CMSize, bb, ee)) = ...
                0.5 * (NegLogLikelihoodOfXandY ...
                - log(K_off) - (DeltaFrame-1)*log(1-K_on) - log(K_on));
            
            % Store the corresponding element of the "auxillary" block.
            CMElements(sub2ind(CMSize, ee+NTraj, bb+NTraj)) = ...
                CMElements(sub2ind(CMSize, bb, ee));
        end
    end
end

% Fill in the "birth" block (lower left) of the cost matrix.
% NOTE: We set FirstFrame = 1 instead of min(SMD.FrameNum) based on the
%       assumption that our experiment started at time corresponding to
%       FrameNum = 1, even if we don't have any localizations in that
%       frame.
FirstFrame = 1;
OnesArray = ones(1, NTraj);
for bb = 1:NTraj
    % Find the starting frame for the bb-th trajectory.
    StartIndexCurrent = StartEndIndices(bb, 1);
    StartFrameCurrent = FrameNum(StartIndexCurrent);
    
    % Compute the cost of "birth" for this trajectory (i.e., this is
    % -log(P(traj bb turning on)
    %   * P(traj bb not turning on for StartFrame-FirstFrame frames))
    % NOTE: There are implicit factors of DeltaT = 1 frame multiplying all
    %       transition rates.
    % NOTE: The usage of sub2ind() has the rows and columns reversed from
    %       normal usage.  This is an artifact resulting from differing
    %       conventions used in MATLAB: usually, MATLAB is column oriented,
    %       but our usage of CMIndices here (and sparse()) is row oriented.
    CMElements(sub2ind(CMSize, bb*OnesArray, (NTraj+1):(2*NTraj))) = ...
        -(log(Rho_off*K_on) + (StartFrameCurrent-FirstFrame)*log(1-K_on));
end

% Fill in the "death" block (upper-right) of the cost matrix.
LastFrame = SMD.NFrames;
for ee = 1:NTraj
    % Find the ending frame for the ee-th trajectory.
    EndIndexCurrent = StartEndIndices(ee, 2);
    EndFrameCurrent = FrameNum(EndIndexCurrent);
    
    % Compute the cost of "death" for this trajectory (i.e., this is
    % -log(P(traj ee turning off)
    %   * P(traj ee not turning back on for LastFrame-EndFrame frames)
    % NOTE: There are implicit factors of DeltaT = 1 frame multiplying all
    %       transition rates.
    % NOTE: The usage of sub2ind() has the rows and columns reversed from
    %       normal usage.  This is an artifact resulting from differing
    %       conventions used in MATLAB: usually, MATLAB is column oriented,
    %       but our usage of CMIndices here (and sparse()) is row oriented.
    CMElements(sub2ind(CMSize, (NTraj+1):(2*NTraj), ee*OnesArray)) = ...
        -(log(K_off) + (LastFrame-EndFrameCurrent) * log(1-K_on));
end

% Remove all of the remaining NonLinkMarker elements.
KeepBool = (CMElements ~= NonLinkMarker);
CMElements = CMElements(KeepBool);
CMIndices = CMIndices(KeepBool, :);

% Create the cost matrix.
if CreateSparseMatrix
    % Create the output CostMatrix as a sparse type.
    CostMatrix = sparse(CMIndices(:, 1), CMIndices(:, 2), ...
        double(CMElements));
else
    % Create the output CostMatrix as a full sized matrix.
    CostMatrix = NonLinkMarker * ones(CMSize);
    CostMatrix(sub2ind(CMSize, CMIndices(:, 1), CMIndices(:, 2))) = ...
        CMElements;
end


end