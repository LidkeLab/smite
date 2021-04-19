function [Assign12, Cost12] = solveLAP(CostMatrix, NonlinkMarker)
%solveLAP solves the linear assignment problem specified by CostMatrix.
% This method will treat the costs in the input matrix CostMatrix as the
% costs associated with linking rows and columns together with only one
% link per row/column allowed (linear assignment problem, or LAP).
%
% NOTE: The input CostMatrix must not have any NaN or inf elements, and it
%       must allow for at least one assignment per row/column (i.e., there
%       can't be an entire row of NonlinkMarker's, or for a sparse matrix,
%       there can't be an entire row of 0's).
%
% NOTE: This method is just a wrapper around the C++ code c_lap written by
%       R. Jonker and A. Volgenant, University of Amsterdam (see
%       smite/MATLAB/source/c/c_lap.cpp).
%
% INPUTS:
%   CostMatrix: This is a square matrix whose elements CostMatrix(ii, jj)
%               correspond to the cost associated with linking item ii with
%               jj (e.g., in frame connection, the cost of linking
%               localization ii in frame n with localization jj in frame
%               n+1). (NxN numeric array with no NaN or inf elements)
%   NonlinkMarker: A marker that may or may not be present in CostMatrix
%                  which indicate that we strictly will not link the items
%                  associated with that element of CostMatrix.
%                  (1x1 scalar that is not NaN or inf)(Default = -1)
%
% OUTPUTS:
%   Assign12: An array specifying which assignments were made by solving
%             the LAP.  For example, if Link12(nn) = mm, CostMatrix(mm, nn)
%             would be the cost associated with the optimal assignment.  In
%             the context of frame connection, this could mean, e.g., that
%             localization mm in frame 1 would link to localization nn in
%             frame 2. (Nx1 int32)
%   Cost12: An array providing the costs associated with the assignments in
%           Assign12 (these are just the elements of CostMatrix associated
%           with the assignments in Assign12). (Nx1 double)
%
% REQUIRES:
%   compiled c_lap.mex64

% Created by:
%   David J. Schodt, rewrote old code of unknown author, only minor
%       modifications were made (Lidke Lab, 2020)


% Set defaults if needed.
if (~exist('NonlinkMarker', 'var') || isempty(NonlinkMarker))
    NonlinkMarker = -1;
end

% If the input cost matrix is empty, we can just set the outputs to empty
% and stop here.
CostMatrixSize = size(CostMatrix, 1); % we need this later anyways
if ~CostMatrixSize
    Assign12 = [];
    Cost12 = [];
    return
end

% Validate the input CostMatrix.
CostMatrix = double(CostMatrix);
if (any(isnan(CostMatrix(:))) || any(isinf(CostMatrix(:))))
    error(['SMA_SPT.solveLAP(): The input cost matrix cannot contain ', ...
        'NaNs or infs'])
end
if issparse(CostMatrix)
    % Ensure that the structural rank of the CostMatrix is equal to the
    % number of rows/columns (meaning that a unique assignment can be made
    % for every row/column).
    if (sprank(CostMatrix) < CostMatrixSize)
        error(['SMA_SPT.solveLAP(): The input cost matrix must allow ', ...
            'for at least one assignment per row/column'])
    end
else
    % For non-sparse matrices we can still use the structural rank check,
    % we just want to ignore the NonlinkMarker elements.
    ValidCosts = (CostMatrix ~= NonlinkMarker);
    if (sprank(ValidCosts) < CostMatrixSize)
        error(['SMA_SPT.solveLAP(): The input cost matrix must allow ', ...
            'for at least one assignment per row/column'])
    end
end

% Convert the CostMatrix into a column vector for use in c_lap.mex64
if issparse(CostMatrix)
    [RowIndices, ColumnIndices, CostMatrixStacked] = find(CostMatrix);
else
    [RowIndices, ColumnIndices] = find(ValidCosts);
    CostMatrixStacked = CostMatrix(ValidCosts);
end

% Force appropriate variable types for the c code/for improved memory usage
% NOTE: I don't know why we add the extra 0 at the beginning, so for now
%       I'm just leaving it.
RowIndices = int32([0; RowIndices]);
ColumnIndices = int32([0; ColumnIndices]);
CostMatrixStacked = double([0; CostMatrixStacked]);

% Solve the LAP with the mex code.
StackSize = int32(numel(CostMatrixStacked));
[Assign12, ~, ~, ~] = ...
    c_lap(CostMatrixSize, StackSize, CostMatrixStacked, ...
    RowIndices, [0; find(diff(ColumnIndices)); StackSize]);
Assign12 = Assign12(2:end);
Cost12 = full(CostMatrix(sub2ind([CostMatrixSize, CostMatrixSize], ...
    Assign12, (1:CostMatrixSize).')));


end