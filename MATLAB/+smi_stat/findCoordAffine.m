function [AffineTransform] = findCoordAffine(Coords1, Coords2, MaxDist)
%findCoordAffine finds an affine transform to transform Coords2 to Coords1.
% This method will compute an affine transform which can transform Coords2
% coordinates into alignment with Coords1 coordinates. The cost function for
% the affine transform computation is given as the sum of the thresholded
% nearest-neighbor distances, thresholded to MaxDist.
%
% INPUTS:
%   Coords1: [X, Y] coordinates. (Nx2 numeric array)
%   Coords2: [X, Y] coordinates that are to be aligned to Coords1 via the 
%            output AffineTransform. (Nx2 numeric array)
%   MaxDist: Nearest-neighbor distance threshold between Coords1 and 
%            Coords2 localizations. (Default = 1 pixel)
%
% OUTPUTS:
%   AffineTransform: An affine transform that can be used to transform
%                    Coords1 coordinates into alignment with Coords2
%                    localizations.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Set defaults if needed.
if (~exist('MaxDist', 'var') || isempty(MaxDist))
    MaxDist = 1;
end

% Convert input coordinates to type double.
Coords1 = double(Coords1);
Coords2 = double(Coords2);

% Determine the center of mass shift between the coordinates (will be used
% as an initial guess for the affine transform).
COMShift = mean(Coords1, 1) - mean(Coords2, 1);
XInitial = [1, 0, 0, 1, COMShift];

% Perform an unconstrained fit to estimate the affine transform values.
CostFunction = @(X) nndSum(Coords1, Coords2, affine2d(...
    [X(1), X(2), 0; ...
    X(3), X(4), 0; ...
    X(5), X(6), 1]));
XHat = fminsearch(CostFunction, XInitial);
AffineTransform = affine2d(...
    [XHat(1), XHat(2), 0; ...
    XHat(3), XHat(4), 0; ...
    XHat(5), XHat(6), 1]);


    function Cost = nndSum(Coords1, Coords2, AffineTransform)
        % Transform Coords2 with the affine transform.
        Coords2Transformed = ...
            smi_core.ChannelRegistration.transformCoords(...
            AffineTransform, Coords2);
        
        % Compute the nearest-neighbor distances between Coords1 and the
        % transformed Coords2.
        [~, NNDistances] = knnsearch(Coords2Transformed, Coords1);
        
        % Compute the thresholded sum of the nearest-neighbor distances.
        NNDistances(NNDistances > MaxDist) = MaxDist;
        Cost = sum(NNDistances);
    end


end