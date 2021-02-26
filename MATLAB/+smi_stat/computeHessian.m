function [HessianMatrix] = computeHessian(FunctionHandle, Location, ...
    DeltaHFraction, DeltaHBound)
%computeHessian computes the Hessian of FunctionHandle around ParamsHat.
% This method computes the Hessian matrix of the function handle
% FunctionHandle evaluated at Location.
%
% INPUTS:
%   FunctionHandle: A function handle which expects numel(Location)
%                   parameters, with the indices corresponding to those of
%                   Location.  FunctionHandle must output a single value.
%   Location: The arguments of FunctionHandle at which we'll estimate the
%             Hessian.
%   DeltaHFraction: Fraction of each element of Location, where the step 
%                   size (i.e., Delta h) along each dimension will be set
%                   to, e.g., h_n = DeltaHFraction*Location(n) subject to
%                   the bounds given in DeltaHBound.
%                   (Default = 1e-2 chosen arbitrarily)
%   DeltaHBound: Bounds to the step size Delta h.
%                ([min. value, max. value])
%                (Default = min(Location)*[1e-9; 1e-1] chosen arbitrarily)
%
% OUTPUTS:
%   HessianMatrix: The matrix containing an estimate of the Hessian of
%                  FunctionHandle evaluated at Location.
%                  (numel(Location) x numel(Location) numeric array)

% Created by:
%   David J. Schodt (Lidke Lab, 2021) based on an earlier code from
%       someone (unknown) in the Lidke lab


% Set defaults if needed.
if (~exist('DeltaHFraction', 'var') || isempty(DeltaHFraction))
    DeltaHFraction = 1e-2;
end
if (~exist('DeltaHBound', 'var') || isempty(DeltaHBound))
    DeltaHBound = min(Location)  * [1e-9; 1e-1];
end

% Loop through the parameters and construct the Hessian element by element.
% NOTE: We do this by finding the numerical derivative of the function
%       handle w.r.t. the parameter Location(ii) at two different 
%       positions of Location(jj), and then finding the (1/DeltaH) scaled
%       difference between those two derivatives.
NParams = numel(Location);
HessianMatrix = zeros(NParams);
for ii = 1:NParams
    for jj = 1:NParams
        % Define the step sizes along each dimension for the derivative.
        DeltaI = min(DeltaHBound(2), ...
            max(DeltaHBound(1), DeltaHFraction*Location(ii)));
        DeltaJ = min(DeltaHBound(2), ...
            max(DeltaHBound(1), DeltaHFraction*Location(jj)));
        if (ii == jj)
            % Estimate the derivative of the function handle w.r.t.
            % parameter ii to the left of Location.
            InputParams1 = Location; % initialize
            InputParams2 = Location; % initialize
            InputParams1(ii) = Location(ii) - DeltaI;
            Derivative1 = (1/DeltaI) ...
                * (FunctionHandle(InputParams2) ...
                - FunctionHandle(InputParams1));
            
            % Estimate the derivative of the function handle w.r.t. 
            % parameter ii to the right of Location.
            InputParams1 = Location; % initialize
            InputParams2 = Location; % initialize
            InputParams2(ii) = Location(ii) + DeltaI;
            Derivative2 = (1/DeltaI) ...
                * (FunctionHandle(InputParams2) ...
                - FunctionHandle(InputParams1));
        else
            % Estimate the derivative of the function handle w.r.t. 
            % parameter ii at a set value of parameter jj (set jj to the 
            % left parameter value).
            InputParams1 = Location;
            InputParams2 = Location;
            InputParams1(ii) = Location(ii) - DeltaI/2;
            InputParams2(ii) = Location(ii) + DeltaI/2;
            InputParams1(jj) = Location(jj) - DeltaJ/2;
            InputParams2(jj) = Location(jj) - DeltaJ/2;
            Derivative1 = (1/DeltaI) ...
                * (FunctionHandle(InputParams2) ...
                - FunctionHandle(InputParams1));
            
            % Estimate the derivative of the function handle w.r.t. 
            % parameter ii at a different set value of parameter jj (set 
            % jj to the right parameter value).
            InputParams1 = Location;
            InputParams2 = Location; %
            InputParams1(ii) = Location(ii) - DeltaI/2;
            InputParams2(ii) = Location(ii) + DeltaI/2;
            InputParams1(jj) = Location(jj) + DeltaJ/2;
            InputParams2(jj) = Location(jj) + DeltaJ/2;
            Derivative2 = (1/DeltaI) ...
                * (FunctionHandle(InputParams2) ...
                - FunctionHandle(InputParams1));
        end
        
        % Compute Hessian(ii, jj) based on the difference between the above
        % two partial derivatives.
        HessianMatrix(ii, jj) = (1/DeltaJ) * (Derivative2-Derivative1);
    end
end


end