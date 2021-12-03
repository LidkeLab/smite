function [XY, sigma, combined, SRsave] = ...
   clusterSR(obj, XY_orig, sigma_orig, Sigma_Reg)
% Combine multiple clustered points into single localizations when appropriate
% via a top-down descent through a hierarchal dendrogram relationship between
% points.
% n is the original number of points and m is the dimension.
% n' is the final number of points after combinations have occurred.
%
% INPUTS:
%    obj          various properties used by the algorithms
%       LoS                level of significance [0 <= LoS <= 1; 0.01 typical]
%       Method             H-SET collapse method:
%                             'trivial' (for testing), 'hierarchal_singlelabel'
%       Timing             produce timings
%    XY_orig      n x m matrix of coordinates (nm)
%    sigma_orig   n x m matrix of position uncertainties (1 std deviation) (nm)
%    Sigma_Reg    1 x m array of registration error (1 standard deviation) (nm)
%
% OUTPUTS:
%    XY           n' x m final coordinate matrix (nm)
%    sigma        n' x m final position uncertainty matrix (nm)
%    combined     cell array of indices of combined points per cluster
%    SRsave       temporary storage and final results that should NOT be
%                 modified by the user
%       Sigma_Reg        sigma registration (nm)
%       XY_orig          input (x, y) (nm)
%       Sigma_orig       input sigma (nm)
%       XY               collapsed (x, y) (nm)
%       Sigma            collapsed sigma (nm)
%       Nodes_combined   indices of multiple node collapsed into single nodes

% Created by
%    Michael Wester (2019)

   SRsave.XY_orig = XY_orig;
   SRsave.Sigma_orig = sigma_orig;
   SRsave.Sigma_Reg = Sigma_Reg;

   % Find clusters in the data.
   if obj.Timing
      tic
   end
   switch obj.Method
   case 'trivial'
      % For testing purposes only.
      SRsave.XY = XY_orig;
      SRsave.Sigma = sigma_orig;

   case 'hierarchal_singlelabel'
      % Collapse multiple emitters into single emitters.
      [SRsave.XY, SRsave.Sigma, SRsave.Nodes_combined] = ...
         obj.hierarchalSingleLabel(XY_orig, sigma_orig, Sigma_Reg, obj.LoS);

   otherwise
      error('Unknown method: %s\n', obj.Method);
   end
   if obj.Timing
      toc
   end

   XY = SRsave.XY;
   sigma = SRsave.Sigma;
   combined = SRsave.Nodes_combined;

end
