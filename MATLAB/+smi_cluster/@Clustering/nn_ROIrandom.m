function hfig = nn_ROIrandom(obj, SMD, A_ROI, desc)
%nn_ROIrandom plots the PDF of nearest neighbor distances (NND) for points in a
% ROI vs a theoretical curve based on the same point density.
%
% INPUTS:
%    SMD         (x, y) coordinates in the format
%                    (1) SMD structure: SMD.X, SMD.Y (pixel)
%                    (2) N x 2 array of coordinates (nm)
%    A_ROI    [OPTIONAL] area of the ROI (nm^2)
%             [default is to compute it from the coordinate extremes]
%    desc     [OPTIONAL] description used for the plot's title [default = '']
%
% OUTPUT:
%    hfig     figure handle

% Created by
%    Michael Wester (2020)

   % Check for SMD structure; in this situation, assume the units are pixels.
   SMDstruct = false;
   if isfield(SMD, 'X') & isfield(SMD, 'Y')
      SMDstruct = true;
      XY = [SMD.X, SMD.Y] .* obj.PixelSize;
   % Check for N x 2 matrix of coordinates.
   elseif ismatrix(SMD) && size(SMD, 2) == 2
      XY = SMD;
   else
      error('Unrecognized coordinate input!');
   end

   if ~exist('desc', 'var')
      desc = '';
   else
      desc = regexprep(desc, '_', '\\_');
   end
   % Compute A_ROI from the coordinate extremes.
   if ~exist('A_ROI', 'var')
      A_ROI = (max(X) - min(X)) * (max(Y) - min(Y));
   end

   % Nearest neighbor distances (NND).
   [~, D] = knnsearch(XY, XY, 'K', 2);
   if ~isempty(D) && ~isscalar(D)
      D = D(:, 2);
   end

   P = prctile(D, 99);
   hfig = figure;
   hold on
   h = histogram(D(D < P), 30);
   maxD = max(D(D < P));

   % PDF for a random NN distribution.  Details at
   %    https://en.wikipedia.org/wiki/Mean_inter-particle_distance
   % density (#/nm^2)
   rho = size(XY, 1) / A_ROI;   % (#/nm^2)
   a = 1 / sqrt(pi * rho);   % (nm)
   r = 0 : maxD/1000 : maxD;   % (nm)
   PDF = 2/a * (r/a) .* exp(-(r/a).^2);

   %h.Normalization = 'probability';
   h.Normalization = 'PDF';
   xlabel('NN localization distances (nm)');
   ylabel('PDF');
   title(desc);
   plot(r, PDF, 'r-', 'LineWidth', 2);
   legend('data', 'random');
   hold off

end
