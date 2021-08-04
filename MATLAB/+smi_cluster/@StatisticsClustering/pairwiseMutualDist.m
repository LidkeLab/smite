function results = pairwiseMutualDist(obj, particle_types, SMD1, SMD2)
% Plot pairwise distances and CDFs between two populations of particles.
%
% INPUTS:
%    obj             various properties used by the algorithms
%       Font_props      [{'FontSize', 15, 'FontWeight', 'bold'}]
%       Line_props      [{'LineWidth', 3}]
%       Fig_ext         figure extension
%       ResultsDir      directory to store results
%       BaseName        name to identify saved plots, which will have various
%                       descriptive suffixes attached
%       Nsims           number of random simulations to perform
%    particle_types   cell array of particle types (names)
%                     e.g., particle_types = {'5', '10'};
%    SMD1 and SMD2    (x, y) coordinates of the dataset (nm) in the format
%                     (1) SMD structure: SMD1.X, SMD1.Y, SMD2.X and SMD2.Y,
%                     (2) N x 2 array of coordinates,
%                     (3) RoI struct with fields (nm):
%       ROI               [x_min, x_max, y_min, y_max, {z_min, z_max}]] of ROI
%       X, Y{, Z}         (x, y, z) coordinates of points inside where X, Y, Z
%                         are of the form {[n1 x 1]}
%
% OUTPUTS:
%    results          results structure:
%       bins             number of histogram bins
%       sims             number of random simulations performed
%       x, y             cell array of point coordinates for random simulations
%       xx               
%       f                

% Created by
%    Michael J. Wester (2021)

   base_name = obj.BaseName;

   % Dimension (2D or 3D)
   dim = 2;
   if ~exist('SMD2', 'var')
      if iscell(SMD1) && ismatrix(SMD1) && isfield(SMD1{1}, 'ROI')
         RoI = SMD1;
         if isempty(obj.ROI)
            obj.ROI = RoI{1}.ROI;
         end
         if isfield(RoI{1}, 'Z') && numel(RoI{1}.Z) > 0
            dim = 3;
         end
      else
         error('Unrecognized format for SMD1!');
      end
   elseif ismatrix(SMD1) && ismatrix(SMD2) ...
      && size(SMD1, 2) > 1 && size(SMD2, 2) > 1
      RoI{1}.X = {SMD1(:, 1), SMD2(:, 1)};
      RoI{1}.Y = {SMD1(:, 2), SMD2(:, 2)};
      if size(SMD1, 2) == 3 && size(SMD2, 2) == 3
         dim == 3;
         RoI{1}.Z = {SMD1(:, 3), SMD2(:, 3)};
      end
   elseif isstruct(SMD1) && isstruct(SMD2)
      RoI{1}.X = {SMD1.X, SMD2.X};
      RoI{1}.Y = {SMD1.Y, SMD2.Y};
      if isfield(SMD1, 'Z') && isfield(SMD2, 'Z') ...
         && numel(SMD1.Z) > 0 && numel(SMD2.Z) > 0
         dim == 3;
         RoI{1}.Z = {SMD1.Z, SMD2.Z};
      end
   else
      error('Unrecognized format for SMD1/SMD2!');
   end

   bins = 100;
   sims = obj.Nsims;
   %sims = 20;

   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on

   if length(particle_types) ~= 2
      error('Mutual distances can only be computed between 2 particle types!');
   end

% P   cell array of point coordinates (n x 2 per cell element)
%     e.g., P{1} = [x1, y1];   P{2} = [x2, y2];
%     where P{1} is (n1 x 2), P{2} is (n2 x 2)

   if dim == 2
      P{1} = [RoI{1}.X{1}, RoI{1}.Y{1}];
      P{2} = [RoI{1}.X{2}, RoI{1}.Y{2}];
   else
      P{1} = [RoI{1}.X{1}, RoI{1}.Y{1}, RoI{1}.Z{1}];
      P{2} = [RoI{1}.X{2}, RoI{1}.Y{2}, RoI{1}.Z{2}];
   end

   X = P{1};
   Y = P{2};
   D = pdist2(X, Y);
   D = D(:);
%  M = mean(D);
%  S = std(D);

   % Establish the positions of the bin centers (x) based on a normal random
   % distribution that inherits its characteristics from the coordinates in
   % X given the specified number of bins (also, see comment below).
%  Xr = M + S*randn(size(X));
%  Yr = M + S*randn(size(Y));
%  Dr = pdist2(Xr, Yr);
%  [~, x] = hist(Dr, bins);

%  y = hist(D, x);
   [y, x] = hist(D, bins);
   plot(x, y, 'k-', obj.Line_props{:});

   % Create a random distribution with the same characteristics as X:
   % based on the same number of points, mean, standard deviation.  Average
   % over the specified number of simulations (sims).
%  YR = 0;
%  for i = 1 : sims
%     Xr = M + S*randn(size(X));
%     Dr = pdist(Xr);
%     yr = hist(Dr, x);
%     YR = YR + yr;
%  end
%  yr = YR / sims;
%  plot(x, yr, 'r--', obj.Line_props{:});

%  legend('data', 'random', 'Location', 'NorthEast');
   title(['Pairwise Distance PDF for ', base_name, '\_', ...
          particle_types{1}, ',', particle_types{2}]);
   xlabel('distance (nm)');
   ylabel('frequency');
   hold off
   name = fullfile(obj.ResultsDir, [base_name, '_', particle_types{1}, ...
                                    ',', particle_types{2},            ...
                                    '_pairwisePDF2']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on

   [f, xx] = ecdf(sort(D));
   plot(xx, f, 'k-', obj.Line_props{:});

   %f = cumsum(1 : numel(D));
   %f = f / f(end);
   %xx = sort(D);
   %plot(xx, f, 'b-', obj.Line_props{:});

%  [fr, xxr] = ecdf(sort(Dr));
%  plot(xxr, fr, 'r--', obj.Line_props{:});

   %fr = cumsum(1 : numel(Dr));
   %fr = fr / fr(end);
   %xxr = sort(Dr);
   %plot(xxr, fr, 'm--', obj.Line_props{:});

%  legend('data', 'random', 'Location', 'SouthEast');
   title(['Pairwise Distance CDF for ', base_name, '\_', ...
          particle_types{1}, ',', particle_types{2}]);
   xlabel('distance (nm)');
   ylabel('frequency');
   hold off
   name = fullfile(obj.ResultsDir, [base_name, '_', particle_types{1}, ...
                                    ',', particle_types{2},            ...
                                    '_pairwiseCDF2']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end

results.bins = bins;
results.sims = sims;
results.x = x;
results.y = y;
results.xx = xx;
results.f = f;

end
