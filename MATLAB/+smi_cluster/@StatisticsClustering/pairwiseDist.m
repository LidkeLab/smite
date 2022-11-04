function results = pairwiseDist(obj, particle_types, SMD)
%pairwiseDist plots pairwise distances and CDFs.
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
%    particle_types  cell array of particle types (names)
%                    e.g., particle_types = {'5', '10'};
%    SMD              (x, y) coordinates of the dataset (nm) in the format
%                     (1) SMD structure: SMD.X and SMD.Y,
%                     (2) N x 2 array of coordinates,
%                     (3) RoI struct with fields (nm):
%       ROI               [x_min, x_max, y_min, y_max, {z_min, z_max}]] of ROI
%       X, Y{, Z}         (x, y, z) coordinates of points inside where X, Y, Z
%                         are of the form {[n1 x 1]}
%
% OUTPUTS:
%    results         results structure:
%       bins            number of histogram bins
%       sims            number of random simulations performed
%       x_, y_          cell array of point coordinates for random simulations
%       y_r             mean y-random coordinates

% Created by
%    Michael J. Wester (2020)

   base_name = obj.BaseName;
   base_text = regexprep(base_name, '_', '\\_');

   % Dimension (2D or 3D)
   dim = 2;
   if iscell(SMD) && ismatrix(SMD) && isfield(SMD{1}, 'ROI')
      RoI = SMD;
      if isfield(RoI{1}, 'Z') && numel(RoI{1}.Z) > 0
         dim = 3;
      end
   elseif ismatrix(SMD) && size(SMD, 2) > 1
      RoI{1}.X = {SMD(:, 1)};
      RoI{1}.Y = {SMD(:, 2)};
      if size(SMD, 2) == 3
         dim == 3;
         RoI{1}.Z = {SMD1(:, 3)};
      end
   elseif isstruct(SMD)
      RoI{1}.X = {SMD.X};
      RoI{1}.Y = {SMD.Y};
      if isfield(SMD, 'Z') && numel(SMD.Z) > 0
         dim == 3;
         RoI{1}.Z = {SMD.Z};
      end
   else
      error('Unrecognized format for SMD1/SMD2!');
   end

   bins = 100;
   sims = obj.Nsims;
   %sims = 20;

   for m = 1 : length(particle_types)
      if ~isempty(obj.Fig_ext)
         figure('Visible', 'off');
      else
         figure;
      end
      axes(obj.Font_props{:})
      hold on

      if dim == 2
         X = [RoI{1}.X{m}, RoI{1}.Y{m}];
      else
         X = [RoI{1}.X{m}, RoI{1}.Y{m}, RoI{1}.Z{m}];
      end
      D = pdist(X);
      M = mean(D);
      S = std(D);

      % Establish the positions of the bin centers (x) based on a normal random
      % distribution that inherits its characteristics from the coordinates in
      % X given the specified number of bins (also, see comment below).
      Xr = M + S*randn(size(X));
      Dr = pdist(Xr);
      [~, x] = hist(Dr, bins);

      y = hist(D, x);
      plot(x, y, 'k-', obj.Line_props{:});

      % Create a random distribution with the same characteristics as X:
      % based on the same number of points, mean, standard deviation.  Average
      % over the specified number of simulations (sims).
      YR = 0;
      for i = 1 : sims
         Xr = M + S*randn(size(X));
         Dr = pdist(Xr);
         yr = hist(Dr, x);
         YR = YR + yr;
      end
      yr = YR / sims;
      plot(x, yr, 'r--', obj.Line_props{:});

      x_{m} = x;
      y_{m} = y;
      y_r{m} = yr;

      legend('data', 'random', 'Location', 'NorthEast');
      title(['Pairwise Distance PDF for ', base_text, '\_', ...
             particle_types{m}]);
      xlabel('distance (nm)');
      ylabel('frequency');
      hold off
      name = fullfile(obj.ResultsDir, ...
                      [base_name, '_', particle_types{m}, '_pairwisePDF']);
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

      [fr, xxr] = ecdf(sort(Dr));
      plot(xxr, fr, 'r--', obj.Line_props{:});

      %fr = cumsum(1 : numel(Dr));
      %fr = fr / fr(end);
      %xxr = sort(Dr);
      %plot(xxr, fr, 'm--', obj.Line_props{:});

      legend('data', 'random', 'Location', 'SouthEast');
      title(['Pairwise Distance CDF for ', base_text, '\_', ...
             particle_types{m}]);
      xlabel('distance (nm)');
      ylabel('frequency');
      hold off
      name = fullfile(obj.ResultsDir, ...
                      [base_name, '_', particle_types{m}, '_pairwiseCDF']);
      if ~isempty(obj.Fig_ext)
         print(['-d', obj.Fig_ext], name);
      else
         saveas(gcf, name);
         delete(gcf);
      end
   end

results.bins = bins;
results.sims = sims;
results.x_ = x_;
results.y_ = y_;
results.y_r = y_r;

end
