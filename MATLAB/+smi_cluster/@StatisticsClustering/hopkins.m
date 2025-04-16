function results = hopkins(obj, particle_types, SMD)
% Use the Hopkins' statistic to test the clustering of the points in SMD.
%
% INPUTS:
%    obj             various properties used by the algorithms
%       Font_props      [{'FontSize', 15, 'FontWeight', 'bold'}]
%       Line_props      [{'LineWidth', 3}]
%       Fig_ext         figure extension
%       ResultsDir      directory to store results
%       BaseName        name to identify saved plots, which will have various
%                       descriptive suffixes attached
%       Rate            [ 20]    sampling rate for statistical functions
%       HopTestPts      number of Hopkins' statistic test points
%       HopTests        number of tests to run to produce a Hopkins' statistic
%    particle_types   cell array of particle types (names)
%                     e.g., particle_types = {'5', '10'};
%                         are of the form {[n1 x 1], [n2 x 1]}
%    SMD              (x, y) coordinates of the dataset (nm) in the format
%                     (1) SMD structure: SMD.X and SMD.Y,
%                     (2) N x 2 array of coordinates,
%                     (3) RoI struct with fields (nm):
%       ROI               [x_min, x_max, y_min, y_max, {z_min, z_max}]] of ROI
%       X, Y{, Z}         (x, y, z) coordinates of points inside where X, Y, Z
%                         are of the form {[n1 x 1]}
%
% OUTPUTS:
%    results          results structure:
%       test             number of test points
%       ntests           number of Hopkins statistics to use
%       bins             number of bins for the probability graphs
%       fitting          Gaussian fitting of the data? (boolean)
%       H                Hopkins' statistic
%       xa               analytic distribution for the Hopkins test (x-axis)
%       ya               analytic distribution for the Hopkins test (y-axis)

% Created by
%    Michael Wester and Stanly Steinberg (2008)

base_name = obj.BaseName;
base_text = regexprep(base_name, '_', '\\_');

% Dimension (2D or 3D)
dim = 2;
if iscell(SMD) && ismatrix(SMD) && isfield(SMD{1}, 'ROI')
   RoI = SMD;
   if isempty(obj.ROI)
      obj.ROI = RoI{1}.ROI;
   end
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

% H_nm             horizontal size of the ROI (nm)
% V_nm             vertical   size of the ROI (nm)
% D_nm             depth      size of the ROI (nm) for 3D stats [OPTIONAL]
H_nm = obj.ROI(2) - obj.ROI(1);
V_nm = obj.ROI(4) - obj.ROI(3);
if dim == 3
   D_nm = obj.ROI(6) - obj.ROI(5);
else
   D_nm = -1;
end

%test = 5;        % The number of test points.
%ntests = 1000;   % The number of Hopkins statistics to use.
test = obj.HopTestPts;   % The number of test points.
ntests = obj.HopTests;   % The number of Hopkins statistics to use.
bins = 100;      % The number of bins for the probability graphs.
fitting = false; % Gaussian fitting of the data.

% Compute the Hopkins' statistic for each negative.
for m = 1 : length(particle_types)
   if obj.Verbose >= 2
      fprintf('Compute the Hopkins statistic for %s %s (%dD).\n', ...
              base_name, particle_types{m}, dim);
   end
   if dim == 2
      P = [RoI{1}.X{m}, RoI{1}.Y{m}];
   else
      P = [RoI{1}.X{m}, RoI{1}.Y{m}, RoI{1}.Z{m}];
   end
   X = [];
   if D_nm <= 0   % 2D
      for i = 1 : ntests
         X = [X, smi_cluster.StatisticsClustering.hopkinstat( ...
                    P, H_nm, V_nm, test)];
      end
   else   % 3D
      for i = 1 : ntests
         X = [X, smi_cluster.StatisticsClustering.hopkinstat3( ...
                    P, H_nm, V_nm, D_nm, test)];
      end
   end
   H{m} = X;
end

% The analytic distribution for the Hopkins test:
xa = linspace(0, 1, obj.Rate);
ya = ((xa.^(test-1)).*((1-xa).^(test-1)))/(gamma(test)^2/gamma(2*test));

% Plot the summary statistics.
% The PDF for each negative in the experiment.
if obj.Verbose >= 2
   fprintf('   Plot the Hopkins statistics.\n');
end
for m = 1 : length(particle_types)
   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:})
   hold on
   [X, V] = smi_cluster.StatisticsClustering.histogram(H{m}, [0, 1], bins);
   Pr = V*bins/ntests;
   bar(X, Pr, 1)
   hold on
   if fitting
      % Gaussian fit to the data.
      f = fit(X', Pr', 'gauss2');
      [fit_max, i_max] = max(f(X));
      x_max = X(i_max);
      h = plot(f, 'g');
      set(h, obj.Line_props{1}, obj.Line_props{2});
   end
   plot(xa, ya, 'r', obj.Line_props{:})
   axis([0 1 0 ceil(max(max(Pr), max(ya)))]);
   xlabel('H -- The Hopkins Statistic');
   ylabel('PDF');
   if fitting
      title({['Hopkins PDF for ', base_text, '\_', particle_types{m}], ...
             sprintf('(x, fit)_{max} = (%5.3f, %.3f)', x_max, fit_max)});
   else
      title(['Hopkins PDF for ', base_text, '\_', particle_types{m}]);
   end
   if fitting
      legend('data', 'Gaussian fit', 'random', 'Location', 'NorthWest');
   else
      legend('data', 'random', 'Location', 'NorthWest');
   end
   name = fullfile(obj.ResultsDir, ...
                   [base_name, '_', particle_types{m}, '_hopkinspdf']);
   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
   else
      saveas(gcf, name);
      delete(gcf);
   end
end

hold off

results.test = test;
results.ntests = ntests;
results.bins = bins;
results.fitting = fitting;
results.H = H;
results.xa = xa;
results.ya = ya;

end
