function results = ripley(obj, particle_types, SMD)
% Use Ripley's statistics to test the clustering of the points in P.
%
% INPUTS:
%    obj              various properties used by the algorithms
%       Font_props      [{'FontSize', 15, 'FontWeight', 'bold'}]
%       Line_props      [{'LineWidth', 3}]
%       Fig_ext         figure extension
%       ResultsDir      directory to store results
%       BaseName        name to identify saved plots, which will have various
%                       descriptive suffixes attached
%       Rate            [ 20]    sampling rate for statistical functions
%       Ripley_cutoff   [200]    Ripley distance cutoff (nm)
%    particle_types   cell array of particle types (names)
%                     e.g., particle_types = {'5', '10'};
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
%       r                r-axis: [0 : obj.Rate] * dr
%       K                K(r)

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

M = length(particle_types);
for i = 1 : M
   probes{i} = particle_types{i};
end

colors = {'r', 'g', 'b'};

% Dimension (2D or 3D)
dim = 2;
if D_nm > 0
   dim = 3;
end

% Compute the Ripley statistcs.

% Divide the interval [0, cutoff] into Rate number of pieces of length dr.
% The look at the distance d between two points and add 1 to R(d/dr).
% This provides a PDF for the distances.
dr=obj.Ripley_cutoff/obj.Rate;
for m = 1 : length(particle_types)
   if obj.Verbose >= 2
      fprintf('Compute Ripley statistic for %s %s (%dD).\n', ...
              base_name, particle_types{m}, dim);
   end
   if dim == 2
      P = [RoI{1}.X{m}, RoI{1}.Y{m}];
   else
      P = [RoI{1}.X{m}, RoI{1}.Y{m}, RoI{1}.Z{m}];
   end
   X = P;
   nX = length(X);
   R = zeros(1, obj.Rate);
   if dim == 2
      %for i = 1 : length(X)
      for i = 1 : nX
         for j = 1 : i - 1
            p = ceil(sqrt((X(j,1) - X(i,1))^2 + (X(j,2) - X(i,2))^2)/dr);
            if p > 0 & p <= obj.Rate
               R(p) = R(p) + 1;
            end
         end
      end
   else
      for i = 1 : nX
         for j = 1 : i - 1
            p = ceil(sqrt((X(j,1) - X(i,1))^2 + (X(j,2) - X(i,2))^2 + ...
                          (X(j,3) - X(i,3))^2)/dr);
            if p > 0 & p <= obj.Rate
               R(p) = R(p) + 1;
            end
         end
      end
   end
   % Add a zero to the beginning of R for nicer plotting.
   R = [0, R];
   % Compensate for the j=1:i-1 savings.
   R = 2*R;
   % Convert the PDF to a CDF.
   for i = 2 : length(R)
      R(i) = R(i) + R(i-1);
   end
   % Average over the number of particles.
   %R = R/length(X);
   R = R/nX;
   % Compute the intensity and normalize R.
   if dim == 2
      %lambda = length(X)/(H_nm*V_nm);
      lambda = nX/(H_nm*V_nm);
   else
      lambda = nX/(H_nm*V_nm*D_nm);
   end
   R = R/lambda;
   K(m) = {R};
end

% Plot Ripley
if obj.Verbose >= 2
   fprintf('   Plot Ripley statistic.\n');
end
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
grid on
hold on
% The analytic expected value.
plot([0:obj.Rate]*dr, pi*([0:obj.Rate]*dr).^2, 'k', ...
     obj.Line_props{:})
mx = 0;
for m = 1 : length(particle_types)
   N = K{m};
   mx = max(mx, max(N));
   plot([0:obj.Rate]*dr, N, colors{m}, obj.Line_props{:})
end
legend({'random', probes{:}}, 'Location', 'NorthWest');
axis([0 obj.Ripley_cutoff 0 mx]);
xlabel('r (nm)');
ylabel('K(r)');
title(['Ripley for ', base_text]);
hold off
name = fullfile(obj.ResultsDir, [base_name, '_ripley']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% Plot sqrt of Ripley
if obj.Verbose >= 2
   fprintf('   Plot the L(r) statistic.\n');
end
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
grid on
hold on
mx = 0;
for m = 1 : length(particle_types)
   N = sqrt(K{m}/pi) - [0 : obj.Rate]*dr;
%  mx = max(mx, max(N));
   plot([0:obj.Rate]*dr, N, colors{m}, obj.Line_props{:})
end
legend(probes, 'Location', 'NorthWest');
xlabel('r (nm)');
ylabel('L(r) - r ');
title(['Ripley L(r) - r for ', base_text]);
hold off
name = fullfile(obj.ResultsDir, [base_name, '_l_r_ripley']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% Plot Ripley/(pi*r^2) (in microns)
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
grid on
hold on
mx = 0;
if obj.Verbose >= 2
   fprintf('   Plot the K(r)/(pi r^2) statistic.\n');
end
for m = 1 : length(particle_types)
   % Avoid dividing by 0.
   N = K{m}(2:length(K{m}));
   N = N./(pi*(dr*[1:obj.Rate]).^2);
   mx = max(mx, max(N));
   plot([0:obj.Rate]*dr, [0, N], colors{m}, obj.Line_props{:});
end
legend(probes, 'Location', 'NorthEast');
axis([0 obj.Ripley_cutoff 0 mx]);
xlabel('r (nm)');
ylabel('K(r)/(\pi r^2)');
title(['Normalized Ripley for ', base_text]);
hold off
name = fullfile(obj.ResultsDir, [base_name, '_norm_ripley']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% m       is the number of different types of particles.
% n       is the number of images in the experiment.
% Rate    is the sampling rate.

results.r = [0:obj.Rate]*dr;
results.K = K;

end
