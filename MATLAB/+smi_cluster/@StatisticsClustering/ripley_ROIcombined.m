function results = ripley_ROIcombined(obj, n_ROIs, RoI)
% Use Ripley's statistics to test the clustering of a series of ROIs
% all of the same size.
%
% INPUTS:
%    obj         various properties used by the algorithms
%       Font_props      [{'FontSize', 15, 'FontWeight', 'bold'}]
%       Line_props      [{'LineWidth', 3}]
%       Fig_ext         figure extension
%       ResultsDir      directory to store results
%       BaseName        name to identify saved plots, which will have various
%                       descriptive suffixes attached
%       Rate            [ 20]    sampling rate for statistical functions
%       Ripley_cutoff   [200]    Ripley distance cutoff (nm)
%    n_ROIs      number of ROIs to combine
%    RoI         n_ROIs cell array containing the following fields (nm):
%       ROI      ROI limits in the form
%                   [x_min, x_max, y_min, y_max {, z_min, z_max}]
%       X,Y{,Z}  (x, y, z) coordinates of points inside where X, Y, Z are of
%                the form {[n1 x 1]}, that is, referenced by X{1}
%
% OUTPUTS:
%    results     results structure:
%       r           r-axis: [0 : obj.Rate] * dr
%       K           K(r)

% Originally written by Michael Wester and Stanly Steinberg in 2008; extended
% to 3D and combined ROIs in 2017--2018.

base_name = obj.BaseName;
b_n = regexprep(base_name, '_', '\\_');

% Dimension (2D or 3D)
dim = 2;
if isfield(RoI{1}, 'Z') && numel(RoI{1}.Z) > 0
   dim = 3;
end

% Divide the interval [0, cutoff] into Rate number of pieces of length dr.
% The look at the distance d between two points and add 1 to R(d/dr).
% This provides a PDF for the distances.
dr=obj.Ripley_cutoff/obj.Rate;
K = zeros(1, obj.Rate + 1);
for m = 1 : n_ROIs
   x_min = RoI{m}.ROI(1);
   x_max = RoI{m}.ROI(2);
   y_min = RoI{m}.ROI(3);
   y_max = RoI{m}.ROI(4);
   if dim == 3
      z_min = RoI{m}.ROI(5);
      z_max = RoI{m}.ROI(6);
   end

   if dim == 2
      X = [RoI{m}.X{1}, RoI{m}.Y{1}];
   else
      X = [RoI{m}.X{1}, RoI{m}.Y{1}, RoI{m}.Z{1}];
   end
   nX = length(X);
   R = zeros(1, obj.Rate);
   if dim == 2
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
   R = R/nX;
   % Compute the intensity and normalize R.
   if dim == 2
      lambda = nX/((x_max - x_min)*(y_max - y_min));
   else
      lambda = nX/((x_max - x_min)*(y_max - y_min)*(z_max - z_min));
   end
   R = R/lambda;
   % Accumulate over each ROI.
   K = K + R;
end
% Normalize over all ROIs.
K = K / n_ROIs;

% Plot Ripley
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
grid on
hold on
% The analytic expected value.
A = pi*([0:obj.Rate]*dr).^2;
plot([0:obj.Rate]*dr, A, 'k', obj.Line_props{:})
plot([0:obj.Rate]*dr, K, 'r', obj.Line_props{:})
legend({'random', 'probes'}, 'Location', 'NorthWest');
axis([0 obj.Ripley_cutoff 0 max(K)]);
xlabel('r (nm)');
ylabel('K(r)');
title(['Ripley K [RC] for ', b_n]);
hold off
name = fullfile(obj.ResultsDir, [base_name, '_ripleyK_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% Plot sqrt of Ripley
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
hold on
grid on
M = sqrt(A/pi) - [0 : obj.Rate]*dr;
N = sqrt(K/pi) - [0 : obj.Rate]*dr;
plot([0:obj.Rate]*dr, M, 'k', obj.Line_props{:})
plot([0:obj.Rate]*dr, N, 'r', obj.Line_props{:})
legend({'random', 'probes'}, 'Location', 'Best');
xlabel('r (nm)');
ylabel('L(r) - r ');
title(['Ripley L(r) - r [RC] for ', b_n]);
hold off
name = fullfile(obj.ResultsDir, [base_name, '_ripleyL_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% Plot Ripley/(pi*r^2) (in nm)
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
hold on
grid on
% Avoid dividing by 0.
M = A(2:length(A));
M = M./(pi*(dr*[1:obj.Rate]).^2);
N = K(2:length(K));
N = N./(pi*(dr*[1:obj.Rate]).^2);
plot([0:obj.Rate]*dr, [1, M], 'k', obj.Line_props{:});
plot([0:obj.Rate]*dr, [1, N], 'r', obj.Line_props{:});
legend({'random', 'probes'}, 'Location', 'Best');
axis([0 obj.Ripley_cutoff 0 max(N)]);
xlabel('r (nm)');
ylabel('K(r)/(\pi r^2)');
title(['Normalized Ripley [RC] for ', b_n]);
hold off
name = fullfile(obj.ResultsDir, [base_name, '_ripleyNormalized_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

results.r = [0:obj.Rate]*dr;
results.K = K;

end
