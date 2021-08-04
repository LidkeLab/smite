function results = ...
   bivariateRipley_ROIcombined(obj, particle_types, n_ROIs, RoI)
% Use bivariate Ripley statistics to test the clustering of a series of ROIs
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
%
%       Rate            [ 20]    sampling rate for statistical functions
%       Ripley_cutoff   [200]    Ripley distance cutoff (nm)
%       Confidence      [2.576]  bivariate Ripley confidence interval:
%                                   95% confidence interval: Confidence = 1.96
%                                   99% confidence interval: Confidence = 2.576
%       Nsims         = [100]    bivariate Ripley simulations to run
%    particle_types   cell array of particle types (names)
%                     e.g., particle_types = {'5', '10'};
%    n_ROIs           number of ROIs to combine
%    RoI              n_ROIs cell array containing the following fields (nm):
%       ROI              ROI limits in the form
%                           [x_min, x_max, y_min, y_max {, z_min, z_max}]
%       X, Y, Z          cell arrays of (x, y {, z}) coordinates
%                           e.g., X{1} = 1st label x-coordinates, X{2} = 2nd
%                           label where X{1} is (n1 x 1), X{2} is (n2 x 1)
%
% OUTPUTS:
%    results          results structure:
%       R               L(R) statistic
%       E               simulation results
%       A               simulation means
%       High            A + stderr
%       Low             A - stderr

% Created by
%    Michael Wester and Stanly Steinberg (2008)

base_name = obj.BaseName;

% Dimension (2D or 3D)
dim = 2;
if isfield(RoI{1}, 'Z') && numel(RoI{1}.Z) > 0
   dim = 3;
end

% Compute the bivariate Ripley statistics.
% Add Q for second species
if obj.Verbose >= 2
   fprintf('Compute the bivariate Ripley statistic for %s.\n', base_name);
end

% Divide the interval [0, cutoff] into Rate number of pieces of length dr.
% The look at the distance d between two points and add 1 to R(d/dr).
% This provides a PDF for the distances.
dr = obj.Ripley_cutoff/obj.Rate;
RipleyK = zeros(1, obj.Rate + 1);;
% H_nm   horizontal size of the ROI (nm)
% V_nm   vertical   size of the ROI (nm)
% D_nm   depth      size of the ROI (nm)
for m = 1 : n_ROIs
   x_min = RoI{m}.ROI(1);
   x_max = RoI{m}.ROI(2);
   y_min = RoI{m}.ROI(3);
   y_max = RoI{m}.ROI(4);
   H_nm = x_max - x_min;
   V_nm = y_max - y_min;
   if dim == 3
      z_min = RoI{m}.ROI(5);
      z_max = RoI{m}.ROI(6);
      D_nm = z_max - z_min;
   end

   X = [RoI{m}.X{1}, RoI{m}.Y{1}];
   Y = [RoI{m}.X{2}, RoI{m}.Y{2}];
   if dim == 3
      X = [RoI{m}.X{1}, RoI{m}.Y{1}, RoI{m}.Z{1}];
      Y = [RoI{m}.X{2}, RoI{m}.Y{2}, RoI{m}.Z{2}];
   end
   nX = length(X);
   nY = length(Y);
   R = zeros(1, obj.Rate);
   if dim == 2
      %for i = 1 : length(X)
      %   for j = 1 : length(Y)
      for i = 1 : nX
         for j = 1 : nY
            p = ceil(sqrt((X(i,1) - Y(j,1))^2 + (X(i,2) - Y(j,2))^2)/dr);
            if p > 0 & p <= obj.Rate
               R(p) = R(p) + 1;
            end
         end
      end
   else   % dim == 3
      for i = 1 : nX
         for j = 1 : nY
            p = ceil(sqrt((X(i,1) - Y(j,1))^2 + (X(i,2) - Y(j,2))^2 + ...
                          (X(i,3) - Y(j,3))^2)/dr);
            if p > 0 & p <= obj.Rate
               R(p) = R(p) + 1;
            end
         end
      end
   end

   % Add a zero to the beginning of R for nicer plotting.
   R = [0, R];
   % Convert the PDF to a CDF.
   for i = 2 : length(R)
      R(i) = R(i) + R(i-1);
   end
   
   % Compute the intensity and normalize R.
   % This needs a fix. xxx
   % xxx should use the area of the crop_area (err crop region).
   % A = H_nm*V_nm = (x_max - x_min)*(y_max - y_min)
   if dim == 2
      lambda1 = nX/(H_nm*V_nm);
      lambda2 = nY/(H_nm*V_nm);
      R = R/(lambda1*lambda2*(H_nm*V_nm));
   else
      lambda1 = nX/(H_nm*V_nm*D_nm);
      lambda2 = nY/(H_nm*V_nm*D_nm);
      R = R/(lambda1*lambda2*(H_nm*V_nm*D_nm));
   end
   K(1, 2) = {R};
   % Accumulate over each ROI.
   RipleyK = RipleyK + R;
end
% Normalize over all ROIs.
RipleyK = RipleyK / n_ROIs;
K(1, 2) = {RipleyK};

% Plot sqrt of Ripley
if obj.Verbose >= 2
   fprintf('   Plotting the L(r) statistic for %s.\n', base_name);
end
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
grid on
hold on
N = sqrt(K{1, 2}/pi) - [0:obj.Rate]*dr;
plot([0:obj.Rate]*dr, N, '-r', obj.Line_props{:})
xlabel('Distance (nm)');
ylabel('L(r)-r ');
NN = N;

for s = 1 : obj.Nsims
   % Create random point vectors for simulations
   % and calculate Ripley statistic
   if dim == 2
      %Xr = [rand(length(X),1)*H_nm rand(length(X),1)*V_nm];
      %Yr = [rand(length(Y),1)*H_nm rand(length(Y),1)*V_nm];
      Xr = [rand(nX,1)*H_nm rand(nX,1)*V_nm];
      Yr = [rand(nY,1)*H_nm rand(nY,1)*V_nm];
      Rr = zeros(1, obj.Rate);
      %for i = 1 : length(Xr)
      %   for j = 1 : length(Yr)
      for i = 1 : nX
         for j = 1 : nY
            p=ceil(sqrt((Xr(i,1) - Yr(j,1))^2 + (Xr(i,2) - Yr(j,2))^2)/dr);
            if p > 0 & p <= obj.Rate
               Rr(p) = Rr(p) + 1;
            end
         end
      end
   else
      Xr = [rand(nX,1)*H_nm rand(nX,1)*V_nm rand(nX,1)*D_nm];
      Yr = [rand(nY,1)*H_nm rand(nY,1)*V_nm rand(nY,1)*D_nm];
      Rr = zeros(1, obj.Rate);
      for i = 1 : nX
         for j = 1 : nY
            p=ceil(sqrt((Xr(i,1) - Yr(j,1))^2 + (Xr(i,2) - Yr(j,2) + ...
                         Xr(i,3) - Yr(j,3))^2)/dr);
            if p > 0 & p <= obj.Rate
               Rr(p) = Rr(p) + 1;
            end
         end
      end
   end
   % Add a zero to the beginning of R for nicer plotting.
   Rr = [0, Rr];
   % Convert the PDF to a CDF.
   for i = 2 : length(Rr)
      Rr(i) = Rr(i) + Rr(i-1);
   end
   % Compute the intensity and normalize R.
% Fix area here also xxx.
   if dim == 2
      %lambda1r = length(Xr)/(H_nm*V_nm);
      %lambda2r = length(Yr)/(H_nm*V_nm);
      lambda1r = nX/(H_nm*V_nm);
      lambda2r = nY/(H_nm*V_nm);
      Rr = Rr/(lambda1*lambda2*(H_nm*V_nm));
   else
      lambda1r = nX/(H_nm*V_nm*D_nm);
      lambda2r = nY/(H_nm*V_nm*D_nm);
      Rr = Rr/(lambda1*lambda2*(H_nm*V_nm*D_nm));
   end
   N = sqrt(Rr/pi) - [0:obj.Rate]*dr;
   % Plot all random simulations
%  plot([0:obj.Rate]*dr, N, obj.Line_props{:})
   hold on
   % Put Ripley data into a matrix
   E(s, 1:length(Rr))=N;
end

for ss = 1 : length(Rr)
   % Calculate standard deviation of all simulations at particular dr
   S(ss, 1) = std(E(:, ss));
   % Calculate mean of all simulations at particular dr 
   A(ss) = mean(E(:, ss));
   High(ss) = A(ss)+obj.Confidence*S(ss);
   Low(ss)  = A(ss)-obj.Confidence*S(ss);
end

% Plot data
% grid off
% set(gca, 'Box', 'on', 'LineWidth', 3, 'FontSize', 18, 'FontWeight', 'bold')
plot([0:obj.Rate]*dr, High, '--k', obj.Line_props{:});
plot([0:obj.Rate]*dr, Low,  '--k', obj.Line_props{:});
plot([0:obj.Rate]*dr, H_nm, '*k', obj.Line_props{:});

% Fit the figure
b_n = regexprep(base_name, '_', '\\_');
p_t = particle_types{1};
for i = 2 : length(particle_types)
   p_t = [p_t ',' particle_types{i}];
end
p_t = regexprep(p_t, '_', '\\_');

%axis([0 200 min(Low)-5 max(NN)+5])
axis([0 obj.Ripley_cutoff min(Low)-5 max(NN)+5])
title(['Bivariate Ripley Analysis for ', b_n, ' (', p_t, ')'])
legend('data', 'confidence', 'Location', 'Best')
name = fullfile(obj.ResultsDir, [base_name, '_bivripley_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% m       is the number of different types of particles.
% n       is the number of images in the experiment.
% Rate    is the sampling rate.

results.R = R;
results.E = E;
results.A = A;
results.High = High;
results.Low  = Low;

end
