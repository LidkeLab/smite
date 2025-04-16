function H = hopkins_ROIcombined(obj, n_ROIs, RoI)
% Use Hopkins' statistics to test the clustering of a series of ROIs.
%
% INPUTS:
%    obj         various properties used by the algorithms
%       Font_props      [{'FontSize', 15, 'FontWeight', 'bold'}]
%       Fig_ext         figure extension
%       ResultsDir      directory to store results
%       BaseName        name to identify saved plots, which will have various
%                       descriptive suffixes attached
%       Xlim            []       x-axis limits if defined
%       Ylim            []       y-axis limits if defined
%       HopTestPts      number of Hopkins' statistic test points
%       HopTests        number of tests to run to produce a Hopkins' statistic
%    n_ROIs      number of ROIs to combine
%    ROIs        n_ROIs cell array containing the following fields (nm):
%       ROI      ROI limits in the form
%                   [x_min, x_max, y_min, y_max {, z_min, z_max}]
%       X,Y{,Z}  (x, y, z) coordinates of points inside where X, Y, Z are of
%                the form {[n1 x 1]}, that is, referenced by X{1}
% OUTPUTS:
%    H           mean Hopkin's statistic (over ntests) for each ROI

% Originally written by Michael Wester and Stanly Steinberg in 2008; extended
% to 3D and combined ROIs in 2017--2018.

base_name = obj.BaseName;

%test = 5;        % The number of test points.
%ntests = 1000;   % The number of Hopkins statistics to use.
test = obj.HopTestPts;   % The number of test points.
ntests = obj.HopTests;   % The number of Hopkins statistics to use.

% Dimension (2D or 3D)
dim = 2;
if isfield(RoI{1}, 'Z') && numel(RoI{1}.Z) > 0
   dim = 3;
end

% Compute ntests Hopkins' statistics (using test probes) for each ROI.
H  = zeros(1, n_ROIs);
HS = zeros(1, ntests);
for i = 1 : n_ROIs
   if dim == 2
      X = [RoI{i}.X{1}, RoI{i}.Y{1}];
   else
      X = [RoI{i}.X{1}, RoI{i}.Y{1}, RoI{i}.Z{1}];
   end
   ROI = RoI{i}.ROI;
   x_min = ROI(1);
   x_max = ROI(2);
   y_min = ROI(3);
   y_max = ROI(4);
   % hopkinstat assumes the region is [0, x_max - x_min] x [0, y_max - y_min],
   % so shift the coordinates appropriately.
   X(:, 1) = X(:, 1) - x_min;
   X(:, 2) = X(:, 2) - y_min;
   if dim == 2
      for j = 1 : ntests
         HS(j) = smi_cluster.StatisticsClustering.hopkinstat( ...
                    X, x_max - x_min, y_max - y_min, test);
      end
   else
      z_min = ROI(5);
      z_max = ROI(6);
      X(:, 3) = X(:, 3) - z_min;
      for j = 1 : ntests
         HS(j) = smi_cluster.StatisticsClustering.hopkinstat3( ...
                    X, x_max - x_min, y_max - y_min, ...
                                            z_max - z_min, test);
      end
   end
   % H(i) will be the mean of the ntests Hopkins' statistics for ROI i.
   H(i) = mean(HS);
end

if isempty(base_name)
   return;
end

base_text = regexprep(base_name, '_', '\\_');

% Plot a histogram of the findings.
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
hold on
histogram(H, 25);
xlim([0, 1]);
if ~isempty(obj.Xlim)
    xlim(obj.Xlim);
end
title([base_text, ' [all ROIs]']);
xlabel('H (Hopkins statistic)');
ylabel('frequency');
hold off
name = fullfile(obj.ResultsDir, [base_name, '_hopkins_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

% Plot a PDF of the findings.
if ~isempty(obj.Fig_ext)
   figure('Visible', 'off');
else
   figure;
end
axes(obj.Font_props{:})
hold on
% For some reason, Normalization = PDF does not work correctly, but probability
% does, so ...
%histogram(H, 25, 'Normalization', 'PDF');
histogram(H, 25, 'Normalization', 'probability');
xlim([0, 1]);
ylim([0, 1]);
if ~isempty(obj.Xlim)
    xlim(obj.Xlim);
end
if ~isempty(obj.Ylim)
    ylim(obj.Ylim);
end
title([base_text, ' [all ROIs]']);
xlabel('H (Hopkins statistic)');
ylabel('probability');
hold off
name = fullfile(obj.ResultsDir, [base_name, '_hopkins_PDF_RC']);
if ~isempty(obj.Fig_ext)
   print(['-d', obj.Fig_ext], name);
else
   saveas(gcf, name);
   delete(gcf);
end

end
