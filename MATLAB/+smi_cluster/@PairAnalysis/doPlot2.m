function doPlot2(n_ROIs, RoI, desc, particles, results_dir, Color, plotting)
% Plot 2D ROIs.
%
% INPUTS:
%    n_ROIs        number of ROIs
%    RoI           data structure containing x and y coordinates per ROI
%    desc          string identifying the analysis
%    particles     string array describing the two particles
%    results_dir   output directory
%    Color         label 1 and label 2 colors on display
%    plotting      true if producing plots
%
% OUTPUTS:
%    Figures *_ROI*_L1+L2   Color(1) L1 and Color(2) L2 localizations

% Created by
%    Michael J. Wester (2022)

   if ~plotting
      return;
   end
   for i = 1 : n_ROIs
      h = figure;
      axes('FontSize', 15, 'FontWeight', 'bold');
      hold on
      plot(RoI{i}.X{1}, RoI{i}.Y{1}, [Color(1), '.']);
      plot(RoI{i}.X{2}, RoI{i}.Y{2}, [Color(2), '.']);
      txt = sprintf('%s ROI%d %s+%s', regexprep(desc, '_', '\\_'), i, ...
                    particles{1}, particles{2});
      title(txt);
      xlabel('x (xm)');
      ylabel('y (xm)');
      hold off
      txt = sprintf('%s_ROI%d_%s+%s', desc, i, particles{1}, particles{2});
      saveas(h, fullfile(results_dir, sprintf('%s.png', txt)));
      close
   end

end
