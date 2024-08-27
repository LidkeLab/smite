function P = plotCombined(obj, y, bin_width, x_label, ...
                          legend_labels, x_abbrev, colors, line_type)
% Combined frequency, CDF, PDF, plotSpread, box and bar plots of the arrays y.
%
% INPUTS:
%    obj         various properties used by the algorithms
%       Font_props      [{'FontSize', 15, 'FontWeight', 'bold'}]
%       Line_props      [{'LineWidth', 3}]
%       Fig_ext         figure extension
%       BaseName        name to identify saved plots, which will have various
%                       descriptive suffixes attached 
%       PlotDo          ['fnpcCsSxb'] plots to do:
%                                   'f'   frequency
%                                   'n'   normalized
%                                   'p'   PDF
%                                   'c'   CDF
%                                   'C'   CDF (alternative)
%                                   's'   PlotSpread
%                                   'S'   PlotSpread (bars for mean & median)
%                                   'x'   box
%                                   'b'   bar
%       ShowMM          [1]      red mean and green median (2 only mean, 3 only
%                                median) for PlotSpread plots
%       LinLog          ['plot'] options for CDF2 plots are:
%                                'plot', 'semilogx', 'semilogy', 'loglog'
%       ResultsDir      directory to store results
%       Xlim            []       x-axis limits if defined
%       Ylim            []       y-axis limits if defined
%       LegendTitle     ['']     legend title if specified
%       CSV             [false]  produce a CSV file of the data if true
%    y           cell array of data arrays (need not be the same length)
%    bin_width   bin width for histogram plots
%    x_label     text for the x-label
%    legend_labels   {legend_label1, legend_label2, ...}
%    x_abbrev    abbreviated identifier for plots used in constructing the
%                filename
%    colors      [OPTIONAL] colors for plots     (CDF alternative)
%    line_type   [OPTIONAL] line types for plots (CDF alternative)
%
% OUTPUTS:
%    P           P-value for KS test on y{i} and y{j}
%
% REQUIRES:
%    PlotSpread

% Created by
%    Michael Wester (2019)

   base_name = obj.BaseName;

   n = numel(y);

   % If the data contains an Inf, this will cause havoc, so skip plotting
   if any(cellfun(@(yy) any(isinf(yy)), y))
      warning('%s has Infs in its data---skipping plots!', x_abbrev);
      P = zeros(n);
      return;
   end

   if ~exist('line_type', 'var')
      line_type = cell(1, n);
      for i = 1 : n
         line_type{i} = '-';
      end
   end
   if ~exist('colors', 'var')
      colors = ['b', 'r', 'g', 'k', 'c', 'm'];
      if n > 6
         colors = [colors, colors, colors];
         line_type = {'-', '-', '-', '-', '-', '-',       ...
                      '--', '--', '--', '--', '--', '--', ...
                      ':', ':', ':', ':', ':', ':'};
      end
   end

   base_text = regexprep(base_name, '_', '\\_');
   x_text = regexprep(x_label, '_', '\\_');
   legend_text = {};
   for i = 1 : numel(legend_labels)
      legend_text{i} = regexprep(legend_labels{i}, '_', '\\_');
   end

   if obj.CSV
      produceCSVfile(y, obj.ResultsDir, base_name, x_label, legend_labels, ...
                        x_abbrev);
   end

   % frequency
   if any(obj.PlotDo == 'f')
      figure;
      %axes(obj.Font_props{:})
      hold on
      for i = 1 : n
         h(i) = histogram(y{i});
         if ~isempty(bin_width)
            h(i).BinWidth = bin_width;
         end
      end
      if ~isempty(obj.Xlim)
         xlim(obj.Xlim);
      end
      title(base_text);
      xlabel(x_text);
      ylabel('frequency');
      if ~isempty(legend_text)
         lgd = legend(legend_text, 'Location', 'best');
         if ~isempty(obj.LegendTitle)
            title(lgd, obj.LegendTitle);
         end
      end
      hold off
      name = fullfile(obj.ResultsDir, [base_name, x_abbrev, '_freq']);
      if ~isempty(obj.Fig_ext)
         %saveas(gcf, name, obj.Fig_ext)
         print(gcf, name, ['-d', obj.Fig_ext], '-noui')
      end
      saveas(gcf, name, 'fig');
      close
   end

   % normalized
   if any(obj.PlotDo == 'n')
      figure;
      %axes(obj.Font_props{:})
      hold on
      hedge_min = 0;
      hedge_max = -1;
      for i = 1 : n
         [hcounts, hedges] = histcounts(y{i});
         if hedges(end) > hedge_max;
            hedge_min = hedges(1);
            hedge_max = hedges(end);
            delta = hedges(2) - hedges(1);
         end
      end
      edges = hedge_min : delta : hedge_max;

      for i = 1 : n
         [hcounts, hedges] = histcounts(y{i}, edges);
         %if ~isempty(bin_width)
         %   h(i).BinWidth = bin_width;
         %end
         if max(hcounts) == 0
            counts = hcounts;
         else
            counts = hcounts / max(hcounts);
         end
         h = histogram('BinEdges', edges, 'BinCounts', counts);
      end
      if ~isempty(obj.Xlim)
         xlim(obj.Xlim);
      end
      title(base_text);
      xlabel(x_text);
      ylabel('normalized frequency');
      if ~isempty(legend_text)
         lgd = legend(legend_text, 'Location', 'best');
         if ~isempty(obj.LegendTitle)
            title(lgd, obj.LegendTitle);
         end
      end
      hold off
      name = fullfile(obj.ResultsDir, [base_name, x_abbrev, '_normalized']);
      if ~isempty(obj.Fig_ext)
         print(gcf, name, ['-d', obj.Fig_ext], '-noui')
      end
      saveas(gcf, name, 'fig');
      close
   end

   % PDF
   if any(obj.PlotDo == 'p')
      figure;
      %axes(obj.Font_props{:})
      hold on
      for i = 1 : n
         h(i) = histogram(y{i});
         h(i).Normalization = 'probability';
         if ~isempty(bin_width)
            h(i).BinWidth = bin_width;
         end
      end
      if ~isempty(obj.Xlim)
         xlim(obj.Xlim);
      end
      title(base_text);
      xlabel(x_text);
      ylabel('PDF');
      if ~isempty(legend_text)
         lgd = legend(legend_text, 'Location', 'best');
         if ~isempty(obj.LegendTitle)
            title(lgd, obj.LegendTitle);
         end
      end
      hold off
      name = fullfile(obj.ResultsDir, [base_name, x_abbrev, '_PDF']);
      if ~isempty(obj.Fig_ext)
         print(gcf, name, ['-d', obj.Fig_ext], '-noui')
      end
      try
         saveas(gcf, name, 'fig');
      catch ME
         fprintf('### plotCombined: PROBLEM saving %s ###\n', name);
         fprintf('%s\n', ME.identifier);
         fprintf('%s\n', ME.message);
         return;
      end
      close
   end

   % CDF
   if any(obj.PlotDo == 'c')
      figure;
      %axes(obj.Font_props{:})
      hold on
      for i = 1 : n
         h(i) = histogram(y{i});
         h(i).Normalization = 'cdf';
         if ~isempty(bin_width)
            h(i).BinWidth = bin_width;
         end
      end
      BinLimits(1) = min([h.BinLimits]);
      BinLimits(2) = max([h.BinLimits]);
      for i = 1 : n
         h(i).BinLimits = BinLimits;
      end
      if ~isempty(obj.Xlim)
         xlim(obj.Xlim);
      end
      title(base_text);
      xlabel(x_text);
      ylabel('CDF');
      if ~isempty(legend_text)
         lgd = legend(legend_text, 'Location', 'best');
         if ~isempty(obj.LegendTitle)
            title(lgd, obj.LegendTitle);
         end
      end
      hold off
      %X1 = h1.BinEdges;
      %Y1 = h1.Values;
      name = fullfile(obj.ResultsDir, [base_name, x_abbrev, '_CDF']);
      if ~isempty(obj.Fig_ext)
         print(gcf, name, ['-d', obj.Fig_ext], '-noui')
      end
      saveas(gcf, name, 'fig');
      close
   end

   % CDF (alternative)
   if any(obj.PlotDo == 'C')
      % Eliminate plots when confronted with completely empty data.
      isnull = all(cellfun(@isempty, y));
      figure;
      %axes(obj.Font_props{:})
      %hold on
      if ~isempty(obj.Xlim)
         x_max = obj.Xlim(2);
      elseif ~isnull
         %x_max = max(cellfun(@max, y));
         x_max = 0;
         for i = 1 : n
            if ~isempty(y{i})
               x_max = max(x_max, max(y{i}));
            end
         end
      else
         x_max = 1;
      end
      k = 0;
      if ~isnull
         T = cell(1, n);
         X = cell(1, n);
         Y = cell(1, n);
         for i = 1 : n
            if ~isempty(y{i})
               y_i = y{i};
               if any(isnan(y_i))
                  y_i(isnan(y_i)) = [];
               end
               if ~isempty(y_i)
                  k = k + 1;
                  T{i} = tabulate(y_i);
                  X{i} = [0; T{i}(:, 1); x_max];
                  Y{i} = [0; cumsum(T{i}(:, 2)) / numel(y_i); 1];
                  switch obj.LinLog
                  case 'semilogx'
                     semilogx(X{i}, Y{i}, [colors(i), line_type{i}], ...
                              obj.Line_props{:});
                  case 'semilogy'
                     semilogy(X{i}, Y{i}, [colors(i), line_type{i}], ...
                              obj.Line_props{:});
                  case 'loglog'
                     loglog(X{i}, Y{i}, [colors(i), line_type{i}], ...
                            obj.Line_props{:});
                  otherwise
                     plot(X{i}, Y{i}, [colors(i), line_type{i}], ...
                          obj.Line_props{:});
                  end
                  if k == 1
                     hold on
                  end
               end
            end
         end
      end
      if ~isempty(obj.Xlim)
         xlim(obj.Xlim);
      end
      if ~isempty(obj.Ylim)
         ylim(obj.Ylim);
      end
      title(base_text);
      xlabel(x_text);
      ylabel('CDF');
      if ~isempty(legend_text)
         j = 0;
         leg_text = {};
         for i = 1 : n
            if ~isempty(y{i}) && ~all(isnan(y{i}))
               j = j + 1;
               leg_text{j} = legend_text{i};
            end
         end
         lgd = legend(leg_text, 'Location', 'best');
         if ~isempty(obj.LegendTitle)
            title(lgd, obj.LegendTitle);
         end
      end
      hold off
      name = fullfile(obj.ResultsDir, [base_name, x_abbrev, '_CDF2']);
      if ~isempty(obj.Fig_ext)
         print(gcf, name, ['-d', obj.Fig_ext], '-noui')
      end
      saveas(gcf, name, 'fig');
      close

      if obj.CSV
         y_CDF      = cell(2*n, 1);
         labels_CDF = cell(2*n, 1);
         for i = 1 : n
            j = 2*i - 1;
            y_CDF{j}     = X{i};
            y_CDF{j + 1} = Y{i};
            labels_CDF{j}     = [legend_labels{i}, ':x'];
            labels_CDF{j + 1} = [legend_labels{i}, ':y'];
         end
         produceCSVfile(y_CDF, obj.ResultsDir, base_name, x_label, ...
                        labels_CDF, [x_abbrev, '_CDF2']);
      end
   end

   % PlotSpread
   if any(obj.PlotDo == 's')
      figure;
      %axes(obj.Font_props{:})
      if ~isempty(obj.Xlim)
         x_max = obj.Xlim(2);
      elseif ~all(cellfun(@isempty, y))
         %x_max = max(cellfun(@max, y));
         x_max = 0;
         for i = 1 : n
            if ~isempty(y{i})
               x_max = max(x_max, max(y{i}));
            end
         end
      else
         x_max = 1;
      end
      plotSpread(y, 'showMM', obj.ShowMM, 'xNames', legend_text);
      hold on
      if ~isempty(obj.Xlim)
         ylim(obj.Xlim);
      end
      title(base_text);
      ylabel(x_text);
      if n > 1
         set(gca, 'XTickLabelRotation', 45);
      end
      hold off
      name = fullfile(obj.ResultsDir, [base_name, x_abbrev, '_PS']);
      if ~isempty(obj.Fig_ext)
         print(gcf, name, ['-d', obj.Fig_ext], '-noui')
      end
      saveas(gcf, name, 'fig');
      close
   end

   % PlotSpread (bars for mean & median)
   if any(obj.PlotDo == 'S')
      figure;
      %axes(obj.Font_props{:})
      if ~isempty(obj.Xlim)
         x_max = obj.Xlim(2);
      elseif ~all(cellfun(@isempty, y))
         %x_max = max(cellfun(@max, y));
         x_max = 0;
         for i = 1 : n
            if ~isempty(y{i})
               x_max = max(x_max, max(y{i}));
            end
         end
      else
         x_max = 1;
      end
      plotSpreadMM(y, 'showMM', obj.ShowMM + 10, 'xNames', legend_text);
      hold on
      if ~isempty(obj.Xlim)
         ylim(obj.Xlim);
      end
      if n == 1
         title(sprintf('%s\nmean = %g, median = %g', ...
                       base_text, mean(y{1}), median(y{1})));
      else
         title(base_text);
      end
      ylabel(x_text);
      if n > 1
         set(gca, 'XTickLabelRotation', 45);
      end
      hold off
      name = fullfile(obj.ResultsDir, [base_name, x_abbrev, '_PSMM']);
      if ~isempty(obj.Fig_ext)
         print(gcf, name, ['-d', obj.Fig_ext], '-noui')
      end
      saveas(gcf, name, 'fig');
      close
   end

   % box plot
   if any(obj.PlotDo == 'x')
      Y  = [];
      ID = [];
      for i = 1 : n
         if isrow(y{i})
            Y = [Y, y{i}];
         else
            Y = [Y, y{i}'];
         end
         ID = [ID, repmat(i, 1, numel(y{i}))];
      end
      if ~isempty(legend_text)
         boxplot(Y, char(legend_labels{ID}));
      else
         boxplot(Y);
      end
      hold on
      title(base_text);
      ylabel(x_text);
      if n > 1
         set(gca, 'XTickLabelRotation', 45);
      end
      hold off
      name = fullfile(obj.ResultsDir, [base_name, x_abbrev, '_box']);
      if ~isempty(obj.Fig_ext)
         print(gcf, name, ['-d', obj.Fig_ext], '-noui')
      end
      saveas(gcf, name, 'fig');
      close
   end

   % bar graph (with error spread)
   if any(obj.PlotDo == 'b')
      m = arrayfun(@(i) mean(y{i}), 1:n);
      s = arrayfun(@(i)  std(y{i}), 1:n);
      if ~isempty(legend_text)
         bar(categorical(legend_text, legend_text), m);
      else
         bar( m);
      end
      hold on
      errorbar(1:n, m, s, '.', 'CapSize', 25);
      title(base_text);
      ylabel(['mean ', x_text]);
      if n > 1
         set(gca, 'XTickLabelRotation', 45);
      end
      hold off
      name = fullfile(obj.ResultsDir, [base_name, x_abbrev, '_bar']);
      if ~isempty(obj.Fig_ext)
         print(gcf, name, ['-d', obj.Fig_ext], '-noui')
      end
      saveas(gcf, name, 'fig');
      close
   end

   % KS test
   P = zeros(n);
   for i = 1 : n
      P(i, i) = 1;
      for j = i+1 : n
         if ~isempty(y{i}) && ~isempty(y{j}) ...
            && ~all(isnan(y{i})) && ~all(isnan(y{j}))
            [~, P(i, j)] = kstest2(y{i}, y{j});
         else
            P(i, j) = NaN;
         end
         P(j, i) = P(i, j);
      end
   end

   % Various statistics like mean, stdev, median, mode
   printStats(y, x_label, legend_labels, obj.ResultsDir);
end

function produceCSVfile(y, ResultsDir, base_name, x_label, legend_labels, ...
                           x_abbrev)
% Save the data in .csv format as a standard feature.

   n = numel(y);

   name = fullfile(ResultsDir, [base_name, x_abbrev, '.csv']);
   out = fopen(name, 'w');
   Ny = cellfun(@numel, y);   % number of entries per column
   maxNy = max(Ny);   % max number of entries per column
   % Label each column
   for j = 1 : n
      fprintf(out, '%s', legend_labels{j});
      if j < n
         fprintf(out, ',');
      else
         fprintf(out, '\n');
      end
   end
   % Make a header listing the number of data entries per column.
   for j = 1 : n
      fprintf(out, '%d', Ny(j));
      if j < n
         fprintf(out, ',');
      else
         fprintf(out, '\n');
      end
   end
   for i = 1 : maxNy   % i indexes rows
      for j = 1 : n    % j indexes columns
         if i <= Ny(j)
            fprintf(out, '%f', y{j}(i));
         end
         if j < n
            fprintf(out, ',');
         else
            fprintf(out, '\n');
         end
      end
   end
   fclose(out);

end

function printStats(y, x_label, legend_labels, ResultsDir)
%
% Print various statistics (mean, stdev, median, mode) for each array of y.
%
% INPUTS:
%    y               cell array of data arrays (need not be the same length)
%    x_label         text for the x-label
%    legend_labels   {legend_label1, legend_label2, ...}
%    ResultsDir      directory to store results

   out_s = fopen(fullfile(ResultsDir, 'stats.txt'), 'a');

   n = numel(y);
   fprintf(out_s, '%s: mean +/- stdev, median, mode [interval] =\n', x_label);
   for i = 1 : n
      hedge_min = 0;
      hedge_max = -1;
      [hcounts, hedges] = histcounts(y{i});
      if hedges(end) > hedge_max;
         hedge_min = hedges(1);
         hedge_max = hedges(end);
         delta = hedges(2) - hedges(1);
      end
      edges = hedge_min : delta : hedge_max;
      bin = discretize(y{i}, edges);
      m = mode(bin);
      e = edges([m, m + 1]);

      if ~isempty(legend_labels)
         fprintf(out_s, '%2d: %f +/- %f, %f, [%f, %f] %s\n', ...
                 i, mean(y{i}), std(y{i}), median(y{i}), e, legend_labels{i});
      else
         fprintf(out_s, '%2d: %f +/- %f, %f, [%f, %f]\n', ...
                 i, mean(y{i}), std(y{i}), median(y{i}), e);
      end
   end
   fprintf(out_s, '\n');
   fclose(out_s);

end
