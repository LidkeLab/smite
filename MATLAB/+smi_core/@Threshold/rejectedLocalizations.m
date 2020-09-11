function rejectedLocalizations(obj, SMD, options, SaveDir)
% Produce plots of accepted and rejected localization fits,
% individually by % reason rejected and combined by number of reasons rejected
% or by major reason rejected.
%
% INPUT:
%    SMD       SMD (Single Molecule Data) structure.  Fields required:
%          X, Y         x and y coordinates of fitted localizations
%          ThreshFlag   thresholding flag.  Accepted fits have flag = 0, while
%                       rejected fits have a nonzero flag, the value indicating
%                       the reason for rejection.  See setThreshFlag for
%                       further details.
%    options   [OPTIONAL] string indicating which plots to produce or all if
%                         omitted or empty:
%              'R'  rejected due to individual reasons (X,Y,Z; Photons; ...
%              'N'  rejected for the indicated number of reasons
%              'M'  rejected for the indicated major reasons
%    SaveDir   [OPTIONAL] directory in which to save plots produced

%Created by
%    Michael Wester (Lidkelab 2020)

   if ~exist('options', 'var') | isempty(options)
      options = 'RNM';
   end
   Rplot = contains(options, 'R');

   if ~exist('SaveDir', 'var')
       SaveDir = [];
   end

   fprintf('Fits = %d\n', numel(SMD.X));

   % Threshold Flag for the fields or groups of fields named below.
   TF = false(numel(SMD.X), numel(obj.Fields) - 2 - 2);

   % field_bits are the bit numbers corresponding to Fields (above) for the
   % given named fields.

   % Rejection due to X, Y, Z
   field_bits = [strmatch('X', obj.Fields, 'exact'), ...
                 strmatch('Y', obj.Fields, 'exact'), ...
                 strmatch('Z', obj.Fields, 'exact')];
   TF(:, 1) = rejected2DPlot(SMD, field_bits, Rplot, SaveDir, 'X,Y,Z');

   % Rejection due to Photons
   field_bits = strmatch('Photons', obj.Fields, 'exact');
   TF(:, 2) = rejected2DPlot(SMD, field_bits, Rplot, SaveDir, 'Photons');

   % Rejection due to Bg
   field_bits = strmatch('Bg', obj.Fields, 'exact');
   TF(:, 3) = rejected2DPlot(SMD, field_bits, Rplot, SaveDir, 'Bg');

   % Rejection due to PSFSigma
   field_bits = strmatch('PSFSigma',  obj.Fields, 'exact');
   TF(:, 4) = rejected2DPlot(SMD, field_bits, Rplot, SaveDir, ...
                                                     'PSFSigma{,X,Y}');

   % Rejection due to {X,Y,Z}_SE
   field_bits = [strmatch('X_SE', obj.Fields, 'exact'), ...
                 strmatch('Y_SE', obj.Fields, 'exact'), ...
                 strmatch('Z_SE', obj.Fields, 'exact')];
   TF(:, 5) = rejected2DPlot(SMD, field_bits, Rplot, SaveDir, '{X,Y,Z}_SE');

   % Rejection due to PValue
   field_bits = strmatch('PValue', obj.Fields, 'exact');
   TF(:, 6) = rejected2DPlot(SMD, field_bits, Rplot, SaveDir, 'PValue');

   % Combined plot for number of reasons a fit was rejected.
   color = ['b', 'c', 'g', 'y', 'm', 'r'];
   n_colors = numel(color);

   n_reasons = 6;

   fprintf('\n');
   STF = sum(TF, 2);
   if contains(options, 'N')
      figure;
      hold on
      k = 0;
      indx = ~STF;
      % Only include legend entries if the indexing array (indx) has some
      % nonzero entries [tested with sum(indx)].
      if sum(indx) > 0
         fprintf('Accepted = %d\n', sum(indx));
         plot(SMD.X(indx), SMD.Y(indx), 'k.');
         k = k + 1;
         lgd{k} = 'accepted';
      end
      fprintf('Rejected = %d\n', sum(~indx));
      for i = 1 : n_colors - 1
         indx = STF == i;
         if sum(indx) > 0
            fprintf('Rejected for %d reason(s) = %d\n', i, sum(indx));
            plot(SMD.X(indx), SMD.Y(indx), [color(i), '.']);
            k = k + 1;
            lgd{k} = sprintf('%d', i);
         end
      end
      indx = STF >= n_colors;
      if sum(indx) > 0
         fprintf('Rejected for %d+ reasons = %d\n', n_colors, sum(indx));
         plot(SMD.X(indx), SMD.Y(indx), [color(n_colors), '.']);
         k = k + 1;
         lgd{k} = sprintf('%d+', n_colors);
      end
      xlabel('x (pixels)');
      ylabel('y (pixels)');
      reason = 'NumberReasons';
      title(sprintf('Rejected: %s', reason));
      legend(lgd, 'Location', 'Best');
      hold off
      clear lgd
      if ~isempty(SaveDir)
         saveas(gcf, fullfile(SaveDir, reason), 'png');
      end
   end

   if contains(options, 'M')
      % Combined plot for major reasons a fit was rejected.
      figure;
      hold on
      color = [color(n_colors), color(1 : n_colors - 1)];
      k = 0;
      indx = ~STF;
      if sum(indx) > 0
         plot(SMD.X(indx), SMD.Y(indx), 'k.');
         k = k + 1;
         lgd{k} = 'accepted';
      end
      for i = 1 : n_colors
         switch i
         case 1
            indices = [];
            lgd_text = 'multiple';
         case 2
            indices = [1, 5];
            lgd_text = 'X,Y,Z/\_SE';
         case 3
            indices = 2;
            lgd_text = 'Photons';
         case 4
            indices = 3;
            lgd_text = 'Bg';
         case 5
            indices = 4;
            lgd_text = 'PSFSigma';
         case 6
            indices = 6;
            lgd_text = 'PValue';
         %case 7
         %   indices = [3, 7, 9, 10, 12, 13];
         %   lgd_text = 'otherwise';
         end
         % Rejected if any of the specified bits is nonzero for a single
         % (combined) reason or for multiple reasons (STF > 1).  This is not
         % quite consistent as a combined reason may be taken as a multiple
         % reason in some situations, so plot multiple reasons first and let
         % combined reasons overwrite.  Of course, this may completely wipe out
         % multiple reasons while still leaving an entry in the legend.
         if i == 1
            indx = STF > 1;
         else
            indx =  any(TF(:, indices), 2) & ...
                   ~any(TF(:, setdiff(1:n_reasons, indices)), 2);
         end
         if sum(indx > 0)
            plot(SMD.X(indx), SMD.Y(indx), [color(i), '.']);
            k = k + 1;
            lgd{k} = lgd_text;
         end
      end
      xlabel('x (pixels)');
      ylabel('y (pixels)');
      reason = 'MajorReasons';
      title(sprintf('Rejected: %s', reason));
      legend(lgd, 'Location', 'Best');
      hold off
      if ~isempty(SaveDir)
         saveas(gcf, fullfile(SaveDir, reason), 'png');
      end
   end

end % rejectedLocalizations

% -----------------------------------------------------------------------------

function TF = rejected2DPlot(SMD, field_bits, Rplot, SaveDir, reason)
% Make an accept/reject plot for a particular specified reason,
% but only if the number of rejections is nonzero.
%
% INPUT:
%    SMD          SMD (Single Molecule Data) structure.  Fields required:
%          X, Y         x and y coordinates of fitted localizations
%          ThreshFlag   thresholding flag.  Accepted fits have flag = 0, while
%                       rejected fits have a nonzero flag, the value indicating
%                       the reason for rejection.  See setThreshFlag for
%                       further details.
%    field_bits   bit numbers corresponding to combined fields
%    Rplot        if true, produce a plot
%    SaveDir      directory in which to save the plot produced
%    reason       string describing the reason that the fits were rejected for
%                 the given set of field_bits
%
% OUTPUT:
%    TF           combined threshold flag for the set of field_bits.  For
%                 example, if field_bits = [1, 2, 3], corresponding, in turn,
%                 to rejection by X, Y and Z , then TF will be fits
%                 rejected due to the combination of X, Y _or_ Z.

   TFs = arrayfun(@(i) bitget(SMD.ThreshFlag, i), field_bits, ...
                  'UniformOutput', false);
   TF = logical(TFs{1});
   for i = 2 : numel(field_bits)
      TF = TF | TFs{i};
   end

   STF = sum(TF);
   if STF > 0
      fprintf('Rejected due to %s = %d\n', reason, STF);

      if Rplot
         figure;
         hold on
         k = 0;
         if sum(~TF) > 0
            plot(SMD.X(~TF), SMD.Y(~TF), 'k.');
            k = k + 1;
            lgd{k} = 'accepted';
         end
         plot(SMD.X(TF), SMD.Y(TF), 'r.');
         k = k + 1;
         lgd{k} = 'rejected';
         legend(lgd, 'Location', 'Best');
         xlabel('x (pixels)');
         ylabel('y (pixels)');
         treason = regexprep(reason, '([{}_])', '\\$1');
         title(sprintf('Rejected: %s', treason));
         hold off

         if ~isempty(SaveDir)
            saveas(gcf, fullfile(SaveDir, reason), 'png');
         end
      end
   end

end % rejected2DPlot
