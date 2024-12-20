function results = pair_correlation(obj, SMD1, SMD2)
% Perform pair correlation.
%
% INPUTS:
%    obj             various properties used by the algorithms
%       BaseName     descriptive name for the results files
%       Fig_ext      figure extension
%       Font_props   [{'FontSize', 15, 'FontWeight', 'bold'}]
%       HistBinSize  pixel size (nm)
%       Lines        if true, plot lines rather than points for g(r) vs. r
%       ResultsDir   directory to store results
%       ROI          [x_min, x_max, y_min, y_max] (nm)
%       Rmax_axis    sets plotting limit if > 0
%    SMD1 and SMD2   (x, y) coordinates of the two datasets (nm) in the format
%                    (1) SMD structures: SMD1.X, SMD1.Y, SMD2.X and SMD2.Y,
%                    (2) N x 2 array of coordinates.
%                    SMD2 is optional and if omitted or empty requests
%                    auto-correlation rather than cross-correlation.  Units are
%                    assumed to be nm.
%
% OUTPUTS:
%    results         structure containing various results from the algorithm
%       G               2D pair-correlation function values
%       r               radius values (nm)
%       g               angularly averaged pair correlation function
%       dg              errors on angularly averaged g
%       model           model calculated at estimated value
%                    model results
%       objs_per_domain
%       objs_per_domain_SE
%       sigma_domain
%       sigma_domain_SE
%       localizations_per_obj
%       localizations_per_obj_SE
%       sigma_localization
%       sigma_localization_SE

% Modified from code originally written by Carolyn Pehlke.

   % Determine correlation type.
   if ~exist('SMD2', 'var')
      corr_type = 'A';   % auto-correlation
   else
      if isempty(SMD2)
         corr_type = 'A';   % auto-correlation
      else
         corr_type = 'C';   % cross-correlation
      end
   end

   % Allow multiple input formats:
   if isstruct(SMD1) && isfield(SMD1, 'X') && isfield(SMD1, 'Y')
      XY1 = [SMD1.X, SMD1.Y];
   elseif ismatrix(SMD1) && size(SMD1, 2) == 2
      XY1 = SMD1;
   else
      error('SMD1 is not an SMD structure or an N x 2 matrix!');
   end

   if exist('SMD2', 'var')
      if isstruct(SMD2) && isfield(SMD2, 'X') && isfield(SMD2, 'Y')
         XY2 = [SMD2.X, SMD2.Y];
      elseif ismatrix(SMD2) && size(SMD2, 2) == 2
         XY2 = SMD2;
      else
         error('SMD2 is not an SMD structure or an N x 2 matrix!');
      end
   end

   base_name     = obj.BaseName;
   hist_bin_size = obj.HistBinSize;

   if ~isempty(obj.ROI)
      x_min = obj.ROI(1);
      x_max = obj.ROI(2);
      y_min = obj.ROI(3);
      y_max = obj.ROI(4);
      ROI_size = min(obj.ROI([2, 4]) - obj.ROI([1, 3]));
   else
      if corr_type == 'C'
         x_min = min([XY1(:, 1); XY2(:, 1)]);
         x_max = max([XY1(:, 1); XY2(:, 1)]);
         y_min = min([XY1(:, 2); XY2(:, 2)]);
         y_max = max([XY1(:, 2); XY2(:, 2)]);
      else
         x_min = min(XY1(:, 1));
         x_max = max(XY1(:, 1));
         y_min = min(XY1(:, 2));
         y_max = max(XY1(:, 2));
      end
      ROI_size = min(x_max - x_min, y_max - y_min);
   end

   % Set to 1 to get a figure of g(r) with error bars.
   flag = 0;

   % Compute the number of pixels in x and y.
   imszX = round((x_max - x_min) / hist_bin_size);
   imszY = round((y_max - y_min) / hist_bin_size);
   % Create a blank image.
   im1 = zeros(imszX, imszY);
   if corr_type == 'C'
      im2 = zeros(imszX, imszY);
   end
   % Convert (x, y) coordinates into pixel units.
   x1 = round((XY1(:, 1) - x_min) / hist_bin_size) + 1;
   y1 = round((XY1(:, 2) - y_min) / hist_bin_size) + 1;
   if corr_type == 'C'
      x2 = round((XY2(:, 1) - x_min) / hist_bin_size) + 1;
      y2 = round((XY2(:, 2) - y_min) / hist_bin_size) + 1;
   end
   % Get the pixels within the image size.
   mask1 = (x1 > 0) & (y1 > 0) & (x1 <= imszX) & (y1 <= imszY);
   x1 = x1(mask1);
   y1 = y1(mask1);
   if corr_type == 'C'
      mask2 = (x2 > 0) & (y2 > 0) & (x2 <= imszX) & (y2 <= imszY);
      x2 = x2(mask2);
      y2 = y2(mask2);
   end
   % Make a histogram image.
   for i = 1 : size(x1, 1)
      im1(x1(i), y1(i)) = im1(x1(i), y1(i)) + 1;
   end
   if corr_type == 'C'
      for i = 1 : size(x2, 1)
         im2(x2(i), y2(i)) = im2(x2(i), y2(i)) + 1;
      end
   end
   % Establish rmax as half the size of the ROI in pixels.
   rmax = round(ROI_size / (2 * hist_bin_size));
   if corr_type == 'C'
      % Pair crosscorrelation using Veatch method.
      [G, r, g, dg, ~, rmax] = ...
         smi_cluster.PairCorrelation.get_crosscorr(im1, im2, ...
            ones(size(im1)), rmax, flag);
   else
      % Pair autocorrelation using Veatch method.
      [G, r, g, dg, ~, rmax] = ...
         smi_cluster.PairCorrelation.get_autocorr(im1, ones(size(im1)), ...
            rmax, flag);
   end
   % Convert back to nm
   r_nm = r * hist_bin_size;

   if corr_type == 'C'
      rho1 = mean(mean(im1));
      rho2 = mean(mean(im2));
      rho = sqrt(rho1 * rho2);
      paircorr = [{[im1, im2]}, {[r', g', dg']}, {rho}, {G}];
   else
      rho = mean(mean(im1));
      paircorr = [{im1}, {[r', g', dg']}, {rho}, {G}];
   end

   [estimates, errors, model] = ...
      smi_cluster.PairCorrelation.pc_GaussFit(paircorr{2}(:,1), ...
         paircorr{2}(:,2), rmax*obj.RmaxFitFactor, paircorr{3});
%        paircorr{2}(:,2), rmax, paircorr{3});
   estimates = abs(estimates);

   if ~isempty(obj.Fig_ext)
      figure('Visible', 'off');
   else
      figure;
   end
   axes(obj.Font_props{:});
   hold on
   if ~isdir(obj.ResultsDir)
      mkdir(obj.ResultsDir);
   end
   if corr_type == 'C'
      name = fullfile(obj.ResultsDir, [base_name, '_crosscorr']);
   else
      name = fullfile(obj.ResultsDir, [base_name, '_autocorr']);
   end
   if obj.Lines
      plot(r_nm(2:end),paircorr{2}(2:end,2),'k.-', ...
           'MarkerSize',20,'LineWidth',2)
   else
      plot(r_nm(2:end),paircorr{2}(2:end,2),'k.','MarkerSize',10,'LineWidth',2)
   end
   plot(r_nm(2:end),ones(1, numel(r) - 1),'b:','LineWidth',3)
   plot(r_nm(2:end),model(2:end),'--r','LineWidth',3)
   if corr_type == 'C'
      flegend{1} = 'Cross-Correlation';
   else
      flegend{1} = 'Auto-Correlation';
   end
   flegend{2} = 'g(r) Random';
   flegend{3} = 'Fit';
   axis tight
   if obj.Rmax_axis > 0
      xlim([0, obj.Rmax_axis]);
   end
   title(regexprep(base_name, '_', '\\_'));
   xlabel('r (nm)');
   ylabel('g(r)');
   legend(flegend, 'Location', 'Best');
   hold off

   if ~isempty(obj.Fig_ext)
      print(['-d', obj.Fig_ext], name);
      saveas(gcf, name);
   %else
   %   saveas(gcf, name);
   %   delete(gcf);
   end

   Afound=estimates(1);
   s_d_found_pixels=estimates(2);
   Bfound=estimates(3);
   s_l_found_pixels=estimates(4);

   Afound_SE=errors(1);
   s_d_found_pixels_SE=errors(2);
   Bfound_SE=errors(3);
   s_l_found_pixels_SE=errors(4);

   s_d_found=s_d_found_pixels*hist_bin_size;
   s_l_found=s_l_found_pixels*hist_bin_size;
   s_d_found_SE=s_d_found_pixels_SE*hist_bin_size;
   s_l_found_SE=s_l_found_pixels_SE*hist_bin_size;

   % display results
   if corr_type == 'C'
      fprintf('Pair Cross-correlation for %s:\n\n', base_name);
   else
      fprintf('Pair Auto-correlation for %s:\n\n', base_name);
   end

   fprintf('Objects per domain ');
   fprintf('found:\t%5.3g +/- %5.3g\n',Afound,Afound_SE);

   fprintf('Domain size sigma ');
   fprintf('found:\t%5.3g +/- %5.3g\n',s_d_found,s_d_found_SE);

   fprintf('Localizations per object ');
   fprintf('found:\t%5.3g +/- %5.3g\n',Bfound,Bfound_SE);

   fprintf('Localization precision sigma ');
   fprintf('found:\t%5.3g +/- %5.3g\n',s_l_found,s_l_found_SE);

   results.G  = G;
   results.r  = r_nm;
   results.g  = g;
   results.dg = dg;
   results.objs_per_domain          = Afound;
   results.objs_per_domain_SE       = Afound_SE;
   results.sigma_domain             = s_d_found;
   results.sigma_domain_SE          = s_d_found_SE;
   results.localizations_per_obj    = Bfound;
   results.localizations_per_obj_SE = Bfound_SE;
   results.sigma_localization       = s_l_found;
   results.sigma_localization_SE    = s_l_found_SE;
   results.model = model;

end
