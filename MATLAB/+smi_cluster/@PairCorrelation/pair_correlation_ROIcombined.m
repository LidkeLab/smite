function results = pair_correlation_ROIcombined(obj, n_labels, n_ROIs, ...
                                                ROIs, label_num)
% Combine ROIs while performing pair correlation.
% Modified from code originally written by Carolyn Pehlke.
%
% INPUTS:
%    obj           various properties used by the algorithms
%       BaseName      descriptive name for the results files
%       Fig_ext       figure extension
%       Font_props    [{'FontSize', 15, 'FontWeight', 'bold'}]
%       HistBinSize   histogram bin size sometims known as pixel size (nm)
%       Lines         if true, plot lines rather than points for g(r) vs. r
%       ResultsDir    directory to store results
%       Rmax_axis     sets plotting limit if > 0
%    n_labels      number of different labeled particles:
%                  1 -> auto-correlation, 2 -> cross_correlation
%    n_ROIs        number of ROIs to combine
%    ROIs          contains (x, y) coordinates of the two datasets (nm)---see
%                  smi_helpers.ROITools
%    label_num     [OPTIONAL] if n_labels = 1, then whether x1, y1 or x2, y2
%                  should be used for the coordinates contained in ROIs
%
% OUTPUTS:
%    results       structure containing various results from the algorithm

   if n_labels == 2
      corr_type = 'C';   % cross-correlation
   else
      corr_type = 'A';   % auto-correlation
   end

   base_name     = obj.BaseName;
   hist_bin_size = obj.HistBinSize;

   ROI_size = zeros(n_ROIs, 1);
   IM1 = cell(n_ROIs, 1);
   if corr_type == 'C'
      IM2 = cell(n_ROIs, 1);
   end
   for j = 1 : n_ROIs
      x_min = ROIs{j}.ROI(1);
      x_max = ROIs{j}.ROI(2);
      y_min = ROIs{j}.ROI(3);
      y_max = ROIs{j}.ROI(4);
      ROI_size(j) = min(x_max - x_min, y_max - y_min);

      % Compute the number of pixels in x and y.
      imszX = round((x_max - x_min) / hist_bin_size);
      imszY = round((y_max - y_min) / hist_bin_size);
      % Create a blank image.
      im1 = zeros(imszX, imszY);
      if corr_type == 'C'
         im2 = zeros(imszX, imszY);
      end
      % Convert (x, y) coordinates into pixel units.
      if ~exist('label_num', 'var')
         label_num = 1;
      end
      if label_num == 1
         x1 = round((ROIs{j}.X{1} - x_min) / hist_bin_size) + 1;
         y1 = round((ROIs{j}.Y{1} - y_min) / hist_bin_size) + 1;
      elseif label_num == 2
         x1 = round((ROIs{j}.X{2} - x_min) / hist_bin_size) + 1;
         y1 = round((ROIs{j}.Y{2} - y_min) / hist_bin_size) + 1;
      else
         error('Invalid label_num: %d!', label_num);
      end
      if corr_type == 'C'
         x2 = round((ROIs{j}.X{2} - x_min) / hist_bin_size) + 1;
         y2 = round((ROIs{j}.Y{2} - y_min) / hist_bin_size) + 1;
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
      IM1{j} = im1;
      if corr_type == 'C'
         IM2{j} = im2;
      end
   end
   % Establish rmax as half the size of the ROI in pixels.
   rmax = round(min(ROI_size) / (2 * hist_bin_size));

   if corr_type == 'C'
      % Pair crosscorrelation using Veatch method.
      [G, r, g, dg, rmax] = ...
         smi_cluster.PairCorrelation.get_corr(n_ROIs, rmax, IM1, IM2);
   else
      % Pair autocorrelation using Veatch method.
      [G, r, g, dg, rmax] = ...
         smi_cluster.PairCorrelation.get_corr(n_ROIs, rmax, IM1);
   end

   % Convert back to nm
   r_nm = r * hist_bin_size;

   if corr_type == 'C'
      i1 = 0;   i2 = 0;   % intensity sums
      p1 = 0;   p2 = 0;   % pixel sums
      for i = 1 : n_ROIs
         i1 = i1 + sum(sum(IM1{i}));
         i2 = i2 + sum(sum(IM2{i}));
         p1 = p1 + prod(size(IM1{i}));
         p2 = p2 + prod(size(IM2{i}));
      end
      rho1 = i1 / p1;
      rho2 = i2 / p2;
      %rho1 = mean(mean(im1));
      %rho2 = mean(mean(im2));
      rho = sqrt(rho1 * rho2);
      %paircorr = [{[im1, im2]}, {[r', g', dg']}, {rho}, {G}];
   else
      i1 = 0;   % intensity sum
      p1 = 0;   % pixel sum
      for i = 1 : n_ROIs
         i1 = i1 + sum(sum(IM1{i}));
         p1 = p1 + prod(size(IM1{i}));
      end
      rho = i1 / p1;
      %rho = mean(mean(im1));
      %paircorr = [{im1}, {[r', g', dg']}, {rho}, {G}];
   end

   [estimates, errors, model] = ...
      smi_cluster.PairCorrelation.pc_GaussFit(r', g', rmax, rho);
   %  SMA_Cluster.pc_GaussFit(paircorr{2}(:,1), paircorr{2}(:,2), rmax, ...
   %                          paircorr{3});
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
      name = fullfile(obj.ResultsDir, [base_name, '_crosscorrR']);
   else
      name = fullfile(obj.ResultsDir, [base_name, '_autocorrR']);
   end
   if obj.Lines
      plot(r_nm(2:end),g(2:end),'k.-','MarkerSize',20,'LineWidth',2)
      %plot(r_nm(2:end),paircorr{2}(2:end,2),'k.-', ...
      %     'MarkerSize',20,'LineWidth',2)
   else
      plot(r_nm(2:end),g(2:end),'k.','MarkerSize',10,'LineWidth',2)
      %plot(r_nm(2:end),paircorr{2}(2:end,2),'k.','MarkerSize',10,'LineWidth',2)
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
   else
      saveas(gcf, name);
      delete(gcf);
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
