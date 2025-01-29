function [n_ROIs, ROI, index_ROI] = ...
   get_ROI_GaussIm(obj, X, Y, x_size, y_size, txt, SMD)
% Display data for ROI selection with gaussianImage.
%
% Use the mouse to select ROIs (regions of interest):
%    left click  chooses the center of a fixed size (x_size x y_size) region
%    right click chooses an adjustable rectangular size region
%    key press:
%       backspace or delete   deletes the previous region
%       anything else         terminates selection
%
% INPUTS:
%    X                cell array of the x-coordinates of the labeled points in
%                     the entire image, where X{i} corresponds to the label i
%                     x-coordinates
%    Y                cell array of the y-coordinates of the labeled points in
%                     the entire image, where Y{i} corresponds to the label i
%                     y-coordinates
%    x_size, y_size   box diameters used when clicking the left mouse button
%    txt              text to label the ROI figure
%    SMD              cell array of SMD structures needed for gaussianImage
%                     (GaussIm)
%
% OUTPUTS:
%    n_ROIs           number of ROIs created
%    ROI              cell array of [xmin, xmax, ymin, ymax] for each ROI
%    index_ROI        cell array of labeled point indices in each ROI, where
%                     index_ROI{i}{j} corresponds to label j in ROI i

% Created by
%    Michael J. Wester (2021)

   pixel2nm = obj.Pixel2nm;
   SRzoom   = obj.SRzoom;

   n_labels = numel(X);

   n_ROIs = 0;   ROI = [];   index_ROI = [];

   BS = char(8);   DEL = char(127);

%figure(2); plot(SMD{2}.X * pixel2nm, SMD{2}.Y * pixel2nm, 'k.');
%figure(3); plot(X{1}, Y{1}, 'k.');

   % selected = 0   terminate selection
   %            1   valid button press
   %            2   ignored or region delete
   cm = hot(256);
   cm(1, :) = [0, 0, 0];
   if n_labels == 1
      GaussIm = smi_vis.GenerateImages.gaussianImage(SMD{1}, SRzoom);
      P = prctile(GaussIm(GaussIm > 0), 99.9);
      GaussIm(GaussIm > P) = P;
   else
      GaussIm1 = smi_vis.GenerateImages.gaussianImage(SMD{1}, SRzoom);
      GaussIm2 = smi_vis.GenerateImages.gaussianImage(SMD{2}, SRzoom);
      P = prctile([GaussIm1(GaussIm1 > 0); GaussIm2(GaussIm2 > 0)], 99.9);
      GaussIm1(GaussIm1 > P) = P;
      GaussIm2(GaussIm2 > P) = P;
      colorvec = smi_helpers.colorVector(obj.Color);
      GaussIm = imfuse(GaussIm1, GaussIm2, 'falsecolor', ...
                       'Scaling', 'joint', 'ColorChannels', colorvec);
   end
   h = imshow(GaussIm, cm);

   %for i = 2 : n_labels
   %   j = obj.Order(i);
   %   plot(X{j}, Y{j}, obj.Color(j));
   %end

   hold on
   axis off
   title(txt);
   xlabel('x');
   ylabel('y');

   done = false;
   while ~done
      clickval = waitforbuttonpress;
      if clickval == 0   % if a mouse button was pressed ...
         clickType = get(gcf, 'SelectionType');
         fprintf('%s: ', clickType);
         switch clickType
         case 'normal'   % left button press: draw a fixed rectangle
            selected = 1;
            curpt = get(gca, 'CurrentPoint');
            curpt = curpt ./ SRzoom .* pixel2nm;
            xmin = curpt(1, 1) - x_size/2;
            xmax = curpt(1, 1) + x_size/2;
            ymin = curpt(1, 2) - y_size/2;
            ymax = curpt(1, 2) + y_size/2;
         case 'alt'      % right button press: draw an adjustable rectangle
            selected = 1;
            rect = getrect;
            rect = rect ./ SRzoom .* pixel2nm;
            xmin = rect(1);
            xmax = rect(1) + rect(3);
            ymin = rect(2);
            ymax = rect(2) + rect(4);
         otherwise       % middle button press (extend) and double click (open)
            selected = 2;
         end
      else   % key was pressed
         charChoice = get(gcf, 'CurrentCharacter');
         % If a backspace or a delete, cancel the previous ROI
         if charChoice == BS | charChoice == DEL
            selected = 2;
            if n_ROIs > 0
               delete(r(n_ROIs));
               %delete(p(n_ROIs));
               delete(t(n_ROIs));

               ROI{n_ROIs} = [];
               index_ROI{n_ROIs} = [];

               fprintf('delete ROI %d\n', n_ROIs);
               n_ROIs = n_ROIs - 1;
            end
         else
            selected = 0;
         end
      end
      if selected == 1
         n_ROIs = n_ROIs + 1;
         fprintf('add ROI %d', n_ROIs);

         xmin_p = xmin * SRzoom / pixel2nm;
         xmax_p = xmax * SRzoom / pixel2nm;
         ymin_p = ymin * SRzoom / pixel2nm;
         ymax_p = ymax * SRzoom / pixel2nm;

         r(n_ROIs) = plot([xmin_p, xmin_p, xmax_p, xmax_p, xmin_p], ...
                          [ymin_p, ymax_p, ymax_p, ymin_p, ymin_p], ...
                          'g-', 'LineWidth', 3);

         % Adjust for image display where y increases going downward as opposed
         % to normal plots where y increases going upward.
         tmp  = ymin;
         ymin = ymax;
         ymax = tmp;
         ymin = SMD{1}.YSize * pixel2nm - ymin;
         ymax = SMD{1}.YSize * pixel2nm - ymax;

         ROI{n_ROIs} = [xmin, xmax, ymin, ymax];
         for i = 1 : n_labels
            index_ROI{n_ROIs}{i} = ...
               xmin <= X{i} & X{i} <= xmax & ymin <= Y{i} & Y{i} <= ymax;
         end

         not_empty = true;
         for i = 1 : n_labels
            fprintf(' (label %d: %d points)', ...
                    i, numel(find(index_ROI{n_ROIs}{i} == true)));
            not_empty = not_empty && any(index_ROI{n_ROIs}{i} == true);
         end
         fprintf('\n');
         if not_empty
            %XX = [];   YY = [];
            %for i = 1 : n_labels
            %   XX = [ XX; X{i}(index_ROI{n_ROIs}{i}) ];
            %   YY = [ YY; Y{i}(index_ROI{n_ROIs}{i}) ];
            %end
            %XX = XX * SRzoom / pixel2nm;
            %YY = (SMD.YSize * pixel2nm - YY) * SRzoom / pixel2nm;
            %p(n_ROIs) = plot(XX, YY, 'g.'); % color the selected points green
            t(n_ROIs) = text((xmin_p + xmax_p)/2, (ymin_p + ymax_p)/2, ...
                             int2str(n_ROIs));
            t(n_ROIs).Color = 'green';
            t(n_ROIs).FontWeight = 'bold';
         else
            fprintf('One label has no points in this ROI!  Deleting ...\n');
            delete(r(n_ROIs));

            ROI{n_ROIs} = [];
            index_ROI{n_ROIs} = [];

            fprintf('delete ROI %d\n', n_ROIs);
            n_ROIs = n_ROIs - 1;
         end
      elseif selected == 0
         done = true;
      end
   end
   hold off

   if length(ROI) > n_ROIs
      ROI = ROI(~cellfun('isempty', ROI));
      index_ROI = index_ROI(~cellfun('isempty', index_ROI));
   end

end
