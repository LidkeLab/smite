function [n_ROIs, ROI, index_ROI] = get_ROI(obj, X, Y, x_size, y_size, txt)
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
%
% OUTPUTS:
%    n_ROIs           number of ROIs created
%    ROI              cell array of [xmin, xmax, ymin, ymax] for each ROI
%    index_ROI        cell array of labeled point indices in each ROI, where
%                     index_ROI{i}{j} corresponds to label j in ROI i

% Created by
%    Michael J. Wester (2021)

   n_labels = numel(X);

   n_ROIs = 0;   ROI = [];   index_ROI = [];

   BS = char(8);   DEL = char(127);

   % selected = 0   terminate selection
   %            1   valid button press
   %            2   ignored or region delete
   h = figure();
   % An idea by Samantha Schwartz
   %set(h, 'Position', [200, 70, 900, 900*(x_size/y_size)]);
   hold on
   for i = 1 : n_labels
      j = obj.Order(i);
      plot(X{j}, Y{j}, [obj.Color(j), '.'], 'MarkerSize', obj.Msize);
   end
   xlabel('x (nm)');
   ylabel('y (nm)');
   title(txt);
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
            xmin = curpt(1, 1) - x_size/2;
            xmax = curpt(1, 1) + x_size/2;
            ymin = curpt(1, 2) - y_size/2;
            ymax = curpt(1, 2) + y_size/2;
         case 'alt'      % right button press: draw an adjustable rectangle
            selected = 1;
            rect = getrect;
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
         fprintf('add ROI %d\n', n_ROIs);

         r(n_ROIs) = plot([xmin, xmin, xmax, xmax, xmin], ...
                          [ymin, ymax, ymax, ymin, ymin], ...
                          'r-', 'LineWidth', 3);

         ROI{n_ROIs} = [xmin, xmax, ymin, ymax];
         for i = 1 : n_labels
            index_ROI{n_ROIs}{i} = ...
               xmin <= X{i} & X{i} <= xmax & ymin <= Y{i} & Y{i} <= ymax;
         end

         not_empty = true;
         for i = 1 : n_labels
            not_empty = not_empty && any(index_ROI{n_ROIs}{i} == true);
         end
         if not_empty
            %XX = [];   YY = [];
            %for i = 1 : n_labels
            %   XX = [ XX; X{i}(index_ROI{n_ROIs}{i}) ];
            %   YY = [ YY; Y{i}(index_ROI{n_ROIs}{i}) ];
            %end
            %p(n_ROIs) = plot(XX, YY, 'r.');   % color the selected points red
            t(n_ROIs) = text((xmin + xmax)/2, (ymin + ymax)/2, ...
                             int2str(n_ROIs));
            t(n_ROIs).Color = 'black';
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
