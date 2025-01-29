function color_vector = colorVector(color_string)
% Converts a string of one character color names into an ordered RGB vector.
% For example, smi_helpers.colorVector('gm') => [2, 1, 2].
%
% INPUT:
%    color_string   string of ordered one character color names, e.g., 'gm'
% OUTPUT:
%    color_vector   corresponding ordered RGB vector, e.g., [2, 1, 2]

% Written by:
%    Michael J. Wester (2025)

   color_vector = [0, 0, 0];   % black
   n_colors = numel(color_string);
   for i = 1 : n_colors
      switch color_string(i)
      case 'k'   % black
         color_vector = color_vector + i * [0, 0, 0];
      case 'r'   % red
         color_vector = color_vector + i * [1, 0, 0];
      case 'g'   % green
         color_vector = color_vector + i * [0, 1, 0];
      case 'b'   % blue
         color_vector = color_vector + i * [0, 0, 1];
      case 'c'   % cyan
         color_vector = color_vector + i * [0, 1, 1];
      case 'm'   % magenta
         color_vector = color_vector + i * [1, 0, 1];
      case 'y'   % yellow
         color_vector = color_vector + i * [1, 1, 0];
      otherwise
         error('Unrecognized color %c', color_string(i));
      end
   end
   if max(color_vector) > n_colors
      error('color_string ''%s'' produces [%d, %d, %d]', ...
            color_string, color_vector(1), color_vector(2), color_vector(3));
   end

end
