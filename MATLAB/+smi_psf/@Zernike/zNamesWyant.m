function names = zNamesWyant(ll)
%zNamesWyant Classical names for Zernike indices ll (Wyant ordering).

% Created by
%    Michael Wester, 2017, Lidkelab.

   names = arrayfun(@zNameWyant, ll, 'UniformOutput', false);

end

function name = zNameWyant(l)
% Classical name for Zernike index l (Noll ordering).

   % Wyant ordering aberration names (0-based).  See
   % http://www.telescope-optics.net/zernike_expansion_schemes.htm
   Wnames = {'Piston',                         ...
             'Tilt Horizontal',                ...
             'Tilt Vertical',                  ...
             'Defocus',                        ...
             'Primary Astigmatism Vertical',   ...
             'Primary Astigmatism Oblique',    ...
             'Primary Coma Horizontal',        ...
             'Primary Coma Vertical',          ...
             'Primary Spherical',              ...
             'Trefoil Oblique',                ...
             'Trefoil Vertical',               ...
             'Secondary Astigmatism Vertical', ...
             'Secondary Astigmatism Oblique',  ...
             'Secondary Coma Horizontal',      ...
             'Secondary Coma Vertical',        ...
             'Secondary Spherical Aberration'};

   if l < 0
      error('zNameWyant: Zernike polynomials not defined for l < 0!');
   else
      if l >= numel(Wnames)
         name = sprintf('Zernike_Wyant %d', l);
      else
         name = Wnames{l + 1};
      end
   end

end