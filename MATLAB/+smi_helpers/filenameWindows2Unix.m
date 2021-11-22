function filenameUnix = filenameWindows2Unix(filenameWindows, prefix)
%filenameWindows2Unix converts the Windows filenameWindows to the Unix filename
% filenameUnix, optionally preceded by prefix.
%
% INPUTS:
%    filenameWindows   Windows filename (character string)
%    prefix            [OPTIONAL] prefix to be prepended to the filename
%
% OUTPUT:
%    filenameUnix      Unix filename (character string)

% Created by
%    Michael Wester (lidkelab 2021)

   % Replace \'s by /'s.
   filename = regexprep(filenameWindows, '\\', '/');
   % Remove leading Windows device specifier if present.
   filename = regexprep(filename, '^[^:]*:', '');
   % Prepend prefix to the Unix filename if supplied.
   if exist('prefix', 'var')
      filename = [prefix, filename];
   end
   filenameUnix = filename;

end
