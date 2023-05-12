function SaveDir = mkSMITETmpDir(subdir1, subdir2)
% Create full pathnames and temporary directories for SMITE testing/examples.
% This will produce SaveDir:
%    tempdir/'smite'/subdir1/subdir2
%
% INPUTS:
%    subdir1   subdirectory name, typically 'unitTest' or 'examples'
%    subdir2   lower level subdirectory name that identifies the type of test
%
% OUTPUT:
%    SaveDir   full pathname for temporary save directory
%              (tempdir/'smite'/subdir1/subdir2)

   SaveDir = fullfile(tempdir, 'smite');
   if ~isfolder(SaveDir)
      mkdir(SaveDir)
   end
   SaveDir = fullfile(SaveDir, subdir1);
   if ~isfolder(SaveDir)
      mkdir(SaveDir)
   end
   SaveDir = fullfile(SaveDir, subdir2);
   if ~isfolder(SaveDir)
      mkdir(SaveDir)
   end

end
