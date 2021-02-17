function [SMR] = applyThresh(SMD, Verbose)
%   applyThresh applies ThreshFlag to perform thresholding on SMD,
%   where SMD can be any appropriate smite coordinate containing structure.
%
% INPUT:
%   SMD       smite coordinate containing structure
%   Verbose   [OPTIONAL, Default = 1] verbosity level
%
% OUTPUT:
%   SMR       updated object with ThreshFlag applied

%Created by
%   Michael J. Wester (2020) and Sandeep Pallikkuth, Lidke Lab, 2017.

   if ~exist('Verbose', 'var')
      Verbose = 1;
   end

   SMR = SMD;

   SMD_FNames = fieldnames(SMD);
   sizeX = size(SMD.X, 1);
   for nn = 1:length(SMD_FNames)
      fn = SMD.(SMD_FNames{nn});
      % Apply threshold to compatibly sized arrays ...
      if ~isempty(fn) && size(fn, 1) == sizeX
         SMR.(SMD_FNames{nn}) = ...
            SMD.(SMD_FNames{nn})([SMD.ThreshFlag] == 0);
      end
   end

   if Verbose >= 1
      fprintf('   Thresholding: %d -> %d localizations\n', ...
              sizeX, numel(SMR.X));
   end
   if Verbose >= 3
     THR = smi_core.Threshold;
     THR.rejectedLocalizations(SMD, '');
   end

end % applyThresh
