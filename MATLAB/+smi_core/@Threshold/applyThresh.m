function [SMR] = applyThresh(SMD)
%   applyThresh applies ThreshFlag to perform thresholding on SMD,
%   where SMD can be any appropriate coordinate containing structure of the
%   SMA_SR-class.
%
% INPUT:
%   SMD   object of SMA_SR-class
%
% OUTPUT:
%   SMR   updated object with ThreshFlag applied

%Created by
%   Michael J. Wester (2020) and Sandeep Pallikkuth, Lidke Lab, 2017.

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

   fprintf('   Thresholding: %d -> %d localizations\n', sizeX, numel(SMR.X));

end % applyThresh
