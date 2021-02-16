function SMD_True_Labeled = applyLabelEffic(obj, SMD_True)
%applyLabelEffic applies labeling efficiency to an existing set of fluorophores.
%
% INPUTS:
%    obj                object of class SimSMLM defined for:
%       LabelingEfficiency   fractional labeling efficiency in range [0 - 1]
%       Verbose              verbosity level
%    SMD_True           an SMD dataset of true localizations
% OUTPUT:
%    SMD_True_Labeled   an SMD dataset in which unlabeled localizations have
%                       been removed

   Labeled = ones(size(SMD_True.X));
   Labeled(rand(size(SMD_True.X)) > obj.LabelingEfficiency) = 0;
   SMD_True_Labeled = SMD_True;
   SMD_True_Labeled.X(~Labeled) = [];
   SMD_True_Labeled.Y(~Labeled) = [];

   if obj.Verbose >= 1
      fprintf('%d%% labeling efficiency: %d -> %d localizations\n', ...
              round(100 * obj.LabelingEfficiency), numel(Labeled),  ...
              sum(Labeled));
   end

end
