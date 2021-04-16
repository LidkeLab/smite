function genModel(obj)
%genModel generates SMD_Labeled and SMD_Model from SMD_True.
%
% Typical data flows are
%    produce noisy coordinates:
%       SMD_True -> SMD_True_Labeled -> SMD_Model -> SMD_Model_Noisy
%    produce noisy image stacks
%       SMD_True -> SMD_True_Labeled -> SMD_Model -> Model -> Data
% where
%   SMD_True      true locations of localizations
%   SMD_Labeled   obj.LabelingEfficiency applied to SMD_True localizations,
%                 removing localizations that are not labeled
%   SMD_Model     blinks generated for SMD_True_Labeled localizations applied
%
% INPUTS:
%    obj          smi_sim.SimSMLM object (see SimSMLM for properties)
   
% Created by
%    Sajjad Khan and Michael J. Wester (Lidkelab 2021)

   % Apply labeling efficiency.
   if isempty(obj.LabelingEfficiency)
      error('applyLabelEffic must define LabelingEfficiency!');
   end
   % Transform obj.SMD_True into obj.SMD_Labeled.
   obj.applyLabelEffic();
   
   % Generate blinks (units are pixels).
   if isempty(obj.StartState)
      error('genBlinks must define StartState!');
   end
   % Transform obj.SMD_Labeled into obj.SMD_Model.
   obj.genBlinks(obj.StartState);

end
