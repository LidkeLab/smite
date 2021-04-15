function genData(obj)
%simData generates output given SMD_True as input and output selection
%
% Generic note: SMD* below are SMD structures with various fields filled in as
% appropriate at that stage.  Model and Data are image stacks (n x n x f),
% where n is the linear size of the image in pixels and f is the total number
% of frames to be generated (f = obj.NDatasets * obj.NFrames).
%
% Typical data flows are
%    produce noisy coordinates:
%       SMD_True -> SMD_True_Labeled -> SMD_Model -> SMD_Model_Noisy
%    produce noisy image stacks
%       SMD_True -> SMD_True_Labeled -> SMD_Model -> Model -> Data
%
% INPUTS:
%    obj               smi_sim.SimSMLM object (see SimSMLM for properties)
%    pattern           one of the patterns below with the given dependencies:
%                         SiemensStar   obj.NWings
%                         kTets         obj.OrderkTet, obj.RadiuskTet
%
% OUTPUTS:
%   SMD_True           true locations of localizations
%   SMD_True_Labeled   obj.LabelingEfficiency applied to SMD_True
%                      localizations, removing localizations that are not
%                      labeled
%   SMD_Model          blinks generated for SMD_True_Labeled localizations
%   SMD_Model_Noisy    SMD_Model with positional and intensity noise added
%   Model              Gaussian blob image stack produced from SMD_Model
%   Data               Model image stack to which Poisson noise has been
%                      applied
   
% Created by
%    Sajjad Khan and Michael J. Wester (Lidkelab 2021)


   % Apply labeling efficiency.
   obj.SMD_Labeled = obj.applyLabelEffic(obj.SMD_True);
   
   % Generate blinks (units are pixels).
   if isempty(obj.StartState)
      error('genBlinks must define StartState!');
   end
   obj.SMD_Model = obj.genBlinks(obj.SMD_Labeled, obj.StartState);
   
     
end
