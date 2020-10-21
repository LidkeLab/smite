function SMF = setSMFDatasetList(SMF)
% setSMFDatasetList sets SMF.Data.DatasetList from DatasetMods and NDatasets.
%
% INPUT:
%
% SMF:  Fields that are structures: 
%    Data:
%       FileName:    File name (cell array of char array)
%       FileDir:     File directory (char array)
%       DatasetMods: Cell array containing datasets to be used/excluded from
%                    analysis (Mods <-> modifiers). This is meant to be the
%                    user-facing lists which define DatasetList, meaning that 
%                    this is what would be set in the GUI. DatasetMods{1} will 
%                    contain an array of the "inclusion" dataset numbers and 
%                    DatasetMods{2} will contain an array of the "exclusion" 
%                    datasets. DatasetList will be set elsewhere (e.g., 
%                    smi_core.LoadData) to include the set 
%                       intersect(intersect(1:NDatasets, DatasetMods{1}), ...
%                                 setdiff(1:NDatasets, DatasetMods{2})) 
%                    unless DatasetMods{1} is empty, in which case the first
%                    parantheses term is dropped. For example, if 
%                    NDatasets = 20, and you only want to analyze datasets 1:5, 
%                    you can set DatasetMods{1} = 1:5. If you further decide to
%                    exclude datsaets 2 and 4, you could set 
%                    DatasetMods{2} = [2, 4]. 
%                    (cell array of int32 arrays)(Default={[]; []})
%
% OUTPUT:
%
% SMF:  Fields that are structures: 
%    Data:
%       DatasetList: List of datasets of the raw data to be analyzed.
%                    (array of int32)(Default=int32([]))

   NDatasets = smi_core.LoadData.countNDatasets(SMF);
   DatasetList = intersect(intersect(1:NDatasets, SMF.DatasetMods{1}), ...
                           setdiff(1:NDatasets, SMF.DatasetMods{2}));
   SMF.DatasetList = DatasetList;

end
