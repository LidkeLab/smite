function [Channel01, Version, NDatasets, TopList, InstrList] = ...
         seqH5Data(FilePath)
%seqH5Data parses the HDF5 file structure located in FilePath and extracts
% useful quantities.
%
% INPUT:
%    FilePath       full file path of HDF5 file
%
% OUTPUT:
%    Channel01      top-level HDF5 group under which the data is located
%    Version        HDF5 file version readable by SMITE as deduced by the HDF5
%                   file structure
%    NDatasets      number of datasets
%    TopList        top-level groups, e.g., '/Channel01' for all versions,
%                   '/Calibration' and '/Metadata' for versions >= SEQv2
%    InstrList      list of instruments used in the data collection

% Created by
%    Michael J. Wester (2024)

   try
      h = h5info(FilePath);
   catch
      error('seqH5Data(): Invalid H5 format: %s', FilePath)
   end
   NumTopGroups = numel(h.Groups);
   if NumTopGroups == 0
      error('seqH5Data(): No top-level groups: %s', FilePath)
   elseif NumTopGroups == 1
      TopList = {h.Groups.Name};
      Channel01 = h.Groups(1);
      %if any(contains({ Channel01.Groups.Groups.Name }, 'Data0001'))
      try
         tf = any(contains({ Channel01.Groups.Groups.Name }, 'Data0001'));
      catch ME
         tf = false;
      end
      if tf
         Version = 'SEQv1';
         NDatasets = numel(Channel01.Groups.Groups);
         InstrList = { Channel01.Groups.Groups(1).Groups.Name };
      else
         Version = 'SEQv0';
         NDatasets = numel(Channel01.Groups.Datasets);
         InstrList = { Channel01.Groups(1).Groups.Name };
      end
   else
      TopList = {h.Groups.Name};
      if any(contains(TopList, '/Channel01'))
         for i = 1 : NumTopGroups
            if strcmp(TopList{i}, '/Channel01')
               Channel01 = h.Groups(i);
               break;
            end
         end
         NDatasets = numel(Channel01.Groups.Groups);
         InstrList = { Channel01.Groups.Groups(1).Groups.Name };
      else
         error('seqH5Data(): No Channel01 found: %s', FilePath)
      end
      if any(contains(TopList, '/Metadata'))
         for i = 1 : NumTopGroups
            if strcmp(TopList{i}, '/Metadata')
               Metadata = h.Groups(i);
               Version = Metadata.Attributes.Value;
               break;
            end
         end
      else
         error('seqH5Data(): No Metadata found: %s', FilePath)
      end
   end

end
