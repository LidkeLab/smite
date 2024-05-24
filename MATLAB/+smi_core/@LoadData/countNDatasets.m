function NDatasets = countNDatasets(SMF)
% countNDatasets counts the number of datasets in SMF.Data.FileName .
%
% INPUT:
%
% SMF:  Fields that are structures: 
%    Data:
%       FileName:     File name (cell array of char array)
%       FileDir:      File directory (char array)
%       FileType:     Type of data specified by FileName. If using a custom 
%                     extension, you must set this field manually to the true 
%                     underlying file type (e.g., if using a .mat file saved 
%                     as exFile.spt, set obj.Data.FileType = 'mat')
%                     (char array)(Default set to extension of FileName{1})
%       DataVariable: Name of variable saved in FileName which contains the
%                     raw data. (char array)(Default='sequence')
%
% OUTPUT:
%
% NDatasets           Number of datasets counted

% Created by:
%    Michael J. Wester and David J. Schodt (Lidkelab, 2020)

    switch SMF.Data.FileType
        case 'mat'
            if iscell(SMF.Data.FileName)
                NDatasets = numel(SMF.Data.FileName);
                if NDatasets>1
                    return
                end
            end
            
            if iscell(SMF.Data.DataVariable)
                NDatasets = numel(SMF.Data.DataVariable);
                return
            end
            NDatasets = 1;

        case 'ics'
            NDatasets = numel(SMF.Data.FileName);

        case 'h5'
            FilePath = fullfile(SMF.Data.FileDir,SMF.Data.FileName{1});
            [~, H5Version, NDatasets, ~, ~] = ...
                smi_core.LoadData.seqH5Data(FilePath);
            % NDatasets computed above is sufficient for datasets produced in
            % the 2020s, but leave old code in for now for compatability with
            % older H5 formats possibility not recognized by seqH5Data.
            %return;

            % find number of datasets in h5 file
            HD5Info = h5info(FilePath);
            
            % Define a flag to indicate the .h5 file structure: 
            % 0 indicates that all of the data exists in a single 
            % group, 1 indicates each dataset exists in its own 
            % group.
            if any(strcmp(H5Version, {'SEQv0', 'SEQv1'}))
                DataStructFlag = isempty(HD5Info.Groups.Groups.Datasets);
            end
            % Index of h5 file channel containing data.
            ChannelIdx = 1;

            for ii = 1 : numel(HD5Info.Groups)
                if strcmp(HD5Info.Groups(ii).Name,'/Data')
                    DataGroup = HD5Info.Groups(ii);
                    ChannelName = sprintf('Channel%02i',ChannelIdx);
                    break
                elseif strcmp(HD5Info.Groups(ii).Name,'/Channel01')
                    DataGroup = HD5Info.Groups(ii);
                    ChannelName = sprintf('Zposition%03i',ChannelIdx);
                    break
               end
            end
            if strcmp(H5Version, 'SEQv2') % NEW
                DataStructFlag = isempty(DataGroup.Datasets);
            end
            for ii = 1 : numel(DataGroup.Groups)
                if strcmp(DataGroup.Groups(ii).Name,['/Data/' ChannelName])
                    ChannelGroup = DataGroup.Groups(ii);
                    break
                elseif strcmp(DataGroup.Groups(ii).Name, ...
                             ['/Channel01/' ChannelName])
                    if DataStructFlag 
                        % Each dataset exists in its own group.
                        ChannelGroup = DataGroup.Groups(ii).Groups;
                    else
                        % All of the datasets exist in a single
                        % group.
                        ChannelGroup = DataGroup.Groups(ii);
                    end
                    break
                end
            end
            if DataStructFlag 
                % Each dataset exists in its own group.
                NDatasets = numel(ChannelGroup);
            else
                % All of the datasets exist in a single
                % group.
                NDatasets = numel(ChannelGroup.Datasets);
            end
    end
end
