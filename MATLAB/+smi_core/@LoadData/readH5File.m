function [H5Structure] = readH5File(FilePath, GroupName)
%Extracts contents of an h5 file into H5Structure.
% This method will extract the Data and Attributes from a group
% named GroupName in the .h5 file specified by FilePath.
%
% Examples:
%   H5Structure = readH5File('C:\file.h5') will extract all
%       contents of file.h5 and store them in H5Structure.
%   H5Structure = readH5File('C:\file.h5', 'Laser647') will extract
%       contents of the group 'Laser647' from file.h5 given only
%       the group name.
%   H5Structure = readH5File('C:\file.h5', ...
%       '/Channel01/Zposition001/Laser647') will extract contents
%       of the group 'Laser647' from file.h5 given a full group
%       path.
%
% INPUTS: 
%   FilePath: String containing the path to the .h5 file of interest.
%   GroupName: (optional) Name of a specific group in the .h5 file to be
%              extracted.
%
% OUTPUTS:
%   H5Structure: Structured array containing the information extracted from
%                the .h5 file at FilePath.
%
% REQUIRES:
%   MATLAB 2016b or later
%
% CITATION:

% Created by:
%   David J. Schodt (LidkeLab, 2018)


% Ensure that FilePath points to a valid file.
% NOTE: == 2 means a file indeed exists at FilePath
if (exist(FilePath, 'file') ~= 2)
    error(['File specified by FilePath = ', ...
        '''%s'' does not exist.'], FilePath)
end

% If GroupName was not specified, set a flag to indicate we want to extract 
% all contents from the .h5 file.
SaveAll = ~exist('GroupName', 'var');

% Read in all of the information available about the structure of the h5
% file using MATLAB's built-in h5info() method.
CurrentGroups = h5info(FilePath);

% Determine the .h5 file structure being used (i.e. is each
% dataset in its own group or does one group contain all of
% the datasets).
% NOTE: DataFormat==1 means each dataset is in its own group,
%       DataFormat==0 means each dataset is contained in one
%       "supergroup" of all datasets.
DataFormat = contains(CurrentGroups.Groups.Groups.Groups(1).Name, 'Data');

% Search the .h5 file for the desired groups and store their location
% within the file for later extraction.
DesiredGroups = [];
if DataFormat
    % For .h5 files in the DataFormat, we need to continue searching
    % through the data groups so that we find all instances of the
    % requested group.
    DataGroupsSearched = false;
else
    % For .h5 files that do not have the DataFormat, there is only one data
    % group to search.
    DataGroupsSearched = true;
end
while (isempty(DesiredGroups) || ~DataGroupsSearched)
    % If we are extracting all contents in the .h5 file at FilePath, we 
    % don't want to perform the search procedure in this while loop.  Set 
    % DesiredGroups to CurrentGroups (all of the groups found in the .h5 
    % file) and break out of the loop.
    if SaveAll
        DesiredGroups = CurrentGroups;
        break
    end
    
    % Search the current groups for the group of interest.
    DesiredGroups = [DesiredGroups, ...
        findGroupPaths(CurrentGroups, GroupName)];
    
    % If using the data group format and we are at the .h5 structure data
    % group level, we need to search for the group of interest within the
    % sub-groups of each data group.  Otherwise, we need to move one level
    % deeper into the structure to continue the search.
    CurrentGroupNameFormat = extractGroupName(CurrentGroups(1).Name);
    if contains(CurrentGroupNameFormat, 'Data')
        % We are at the data group level and must search each data 
        % group individually.
        for ii = 1:numel(CurrentGroups)
            DesiredGroups = [DesiredGroups, ...
                findGroupPaths(CurrentGroups(ii).Groups, GroupName)];
        end
        
        % Set the DataGroupsSearched flag to indicate we've searched all of
        % the data groups.
        DataGroupsSearched = true;
    else
        % Move one level deeper into the structure and continue searching
        % for the group of interest.
        try
            % Attempt to move one level deeper into the structure.
            CurrentGroups = CurrentGroups.Groups;
        catch
            % We've reached the end of the structure, there are no more 
            % groups to explore.
            error(['No group named ''%s'' was found in ', ...
                'the file specified by FilePath = ''%s'''], ...
                GroupName, FilePath)
        end
    end
end

% Now that we've found the desired group(s), store the
% Attributes, Data, and Children in a more useable format.
for ii = 1:numel(DesiredGroups)
    % Set DesiredGroup = DesiredGroups(ii) to reduce indexing
    % clutter/improve code readability.
    DesiredGroup = DesiredGroups(ii);
    
    % Determine if the current group is the child of a data group (this 
    % will be needed later).
    DesiredGroupParentName = extractGroupName(DesiredGroup.Name, 2);
    IsDataGroupChild = (DataFormat ...
        && contains(DesiredGroupParentName, 'Data'));
    
    % If the desired group has attributes, store them in the
    % output structure.
    if ~isempty(DesiredGroup.Attributes)
        for jj = 1:numel(DesiredGroup.Attributes)
            % Store each attribute as a field in the output structure 
            % accesible directly by that attributes name.  If the
            % DesiredGroup is a data group, we'll place the attribute
            % information one level deeper into the output structure.
            AttributeName = DesiredGroup.Attributes(jj).Name;
            if IsDataGroupChild
                % For children of a datagroup, we need a different path
                % format within the structure (for consistency with the
                % structure produced for non-datagroup files).
                H5Structure(ii).Attributes.(AttributeName) = ...
                    DesiredGroup.Attributes(jj).Value;
            else
                H5Structure.Attributes.(AttributeName) = ...
                    DesiredGroup.Attributes(jj).Value;
            end
        end
    else
        % Create an empty field for the Attributes to prevent
        % issues with functions that may be using this method, proceeding
        % depending on whether or not we are currently storing a data group
        % (for files that contain separate groups for each dataset).
        if IsDataGroupChild
            % For children of a datagroup, we need a different path
            % format within the structure (for consistency with the
            % structure produced for non-datagroup files).
            H5Structure(ii).Attributes = [];
        else
            H5Structure.Attributes = [];
        end
    end
    
    % If the desired group has data, store the data in the output OldStruct
    % structure.
    if ~isempty(DesiredGroup.Datasets)
        for jj = 1:numel(DesiredGroup.Datasets)
            % Read the ii-th dataset from the h5 file and store the data 
            % in a field of the output structure, matching the datasets 
            % name in the h5 file to the name of the field in the output 
            % structure. If the DesiredGroup is a data group, we'll place 
            % the dataset information one level deeper into the output 
            % structure.
            DatasetName = DesiredGroup.Datasets(jj).Name;
            if IsDataGroupChild
                % For children of a datagroup, we need a different path
                % format within the structure (for consistency with the
                % structure produced for non-datagroup files).
                H5Structure(ii).Data.(DatasetName) = h5read(...
                    FilePath, [DesiredGroup.Name, '/', ...
                    DesiredGroup.Datasets(jj).Name]);
            else
                H5Structure.Data.(DatasetName) = ...
                    h5read(FilePath, [DesiredGroup.Name, '/', ...
                    DesiredGroup.Datasets(jj).Name]);
            end
        end
    else
        % Create an empty field for the Data to prevent issues
        % with functions that may be using this method, proceeding
        % depending on whether or not we are currently storing a data group
        % (for files that contain separate groups for each dataset).
        if IsDataGroupChild
            % For children of a datagroup, we need a different path
            % format within the structure (for consistency with the
            % structure produced for non-datagroup files).
            H5Structure(ii).Data = [];
        else
            H5Structure.Data = [];
        end
    end
    
    % If the desired group has Children (subgroups), recurse
    % through those to store their Attributes, Data, and
    % Children.
    if ~isempty(DesiredGroup.Groups)
        % Create a cell array of subgroup names (if subgroups
        % exist).
        % NOTE: If you are comparing the output structure to
        % the Data, Attributes, and Children format used in the
        % exportState() method, the subgroups are assumed to be
        % the Children.
        SubgroupNames = {DesiredGroup.Groups.Name};
        for jj = 1:numel(SubgroupNames)
            % Iteratively explore subgroups of the desired
            % group to store their attributes and data.
            SubgroupStructure = smi_core.LoadData.readH5File(...
                FilePath, SubgroupNames{jj});
            
            % Remove the path information from the subgroup name, e.g. 
            % /Channel01/Zposition001 will become Zposition001.
            SubgroupName = extractGroupName(SubgroupNames{jj});

            % Store the subgroup structure into the output H5Structure.
            H5Structure.Children.(SubgroupName) = SubgroupStructure;
        end
    else
        % Create an empty field for the Children to prevent
        % issues with functions that may be using this method, proceeding
        % depending on whether or not we are currently storing a data group
        % (for files that contain separate groups for each dataset).
        if IsDataGroupChild
            % For children of a datagroup, we need yet another path
            % format within the structure (for consistency with the
            % structure produced for non-datagroup files).
            H5Structure(ii).Children = [];
        else
            H5Structure.Children = [];
        end
    end
end

end

function [DesiredGroups] = findGroupPaths(CurrentGroups, GroupName)
% This function will create a list of paths to a group with name GroupName
% within the set of groups CurrentGroups.

% Iterate through each of the CurrentGroups to search for GroupName.
DesiredGroups = [];
CurrentGroupNames = {CurrentGroups.Name};
for ii = 1:numel(CurrentGroupNames)
    % If the GroupName isn't provided as a full path, simplify the current 
    % group name to remove group structure e.g. if
    % CurrentGroupNames{ii} = '/Channel01/Zposition001', simplify it to 
    % 'ZPosition001' for comparison to GroupName.
    if (GroupName(1) == '/')
        % If the first character of the input GroupName is '/', assume
        % that a full path was provided to the desired group.
        CurrentGroupName = CurrentGroupNames{ii};
    else
        % The full path was not given, just a raw group name was given.
        CurrentGroupName = extractGroupName(CurrentGroupNames{ii});
    end
    
    % Check if the CurrentGroupName matches the desired GroupName, storing
    % the path to the current group if we have a match.
    if strcmp(CurrentGroupName, GroupName)
        DesiredGroups = [DesiredGroups, CurrentGroups(ii)];
    end
end

end

function [GroupName] = extractGroupName(FullGroupName, PathLength)
% This function will extract the group name from a full group
% path within an .h5 file, e.g. input
% FullGroupName = '/Channel01/Zposition001/Data001'
% will yield GroupName = 'Data001'.  The optional input PathLength
% specifies how many levels from the bottom to extract, e.g. in the
% previous example PathLength would be 1, but if PathLength was 2 we would
% have GroupName = 'Zposition001'.

% Set default parameters if needed.
if ~exist('PathLength', 'var')
    PathLength = 1;
end

% If the last character of FullGroupName is a '/', we should remove it.
if ((FullGroupName(end)=='/') && (numel(FullGroupName)>1))
    FullGroupName = FullGroupName(end-1);
end

% Find the last PathLength indices within FullGroupName corresponding to
% the appearance of a '/' character.
LastSlashIndex = find(FullGroupName=='/', PathLength, 'last');

% Return either the GroupName of interest or the name of the lowest level
% group in FullGroupName as appropriate.
if (numel(LastSlashIndex) > 1)
    GroupName = FullGroupName((LastSlashIndex(1)+1):(LastSlashIndex(2)-1));
else
    GroupName = FullGroupName((LastSlashIndex+1):end);
end


end