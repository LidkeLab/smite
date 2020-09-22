classdef LoadData < handle
    % LoadData has methods to load raw microsope dataset from .mat, .ics or
    % .h5 format files.
    %
    % INPUTS:
    %   SMF -
    %   varargin -
    %
    % OUTPUTS:
    %   Data -
    %   SMF -
    %
    % REQUIRES:
    %   MATLAB 2014a or later versions
    %   Dipimage (www.diplib.org)
    %
    % CITATION:
    %   Sandeep Pallikkuth, Lidke Lab, 2020
    
    properties
        SMF
    end
    
    properties (SetAccess = 'protected')
        FileType           % Filetype determined from SMF.FileName
        FullFileName       % Full file path
    end
    
    methods
        function [obj, Data] = loadData(SMF,varargin)
            % INPUT
            %    FileName - full path to datafile, must be mat, ics or h5
            %    varargin - input parameters, different for each file type:
            %       mat: MatVarName, name of matlab variable containing the data
            %           Eg:  [Data] = loadData('FileName.mat',sequence)
            %       ics: no parameters needed
            %       h5: DatasetIdx, index of dataset to be loaded from h5 file
            %           ChannelIdx (optional), index of channel to be loaded, default is 1
            %           Eg:  [Data] = loadData('FileName.h5',1,1)
            % OUTPUT
            %    Data - data loaded from FileName, converted to type single
            
            if isfield(SMF,'FileName')
                obj.FileType = SMF.FileName(end-2:end);
                obj.FullFileName=fullfile(SMF.FileDir,SMF.FileName);
            else
                [filename, pathname]=uigetfile('Y:\*.mat;*.mat;*.ics;*.h5','MultiSelect','on');
                SMF.FileName=filename;
                SMF.FileDir=pathname;
                obj.FullFileName=fullfile(SMF.FileDir,SMF.FileName);
            end
            
            switch obj.FileType
                case 'mat'
                    % loading .mat file/s
                    Data=loadDataMat(obj.FullFileName,varargin);
                case 'ics'
                    % loading images from .ics file
                    Data=loadDataIcs(obj.FullFileName,varargin);
                case '.h5'
                    % loading images from .h5 file
                    Data=loadDataH5(obj.FullFileName,varargin);
            end
        end
        
        function [Data]=loadDataMat(FullFileName,varargin)
            if isempty(varargin)
                error('smi_core.LoadData.loadDataMat:noMatVarName','No Matlab variable name given for mat file, please give as input: Data = smi_core.LoadData.loadDataMat(FileName,MatVarName)')
            end
            MatVarName = varargin{1};
            if iscell(MatVarName)
                MatVarName=MatVarName{1};
            end
            % load data
            tmp = load(FullFileName,MatVarName);
            Data = single(tmp.(MatVarName));
        end
        
        function [Data]=loadDataIcs(FullFileName,varargin)
            tmp = readim(FullFileName);
            Data = single(tmp);
        end
        
        function [Data] = loadDataH5(FullFileName,varargin)
            % INPUT
            %    FullFileName - full path to datafile, must be h5
            %    varargin - input parameters
            %       h5: DatasetIdx, index of dataset to be loaded from h5 file
            %           ChannelIdx (optional), index of channel to be loaded, default is 1
            %           Eg:  [Data] = loadDataH5('FileName.h5',1,1)
            % OUTPUT
            %    Data - data loaded from FileName, converted to type single
            
            if isempty(varargin)
                error('smi_core:loadData:noDatasetIndex','No dataset index given for h5 file, please give as input: Data = smi_core.LoadData.loadData(FileName,DatasetIndex)')
            elseif nargin == 2
                ChannelIdx = 1;
            elseif nargin > 2
                ChannelIdx = varargin{2};
            end
            
            DatasetIdx = varargin{1};
            HD5Info = h5info(FullFileName);
            
            % Define a flag to indicate the .h5 file structure: 0 indicates
            % that all of the data exists in a single group, 1 indicates each
            % dataset exists in its own group.
            DataStructFlag = isempty(HD5Info.Groups.Groups.Datasets);
            
            % setup directory into H5 file
            for ii = 1 : numel(HD5Info.Groups)
                if strcmp(HD5Info.Groups(ii).Name,'/Data')
                    ChannelName = sprintf('Channel%02i',ChannelIdx);
                    DataSetName = sprintf('Data%04i',DatasetIdx);
                    DataName = sprintf('/Data/%s/%s',ChannelName,DataSetName);
                elseif strcmp(HD5Info.Groups(ii).Name,'/Channel01')
                    ChannelName = sprintf('Zposition%03i', ChannelIdx);
                    DataSetName = sprintf('Data%04i', DatasetIdx);
                    
                    % Define the name (including path in .h5 file) to the data.
                    if DataStructFlag % each dataset exists in its own group
                        DataName = sprintf('/Channel01/%s/%s/%s', ...
                            ChannelName, DataSetName, DataSetName);
                    else % all of the datasets exist in a single group
                        DataName = sprintf('/Channel01/%s/%s', ...
                            ChannelName, DataSetName);
                    end
                end
            end
            % check whether channel and dataset exist
            for ii = 1 : numel(HD5Info.Groups)
                if strcmp(HD5Info.Groups(ii).Name,'/Data')
                    DataGroup = HD5Info.Groups(ii);
                    break
                elseif strcmp(HD5Info.Groups(ii).Name,'/Channel01')
                    DataGroup = HD5Info.Groups(ii);
                    break
                end
            end
            % check channel input
            ChannelExists = 0;
            for ii = 1 : numel(DataGroup.Groups)
                if strcmp(DataGroup.Groups(ii).Name,['/Data/' ChannelName])
                    ChannelGroup = DataGroup.Groups(ii);
                    ChannelExists = 1;
                    
                    % Create a list of the datasets in the current group.
                    DatasetList = DataGroup.Groups(ii).Datasets;
                    break
                elseif strcmp(DataGroup.Groups(ii).Name,['/Channel01/' ChannelName])
                    if DataStructFlag % each dataset exists in its own group
                        ChannelGroup = DataGroup.Groups(ii);
                        ChannelExists = 1;
                        
                        % Create a list of the datasets in the current group.
                        for jj = 1:numel(...
                                {DataGroup.Groups(ii).Groups.Datasets})
                            DatasetList(jj) = ...
                                DataGroup.Groups(ii).Groups(jj).Datasets;
                        end
                    else % all of the datasets exist in a single group
                        ChannelGroup = DataGroup.Groups(ii);
                        ChannelExists = 1;
                        
                        % Create a list of the datasets in the current group.
                        DatasetList = DataGroup.Groups(ii).Datasets;
                    end
                    break
                end
            end
            if ~ChannelExists
                error('smi_core:loadData:unknownChannel','h5 file does not contain %s, cannot load data from that channel',ChannelName);
            end
            % check dataset input
            DataSetExists = 0;
            for ii = 1 : numel(DatasetList)
                if strcmp(DatasetList(ii).Name, DataSetName)
                    DataSetExists = 1;
                    break
                end
            end
            if ~DataSetExists
                error('smi_core:loadData:unknownDataset','h5 file does not contain dataset %s in %s, cannot load data',DataSetName,ChannelName);
            end
            % load data
            Data=single(h5read(FileName,DataName));
            otherwise
                % unknown file type
                error('smi_core:loadData:unknownFileType','Unknown type of file, datafile should be of type mat, ics or h5')
        end
    end
end

