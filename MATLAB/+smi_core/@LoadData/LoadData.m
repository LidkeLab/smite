classdef LoadData < handle
    % LoadData has methods to load raw microsope dataset from .mat, .ics or
    % .h5 format files.
    %
    % INPUTS:
    %   SMF - Single Molecule Fitting structure. The name and directory
    %   location of the data file is obtained from fields
    %   (SMF.Data.FileName and SMF.Data.FileDir) of SMF. If the fields are
    %   non-existent or are empty, data file selection will be required
    %   with a pop-up window.
    %   varargin - input arguments helping to locate the 'data' to be
    %   loaded from the data file. This varies for the different file
    %   extensions (.mat, .ics and .h5). For details please see
    %   documentation for the methods.
    %
    % OUTPUTS:
    %   Data - data loaded from FileName, converted to type single
    %   SMF - SMF structure modified by adding fields FileName and FileDir,
    %   if not present.
    %
    % METHODS:
    %   loadDataMat, loadDataIcs, loadDataH5
    %
    % CITATION:
    %   Sandeep Pallikkuth, Lidke Lab, 2020
    
    properties (SetAccess = 'protected')
        FileType           % Filetype determined from SMF.Data.FileName
        FullFileName       % Full file path
    end
    
    methods
        function obj = LoadData()
        % Constructor.
        end

        function [obj,Data,SMF] = loadRawData(obj,SMF,varargin)
            % INPUT
            %    SMF - single molecule fitting structure
            %    varargin - input parameters, different for each file type:
            %       mat: MatVarName, name of matlab variable containing the data
            %            DatasetNum, number indicating the file in
            %            SMF.Data.FileName cell array
            %       ics: no parameters needed
            %       h5: DatasetIdx, index of dataset to be loaded from h5 file
            %           ChannelIdx (optional), index of channel to be loaded, default is 1
            % OUTPUT
            %    Data - data loaded from FileName, converted to type single
            %    SMF - single molecule fitting structure with FileName,
            %    FileDir and ResultsDir added.
            
            if isfield(SMF.Data,'FileName')&&~isempty(SMF.Data.FileName)
                % determining FileTye from SMF entries of FileName
                if iscell(SMF.Data.FileName)
                    obj.FileType=SMF.Data.FileName{1}(end-2:end);
                else
                    obj.FileType = SMF.Data.FileName(end-2:end);
                end
                obj.FullFileName=fullfile(SMF.Data.FileDir,SMF.Data.FileName);
            else
                % If Filename entry not present in SMF, option to choose
                % the file/s
                [filename, pathname]=uigetfile('Y:\*.mat;*.mat;*.ics;*.h5','MultiSelect','on');
                SMF.Data.FileName=filename;
                SMF.Data.FileDir=pathname;
                SMF.Data.ResultsDir=fullfile(pathname,'Results');
                obj.FullFileName=fullfile(SMF.Data.FileDir,SMF.Data.FileName);
                if iscell(SMF.Data.FileName)
                    obj.FileType=SMF.Data.FileName{1}(end-2:end);
                else
                    obj.FileType = SMF.Data.FileName(end-2:end);
                end
            end
            
            switch obj.FileType
                case 'mat'
                    % loading images from .mat file/s
                    Data=smi_core.LoadData.loadDataMat(obj.FullFileName,varargin);
                case 'ics'
                    % loading images from .ics file
                    Data=smi_core.LoadData.loadDataIcs(obj.FullFileName,varargin);
                case '.h5'
                    % loading images from .h5 file
                    Data=smi_core.LoadData.loadDataH5(obj.FullFileName,varargin);
            end
            
            % If SMF.Data.DataROI is empty, set a default.
            if isempty(SMF.Data.DataROI)
                % For now, I'm doing this instead of 
                % [1, 1, size(Data, [1, 2])] because size(Data, [1, 2])
                % didn't work until MATLAB 2019b.
                DataSize = size(Data);
                SMF.Data.DataROI = [1, 1, DataSize(1:2)];
            end
        end
    end
    
    methods (Static)
        function [Data]=loadDataMat(FullFileName,varargin)
            % static method for loading .mat file/s
            % INPUT
            %    FullFileName - full path to datafile, must be mat
            %    varargin - input parameters
            %           MatVarName, name of matlab variable containing the data
            %            DatasetNum, number indicating the file in
            %            SMF.Data.FileName cell array
            %           Eg:  [Data] = smi_core.LoadData.loadDataMat('FullFileName.mat',sequence,DatasetNum)
            % OUTPUT
            %    Data - data loaded from FileName, converted to type single
            if isempty(varargin)
                error('smi_core.LoadData.loadDataMat:','Not enough input arguments: Data = smi_core.LoadData.loadDataMat(FileName,MatVarName,DatasetNumber)');
            elseif nargin<2
                error('smi_core.LoadData.loadDataMat:','Not enough input arguments: Data = smi_core.LoadData.loadDataMat(FileName,MatVarName,DatasetNumber)');
            end
            while any(cellfun(@iscell,varargin))
                varargin = [varargin{cellfun(@iscell,varargin)} varargin(~cellfun(@iscell,varargin))];
            end

            MatVarName = varargin{1};
            if iscell(MatVarName)
                MatVarName=MatVarName{1};
            end
            DatasetNum = varargin{2};
            % load data
            if iscell(FullFileName)
                tmp = load(FullFileName{DatasetNum},MatVarName);
                Data = single(tmp.(MatVarName));
            else
                tmp = load(FullFileName,MatVarName);
                Data = single(tmp.(MatVarName));
            end
                
        end
        
        function [Data]=loadDataIcs(FullFileName,varargin)
            % static method for loading .ics file
            % no input arguments needed
            % INPUT
            %    FullFileName - full path to datafile, must be .ics
            %    Eg:  [Data] = smi_core.LoadData.loadDataIcs('FullFileName.ics')
            % OUTPUT
            %    Data - data loaded from FileName, converted to type single
            tmp = readim(FullFileName);
            Data = single(tmp);
        end
        
        function [Data] = loadDataH5(FullFileName,varargin)
            % static method for loading .h5 file
            % INPUT
            %    FullFileName - full path to datafile, must be .h5
            %    varargin - input parameters
            %       h5: DatasetIdx, index of dataset to be loaded from h5 file
            %           ChannelIdx (optional), index of channel to be loaded, default is 1
            %           Eg:  [Data] = smi_core.LoadData.loadDataH5('FileName.h5',1,1)
            % OUTPUT
            %    Data - data loaded from FileName, converted to type single
            
            if isempty(varargin)
                error('smi_core:LoadData:noDatasetIndex','No dataset index given for h5 file, please give as input: Data = smi_core.LoadData.loadDataH5(FileName,DatasetIndex)')
            elseif nargin == 2
                ChannelIdx = 1;
            elseif nargin > 2
                ChannelIdx = varargin{2};
            end
            while any(cellfun(@iscell,varargin))
                varargin = [varargin{cellfun(@iscell,varargin)} varargin(~cellfun(@iscell,varargin))];
            end
            
            DatasetIdx = varargin{1};
            
            % FullFileName may be a cell array.
            if iscell(FullFileName)
                FullFileName = FullFileName{1};
            end
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
            Data=single(h5read(FullFileName,DataName));
        end

        DatasetList = setSMFDatasetList(SMF);
        NDatasets = countNDatasets(SMF);
        
        function [Gain, Offset, ReadNoise, CalibrationROI] = ...
                loadDataCalibration(SMF)
            %loadDataCalibration loads calibration data from a file.
            % This method attempts to load the gain, offset, and readnoise
            % arrays from a file specified by SMF.Data.CalibrationFilePath. 
            % If these arrays aren't found in the file, a warning will be 
            % issued.
            
            % Load a "Params" struct from the file.
            load(SMF.Data.CalibrationFilePath, 'Params')
            
            % Attempt to extract the appropriate fields from 'Params'.
            Gain = Params.Gain;
            Offset = Params.CCDOffset;
            ReadNoise = Params.CCDVar;
            
            % If the ROI is defined as [YStart, YEnd, XStart, XEnd],
            % we'll want to swap the values to have
            % [YStart, XStart, YEnd, XEnd].
            CalibrationROI = Params.CameraObj.ROI;
            if ((CalibrationROI(2)>CalibrationROI(1)) ...
                    && (CalibrationROI(4)>CalibrationROI(3)))
                % If either of these are met, we'll assume it's defined
                % in the old style. Note that this might fail, but I
                % haven't thought of a better approach -D.J.S. 2/21
                OldROI = CalibrationROI;
                CalibrationROI = ...
                    [CalibrationROI(1), CalibrationROI(3), ...
                    CalibrationROI(2), CalibrationROI(4)];
                warning(['The CalibrationROI in %s was assumed to be ', ...
                    'in the format [YStart, YEnd, XStart, XEnd] and ', ...
                    'thus was re-ordered from ', ...
                    '[%i, %i, %i, %i] to [%i, %i, %i, %i]!'], ...
                    SMF.Data.CalibrationFilePath, ...
                    OldROI(1), OldROI(2), OldROI(3), OldROI(4), ...
                    CalibrationROI(1), CalibrationROI(2), ...
                    CalibrationROI(3), CalibrationROI(4))
            end
        end
    
    end
    
end
