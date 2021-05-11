function [CellLabelStruct] = concatenateResults(PublishedResultsStruct)
%concatenateResults concatenates results across multiple datasets.
% This method will take the input PublishedResultsStruct and concatenate
% the results for each Cell/Label pair.  I.e., for a single label on a
% single cell multiple datasets may have been taken, and this method will
% concatenate the info. in PublishedResultsStruct so that the output struct
% CellLabelStruct will contain, for example, one array with the maximum
% cross-correlation for brightfield registration across all datasets for
% that cell and label.  E.g. PublishedResultsStruct.Cell_01_Label_01 may
% itself have 5 fields, each corresponding to a different .h5 file, with
% each of those 5 fields having a field MaxCorr.  The output
% CellLabelStruct will then have a field
% CellLabelStruct.Cell_01_Label_01.MaxCorr, an Nx1 vector which would be
% the concatenation of the 5 vectors found in
% PublishedResultsStruct.Cell_01_Label_01.Data_*.MaxCorr .
%
% INPUTS:
%   PublishedResultsStruct: Structured array containing results computed in
%                           SMA_Publish for each Cell/Label/Dataset
%                           combination.
%
% OUTPUTS:
%   CellLabelStruct: A structured array which will contain the same
%                    information contained in PublishResultsStruct, but
%                    with that information concatenated across all datasets
%                    for a given Cell/Label combination.

% Created by:
%   David J. Schodt (Lidke Lab, 2018)


% Pre-allocate the CellLabelStruct.
CellLabelStruct = struct();

% Loop through each Cell/Label pair field in PublishedResultsStruct.
CellLabelFields = fieldnames(PublishedResultsStruct);
for ii = 1:numel(CellLabelFields)
    % Isolate the ii-th field name.
    CellLabelField = CellLabelFields{ii};
    
    % Initialize our arrays of concatenated data (we can't specify a size
    % because the size of each field in each dataset might change).
    DiffImages = [];
    OverlayImages = [];
    SSIM = [];
    MaxCorr = [];
    MaxCorrFit = [];
    MaxIterReached = [];
    OffsetFitSuccess = [];
    
    % Loop through each dataset in the current CellLabelField.
    DataFields = fieldnames(PublishedResultsStruct.(CellLabelField));
    for jj = 1:numel(DataFields)
        % Isolate the jj-th field name.
        DataField = DataFields{jj};
        
        % Concatenate the needed fields as appropriate.
        DataStruct = PublishedResultsStruct.(CellLabelField).(DataField);
        DiffImages = cat(3, DiffImages, DataStruct.DiffImages);
        OverlayImages = cat(4, OverlayImages, DataStruct.OverlayImages);
        SSIM = [SSIM; DataStruct.SSIM];
        MaxCorr = [MaxCorr; DataStruct.MaxCorr];
        MaxCorrFit = [MaxCorrFit; DataStruct.MaxCorrFit];
        MaxIterReached = [MaxIterReached; DataStruct.MaxIterReached];
        OffsetFitSuccess = [OffsetFitSuccess; ...
            DataStruct.OffsetFitSuccess];
    end
    
    % Store the concatenated arrays in the output structure.
    CellLabelStruct.(CellLabelField).DiffImages = DiffImages;
    CellLabelStruct.(CellLabelField).OverlayImages = OverlayImages;
    CellLabelStruct.(CellLabelField).SSIM = SSIM;
    CellLabelStruct.(CellLabelField).MaxCorr = MaxCorr;
    CellLabelStruct.(CellLabelField).MaxCorrFit = MaxCorrFit;
    CellLabelStruct.(CellLabelField).MaxIterReached = MaxIterReached;
    CellLabelStruct.(CellLabelField).OffsetFitSuccess = ...
        OffsetFitSuccess;
end


end