function [] = genOverlayResults(obj)
%genOverlayResults generates some two-color overlay results.
% This method will generate several results pertaining to two-color
% overlays of super-resolution data.


% Generate a list of the available overlay data.
ResultsStructDir = fullfile(obj.SaveBaseDir, 'ResultsStructs');
Label1Results = dir(fullfile(ResultsStructDir, '*Label_01*'));
Label1Paths = fullfile(ResultsStructDir, ...
    {Label1Results(~[Label1Results.isdir]).name});
NOverlays = numel(Label1Paths);
if (NOverlays == 0)
    return
end

% Loop through all overlays images and compute the shift between labels.
ImageShift = zeros(NOverlays, 2);
ImageShiftLocal = cell(NOverlays, 3);
AffineTransforms = cell(NOverlays, 1);
SubROIDivisor = 10;
for nn = 1:NOverlays
    % Display a status message in the command line.
    if obj.Verbose
        fprintf(['Publish.performFullAnalysis(): ', ...
            'Computing shift for overlay %i of %i...\n'], ...
            nn, NOverlays);
    end
    
    % Load the results structs into the workspace.
    load(Label1Paths{nn}, 'SMD')
    SMD1 = SMD;
    Label2Path = strrep(Label1Paths{nn}, 'Label_01', 'Label_02');
    if exist(Label2Path, 'file')
        load(Label2Path, 'SMD')
        SMD2 = SMD;
    else
        continue
    end
    
    % Estimate the global shift between the labels.
    ImageShift(nn, :) = smi_core.DriftCorrection.regViaDC(SMD1, SMD2);
    
    % Estimate the local shift bewteen sub-ROIs of the labels.
    [ImageShiftLocal{nn, 1}, ...
        ImageShiftLocal{nn, 2}, ...
        ImageShiftLocal{nn, 3}] = ...
        obj.estimateLocalCoordShifts(SMD1, SMD2, ...
        [1, 1] * SMD1.XSize / SubROIDivisor);
    
    % Compute an affine transform between the coordinates.
    AffineTransforms{nn} = ...
        smi_stat.findCoordAffine([SMD1.X, SMD1.Y], [SMD2.X, SMD2.Y], 50);
end

% Concatenate the max. correlation coefficients from the image in each 
% overlay image (keeping the labels separate) so that we have just one 
% array containing the max. correlation coefficients from all of the
% overlay images produced.  Repeat for the registration errors.
MakeShiftVsCorrPlots = ...
    all(isfield(obj.ResultsStruct, {'RegError', 'MaxCorr'}));
if MakeShiftVsCorrPlots
    ConcatenatedRegError = [{obj.ResultsStruct(:, 1).RegError}.',  ...
        {obj.ResultsStruct(:, 2).RegError}.'];
    ConcatenatedMaxCorr = [{obj.ResultsStruct(:, 1).MaxCorr}.', ...
        {obj.ResultsStruct(:, 2).MaxCorr}.'];
else
    ConcatenatedRegError = [];
    ConcatenatedMaxCorr = [];
end

% Save some of the overlay results.
OverlayInfoStruct.ImageShift = ImageShift;
OverlayInfoStruct.ImageShiftLocal = ImageShiftLocal;
OverlayInfoStruct.RegError = ConcatenatedRegError;
OverlayInfoStruct.MaxCorr = ConcatenatedMaxCorr;
OverlayInfoStruct.PixelSize = obj.SMF.Data.PixelSize;
save(fullfile(obj.SaveBaseDir, 'OverlayInfoStruct.mat'), ...
    'OverlayInfoStruct');
save(fullfile(obj.SaveBaseDir, 'AffineTransforms.mat'), ...
    'AffineTransforms', 'Label1Paths')

% Generate the overlay plots across all cells.
if MakeShiftVsCorrPlots
    obj.makeOverlayPlots(ImageShift, ...
        ConcatenatedRegError, ...
        ConcatenatedMaxCorr, ...
        obj.SMF.Data.PixelSize, ...
        obj.SaveBaseDir)
end


end