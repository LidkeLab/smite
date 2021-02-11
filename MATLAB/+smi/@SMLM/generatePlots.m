function generatePlots(obj, ShowPlots, PlotDo)
%generatePlots creates all histograms and plots for an SMD structure.
%
% INPUT:
%    obj         SMLM object
%       obj.SMD     Single Molecule Data structure
%       obj.SMF     Single Molecule Fitting structure
%    ShowPlots:  Flag for showing plots on the screen (Default = false)
%    PlotDo:     Plots to make chosen from the following list:
%                "Photons"     intensity (estimated photons) histogram
%                "Bg"          background intensity histogram
%                "PSFSigma"    sigma of 2D Gaussian PSF model histogram
%                "PValue"      P-value of fit histogram
%                "X_SE"        standard error in estimated X position histogram
%                "Y_SE"        standard error in estimated Y position histogram
%                "Z_SE"        standard error in estimated Z position histogram
%                "NCombined"   number of connection localizations histogram
%                "DriftX"      cumulative x-drift
%                "DriftY"      cumulative y-drift
%                "DriftZ"      cumulative z-drift
%                "Drift"       estimated 2D or 3D drift
%                "FitFrame"    number of fits per frame
%                "DriftIm"     2D drift image from SR data
%                "GaussIm"     2D Gaussian blob image from SR data
%                "HistIm"      2D histogram image from SR data
%                "CircleIm"    2D Circle image from SR data
%                "CircleImDrift" 2D circle image color coded by time
%                (Default is to make all plots)
%                For example, PlotDo = ["PValue", "FitFrame", DriftIm"]
%                NOTE: plots will only be produced if there is corresponding
%                      data in the SMD structure!
%
% OUTPUT:
%    The figures are saved in .png format in the given SMF.Data.ResultsDir
%
% REQUIRES:
%    Dipimage toolbox (http://www.diplib.org/)

% Created by:
%    Hanieh Mazloom-Farsibaf, Marjolein Meddens Apr 2017 (Keith A. Lidke's lab)
%    Michael J Wester (2020)

if obj.Verbose >= 1
   fprintf('Generating output plots ...\n');
end

if ~exist('ShowPlots', 'var')
   ShowPlots = false;
end

if ~exist('PlotDo', 'var') || isempty(PlotDo)
   PlotDo = ["Photons", "Bg", "PSFSigma", "PValue", "X_SE", "Y_SE", "Z_SE", ...
             "NCombined", "DriftX", "DriftY", "DriftZ", "Drift", "FitFrame",...
             "DriftIm", "GaussIm", "HistIm", "CircleIm", "CircleImDrift"];
end

% PlotSaveDir is the path to the directory for saving plots in .png format
PlotSaveDir = obj.SMF.Data.ResultsDir;
SMD = obj.SMD;

if isempty(SMD.X)
   fprintf('No localization data to plot!\n');
   return;
end

if ismember("Photons", PlotDo)
   %create Photons histogram
   plotAndSaveHist('Photons','Intensity')
end

if ismember("Bg", PlotDo)
   %create Bg histogram
   plotAndSaveHist('Bg','Background')
end

if ismember("PSFSigma", PlotDo)
   plotAndSaveHist('PSFSigma','PSFSigma')
end

if ismember("PValue", PlotDo)
   %create PValue histogram
   plotAndSaveHist('PValue','P value')
end

if ismember("X_SE", PlotDo)
   %create X_SE histogram
   plotAndSaveHist('X_SE','X std error')
end

if ismember("Y_SE", PlotDo)
   %create Y_SE histogram
   plotAndSaveHist('Y_SE','Y std error')
end

if ismember("Z_SE", PlotDo)
   %create Z_SE histogram
   plotAndSaveHist('Z_SE','Z std error')
end

if ismember("NCombined", PlotDo)
   %create Number of Connected localization histogram
   plotAndSaveHist('NCombined','Connected emitters')
end

%cumulative of DriftX, DriftY {, DriftZ} and total drift
if isfield(SMD,'DriftX') && ~isempty(SMD.DriftX) && ...
   isfield(SMD,'DriftY') && ~isempty(SMD.DriftY)
   if ismember("DriftX", PlotDo)
      plotAndSaveCum('DriftX','DriftX')
   end
   if ismember("DriftY", PlotDo)
      plotAndSaveCum('DriftY','DriftY')
   end
   if isfield(SMD,'DriftZ') && ~isempty(SMD.DriftZ) && ismember("DriftZ",PlotDo)
      plotAndSaveCum('DriftZ','DriftZ')
   end 

   % Estimated 2D or 3D drift
   if ismember("Drift", PlotDo)
      DC = smi_core.DriftCorrection;
      DC.PixelSizeZUnit = obj.SMF.DriftCorrection.PixelSizeZUnit;
      figDC = DC.plotDriftCorrection(SMD, 'A');
      FileName = [BaseName '_DriftCorrection'];
      saveas(gcf, fullfile(PlotSaveDir, FileName), 'png');
      if ~ShowPlots; close(gcf); end
   end
end

% BaseName is used to label plot files.
[~,BaseName,~] = fileparts(obj.SMF.Data.FileName{1});

if ismember("FitFrame", PlotDo)
   Nloc_frame = [];
   FitFrame = [];
   % Number of localizations per frame
   for jj=1:max(SMD.DatasetNum)
       for ii=1:max(SMD.FrameNum)
           idx=find(SMD.FrameNum==ii & SMD.DatasetNum==jj);
           Nloc_frame{jj}(ii)=length(idx);
       end
   end
   if length(Nloc_frame)==1
       FitFrame = Nloc_frame{1};
   end
   for ii=1:length(Nloc_frame)-1
       FitFrame=cat(2,Nloc_frame{ii},Nloc_frame{ii+1});
   end

   % plot fits per frame
   figure;
   Frames=1:length(FitFrame);
   plot(Frames,FitFrame);
   xlabel('Frames');
   ylabel('Number of Fits');
   title('Number of fits per frame');
   FileName = [BaseName '_FitsPerFrame.png'];
   saveas(gcf, fullfile(PlotSaveDir, FileName), 'png');

   if ~ShowPlots; close(gcf); end
end

SRImageZoom = 10;

if ismember("DriftIm", PlotDo)
   % Drift image
   [~, DriftImRGB] = smi_vis.GenerateImages.driftImage(SMD, SRImageZoom);
   dipshow(DriftImRGB);
   FileName = [BaseName '_DriftImage'];
   saveas(gcf, fullfile(PlotSaveDir, FileName), 'png');
   if ~ShowPlots; close(gcf); end
end

if ismember("GaussIm", PlotDo)
   % Gaussian image
   [GaussIm] = smi_vis.GenerateImages.gaussianImage(SMD, SRImageZoom);
   dipshow(GaussIm);
   FileName = [BaseName '_GaussImage'];
   saveas(gcf, fullfile(PlotSaveDir, FileName), 'png');
   if ~ShowPlots; close(gcf); end
end

if ismember("HistIm", PlotDo)
   % Histogram image
   [~, HistImRGB] = smi_vis.GenerateImages.histogramImage(SMD, SRImageZoom);
   dipshow(HistImRGB);
   FileName = [BaseName '_HistImage'];
   saveas(gcf, fullfile(PlotSaveDir, FileName), 'png');
   if ~ShowPlots; close(gcf); end
end

if ismember("CircleIm", PlotDo)
    % Generate a circle image.
    [~, CircleImageRGB] = smi_vis.GenerateImages.circleImage(...
        SMD, [], SRImageZoom);
    FileName = [BaseName, '_CircleImage.png'];
    imwrite(CircleImageRGB, fullfile(PlotSaveDir, FileName))
    if ShowPlots
        CircleImFigure = figure();
        CircleImAxes = axes(CircleImFigure);
        imshow(CircleImageRGB, [], 'Parent', CircleImAxes)
    end
end

if ismember("CircleImDrift", PlotDo)
    % Generate a circle image.
    ColorMap = parula(numel(SMD.FrameNum));
    [~, SortIndices] = sort(SMD.FrameNum);
    ColorMap = ColorMap(SortIndices, :);
    [~, CircleImageRGB] = smi_vis.GenerateImages.circleImage(...
        SMD, ColorMap, SRImageZoom);
    FileName = [BaseName, '_CircleImageDrift.png'];
    imwrite(CircleImageRGB, fullfile(PlotSaveDir, FileName))
    if ShowPlots
        CircleImDriftFigure = figure();
        CircleImDriftAxes = axes(CircleImDriftFigure);
        imshow(CircleImageRGB, [], 'Parent', CircleImDriftAxes)
    end
end

    % nested function for plotting and saving histograms
    function plotAndSaveHist(FieldName,HistName)
        if isfield(SMD,FieldName) && ~isempty(SMD.(FieldName))
            Vector_in=SMD.(FieldName);
            FigH = smi_vis.GenerateImages.plotHistogram(Vector_in,HistName);
            [~,BaseName,~] = fileparts(obj.SMF.Data.FileName{1});
            FileName = [BaseName '_' regexprep(HistName,' ','_') '_Hist.png'];
            saveas(FigH,(fullfile(PlotSaveDir,FileName)),'png');
            if ~ShowPlots; close(gcf); end
        end
    end

    % nested function for plotting and saving cumulative Drift
    function plotAndSaveCum(FieldName,CumName)
        if isfield(SMD,FieldName) && ~isempty(SMD.(FieldName))
            FigH = obj.plotCumDrift(SMD,FieldName);
            [~,BaseName,~] = fileparts(obj.SMF.Data.FileName{1});
            FileName = [BaseName '_' CumName '_Cum.png'];
            saveas(FigH,(fullfile(PlotSaveDir,FileName)),'png');
            if ~ShowPlots; close(gcf); end
        end
    end

end
