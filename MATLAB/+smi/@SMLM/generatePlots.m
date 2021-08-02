function generatePlots(obj, PlotSaveDir1, PlotSaveDir2, ShowPlots, PlotDo)
%generatePlots creates all histograms and plots for an SMD structure.
%
% INPUT:
%    obj          SMLM object
%       obj.SMD      Single Molecule Data structure
%       obj.SMF      Single Molecule Fitting structure
%       obj.SRImageZoom    magnification factor for SR     images
%       obj.SRCircImZoom   magnification factor for circle images
%    PlotSaveDir1 Directory in which to save especially useful (priority 1)
%                 plots, like GaussIm
%    PlotSaveDir2 Directory in which to save all the other (priority 2) plots
%                 (typically, a subdirectory of PlotSaveDir1)
%    ShowPlots:   Flag for showing plots on the screen (Default = false)
%    PlotDo:      Plots to make chosen from the following list:
%                 "Photons"    intensity (estimated photons) histogram
%                 "Bg"         background intensity histogram
%                 "PSFSigma"   sigma of 2D Gaussian PSF model histogram
%                 "PValue"     P-value of fit histogram
%                 "X_SE"       standard error in estimated X position histogram
%                 "Y_SE"       standard error in estimated Y position histogram
%                 "Z_SE"       standard error in estimated Z position histogram
%                 "NCombined"  number of connection localizations histogram
%                 "DriftX"     cumulative x-drift
%                 "DriftY"     cumulative y-drift
%                 "DriftZ"     cumulative z-drift
%                 "CumDrift"   estimated 2D or 3D cumulative drift
%                 "Drift"      estimated 2D or 3D absolute drift
%                 "FitFrame"   number of fits per frame
%                 "DriftIm"    2D drift image from SR data
%                 "GaussIm"    2D Gaussian blob image from SR data
%                 "HistIm"     2D histogram image from SR data
%                 "CircleIm"   2D Circle image from SR data
%                 "CircleImDrift" 2D circle image color coded by time
%                 (Default is to make all plots)
%                 For example, PlotDo = ["PValue", "FitFrame", DriftIm"]
%                 NOTE: plots will only be produced if there is corresponding
%                      data in the SMD structure!
%
% OUTPUT:
%    The figures are saved in .png format in PlotSaveDir1/2.
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
   % Removed: DriftX, DriftY, DritfZ, HistIm
   PlotDo = ["Photons", "Bg", "PSFSigma", "PValue", "X_SE", "Y_SE", "Z_SE", ...
             "NCombined", "CumDrift", "FitFrame", "DriftIm",                ...
             "GaussIm", "CircleIm", "CircleImDrift"];
end

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
   %plotAndSaveHist('PValue','P value')
   plotAndSavePValueHist('PValue','P value')
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

   % Estimated 2D or 3D cumulative drift
   if ismember("CumDrift", PlotDo)
      DC = smi_core.DriftCorrection;
      DC.PixelSizeZUnit = obj.SMF.DriftCorrection.PixelSizeZUnit;
      figDC = DC.plotDriftCorrection(SMD, 'R');
      if ~isempty(PlotSaveDir2)
         FileName = [BaseName '_CumDriftCorrection'];
         saveas(gcf, fullfile(PlotSaveDir2, FileName), 'png');
      end
      if ~ShowPlots; close(gcf); end
   end

   % Estimated 2D or 3D drift
   if ismember("Drift", PlotDo)
      DC = smi_core.DriftCorrection;
      DC.PixelSizeZUnit = obj.SMF.DriftCorrection.PixelSizeZUnit;
      figDC = DC.plotDriftCorrection(SMD, 'A');
      if ~isempty(PlotSaveDir2)
         FileName = [BaseName '_DriftCorrection'];
         saveas(gcf, fullfile(PlotSaveDir2, FileName), 'png');
      end
      if ~ShowPlots; close(gcf); end
   end
end

% BaseName is used to label plot files.
[~,BaseName,~] = fileparts(obj.SMF.Data.FileName{1});

if ismember("FitFrame", PlotDo)
   Nloc_frame = [];
   % Number of localizations per frame
   for jj=1:max(SMD.DatasetNum)
       for ii=1:max(SMD.FrameNum)
           idx=find(SMD.FrameNum==ii & SMD.DatasetNum==jj);
           Nloc_frame{jj}(ii)=length(idx);
       end
   end
   FitFrame = Nloc_frame{1};
   for ii=1:length(Nloc_frame)-1
       FitFrame=cat(2, FitFrame, Nloc_frame{ii+1});
   end

   % plot fits per frame
   figure;
   Frames=1:length(FitFrame);
   plot(Frames,FitFrame);
   xlabel('Frames');
   ylabel('Number of Fits');
   title('Number of fits per frame');
   if ~isempty(PlotSaveDir2)
      FileName = [BaseName '_FitsPerFrame.png'];
      saveas(gcf, fullfile(PlotSaveDir2, FileName), 'png');
   end

   if ~ShowPlots; close(gcf); end
end

if ismember("DriftIm", PlotDo)
   % Drift image
   [~, DriftImRGB] = smi_vis.GenerateImages.driftImage(SMD, obj.SRImageZoom);
   if ~isempty(PlotSaveDir2)
      FileName = [BaseName, '_DriftImage.png'];
      imwrite(single(DriftImRGB), fullfile(PlotSaveDir2, FileName))
   end
   if ShowPlots
      DriftImFigure = figure();
      DriftImAxes = axes(DriftImFigure);
      imshow(DriftImRGB, [], 'Parent', DriftImAxes)
   end
end

if ismember("GaussIm", PlotDo)
   % Gaussian image
   [GaussIm] = smi_vis.GenerateImages.gaussianImage(SMD, obj.SRImageZoom);
   if ~isempty(PlotSaveDir1)
      FileName = [BaseName, '_GaussImage.png'];
      imwrite(GaussIm, fullfile(PlotSaveDir1, FileName))
   end
   if ~isempty(PlotSaveDir2)
      FileName = [BaseName, '_GaussImage.png'];
      imwrite(GaussIm, fullfile(PlotSaveDir2, FileName))
   end
   if ShowPlots
      GaussImFigure = figure();
      GaussImAxes = axes(GaussImFigure);
      imshow(GaussIm, [], 'Parent', GaussImAxes)
   end
end

if ismember("HistIm", PlotDo)
   % Histogram image
   [HistIm] = smi_vis.GenerateImages.histogramImage(SMD, obj.SRImageZoom);
   if ~isempty(PlotSaveDir2)
      FileName = [BaseName, '_HistImage.png'];
      imwrite(single(HistIm), fullfile(PlotSaveDir2, FileName))
   end
   if ShowPlots
      HistImFigure = figure();
      HistImAxes = axes(HistImFigure);
      imshow(HistIm, [], 'Parent', HistImAxes)
   end
end

if ismember("CircleIm", PlotDo)
   % Generate a circle image.
   [~, CircleImageRGB] = smi_vis.GenerateImages.circleImage(...
      SMD, [], obj.SRCircImZoom);
   if ~isempty(PlotSaveDir2)
      FileName = [BaseName, '_CircleImage.png'];
      imwrite(CircleImageRGB, fullfile(PlotSaveDir2, FileName))
   end
   if ShowPlots
      CircleImFigure = figure();
      CircleImAxes = axes(CircleImFigure);
      imshow(CircleImageRGB, [], 'Parent', CircleImAxes)
   end
end

if ismember("CircleImDrift", PlotDo)
    % Generate a circle image.
    ColorMap = parula(SMD.NFrames * SMD.NDatasets);
    [~, CircleImageRGB] = smi_vis.GenerateImages.circleImage(...
        SMD, ColorMap(SMD.NFrames*(SMD.DatasetNum-1) + SMD.FrameNum, :), ...
        obj.SRCircImZoom);
    if ~isempty(PlotSaveDir2)
       FileName = [BaseName, '_CircleImageDrift.png'];
       imwrite(CircleImageRGB, fullfile(PlotSaveDir2, FileName))
   end
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
          if ~isempty(PlotSaveDir2)
             FileName = [BaseName '_' regexprep(HistName,' ','_') '_Hist.png'];
             saveas(FigH,(fullfile(PlotSaveDir2,FileName)),'png');
          end
          if ~ShowPlots; close(gcf); end
       end
    end

    % nested function for plotting and saving the PValue histogram, which is
    % plotted specially with log axes.
    function plotAndSavePValueHist(FieldName,HistName)
       if isfield(SMD,FieldName) && ~isempty(SMD.(FieldName))
          Vector_in=SMD.(FieldName);
          % Convert data values near zero to slightly larger values so that
          % they show up in the histogram plot below.
          %Vector_in(Vector_in <= 0.001) = 0.001;
          FigH = smi_vis.GenerateImages.plotHistogram(Vector_in,HistName);
          %set(gca, 'xscale', 'log');
          set(gca, 'yscale', 'log');
          xlim([0, 1]);

          [~,BaseName,~] = fileparts(obj.SMF.Data.FileName{1});
          if ~isempty(PlotSaveDir2)
             FileName = [BaseName '_' regexprep(HistName,' ','_') '_Hist.png'];
             saveas(FigH,(fullfile(PlotSaveDir2,FileName)),'png');
          end
          if ~ShowPlots; close(gcf); end
       end
    end

    % nested function for plotting and saving cumulative Drift
    function plotAndSaveCum(FieldName,CumName)
       if isfield(SMD,FieldName) && ~isempty(SMD.(FieldName))
          FigH = smi_core.DriftCorrection.plotCumDrift(SMD,FieldName);
          [~,BaseName,~] = fileparts(obj.SMF.Data.FileName{1});
          if ~isempty(PlotSaveDir2)
             FileName = [BaseName '_' CumName '_Cum.png'];
             saveas(FigH,(fullfile(PlotSaveDir2,FileName)),'png');
          end
          if ~ShowPlots; close(gcf); end
       end
    end

end
