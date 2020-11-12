function generatePlots(obj, ShowPlots, PlotDo)
%generatePlots creates all histograms and plots for an SMD structure.
%
% INPUT:
%    obj         SMLM object
%       obj.SMD     Single Molecule Data structure
%       obj.SMF     Single Molecule Fitting structure
%    ShowPlots:  Flag for showing plots on the screen (Default = false)
%    PlotDo:     Plots to make chosen from the following list:
%                "Photons", "Bg", "PSFSigma", "Pvalue", "X_SE", "Y_SE", "Z_SE",
%                "NCombined", "DriftX", "DriftY", "DriftZ", "FitFrame",
%                "DriftIm", "GaussIm", "HistIm", "Drift"
%                (Default is to make all plots)
%                For example, PlotDo = ["Pvalue", "FitFrame", DriftIm"]
%
% OUTPUT:
%    The figures are saved in .png format in the given directory
%
% REQUIRES:
%    Dipimage toolbox (http://www.diplib.org/)

% Created by:
%    Hanieh Mazloom-Farsibaf, Marjolein Meddens Apr 2017 (Keith A. Lidke's lab)

fprintf('Generating output plots ...\n');

if ~exist('ShowPlots', 'var')
   ShowPlots = false;
end

if ~exist('PlotDo', 'var') || isempty(PlotDo)
   PlotDo = ["Photons", "Bg", "PSFSigma", "Pvalue", "X_SE", "Y_SE", "Z_SE", ...
             "NCombined", "DriftX", "DriftY", "DriftZ", "FitFrame", ...
             "DriftIm", "GaussIm", "HistIm", "Drift"];
end

% PlotSaveDir is the path to the directory for saving plots in .png format
PlotSaveDir = obj.SMF.Data.ResultsDir;
SMD = obj.SMD;

if matches("Photons", PlotDo)
   %create Photons histogram
   plotAndSaveHist('Photons','Intensity')
end

if matches("Bg", PlotDo)
   %create Bg histogram
   plotAndSaveHist('Bg','Background')
end

if matches("PSFSigma", PlotDo)
   plotAndSaveHist('PSFSigma','PSFSigma')
end

if matches("Pvalue", PlotDo)
   %create Pvalue histogram
   plotAndSaveHist('Pvalue','P value')
end

if matches("X_SE", PlotDo)
   %create X_SE histogram
   plotAndSaveHist('X_SE','X std error')
end

if matches("Y_SE", PlotDo)
   %create Y_SE histogram
   plotAndSaveHist('Y_SE','Y std error')
end

if matches("Z_SE", PlotDo)
   %create Z_SE histogram
   plotAndSaveHist('Z_SE','Z std error')
end

if matches("NCombined", PlotDo)
   %create Number of Connected localization histogram
   plotAndSaveHist('NCombined','Connected emitters')
end

%cumulative of DriftX, DriftY {, DriftZ} and total drift
if isfield(SMD,'DriftX') && ~isempty(SMD.DriftX) && ...
   isfield(SMD,'DriftY') && ~isempty(SMD.DriftY)
   if matches("DriftX", PlotDo)
      plotAndSaveCum('DriftX','DriftX')
   end
   if matches("DriftY", PlotDo)
      plotAndSaveCum('DriftY','DriftY')
   end
   if isfield(SMD,'DriftZ') && ~isempty(SMD.DriftZ) && matches("DriftZ",PlotDo)
      plotAndSaveCum('DriftZ','DriftZ')
   end 

   % Estimated 2D or 3D drift
   if matches("Drift", PlotDo)
      DC = smi_core.DriftCorrection;
      DC.PixelSizeZUnit = obj.SMF.DriftCorrection.PixelSizeZUnit;
      figDC = DC.plotDriftCorrection(SMD, 'A');
      FileName = [BaseName '_DriftCorrection'];
      saveas(gcf, fullfile(PlotSaveDir, FileName), 'png');
      if ~ShowPlots; close(gcf); end
   end
end

end

% BaseName is used to label plot files.
[~,BaseName,~] = fileparts(obj.SMF.Data.FileName{1});

if matches("FitFrame", PlotDo)
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

if matches("DriftIm", PlotDo)
   % Drift image
   [~, DriftImRGB] = smi_vis.GenerateImages.driftImage(SMD, SRImageZoom);
   dipshow(DriftImRGB);
   FileName = [BaseName '_DriftImage'];
   saveas(gcf, fullfile(PlotSaveDir, FileName), 'png');
   if ~ShowPlots; close(gcf); end
end

if matches("GaussIm", PlotDo)
   % Gaussian image
   [GaussIm] = smi_vis.GenerateImages.gaussianImage(SMD, SRImageZoom);
   dipshow(GaussIm);
   FileName = [BaseName '_GaussImage'];
   saveas(gcf, fullfile(PlotSaveDir, FileName), 'png');
   if ~ShowPlots; close(gcf); end
end

if matches("HistIm", PlotDo)
   % Histogram image
   [~, HistImRGB] = smi_vis.GenerateImages.histogramImage(SMD, SRImageZoom);
   dipshow(HistImRGB);
   FileName = [BaseName '_HistImage'];
   saveas(gcf, fullfile(PlotSaveDir, FileName), 'png');
   if ~ShowPlots; close(gcf); end
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
