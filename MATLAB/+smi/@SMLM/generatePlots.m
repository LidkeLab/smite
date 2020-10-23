function generatePlots(obj, ShowPlots)
%generatePlots creates all histograms and plots for an SMD structure.
%
% INPUT:
%    obj         SMLM object
%       obj.SMD     Single Molecule Data structure
%       obj.SMF     Single Molecule Fitting structure
%    ShowPlots:  Flag for showing plots on the screen (Default = false)
%
% OUTPUT:
%    The figures are saved in given directory or only show the histogram
%    The figures are saved in .png format in the given directory
%    PlotSaveDir Path to save plot in .png format (Default=obj.SaveDir)

% Created by:
%    Hanieh Mazloom-Farsibaf, Marjolein Meddens Apr 2017 (Keith A. Lidke's lab)

fprintf('Generating output plots ...\n');

if ~exist('ShowPlots', 'var')
   ShowPlots = false;
end

PlotSaveDir = obj.SMF.Data.ResultsDir;
SMD = obj.SMD;

%create Photons histogram
plotAndSaveHist('Photons','Intensity')

%create Bg histogram
plotAndSaveHist('Bg','Background')

%create PSFSigma histogram
plotAndSaveHist('PSFSigma','PSFSigma')

%create Pvalue histogram
plotAndSaveHist('Pvalue','P value')

%create X_SE histogram
plotAndSaveHist('X_SE','X standard error')

%create Y_SE histogram
plotAndSaveHist('Y_SE','Y standard error')

%create Z_SE histogram
plotAndSaveHist('Z_SE','Z standard error')

%create Number of Connected localization histogram
plotAndSaveHist('NCombined','Connected emitters')

%cumulative of DriftX, DriftY
if isfield(SMD,'DriftX') && ~isempty(SMD.DriftX) && isfield(SMD,'DriftY') && ~isempty(SMD.DriftY)
        plotAndSaveCum('DriftX','Drift in X direction')
        plotAndSaveCum('DriftY','Drift in Y direction')
            if isfield(SMD,'DriftZ') && ~isempty(SMD.DriftZ)

        plotAndSaveCum('DriftZ','Drift in Z direction')
            end 
end

% Number of localizations per frame
for jj=1:max(SMD.DatasetNum)
    for ii=1:max(SMD.FrameNum)
        idx=find(SMD.FrameNum==ii & SMD.DatasetNum==jj);
        Nloc_frame{jj}(ii)=length(idx);
    end
end
if length(Nloc_frame)==1
    Fit_Frame = Nloc_frame{1};
end
for ii=1:length(Nloc_frame)-1
    Fit_Frame=cat(2,Nloc_frame{ii},Nloc_frame{ii+1});
end

% plot fit per frame
figure;
Frames=1:length(Fit_Frame);
plot(Frames,Fit_Frame);
xlabel('Frames');
ylabel('Number of Fits');
title('Number of fits per frame');
            [~,BaseName,~] = fileparts(obj.FileList{1});
FileName = [BaseName '_FitsPerFrame.png'];
if ~isempty(PlotSaveDir)
    saveas(gcf, fullfile(PlotSaveDir, FileName), 'png');
end
if ~ShowPlots; close(gcf); end

% nested function for plotting and saving histograms
    function plotAndSaveHist(FieldName,HistName)
        if isfield(SMD,FieldName) && ~isempty(SMD.(FieldName))
            Vector_in=SMD.(FieldName);
            FigH = obj.plotHistograms(Vector_in,HistName);
            [~,BaseName,~] = fileparts(obj.FileList{1});
            FileName = [BaseName '_' HistName '_Hist.png'];
            saveas(FigH,(fullfile(PlotSaveDir,FileName)),'png');
            if ~ShowPlots; close(gcf); end
        end
    end

% nested function for plotting and saving cumulative Drift
    function plotAndSaveCum(FieldName,CumName)
        if isfield(SMD,FieldName) && ~isempty(SMD.(FieldName))
            FigH = obj.plotCumDrift(SMD,FieldName);
            [~,BaseName,~] = fileparts(obj.FileList{1});
            FileName = [BaseName '_' CumName '_Cum.png'];
            saveas(FigH,(fullfile(PlotSaveDir,FileName)),'png');
            if ~ShowPlots; close(gcf); end
        end
    end
end
