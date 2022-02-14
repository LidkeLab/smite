function [PlotFigure, DisplayParams] = plotSepDistribs(PlotFigure, ...
    TRArray, PDFHandles, PDFInputs, DisplayParams)
%plotSepDistribs makes histograms of the distributions of separations.
% This function makes histograms of the separation distributions for each
% state of the HMM.  This should be run after the Viterbi sequences have
% been stored in TRArray.  A model of the expected separation distributions
% is overlain on the observed histograms.
%
% WARNING: For now, this is somewhat of a "skeleton" in that I'm only
%          plotting the distribution for the dimer state.  DJS 22/1/20
%
% INPUTS:
%   PlotFigure: Figure in which we'll make the plots.
%   TRArray: Array of tracking results pairs that have there state sequence
%            field population.
%   PDFHandles: Array of handles to the state emission PDFs (see
%               smi_stat.HMM.generateEmissionPDFs())
%   PDFInputs: Cell array of inputs needed for PDFHandles. Indices should
%              match the ordering defined in
%              smi_stat.HMM.generateEmissionPDFs().  Note that some of
%              these will be populated internally based on values in
%              TRArray.
%   DisplayParams: A structure of display parameters for the plots.
%
% OUTPUTS:
%   PlotFigure: Figure containing the plots.
%   DisplayParams: A structure of display parameters for the plots.

% Created by:
%   David J. Schodt (Lidke Lab, 2022)


% Isolate the separations for each state.
% NStates = numel(PDFHandles);
NStates = 1; % only focus on dimer state, see WARNING above
NPairs = size(TRArray, 1);
Separation = cell(NStates, 1);
for nn = 1:NStates
    for ii = 1:NPairs
        Separation{nn} = [Separation{nn};
            TRArray(ii, 1).Separation(TRArray(ii, 1).StateSequence == nn)];
    end
end

% Define the model for the expected distribution of separations.
Model = cell(NStates, 1);
for nn = 1:NStates
    % Loop over each trajectory pair and add it's contribution to the model
    % distribution.
    SepModel = linspace(min(Separation{nn}), max(Separation{nn}), ...
        numel(Separation{nn})).';
    PDFInputs{1} = SepModel;
    Model{nn} = zeros(size(SepModel));
    for ii = 1:NPairs
        % Isolate some needed arrays from TRArray.
        % NOTE: We typically want TRArray(ii, 2).RegError because channel 1
        %       is used as the reference channel.
        KeepBool1 = (TRArray(ii, 1).StateSequence == nn);
        KeepBool2 = (TRArray(ii, 2).StateSequence == nn);
        AverageSE = [TRArray(ii, 1).AverageSE(KeepBool1), ...
            TRArray(ii, 2).AverageSE(TRArray(ii, 2).StateSequence == nn)];
%         DeltaT = double(diff(sort(TRArray(ii, 1).FrameNum(KeepBool1))));
        RegError = TRArray(ii, 2).RegError(KeepBool2);
%         D = [TRArray(ii, 1).DiffusionCoefficient(KeepBool1), ...
%             TRArray(ii, 2).DiffusionCoefficient(KeepBool2)];

        % Loop over each observation and add its contribution.
        for jj = 1:numel(RegError)
            PDFInputs{2} = AverageSE(jj, :);
%             PDFInputs{3} = DeltaT(jj); % not used, see WARNING above
            PDFInputs{4} = RegError(jj);
%             PDFInputs{6} = D(jj, :); % not used, see WARNING above
            Model{nn} = Model{nn} ...
                + PDFHandles{nn}(PDFInputs)/numel(RegError);
        end
    end
    Model{nn} = Model{nn} / NPairs;

    % Plot the model over a histogram of the observations.
    PlotAxes = subplot(NStates, 1, nn, 'Parent', PlotFigure);
    histogram(PlotAxes, Separation{nn}, 'Normalization', 'pdf')
    hold(PlotAxes, 'on')
    plot(PlotAxes, SepModel, Model{nn})
    xlabel(PlotAxes, 'Separation (pixels)')
    ylabel(PlotAxes, 'PDF(Separation)')
    legend(PlotAxes, {'Observed', 'Model'}, 'Location', 'best')
    title(PlotAxes, DisplayParams.TitleString)
end


end