function [PlotAxes] = plotEnsembleMSD(PlotAxes, MSDStruct, DiffusionStruct)
%plotEnsembleMSD plots an ensemble MSD and an associated fit.
% This method will plot the MSD data in MSDStruct as well as the fit
% information provided in DiffusionStruct.
%
% INPUTS:
%   PlotAxes: Axes in which the plot will be made. (Default = gca())
%   MSDStruct: Structure array with MSD results (see outputs of
%              computeMSD())
%   DiffusionStruct: Structure array containing MSD fit ressults. 
%                    (Default = [], meaning no fit results are plotted).
%
% OUTPUTS:
%   PlotAxes: Axes in which the plot was made.

% Created by:
%   David J. Schodt (Lidke lab, 2021)


% Set defaults if needed.
if (~exist('PlotAxes', 'var') || isempty(PlotAxes))
    PlotAxes = gca();
end
if (~exist('DiffusionStruct', 'var') || isempty(DiffusionStruct))
    DiffusionStruct = [];
end

% Generate the plot.
plot(PlotAxes, MSDStruct.FrameLags, MSDStruct.MSD, '.')
hold(PlotAxes, 'on')
if ~isempty(DiffusionStruct)
    FitParams = DiffusionStruct(2).FitParams;
    FrameArray = MSDStruct.FrameLags([1, numel(MSDStruct.FrameLags)]);
    plot(PlotAxes, FrameArray, FitParams(1)*FrameArray + FitParams(2))
end


end