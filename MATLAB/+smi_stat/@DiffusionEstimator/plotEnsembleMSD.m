function [PlotAxes] = plotEnsembleMSD(PlotAxes, ...
    MSDStruct, DiffusionStruct, UnitFlag)
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
%   UnitFlag: Flag to specify camera units (0) or physical units (1).
%             (Default = 0)
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
if (~exist('UnitFlag', 'var') || isempty(UnitFlag))
    UnitFlag = 0;
end

% Generate the plot.
FrameConversion = UnitFlag/DiffusionStruct(2).FrameRate + ~UnitFlag;
MSDConversion = ~UnitFlag ...
    + UnitFlag*(DiffusionStruct(2).PixelSize^2);
plot(PlotAxes, MSDStruct.FrameLags*FrameConversion, ...
    MSDStruct.MSD*MSDConversion, '.')
hold(PlotAxes, 'on')
if ~isempty(DiffusionStruct)
    FitParams = DiffusionStruct(2).FitParams;
    FrameArray = MSDStruct.FrameLags([1, numel(MSDStruct.FrameLags)]);
    plot(PlotAxes, FrameArray*FrameConversion, ...
        MSDConversion * (FitParams(2)*FrameArray + FitParams(1)))
end
TimeUnit = smi_helpers.stringMUX({'frames', 'seconds'}, UnitFlag);
MSDUnit = smi_helpers.stringMUX(...
    {'pixels^2', 'micrometers^2'}, UnitFlag);
xlabel(PlotAxes, TimeUnit)
ylabel(PlotAxes, MSDUnit)
legend(PlotAxes, {'MSD', 'Fit'}, 'Location', 'best')

end