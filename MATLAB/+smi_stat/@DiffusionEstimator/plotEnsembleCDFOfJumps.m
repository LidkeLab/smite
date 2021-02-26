function [PlotAxes] = plotEnsembleCDFOfJumps(PlotAxes, ...
    MSDEnsemble, DiffusionStruct, DiffusionModel, UnitFlag)
%plotEnsembleCDFOfJumps plots the CDF of jumps and an associated fit.
% This method will plot the CDF of jumps data provided in MSDEnsemble as
% well as the associated fit provided in DiffusionStruct.
%
% INPUTS:
%   PlotAxes: Axes in which the plot will be made. (Default = gca())
%   MSDEnsemble: Structure array with ensemble MSD results (see outputs of
%                computeMSD())
%   DiffusionStruct: Structure array containing MSD fit ressults. 
%                    (Default = [], meaning no fit results are plotted).
%   DiffusionModel: A string specifying the diffusion model to fit to the
%                   MSD. See options in DiffusionEstimator class property
%                   'DiffusionModel'. (Default = 'Brownian')
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
if (~exist('DiffusionModel', 'var') || isempty(DiffusionModel))
    DiffusionModel = 'Brownian';
end
if (~exist('UnitFlag', 'var') || isempty(UnitFlag))
    UnitFlag = 0;
end

% Plot the CDF of jumps.
JumpsConversion = UnitFlag*DiffusionStruct(2).PixelSize + ~UnitFlag;
Jumps = MSDEnsemble.SortedJumps * JumpsConversion;
stairs(PlotAxes, Jumps, MSDEnsemble.CDFOfJumps)

% If needed, plot the fit results.
if ~isempty(DiffusionStruct)
    hold(PlotAxes, 'on')
    switch DiffusionModel
        case {'Brownian', 'brownian'}
            % Plot the CDF of jumps expected for Brownian diffusion.
            ModelCDF = smi_stat.DiffusionEstimator.brownianJumpCDF(...
                DiffusionStruct(2).FitParams, MSDEnsemble.SortedJumps, ...
                MSDEnsemble.FrameLags, MSDEnsemble.NPoints);
            plot(PlotAxes, Jumps, ModelCDF)
    end
end
JumpsUnit = smi_helpers.stringMUX({'pixels', 'micrometers'}, UnitFlag);
xlabel(PlotAxes, JumpsUnit)
ylabel(PlotAxes, 'CDF of jumps')
legend(PlotAxes, {'CDF of jumps', 'Fit'}, 'Location', 'best')


end