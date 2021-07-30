function [PlotAxes] = plotEnsembleMSD(PlotAxes, ...
    MSDEnsemble, DiffusionStruct, DiffusionModel, UnitFlag)
%plotEnsembleMSD plots an ensemble MSD and an associated fit.
% This method will plot the MSD data in MSDEnsemble as well as the fit
% information provided in DiffusionStruct.
%
% INPUTS:
%   PlotAxes: Axes in which the plot will be made. (Default = gca())
%   MSDEnsemble: Structure array with ensemble MSD results (see outputs of
%                computeMSD())
%   DiffusionStruct: Structure array containing MSD fit ressults. 
%                    (Default = [], meaning no fit results are plotted).
%   DiffusionModel: A string specifying the diffusion model to fit to the
%                   MSD. See options in DiffusionEstimator class property
%                   'DiffusionModel'. (Default = 'brownian1c')
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
    DiffusionModel = 'brownian1c';
end
if (~exist('UnitFlag', 'var') || isempty(UnitFlag))
    UnitFlag = 0;
end

% Plot the MSD.
FrameRate = DiffusionStruct(2).FrameRate;
PixelSize = DiffusionStruct(2).PixelSize;
FrameConversion = ~UnitFlag + UnitFlag/FrameRate;
JumpConversion = ~UnitFlag + UnitFlag*PixelSize;
FrameLags = MSDEnsemble.FrameLags * FrameConversion;
plot(PlotAxes, FrameLags, MSDEnsemble.MSD*(JumpConversion^2), '.')
hold(PlotAxes, 'on')

% If needed, plot the fit results.
if ~isempty(DiffusionStruct)
    % Define some unit conversion parameters. This is needed because the
    % DiffusionStruct allows for either camera units or physical units
    % (unlike the MSD structures).
    IsCameraUnits = strcmpi(DiffusionStruct(2).Units, ...
        {'pixels'; 'frames'});
    JumpUnitConversion = IsCameraUnits(1)*JumpConversion ...
        + ~IsCameraUnits(1)*(UnitFlag + ~UnitFlag/PixelSize);
    TimeUnitConversion = IsCameraUnits(2)*FrameConversion ...
        + ~IsCameraUnits(2)*(UnitFlag + ~UnitFlag*FrameRate);
    
    % Plot the MSD fit.
    switch lower(DiffusionModel)
        case 'brownian1c'
            % The Brownian diffusion model suggests the MSD is linear with
            % time.
            FitParams = DiffusionStruct(2).FitParams ...
                * (JumpUnitConversion^2) ./ [1, TimeUnitConversion];
            plot(PlotAxes, ...
                FrameLags, FitParams(2)*FrameLags + FitParams(1))
    end
end
TimeUnit = smi_helpers.arrayMUX({'frames', 'seconds'}, UnitFlag);
MSDUnit = smi_helpers.arrayMUX(...
    {'pixels^2', 'micrometers^2'}, UnitFlag);
xlabel(PlotAxes, sprintf('Time lag (%s)', TimeUnit))
ylabel(PlotAxes, sprintf('MSD (%s)', MSDUnit))
legend(PlotAxes, {'MSD', 'Fit'}, 'Location', 'best')


end