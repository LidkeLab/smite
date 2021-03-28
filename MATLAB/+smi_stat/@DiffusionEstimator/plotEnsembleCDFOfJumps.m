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
%                   'DiffusionModel'. (Default = 'Brownian1C')
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
    DiffusionModel = 'Brownian1C';
end
if (~exist('UnitFlag', 'var') || isempty(UnitFlag))
    UnitFlag = 0;
end

% Plot the CDF of jumps.
FrameRate = DiffusionStruct(2).FrameRate;
PixelSize = DiffusionStruct(2).PixelSize;
FrameConversion = ~UnitFlag + UnitFlag/FrameRate;
JumpConversion = UnitFlag*PixelSize + ~UnitFlag;
Jumps = MSDEnsemble.SortedJumps * JumpConversion;
stairs(PlotAxes, Jumps, MSDEnsemble.CDFOfJumps)

% If needed, plot the fit results.
if ~isempty(DiffusionStruct)
    % Determine how many components are used in the model.
    switch lower(DiffusionModel)
        case 'brownian1c'
            NComponents = 1;
        case 'brownian2c'
            NComponents = 2;
        otherwise
            error('Unknown ''DiffusionModel'' = %s', DiffusionModel)
    end
    
    % Define some unit conversion parameters. This is needed because the
    % DiffusionStruct allows for either camera units or physical units
    % (unlike the MSD structures).
    IsCameraUnits = strcmpi(DiffusionStruct(2).Units, ...
        {'pixels'; 'frames'});
    JumpUnitConversion = IsCameraUnits(1)*JumpConversion ...
        + ~IsCameraUnits(1)*(UnitFlag + ~UnitFlag/PixelSize);
    TimeUnitConversion = IsCameraUnits(2)*FrameConversion ...
        + ~IsCameraUnits(2)*(UnitFlag + ~UnitFlag*FrameRate);
    
    % Plot the CDF fit.
    hold(PlotAxes, 'on')
    FitParams = DiffusionStruct(2).FitParams;
    FitParams(1:NComponents) = FitParams(1:NComponents) ...
        * (JumpUnitConversion^2) / TimeUnitConversion;
    ModelCDF = smi_stat.DiffusionEstimator.brownianJumpCDF(...
        FitParams, Jumps, ...
        MSDEnsemble.FrameLags * FrameConversion, ...
        MSDEnsemble.NPoints, ...
        MSDEnsemble.LocVarianceSum * (JumpConversion^2));
    plot(PlotAxes, Jumps, ModelCDF)
end
JumpsUnit = smi_helpers.stringMUX({'pixels', 'micrometers'}, UnitFlag);
xlabel(PlotAxes, sprintf('Jumps (%s)', JumpsUnit))
ylabel(PlotAxes, 'CDF of jumps')
legend(PlotAxes, {'CDF of jumps', 'Fit'}, 'Location', 'best')


end