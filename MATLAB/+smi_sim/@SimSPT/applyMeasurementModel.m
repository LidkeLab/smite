function [TrajStruct, TrajStructModel] = ...
    applyMeasurementModel(TrajStruct, SimParams)
%applyMeasurementModel simulates measurement effects in trajectories.
% This method simulates some measurement effects (e.g., motion blur and
% missed localizations of visible emitters) for trajectories in
% 'TrajStructModel'.
%
% INPUTS:
%   TrajStruct: Structure containing trajectory data (see
%               smi_sim.SimSPT.simTrajectories())
%   SimParams: Structure of simulation parameters (see
%              smi_sim.SimSPT.defineDefaultParams())
%
% OUTPUTS:
%   TrajStruct: Input 'TrajStruct' after applying the measurement model.
%   TrajStructModel: Input 'TrajStruct' after motion blurring but before
%                    applying measurement noise.

% Created by:
%   David J. Schodt (Lidke Lab, 2021)


% Simulate motion blur by combining simulated subframes into a single
% frame.
NTraj = size(TrajStruct.Trajectories, 1);
NPadding = mod(SimParams.NFrames, SimParams.SubframeDensity);
if (SimParams.SubframeDensity > 1)
    % Grab a few arrays out of the TrajectoryStruct and pad them to ensure
    % each trajectory length is divisible by SubframeDensity.
    NaNPadArray = NaN(NTraj, NPadding);
    IsOnSub = [TrajStruct.IsOn, zeros(NTraj, NPadding, 'logical')];
    PhotonsSub = [TrajStruct.Photons, zeros(NTraj, NPadding)];
    XSub = [TrajStruct.Trajectories(:, :, 1), NaNPadArray];
    YSub = [TrajStruct.Trajectories(:, :, 2), NaNPadArray];
    
    % Apply the measurement model one trajectory at a time (it's a bit
    % easier this way, although probably slower).
    NaNArray = NaN(NTraj, SimParams.NFrames);
    IsOn = zeros(NTraj, SimParams.NFrames, 'logical');
    Photons = NaNArray;
    X = NaNArray;
    Y = NaNArray;
    XY_SE = NaNArray;
    for nn = 1:NTraj
        % Reshape arrays for this trajectory to assist in conversion from
        % subframes to frames.
        IsOnCurrent = reshape(IsOnSub(nn, :), ...
            SimParams.SubframeDensity, SimParams.NFrames).';
        PhotonsCurrent = reshape(PhotonsSub(nn, :), ...
            SimParams.SubframeDensity, SimParams.NFrames).';
        XYVarCurrent = (SimParams.PSFSigma^2) ./ PhotonsCurrent;
        XCurrent = reshape(XSub(nn, :), ...
            SimParams.SubframeDensity, SimParams.NFrames).';
        YCurrent = reshape(YSub(nn, :), ...
            SimParams.SubframeDensity, SimParams.NFrames).';
        
        % Perform the conversion from subframes to frames (approximate
        % motion blurring). For photons, we can just take the sum. The
        % positions and standard errors will be estimated as though each
        % subframe is another observation of the same Gaussian (like we do
        % in frame connection). Note that this probably isn't the best way
        % to simulate motion blur, as the smearing over subframes will
        % likely give us a higher standard error than we're estimating
        % here.
        IsOn(nn, :) = any(IsOnCurrent, 2);
        Photons(nn, :) = sum(PhotonsCurrent, 2);
        XY_SE(nn, :) = sqrt(1 ./ sum(1./XYVarCurrent, 2));
        X(nn, :) = sum(XCurrent./XYVarCurrent, 2) ...
            ./ sum(1./XYVarCurrent, 2);
        Y(nn, :) = sum(YCurrent./XYVarCurrent, 2) ...
            ./ sum(1./XYVarCurrent, 2);
    end
    TrajStructModel.IsOn = IsOn;
    TrajStructModel.Photons = Photons;
    TrajStructModel.Photons_SE = sqrt(Photons);
    TrajStructModel.Bg = SimParams.Bg * ones(size(Photons));
    TrajStructModel.Bg_SE = sqrt(TrajStructModel.Bg);
    TrajStructModel.Trajectories = cat(3, X, Y);
    TrajStructModel.Trajectories_SE = repmat(XY_SE, 1, 1, 2);
else
    TrajStructModel = TrajStruct;
    TrajStructModel.Photons_SE = sqrt(TrajStructModel.Photons);
    TrajStructModel.Bg_SE = sqrt(TrajStructModel.Bg);
    TrajStructModel.Trajectories_SE = ...
        repmat(SimParams.PSFSigma ./ sqrt(TrajStructModel.Photons), 1, 1, 2);
end
TrajStruct = TrajStructModel;
TrajStruct.Trajectories = TrajStruct.Trajectories ...
    + TrajStructModel.Trajectories_SE.*randn(NTraj, SimParams.NFrames, 2);
TrajStruct.Photons = poissrnd(TrajStruct.Photons);
TrajStruct.Bg = poissrnd(TrajStruct.Bg);


end