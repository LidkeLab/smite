function [EmissionMatrix] = ...
    generateEmissionMatrix(PDFHandles, DataCell)
%generateEmissionMatrix computes emission matrix for an observed separation
% This method will compute the emission probability densities of observing
% a separation given a state model (e.g. Dimer, Domain, Free, ...), where
% all model parameters and required data are contained in the DataCell cell
% array.
%
% INPUTS:
%   DataCell: A cell array containing all data/parameters need to compute
%             the emission probabilities for each state.
%             DataCell{1}: Separation in current frame (pixels)
%             DataCell{2}: Standard error in current frame (pixels)
%             DataCell{3}: Frames elapsed since previous observation
%             DataCell{4}: Sigma overlay, i.e., registration error (pixels)
%             DataCell{5}: Dimer separation (pixels)
%             DataCell{6}: Diffusion coefficient for each trajectory.
%                          (pixel^2 / frame)
%             DataCell{7}: Maximum separation between two trajectories
%                          (used in computing initial probability density
%                          of free state)(pixels)
%             DataCell{8:end}: Additional parameters/data, case dependent,
%                              e.g., DataCell{8} is sigma domain for the
%                              dimer, domain, free model
%   PDFHandles: A cell array of function handles, with each handle defining
%               the probability density of the observed separation for a
%               given model state.
%
% OUTPUTS:
%   EmissionMatrix: The probability densities of observing the input
%                   Seperation given a state model specified by the inputs.
%                   This is column vector, with row ii corresponding to the
%                   probability density of observing the input Separation
%                   given that the underlying state is state ii.
%                   (N x numel(PDFHandles) double)

% Created by:
%   David J. Schodt (Lidke Lab, 2020)


% Compute the emission probability densities based on the specified state
% model, storing each of the densities as elements in a column vector in
% EmissionMatrix.
NObservations = numel(DataCell{1});
NStates = numel(PDFHandles);
EmissionMatrix = zeros(NObservations, NStates, 'double');
for ii = 1:NStates
    EmissionMatrix(:, ii) = PDFHandles{ii}(DataCell);
end
EmissionMatrix(:, 1) = PDFHandles{1}(DataCell);


end