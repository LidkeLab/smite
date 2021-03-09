function [FigHandle] = plotXYDriftParametric(SMD)
%plotXYDriftParametric makes a parametric plot of the x,y drift model.
% This method creates a plot of the inter-dataset linear drift model, 
% plotted as drift vectors for each dataset connected head to tail. 
%
% INPUT:
%    SMD: Single molecule results structure containing that results from
%         a previous single molecule analysis.
%
% OUTPUT:
%    FigHandle: Figure handle for the figure containing the parametric plots.

%Created by:
% David James Schodt (Lidke Lab, 2018)

% Check inputs.
if ~exist('SMD', 'var')
    error('You must enter an SMD structure to plot the drift model.')
end

% Create the plot figure.
FigHandle = figure;
hold('on');

% Extract relevant parameters from the SMD structure.
DriftX = SMD.DriftX; % pixels
DriftY = SMD.DriftY;
NDatasets = SMD.NDatasets;

% Compute the drift vectors for each dataset and plot them head to tail.
TailPosition = [0; 0]; % start the first drift vector at [0; 0]
for ii = 1:NDatasets
    % Compute the drift vector for the ii-th dataset.
    DriftVector = [DriftX(end, ii) - DriftX(1, ii); ...
        DriftY(end, ii) - DriftY(1, ii)];
    
    % Plot the drift vector aligned to the head of the previous dataset.
    quiver(TailPosition(1), TailPosition(2), ...
           DriftVector(1), DriftVector(2), 'b', 'AutoScale', 'off');
    
    % Update the tail position for the next iteration (this is just the
    % head position of the current drift vector).
    TailPosition = TailPosition + DriftVector;
end
title('Drift Model')
xlabel('X Drift Model (pixels)')
ylabel('Y Drift Model (pixels)')

end
