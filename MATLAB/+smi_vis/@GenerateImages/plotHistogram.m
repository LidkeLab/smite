function [FigHandle] = plotHistograms(Vector_in, Hist_Name)
%plotHistogram creates a histogram of a field from an SMD structure.
% 
% INPUT:
%    Vector_in   Vector of values to create histogram. 
%    Hist_Name   Name of histogram graph (Default=[])
%
% OUTPUT:
%    FigHandle:  Figure handle of histogram plot    

% Created by
%    Hanieh Mazloom-Farsibaf   Apr 2017 (Keith A. Lidke's lab)

if nargin<1
    error('Inter a vector format input for histogram')
end

if nargin<2 || isempty(Hist_Name)
    Hist_Name=inputname(1);
end

Vector_in = single(Vector_in);
Nbin=single(50); % number of bins in histogram

% Plot the histogram of Hist_Name
FigHandle=figure;
hist(Vector_in,Nbin);
xlabel(sprintf('%s',Hist_Name));
ylabel('Frequency');
title(sprintf('Histogram of %s',Hist_Name));

end
