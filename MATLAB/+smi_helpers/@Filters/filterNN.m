function SMD = filterNN(SMD, Verbose, n_NN, MedianMultiplier)
%filterNN:  Localizations are filtered based on the NND within MedianMultiplier
% times the median of the localization sigma, that is, localizations are
% eliminated if they do not have n_NN neighbors that are within Medianultiplier
% times the localization sigma median.
%
% Do not use on dSTORM data (set n_NN == 0).
%
% INPUTS:
%    SMD                Single Molecule Data structure
%    Verbose            verbosity flag [DEFAULT = false]
%    n_NN               minimum number of neighbors within
%                       MedianMultiplier*Prec_Median required to retain a
%                       localization [DEFAULT = 0]
%    MedianMultiplier   median multiplier for which localizations satisfying
%                       NND > MedianMultiplier*Prec_Median are removed
%                       [DEFAULT = 3]
%
% OUTPUT:
%    SMD                modified Single Molecule Data structure

% Created by
%    Mohamad Fazel and Michael J. Wester (5/25/2022)

if ~exist('Verbose', 'var')
   Verbose = false;
end

if ~exist('MedianMultiplier', 'var')
   MedianMultiplier = 3;
end

if ~exist('n_NN', 'var')
   n_NN = 0;
end

if n_NN > 0
   Prec_Median = median([SMD.X_SE; SMD.Y_SE]);

   % Original implementation.
%  [~,D]=knnsearch([SMD.X,SMD.Y],[SMD.X,SMD.Y],'K',length(SMD.X));
%  D(:,1)=[];
%  ID = D < 3*Prec_Median;
%  N = sum(ID,2);
%  Ind = N >= n_NN;

   % Less memory intensive and faster implementation.
   Prec_Median = median([SMD.X_SE; SMD.Y_SE]);
   % Find the number of nearest neighbors within MedianMultiplier*Prec_Median
   % for each localization.
   ID = rangesearch([SMD.X, SMD.Y], [SMD.X, SMD.Y], ...
                    MedianMultiplier*Prec_Median);
   % Number of neighbors, removing 1 which self counts the localization.
   N = cellfun(@numel, ID) - 1;
   Ind = N >= n_NN;
 
   % Only retain localizations that satisfy the above criteria.
   SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, Ind);

   if Verbose >= 2
      fprintf('Neighbor filtered localizations kept = %d out of %d\n',...
              sum(Ind), numel(Ind));
   end

   if sum(Ind) == 0
      error('No localizations kept!');
   end
end
