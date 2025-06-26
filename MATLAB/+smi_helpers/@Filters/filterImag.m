function SMD = filterImag(SMD, Verbose)
%filterImag: Filter out _SE SMD fields with non-zero imaginary components.
%
% INPUT:
%    SMD     Single Molecule Data structure
%    Verbose verbosity flag [DEFAULT = false]
%
% OUTPUT:
%    SMD   modified Single Molecule Data structure

% Created by
%   Michael J. Wester (5/23/2025)

if ~exist('Verbose', 'var')
   Verbose = false;
end

n_prefilter = numel(SMD.X);

% Find SMD fields that contain an SE in their name.
f = fields(SMD);
SE = f(contains(f, '_SE'));
% Eliminate fields that do not contain n_prefilter elements.
inuse = arrayfun(@(i) numel(SMD.(SE{i})) == n_prefilter, 1:numel(SE));
SEinuse = SE(inuse);

% Check the fields for nonzero imaginary values.
complex = [];
for i = 1 : numel(SEinuse)
   complex = union(complex, find(imag(SMD.(SE{i})) ~= 0));
end

noncomplex = setdiff(1:n_prefilter, complex);
SMD = smi_core.SingleMoleculeData.isolateSubSMD(SMD, noncomplex);

if Verbose >= 2 && numel(SMD.X) ~= n_prefilter
   fprintf('   Complex standard errors removed = %d out of %d\n', ...
           n_prefilter - numel(SMD.X), n_prefilter);
end

end
