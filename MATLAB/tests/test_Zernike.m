% Run various tests on the core functionality of SMITE.  Much output will be
% saved in tempdir/smite/unitTest/name_of_test.  ExpectedResults are provided
% in the directory in which run_tests.m resides where very large files have
% been deleted so as to not bloat up the the SMITE distribution.

% +smi_psf

fprintf('smi_psf.Zernike.unitTest\n');
smi_psf.Zernike.unitTest()
