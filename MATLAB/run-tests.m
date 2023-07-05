% Run various tests on the core functionality of SMITE.  Much output will be
% saved in tempdir/smite/unitTest/name_of_test.  ExpectedResults are provided
% in the directory in which run_tests.m resides where very large files have
% been deleted so as to not bloat up the the SMITE distribution.

% +smi

fprintf('smi.SMLM.unitTest\n');
try
   smi.SMLM.unitTest()
end

fprintf('smi.SPT.unitTestFFGC (frame-to-frame and gap closing processes)\n');
try
   smi.SPT.unitTestFFGC()
end
