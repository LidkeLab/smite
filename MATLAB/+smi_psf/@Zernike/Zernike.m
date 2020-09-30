classdef Zernike < handle

%Zernike includes low-level Zernike polynomial static functions.

% Created by
%    Michael Wester, 2020, Lidkelab.

% =============================================================================
methods(Static)

   NollIndex = zNM2Noll(N, M)
   l = zNM2Wyant(n, m)
   nZ = zNZNoll(n)
   nZ = zNZWyant(n)
   names = zNamesNoll(ll)
   names = zNamesWyant(ll)
   [N, M] = zNoll2NM(NollIndex)
   l_max = zProperNollIndex(l_max)
   [n, m] = zWyant2NM(l)
   unitTest()

end % methods(Static)
% =============================================================================

end % classdef Zernike
