function success = unitTest()
%unitTest Test functionality of the Zernike polynomial portion of the class.

% Created by
%    Michael Wester, 2017, Lidkelab.

   success = 0;

   Z = smi_psf.Zernike();

   fprintf('\n');
   l = Z.zNM2Wyant(3, 2);
   fprintf('Wyant: l = %d (%s)\n', l, char(Z.zNamesWyant(l)));
   [n, m] = Z.zWyant2NM(11);
   fprintf('Wyant: (n, m) = (%d, %d)\n', n, m);
   l = Z.zNM2Wyant(3, -2);
   fprintf('Wyant: l = %d (%s)\n', l, char(Z.zNamesWyant(l)));
   [n, m] = Z.zWyant2NM(12);
   fprintf('Wyant: (n, m) = (%d, %d)\n', n, m);
   l = Z.zNM2Wyant(3, 0);
   fprintf('Wyant: l = %d (%s)\n', l, char(Z.zNamesWyant(l)));
   [n, m] = Z.zWyant2NM(15);
   fprintf('Wyant: (n, m) = (%d, %d)\n', n, m);

   fprintf('\n');
   l = Z.zNM2Noll(3, 3);
   fprintf('Noll:  l = %d (%s)\n', l, char(Z.zNamesNoll(l)));
   [n, m] = Z.zNoll2NM(10);
   fprintf('Noll:  (n, m) = (%d, %d)\n', n, m);
   l = Z.zNM2Noll(3, -3);
   fprintf('Noll:  l = %d (%s)\n', l, char(Z.zNamesNoll(l)));
   [n, m] = Z.zNoll2NM(9);
   fprintf('Noll:  (n, m) = (%d, %d)\n', n, m);
   l = Z.zNM2Noll(4, 0);
   fprintf('Noll:  l = %d (%s)\n', l, char(Z.zNamesNoll(l)));
   [n, m] = Z.zNoll2NM(11);
   fprintf('Noll:  (n, m) = (%d, %d)\n', n, m);

   nmax = 5;

   nZ_terms = Z.zNZWyant(nmax);
   fprintf('Wyant: %d total terms for nmax = %d\n', nZ_terms, nmax);

   nZ_terms = Z.zNZNoll(nmax);
   fprintf('Noll:  %d total terms for nmax = %d\n', nZ_terms, nmax);

   % Make sure that l_max includes both azimuthal terms (cos and sin) for a
   % particular n and m.
   l_max = 5;
   l_max_new = Z.zProperNollIndex(l_max);
   if l_max_new ~= l_max
      l_max = l_max_new;
      fprintf('\n!!!Resetting l_max to %d!!!\n', l_max);
   end

   success = 1;

   fprintf('Done.\n');

end
