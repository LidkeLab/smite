### +smi_psf/@Zernike

Zernike includes low-level Zernike polynomial static functions.

SEE ALSO:
- smi_psf.PointSpreadFunction

---

static methods:
- NollIndex = **[zNM2Noll](zNM2Noll.m)**(N, M):
  Convert Zernike indices (n, m) into l (Noll ordering)
- l = **[zNM2Wyant](zNM2Wyant.m)**(n, m):
  Convert Zernike indices (n, m) into l (Wyant ordering)
- nZ = **[zNZNoll](zNZNoll.m)**(n):
  Number of Zernike polynomials up through order n (Noll ordering)
- nZ = **[zNZWyant](zNZWyant.m)**(n):
  Number of Zernike polynomials up through order n (Wyant ordering)
- names = **[zNamesNoll](zNamesNoll.m)**(ll):
  Classical names for Zernike indices ll (Noll ordering)
- names = **[zNamesWyant](zNamesWyant.m)**(ll):
  Classical names for Zernike indices ll (Wyant ordering)
- [N, M] = **[zNoll2NM](zNoll2NM.m)**(NollIndex):
  Convert Zernike index l into (n, m) (Noll ordering)
- l_max = **[zProperNollIndex](zProperNollIndex.m)**(l_max):
  Make sure that l_max includes both azimuthal terms (cos + sin) for a
  particular n and m.
- [n, m] = **[zWyant2NM](zWyant2NM.m)**(l):
  Convert Zernike index l into (n, m) (Wyant ordering)

unit test:
- success = **[unitTest](unitTest.m)**():
  Test functionality of the Zernike polynomial portion of the class
