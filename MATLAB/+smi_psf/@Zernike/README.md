### +smi_psf/@Zernike

Zernike includes low-level Zernike polynomial static functions.

SEE ALSO:
- smi_psf.PointSpreadFunction

static methods:
- NollIndex = **zNM2Noll**(N, M):
  Convert Zernike indices (n, m) into l (Noll ordering)
- l = **zNM2Wyant**(n, m):
  Convert Zernike indices (n, m) into l (Wyant ordering)
- nZ = **zNZNoll**(n):
  Number of Zernike polynomials up through order n (Noll ordering)
- nZ = **zNZWyant**(n):
  Number of Zernike polynomials up through order n (Wyant ordering)
- names = **zNamesNoll**(ll):
  Classical names for Zernike indices ll (Noll ordering)
- names = **zNamesWyant**(ll):
  Classical names for Zernike indices ll (Wyant ordering)
- [N, M] = **zNoll2NM**(NollIndex):
  Convert Zernike index l into (n, m) (Noll ordering)
- l_max = **zProperNollIndex**(l_max):
  Make sure that l_max includes both azimuthal terms (cos + sin) for a
  particular n and m.
- [n, m] = **zWyant2NM**(l):
  Convert Zernike index l into (n, m) (Wyant ordering)

unit test:
- success = **unitTest**():
  Test functionality of the Zernike polynomial portion of the class
