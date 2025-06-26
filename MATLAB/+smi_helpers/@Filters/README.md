### +smi_helpers/@Filters

Filters operating on SMD (Single Molecule Data) structures used especially
for BaGoL analyses:

```
SR data -> remove localizations with negative coordinates
        -> intensity filter
        -> inflate standard errors
        -> frame connection, removing connections which involve only 1 frame
        -> NND filter --- Do not use on dSTORM data!
        -> BaGoL
```

---

static methods:
- **[filterNonNeg](filterNonNeg.m)**:
  filter out localizations with negative coordinates
- **[filterIntensity](filterIntensity.m)**:
  filter localizations based on intensity
- **[inflateSE](inflateSE.m)**:
  inflate standard errors
- **[filterFC](filterFC.m)**:
  filter out localizations representing fewer than nFC frame connections
- **[filterNN](filterNN.m)**:
  localizations are filtered based on the nearest neighbor distance
- **[filterImag](filterImag.m)**:
  filter out _SE SMD fields with non-zero imaginary components
