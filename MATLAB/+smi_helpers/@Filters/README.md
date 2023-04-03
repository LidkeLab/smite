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

static methods:
- **filterNonNeg**:
  filter out localizations with negative coordinates
- **filterIntensity**:
  filter localizations based on intensity
- **inflateSE**:
  inflate standard errors
- **filterFC**:
  filter out localizations representing fewer than nFC frame connections
- **filterNN**:
  localizations are filtered based on the nearest neighbor distance
