## BaGoL (Bayesian Grouping of Localizations)

This [example](BaGoL_EGFR_dSTORM.m) is slightly modified (`EGF = BaGoL;`
replaced by `EGF = smi.BaGoL;`) from the BaGoL distribution example of
the same name found at
[https://github.com/LidkeLab/BaGoL](https://github.com/LidkeLab/BaGoL)
(under Software).  Two other examples of running BaGoL and expected results
can be found there as well.  This software can also be found in *Code Ocean*:

> Mohamadreza Fazel, Michael J. Wester, David J. Schodt, Sebastian Restrepo
Cruz, Sebastian Strauss, Florian Schueder, Thomas Schlichthaerle, Jennifer M.
Gillette, Diane S. Lidke, Bernd Rieger, Ralf Jungmann, Keith A. Lidke,
"BaGoL (Bayesian Grouping of Localizations) [Source Code]",
*Code Ocean*, November 10, 2022,
[https://codeocean.com/capsule/de6769fd-f009-492d-9055-4c694848b128/](https://codeocean.com/capsule/de6769fd-f009-492d-9055-4c694848b128/)
(DOI: 10.24433/CO.3605166.v1).

Details on [BaGoL properties](../../MATLAB/+smi/@BaGoL/README.md) used in
the above example are described more fully, along with a list of BaGoL methods.

In addition,
[MATLAB/examples/hierBaGoL_wrapper](../../MATLAB/examples/hierBaGoL_wrapper.m)
is a script for processing multiple hierarchical BaGoL datasets or
splitting a single dataset into multiple previously defined ROIs.
See the [summary](hierBaGoL_wrapperSummary.md) for further details.

---

A full reference to the BaGoL paper:

> Mohamadreza Fazel, Michael J. Wester, David J. Schodt, Sebastian Restrepo
Cruz, Sebastian Strauss, Florian Schueder, Thomas Schlichthaerle, Jennifer M.
Gillette, Diane S. Lidke, Bernd Rieger, Ralf Jungmann and Keith A. Lidke,
"High-Precision Estimation of Emitter Positions using Bayesian Grouping of
Localizations" *Nature Communications*, Volume 13, Number 7152,
November 22, 2022, 1--11,
[https://www.nature.com/articles/s41467-022-34894-2](https://www.nature.com/articles/s41467-022-34894-2)
(DOI: 10.1038/s41467-022-34894-2).

### Abstract

Single-molecule localization microscopy super-resolution methods rely
on stochastic blinking/binding events, which often occur multiple times
from each emitter over the course of data acquisition. Typically, the
blinking/binding events from each emitter are treated as independent
events, without an attempt to assign them to a particular emitter. Here,
we describe a Bayesian method of inferring the positions of the
tagged molecules by exploring the possible grouping and combination
of localizations from multiple blinking/binding events. The results are
position estimates of the tagged molecules that have improved localization
precision and facilitate nanoscale structural insights. The Bayesian
framework uses the localization precisions to learn the statistical
distribution of the number of blinking/binding events per emitter and
infer the number and position of emitters. We demonstrate the method on
a range of synthetic data with various emitter densities, DNA origami
constructs and biological structures using DNA-PAINT and dSTORM data. We
show that under some experimental conditions it is possible to achieve
sub-nanometer precision.
