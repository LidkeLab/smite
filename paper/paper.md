---
title: 'SMITE: Single Molecule Imaging Toolbox Extraordinaire (MATLAB)'
tags:
  - MATLAB
  - single molecule localization microscopy (SMLM)
  - single particle tracking (SPT)
  - super resolution
authors:
  - name: David J. Schodt
    orcid: 0000-0002-8986-2736
    affiliation: 1
    equal-contrib: true
  - name: Michael J. Wester
    orcid: 0000-0002-3520-7605
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
    equal-contrib: true
  - name: Mohamadreza Fazel
    orcid: 0000-0002-6215-1336
    affiliation: 1
  - name: Sajjad Khan
    orcid: 0000-0002-6910-5199
    affiliation: 1
  - name: Hanieh Mazloom-Farsibaf
    orcid: 0000-0002-2571-0418
    affiliation: 1
  - name: Sandeep Pallikkuth
    orcid: 0009-0003-0400-6389
    affiliation: 1
  - name: Marjolein B. M. Meddens
    orcid: 0000-0002-9965-1342
    affiliation: 1
  - name: Farzin Farzam
    orcid: 0000-0002-9939-1923
    affiliation: 1
  - name: Eric A. Burns
    orcid: 0000-0002-1625-2400
    affiliation: 3
  - name: William K. Kanagy
    orcid: 0000-0002-5756-9965
    affiliation: 3
  - name: Derek A. Rinaldi
    orcid: 0009-0000-8394-3626
    affiliation: 3
  - name: Elton Jhamba
    orcid: 0000-0002-5272-6466
    affiliation: 3
  - name: Sheng Liu
    orcid: 0000-0003-1225-0763
    affiliation: 1
  - name: Peter K. Relich
    orcid: 0000-0002-6063-6233
    affiliation: 1
  - name: Mark J. Olah
    affiliation: 1
  - name: Stanly L. Steinberg
    affiliation: 2
  - name: Keith A. Lidke^[corresponding author]
    orcid: 0000-0002-9328-4318
    affiliation: 1
    corresponding: true
affiliations:
 - name: Department of Physics and Astronomy, University of New Mexico
   index: 1
 - name: Department of Mathematics and Statistics, University of New Mexico
   index: 2
 - name: Department of Pathology, University of New Mexico Health Sciences Center
   index: 3
date: 29 September 2022
bibliography: paper.bib

---

# Summary

Fluorescence single molecule imaging comprises a variety of techniques that
involve detecting individual fluorescent molecules.  Many of these techniques
involve localizing individual fluorescent molecules with precisions below the
diffraction limit, which limits the spatial resolution of (visible) light-based
microscopes.  These methodologies are widely used to image biological
structures at the nanometer scale by fluorescently tagging the structures of
interest, elucidating details of the biological behavior observed.

Two common techniques are single-molecule localization microscopy (SMLM),
[@Lidke_article:2005; @Betzig_article:2006; @Rust_article:2006;
@Hell_article:2007; @vandeLinde_article:2011; @Fazel_article:2022]
which is used to produce 2D or 3D super-resolution images of static or nearly
static structures, and single-particle tracking (SPT) [@Shen_article:2017],
which follows the time course of one or a very small number of moving tagged
molecules.  SMLM often involves distributions of particles at medium to high
density, while SPT works in a very low density domain.  These procedures all
require intensive numerical computation, and the methods are tightly
interwoven.

# Statement of need

SMITE is a MATLAB-based toolbox that provides analysis tools for fluorescence
single molecule imaging with an emphasis on single molecule localization
microscopy (SMLM) and single-particle tracking (SPT).  The SMITE toolbox
consists of a MATLAB infrastructure with some C and CUDA code embedded to
provide CPU/GPU speed-ups for particularly expensive computations.
The source code for SMITE has been archived to GitHub:
[https://github.com/LidkeLab/smite](https://github.com/LidkeLab/smite)

SMITE is designed around the concept that a parameter structure, the Single
Molecule Fitting (SMF) structure, uniquely and completely defines the data
analysis.  The results are completely contained in a Single Molecule Data (SMD)
structure.  SMITE is designed to make lowest-level tools just as easy to use as
the higher-level application-specific classes.  All tools make use of the SMF
and SMD structures.  SMITE is organized into a set of namespaces that group
similar tools and concepts.  The namespace  `+smi`  contains the highest level
tools that will be the most common entry point for processing SMLM and SPT data
sets. 

Code coverage includes mature SMLM data analysis techniques (applying gain and 
offset corrections to raw data, finding localizations, thresholding
localizations based on various criteria, frame connection and drift
correction), SMLM/SPT simulations, sophisticated SPT analyses, post-processing
clustering and statistical analyses (e.g., diffusion analysis and hidden Markov
models for characterizing dimers in SPT results), a variety of visualizations,
experimental point spread function creation and characterization, all sprinkled
with various examples of usage.  Interaction with these tools is via GUIs or
scripting.  See \autoref{fig:smite_overview} for several examples of SMITE
GUIs.

SMITE is a tool designed to be used by researchers and upper level students
interested in fluorescence single molecule imaging and applications.
Some of the algorithms have already been published: 2D Gaussian blob maximum
likelihood estimate [@Smith_article:2010], frame connection
[@Schodt_article:2021], drift correction [@Wester_article:2021], Bayesian
grouping of localizations [@Fazel_article:2022a], diffusion estimation
[@Relich_article:2016].  However, this is the first time that they have been
integrated together, sharing common data structures.
Applications are described in [@FrancoNitta_article:2021;
@Mazloom-Farsibaf_article:2021; @Bailey_article:2022; @Schodt_article:2023].
Typical raw image data can be found in [@Pallikkuth_data:2018].
A summary of the namespaces and classes in SMITE can be found in the online
documentation at
[https://github.com/LidkeLab/smite/blob/main/doc/SMITEclasses.md](https://github.com/LidkeLab/smite/blob/main/doc/SMITEclasses.md).

SMAP [@Ries_article:2020], an alternative MATLAB integrated SMLM/SPT code, is
GUI oriented, while SMITE was designed to be more focused on scripting
(although many GUIs are available as well) in order to make batch processing
extremely simple.  SMITE, in addition, is designed to operate with HDF5
(Hierarchical Data Format) files which efficiently store very large datasets,
while SMAP preferentially works with TIFF formatted files.  Both SMITE and SMAP
work with separate software to control instruments, MATLAB Instrument Control
(MIC) [@Pallikkuth_article:2018] and Micro-Manager [@Edelstein_article:2014],
respectively.

![SMITE GUIs for (upper left) making movies from SPT trajectories, (upper
right) SMLM analysis, (lower left) channel registration, and (lower right)
inspection of results contained in SMD
structures.\label{fig:smite_overview}](smite_overview.pdf){ width=100% }

# Author Contributions

KAL conceived and supervised development of SMITE and its predecessors.
DJS, MF, HMF, MBMM and KAL coded SR localization techniques.
SP, MJW, MF and HMF implemented thresholding.
DJS and HMF wrote frame connection.
MJW, KAL and DJS developed drift correction based on an earlier version by FF.
MF and KAL developed BaGoL; MJW wrote the interface to SMITE.
DJS, HMF, WKK, DAR and EJ designed and wrote code for single particle tracking
based on ideas from PKR.
DJS, MBMM, HMF and SP created visualizations.
SK, MJW and DJS added SMLM and SPT simulations.
KAL, SL and MJW developed Zernike polynomial point spread function engineering.
MJW packaged clustering techniques, and
SLS and MJW wrote code for various cluster statistics.
DJS added channel registration, various statistics, dimer hidden Markov
modeling and the batch Publish class.
PKR and MJO developed the diffusion estimator, while MJO also contributed to
the change detector.
MJW, SK and EAB added GitHub documentation.
MJW and DJS wrote the manuscript.
All authors reviewed the manuscript.

# Acknowledgements

This work was supported by NIH grants
K12GM088021 (ASERT-IRACDA program),
NCI P30CA118100,
NCI R01CA248166,
NIGMS R01GM109888,
NIGMS R21GM104691,
NIGMS R21GM132716,
NIGMS R35GM126934,
NCRR RR024438,
NIGMS 1R01GM140284,
NIBIB 1R21EB019589,
NIGMS 5P50GM085273 (New Mexico Spatiotemporal Modeling Center);
NSF grants 0954836 and 2039517; DOE grant DE-SC0019267;
support from the University of New Mexico Office of the Vice President for
Research Program for Enhancing Research Capacity; and
supported by grants from NVIDIA and utilized an NVIDIA A6000 GPU.

# References
