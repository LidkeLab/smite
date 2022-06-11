---
title: 'SMITE: Single Molecule Imaging Toolbox Extraordinaire'
tags:
  - MATLAB
  - single molecule localization microscopy (SMLM)
  - single particle tracking (SPT)
  - super resolution
authors:
  - name: David J. Schodt^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-8986-2736
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Michael J. Wester^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-3520-7605
    affiliation: "1, 2"
  - name: Mohamadreza Fazel
    affiliation: 1
  - name: Hanieh Mazloom-Farsibaf
    orcid: 0000-0002-2571-0418
    affiliation: 1
  - name: Sandeep Pallikkuth
    affiliation: 1
  - name: Sajjad Khan
    affiliation: 1
  - name: Marjolein B. M. Meddens
    affiliation: 1
  - name: Farzin Farzam
    affiliation: 1
  - name: Will K. Kanagy
    affiliation: 3
  - name: Derek A. Rinaldi
    affiliation: 3
  - name: Elton Jhamba
    orcid: 0000-0002-5272-6466
    affiliation: 3
  - name: Keith A. Lidke^[corresponding author]
    orcid: 0000-0002-9328-4318
    affiliation: 1
affiliations:
 - name: Department of Physics and Astronomy, University of New Mexico
   index: 1
 - name: Department of Mathematics and Statistics, University of New Mexico
   index: 2
 - name: Department of Pathology, University of New Mexico Health Sciences Center
   index: 2
date: 20 April 2022
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

Two common techniques are single-molecule localization microscropy (SMLM),
which is used to produce 2D or 3D super-resolution images of static or nearly
static structures, and single-particle tracking (SPT), which follows the time
course of one or a very small number of moving tagged molecules.  SMLM often
involves distributions of particles at medium to high density, while SPT works
in a very low density domain.  These procedures all require intensive numerical
computation, and the methods are tightly interwoven.  The SMITE toolbox
consists of a MATLAB infrastructure with some C and CUDA code embedded to
provide CPU/GPU speed-ups for particularly expensive computations.

# Statement of need

SMITE is a MATLAB-based toolbox that provides analysis tools for fluorescence
single molecule imaging with an emphasis on single molecule localization
microscopy (SMLM) and single-particle tracking (SPT).

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
interested in fluorescence single molecule imaging and applications.  Parts of
it have already been or are in the process of being published, e.g., frame
connection [@Schodt_article:2021], drift correction [@Wester_article:2021],
Bayesian grouping of localizations [@Fazel_unpub:2019].  Applications are
described in [@Mazloom-Farsibaf_article:2021; @Bailey_article:2021].
Typical raw image data can be found in [@Pallikkuth_data:2018].

![SMITE GUIs for making movies from SPT trajectories, SMLM analysis, channel
registration, and inspection of results contained in SMD
structures.\label{fig:smite_overview}](smite_overview.pdf){ width=100% }

# References
