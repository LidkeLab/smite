---
title: "smite LLM Interactive Guide - Complete Documentation Index"
last_updated: "2025-10-11"
version: "1.0.0"
status: "complete"
total_documents: 57
---

# smite LLM Interactive Guide

**Single Molecule Imaging Toolbox Extraordinaire**

Complete LLM-optimized documentation covering SMLM, SPT, clustering, and advanced analysis workflows.

> **About This Documentation**
>
> This comprehensive guide includes 57 documents organized across 8 categories, designed for both human readers and AI-assisted development. Each document includes clear prerequisites, working code examples, and cross-references.
>
> **Coverage:** 100% of planned documentation (Phase 1 + Phase 2 complete)

## Quick Navigation

**Choose your learning path:**

- 🚀 [New to smite?](#getting-started) - Installation and first analysis
- 🔬 [SMLM Analysis](#smlm-workflow) - Super-resolution localization microscopy
- 🎯 [Particle Tracking](#spt-workflow) - Single particle tracking and diffusion
- 📦 [Batch Processing](#batch-workflow) - Multi-dataset analysis
- 🧬 [Advanced Techniques](#advanced-topics) - Clustering, 3D, multi-channel, BaGoL
- 📚 [API Reference](#api-reference) - Complete namespace documentation
- 🔧 [Troubleshooting](#troubleshooting) - Problem diagnosis and solutions

---

## Getting Started

**New users start here!**

| Document | Level | Summary |
|----------|-------|---------|
| [Quick Start](getting-started/quickstart.md) | Beginner | Get up and running in 5 minutes |
| [Installation](getting-started/installation.md) | Beginner | Complete setup guide (MATLAB, GPU, toolboxes) |
| [First Analysis](getting-started/first-analysis.md) | Beginner | Run your first SMLM analysis with GUI |
| [Understanding Results](getting-started/understanding-results.md) | Beginner | Interpret localization quality metrics |

**Time to complete:** 1-2 hours for full setup and first analysis

---

## Core Concepts

**Fundamental knowledge for working with smite**

| Document | Level | Summary |
|----------|-------|---------|
| [Architecture](core-concepts/architecture.md) | Intermediate | SMF/SMD paradigm and namespace organization |
| [SMF Structure](core-concepts/smf-structure.md) | Intermediate | Complete parameter structure reference |
| [SMD Structure](core-concepts/smd-structure.md) | Intermediate | Localization results structure |
| [TR Structure](core-concepts/tr-structure.md) | Intermediate | Trajectory-organized tracking results |
| [Coordinate System](core-concepts/coordinate-system.md) | Beginner | Pixel coordinates and transformations |
| [Data Flow](core-concepts/data-flow.md) | Intermediate | Raw data → localization → analysis pipeline |

---

## Workflows

**Complete analysis pipelines from start to finish**

| Document | Level | Summary |
|----------|-------|---------|
| [SMLM Analysis](workflows/smlm-analysis.md) | Intermediate | Full super-resolution workflow |
| [SPT Tracking](workflows/spt-tracking.md) | Intermediate | Single particle tracking pipeline |
| [BaGoL Clustering](workflows/bagol-clustering.md) | Advanced | Bayesian grouping of localizations |
| [Batch Processing](workflows/batch-processing.md) | Intermediate | Multi-dataset automated analysis |
| [Drift Correction](workflows/drift-correction.md) | Intermediate | Stage drift correction methods |
| [Frame Connection](workflows/frame-connection.md) | Intermediate | Linking blinking emitters (SMLM) |
| [Channel Registration](workflows/channel-registration.md) | Advanced | Multi-color alignment |
| [ROI Analysis](workflows/roi-analysis.md) | Intermediate | Region-based analysis |

---

## How-To Guides

**Task-specific practical guides**

### Data Management
- [Load Data](how-to/load-data.md) - Load .h5/.mat files (Beginner)
- [Save and Load SMD](how-to/save-and-load-smd.md) - Persist results (Beginner)
- [Export Results](how-to/export-results.md) - Export to CSV/JSON/HDF5 (Beginner)
- [Combine Datasets](how-to/combine-datasets.md) - Merge multiple SMD/TR (Intermediate)

### Configuration and Setup
- [Create SMF](how-to/create-smf.md) - Build parameter structures (Beginner)
- [Calibrate Camera](how-to/calibrate-camera.md) - Gain/offset/variance (Intermediate)
- [Use GPU](how-to/use-gpu.md) - GPU acceleration setup (Intermediate)
- [Tune Parameters](how-to/tune-parameters.md) - Optimize for best results (Intermediate)

### Analysis Operations
- [Localize Molecules](how-to/localize-molecules.md) - Run localization (Beginner)
- [Threshold Results](how-to/threshold-results.md) - Quality filtering (Beginner)
- [Visualize Results](how-to/visualize-results.md) - Create SR images (Beginner)
- [Work with ROIs](how-to/work-with-rois.md) - ROI operations (Beginner)
- [Simulate Data](how-to/simulate-data.md) - Generate test data (Beginner)

---

## Examples

**Complete working examples you can run**

| Document | Level | Summary |
|----------|-------|---------|
| [Basic Localization](examples/basic-localization.md) | Beginner | Simple localization with validation |
| [Tracking & Diffusion](examples/tracking-diffusion.md) | Intermediate | SPT with MSD analysis |
| [Clustering Analysis](examples/clustering-analysis.md) | Intermediate | DBSCAN and Voronoi clustering |
| [Multi-Channel](examples/multi-channel.md) | Advanced | Two-color SMLM with colocalization |
| [3D Localization](examples/3d-localization.md) | Advanced | Astigmatism-based 3D SMLM |

---

## API Reference

**Complete namespace and class documentation**

### Main Namespaces
- [+smi](api-reference/smi-namespace.md) - Top-level classes (SMLM, SPT, BaGoL, Publish)
- [+smi_core](api-reference/smi-core.md) - Core processing (LocalizeData, DriftCorrection, FrameConnection)
- [+smi_cluster](api-reference/smi-cluster.md) - Clustering algorithms (DBSCAN, Voronoi, H-SET)
- [+smi_stat](api-reference/smi-stat.md) - Statistical methods (DiffusionEstimator, HMM)
- [+smi_vis](api-reference/smi-vis.md) - Visualization tools (GenerateImages, GenerateMovies)
- [+smi_sim](api-reference/smi-sim.md) - Simulation classes (GaussBlobs, SimSMLM, SimSPT)
- [+smi_psf](api-reference/smi-psf.md) - PSF modeling (PointSpreadFunction, 3D)
- [+smi_helpers](api-reference/smi-helpers.md) - Utility functions

### Navigation
- [How to Read API Docs](api-reference/navigation.md) - Understanding method signatures and parameters

---

## Troubleshooting

**Common problems and solutions**

| Document | Level | Summary |
|----------|-------|---------|
| [Installation Issues](troubleshooting/installation-issues.md) | Beginner | Setup and PATH problems |
| [GPU Problems](troubleshooting/gpu-problems.md) | Intermediate | CUDA and GPU errors |
| [Memory Issues](troubleshooting/memory-issues.md) | Intermediate | Out of memory solutions |
| [Compilation Errors](troubleshooting/compilation-errors.md) | Advanced | Mex and CUDA compilation |
| [Localization Quality](troubleshooting/localization-quality.md) | Intermediate | Poor results diagnosis |
| [Common Mistakes](troubleshooting/common-mistakes.md) | Beginner | Frequently made errors |

---

## Reference

**Technical specifications and requirements**

| Document | Level | Summary |
|----------|-------|---------|
| [File Formats](reference/file-formats.md) | Intermediate | HDF5, MAT, calibration files |
| [Data Structures](reference/data-structures.md) | Intermediate | Complete SMF/SMD/TR field reference |
| [Dependencies](reference/dependencies.md) | Beginner | MATLAB toolboxes and requirements |
| [Compilation](reference/compilation.md) | Advanced | Building mex and CUDA files |
| [Testing](reference/testing.md) | Intermediate | Running unit tests |

---

## Learning Paths by Use Case

### SMLM Workflow
**Goal:** Analyze super-resolution SMLM data

1. [Installation](getting-started/installation.md)
2. [First Analysis](getting-started/first-analysis.md)
3. [SMLM Analysis](workflows/smlm-analysis.md)
4. [Drift Correction](workflows/drift-correction.md)
5. [Frame Connection](workflows/frame-connection.md)
6. [Visualize Results](how-to/visualize-results.md)

### SPT Workflow
**Goal:** Track particles and analyze diffusion

1. [Installation](getting-started/installation.md)
2. [Core Concepts: TR Structure](core-concepts/tr-structure.md)
3. [SPT Tracking](workflows/spt-tracking.md)
4. [Tracking & Diffusion Example](examples/tracking-diffusion.md)
5. [+smi_stat API](api-reference/smi-stat.md)

### Batch Workflow
**Goal:** Process multiple datasets automatically

1. [Batch Processing](workflows/batch-processing.md)
2. [Create SMF](how-to/create-smf.md)
3. [Combine Datasets](how-to/combine-datasets.md)
4. [+smi Publish Class](api-reference/smi-namespace.md#publish)

### Advanced Topics
**Goal:** Clustering, 3D, multi-channel analysis

1. [BaGoL Clustering](workflows/bagol-clustering.md)
2. [Clustering Analysis](examples/clustering-analysis.md)
3. [3D Localization](examples/3d-localization.md)
4. [Multi-Channel](examples/multi-channel.md)
5. [Channel Registration](workflows/channel-registration.md)

---

## Documentation Statistics

- **Total Documents:** 57
- **Categories:** 8
- **Complete Examples:** 5
- **API References:** 9
- **Troubleshooting Guides:** 6
- **Total Word Count:** ~150,000 words
- **Coverage:** 100% (all planned documentation complete)

## Phase Completion

- ✅ **Phase 1** (Essential): 12 documents - Complete
- ✅ **Phase 2** (Comprehensive): 45 documents - Complete
- 🔄 **Phase 3** (Maintenance): Ongoing updates

---

## Using This Documentation

### For LLMs and AI Assistants

This documentation is optimized for LLM consumption with:
- Structured frontmatter with metadata
- Clear prerequisites and related documents
- Complete working code examples
- Systematic organization by category
- Machine-readable manifest.json

**Recommended approach:**
1. Start with [manifest.json](manifest.json) for complete document catalog
2. Follow prerequisite chains for learning paths
3. Reference API docs for method signatures
4. Use examples for complete working code

### For Human Users

- **Browse by category** using the tables above
- **Follow learning paths** for structured learning
- **Search by topic** using tags in manifest.json
- **Jump to examples** for quick-start code

### Getting Help

- **Documentation issues:** [GitHub Issues](https://github.com/LidkeLab/smite/issues)
- **Installation help:** See [Installation Issues](troubleshooting/installation-issues.md)
- **Quality problems:** See [Localization Quality](troubleshooting/localization-quality.md)
- **General questions:** Check [Common Mistakes](troubleshooting/common-mistakes.md) first

---

## Document Access

- **Local:** `doc/llm-guide/` in your smite installation
- **GitHub:** https://github.com/LidkeLab/smite/tree/main/doc/llm-guide
- **Raw files:** https://raw.githubusercontent.com/LidkeLab/smite/main/doc/llm-guide/

---

## Acknowledgments

This LLM-optimized documentation was created to enable AI-assisted development with smite. All code examples are tested on GPU hardware and validated against the source code.

**smite Development Team:**
- Lidke Lab, University of New Mexico
- Contributors: See main repository

**Documentation Version:** 1.0.0 (Complete)
**Last Updated:** 2025-10-11
**smite Version Compatibility:** 1.2.6+
